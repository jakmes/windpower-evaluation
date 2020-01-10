
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
set.seed(1)
library("quantreg")
library("splines")
library("qgam")
library("forecast")
data <- read.csv("Task15_W_Zone1_10/Task15_W_Zone1.csv")
names(data) <- c("zone_id", "timestamp", "p", "u10", "v10", "u100", "v100")
data$ws10 <- with(data, sqrt(u10^2 + v10^2))
data$ws100 <- with(data, sqrt(u100^2 + v100^2))
data$dir100 <- with(data, atan(u100/v100))
data <- na.omit(data)

train_length <- 8760L
train <- data[1:train_length,]
test <- tail(data, -train_length)



## Off-line: k-fold Cross-validation ####

# Number of CV folds and stratification
n_fold <- 6L
stratification <-2L
validation_is <- rep(rep(1:n_fold,each=floor(train_length/(stratification*n_fold))),stratification)

# Initialise container for results
cv_results <- data.frame(Fold=1:n_fold,
                         Benchmark=NA,
                         Proposed=NA)
model1_pred <- model2_pred <- rep(NA,train_length)

# Perform CV
for(fold in 1:n_fold){
  
  validation_i <- which(validation_is==fold)
  
  # Model 1: Wind Speed at 100m only
  model1 <- with(train[-validation_i,], lprq(ws100, p, h = 1, tau=0.5))
  model1_pred[validation_i] <- approx(model1$xx, model1$fv, train[validation_i,]$ws100)$y
  cv_results[fold,2] <- mean(abs(train$p[validation_i]-model1_pred[validation_i]),na.rm = T)
  
  # Model 2: Wind Speed and direction
  model2 <- qgam(p~s(ws100,bs="cr",k=20)+s(dir100,bs="cc",k = 36),
                 qu=0.5,
                 data=train[-validation_i,])
  model2_pred[validation_i] <- predict(model2, newdata = train[validation_i,])
  cv_results[fold,3] <- mean(abs(train$p[validation_i]-model2_pred[validation_i]),na.rm=T)
  
}
cv_results[7,1] <- "All"
cv_results[7,2] <- mean(abs(train$p-model1_pred),na.rm=T)
cv_results[7,3] <- mean(abs(train$p-model2_pred),na.rm=T)

# Bootstrap re-sampling MAE
MAE <- NULL
for(i in 1:250) {
  ind <- sample(train_length, replace = TRUE)
  testshort <- train[ind,]
  MAE <- rbind(MAE,  c(Benchmark = mean(abs(model1_pred[ind]-testshort$p),na.rm=T),
                       Proposed = mean(abs(model2_pred[ind]-testshort$p),na.rm=T)))
}
SS <- 1-MAE/matrix(rep(MAE[,1],2), nrow = nrow(MAE))

par(mfrow=c(1,2))
boxplot(MAE*100, ylab = "MAE")
boxplot(100*SS, ylab = "MAE Skill Score", ylim = c(0,5)); par(mfrow=c(1,1))

## Results for paper
require(xtable)
print(xtable(cv_results,digits = 3,display = rep("g",4)),
      include.rownames=FALSE)

## DM Test
e1 <- train$p-model1_pred
e2 <- train$p-model2_pred
notNA <- which(!is.na(e1) & !is.na(e2))
dm.test(e1[notNA],e2[notNA], h = 30,power=1)$p.value


## On-line: Rolling evaluation in test set ####

# Set window length
window_length <- 24*7

# Fit models
model1 <- with(train, lprq(ws100, p, h = 1, tau=0.5))
model2 <- qgam(p~s(ws100,k=20,bs="cr")+s(dir100,bs="cc",k = 36),
               qu=0.5,
               data=train)

# Initialise containers
roll_results <- data.frame(Window=1:floor(nrow(test)/window_length),
                           Model1_mae=NA,
                           Model2_mae=NA)
model1_pred_test <- model2_pred_test <- rep(NA,nrow(test))

# Run sliding window evaluation
for(window in 1:floor(nrow(test)/window_length)){
  
  test_i <- (window-1)*window_length + 1:window_length
  
  model1_pred_test[test_i] <- approx(model1$xx, model1$fv, test[test_i,]$ws100)$y
  roll_results[window,2] <- mean(abs(test$p[test_i]-model1_pred_test[test_i]),na.rm=T)
  
  model2_pred_test[test_i] <- predict(model2, newdata = test[test_i,])
  roll_results[window,3] <- mean(abs(test$p[test_i]-model2_pred_test[test_i]),na.rm=T)
  
}

plot(roll_results$Window,roll_results$Model1_mae,type="b",ylim = c(0,0.2),
     xlab="Window",ylab="Mean Absolute Error")
lines(roll_results$Window,roll_results$Model2_mae,type="b",col=2)

## Results for Paper
colnames(roll_results) <- c("Week","Benchmark","Proposed")
roll_results[,2:3] <- round(roll_results[,2:3], digits = 3)
# for(i in 1:nrow(roll_results)){
#   roll_results[i,which.min(roll_results[i,2:3])+1] <- 
#     paste0("\\textbf{", roll_results[i,which.min(roll_results[i,2:3])+1], "}")
# }
print(xtable(cbind(roll_results[1:10,],roll_results[11:20,],
                   roll_results[21:30,]),align = "lccc|ccc|ccc",digits = 3),
      include.rownames=FALSE,sanitize.text.function=identity)


# Bootstrap re-sampling MAE
MAE <- NULL
MAE_roll <- NULL
for(window in 1:floor(nrow(test)/window_length)){
  for(i in 1:250) {
    
    # Individual Windows
    test_i <- (window-1)*window_length + 1:window_length
    
    ind <- sample(test_i, replace = TRUE)
    testshort <- test[ind,]
    MAE <- rbind(MAE,  c(Window=window,
                         model1 = mean(abs(model1_pred_test[ind]-testshort$p),na.rm=T),
                         model2 = mean(abs(model2_pred_test[ind]-testshort$p),na.rm=T)))
    
    
    # Rolling Aggregate of Windows
    test_i <-  1:(window_length*window)
    ind <- sample(test_i, replace = TRUE)
    testshort <- test[ind,]
    MAE_roll <- rbind(MAE_roll,  c(Window=window,
                                   model1 = mean(abs(model1_pred_test[ind]-testshort$p),na.rm=T),
                                   model2 = mean(abs(model2_pred_test[ind]-testshort$p),na.rm=T)))
    
    
  }
}


## Plots

# Indiv Windows
SS <- cbind(Window=MAE[,1],1-MAE[,2:3]/matrix(rep(MAE[,2],2), nrow = nrow(MAE)))
# boxplot(model1~Window,MAE,ylab = "MAE")
boxplot(model2~Window,SS,
        main="Skill Score",xlab="Week",ylab = "MAE Skill Score")
abline(h=0,lty=2,lwd=2)

# Rolling Accumulation of Windows
SS_roll <- cbind(Window=MAE_roll[,1],1-MAE_roll[,2:3]/matrix(rep(MAE_roll[,2],2), nrow = nrow(MAE_roll)))
# boxplot(model1~Window,MAE_roll, ylab = "MAE")
boxplot(model2~Window,SS_roll,
        main="Evolution of Skill Score over Duration of Evaluation",xlab="Number of Weeks",ylab = "Cumulative MAE Skill Score",ylim=c(-.08,.08))
abline(h=0,lty=2,lwd=2)

## Plot for paper
par(mfrow=c(1,2))
boxplot(model2~Window,SS[which(SS[,1]<=30),],
        main="Individual Week Skill Score",xlab="Week",ylab = "MAE Skill Score",
        ylim=c(-0.15,0.15),axes=F); axis(1); axis(2); box()
abline(h=0,lty=2,lwd=2)
boxplot(model2~Window,SS_roll[which(SS[,1]<=30),],
        main="Aggregate Skill Score",xlab="Number of Weeks",ylab = "MAE Skill Score",
        ylim=c(-0.15,0.15),axes=F); axis(1); axis(2); box()
abline(h=0,lty=2,lwd=2); par(mfrow=c(1,1))

## DM Test
MAE_roll <- data.frame(Week=1:45,
                       Benchmark=NA,
                       Proposed=NA,
                      Difference=NA,
                      `DM p-value`=NA)
for(window in 1:nrow(MAE_roll)){
  
    # test_i <- (window-1)*window_length + 1:window_length
  test_i <- 1:(window_length*window)

  e1 <- model1_pred_test[test_i]-test$p[test_i]
  e2 <- model2_pred_test[test_i]-test$p[test_i]
  MAE_roll$Benchmark[window] <- mean(abs(e1),na.rm = T)
  MAE_roll$Proposed[window] <- mean(abs(e2),na.rm = T)
  MAE_roll$Difference[window] <- 100*(1-mean(abs(e2),na.rm = T)/mean(abs(e1),na.rm = T))
  MAE_roll$DM.p.value[window] <- dm.test(e1,e2, h = 30,power=1)$p.value  
}

MAE_roll[,2:5] <- round(MAE_roll[,2:5], digits = 3)
MAE_roll$Difference <- paste0(round(MAE_roll$Difference,digits = 1),"%")
MAE_roll$DM.p.value[MAE_roll$DM.p.value<0.001] <- "<0.001"
MAE_roll[,2] <- as.character(MAE_roll[,2])
MAE_roll[,3] <- as.character(MAE_roll[,3])
print(xtable(MAE_roll[c(1,5,10,15,30,45),],
             align = "lccccc"),
      include.rownames=FALSE,sanitize.text.function=identity)

# plot(MAE_roll$Week,MAE_roll$DM.p.value)
# abline(h=0.05,lty=2,col=2)
