
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
set.seed(1)
library("quantreg")
library("splines")
data <- read.csv("Task15_W_Zone1_10/Task15_W_Zone1.csv")
names(data) <- c("zone_id", "timestamp", "p", "u10", "v10", "u100", "v100")
data$ws10 <- with(data, sqrt(u10^2 + v10^2))
data$ws100 <- with(data, sqrt(u100^2 + v100^2))
data$dir100 <- with(data, atan(u100/v100))
data <- na.omit(data)

train_length <- 10000L
train <- data[1:train_length,]
test <- tail(data, -train_length)


## Off-line vs on-line testing ####


## Off-line: k-fold Cross-validation

n_fold <- 5
cv_results <- data.frame(Fold=1:n_fold,
                         Model1_mae=NA,
                         Model2_mae=NA)
model1_pred <- model2_pred <- rep(NA,train_length)
for(fold in 1:n_fold){
  
  validation_i <- (fold-1)*floor(train_length/n_fold) + 1:floor(train_length/n_fold)
  
  # Model 1: Wind Speed at 100m only
  model1 <- with(train[-validation_i,], lprq(ws100, p, h = 1, tau=0.5))
  model1_pred[validation_i] <- approx(model1$xx, model1$fv, train[validation_i,]$ws100)$y
  cv_results[fold,2] <- mean(abs(train$p[validation_i]-model1_pred[validation_i]),na.rm = T)
  
  # Model 2: Wind Speed and direction  
  model2 <- with(train[-validation_i,], rq(p~bs(ws100,Boundary.knots = c(0,20),df=10)+
                                             bs(dir100,Boundary.knots = pi/c(-2,2),df=18),
                                           tau=0.5))
  model2_pred[validation_i] <- predict(model2, newdata = train[validation_i,])
  cv_results[fold,3] <- mean(abs(train$p[validation_i]-model2_pred[validation_i]),na.rm=T)
  
}

# cv_results

MAE <- NULL
for(i in 1:250) {
  ind <- sample(train_length, replace = TRUE)
  testshort <- train[ind,]
  MAE <- rbind(MAE,  c(model1 = mean(abs(model1_pred[ind]-testshort$p),na.rm=T),
                       model2 = mean(abs(model2_pred[ind]-testshort$p),na.rm=T)))
}
SS <- 1-MAE/matrix(rep(MAE[,1],2), nrow = nrow(MAE))
boxplot(MAE, ylab = "MAE")
boxplot(SS, ylab = "MAE Skill Score", ylim = c(0, 0.05))




## On-line: Rolling evaluation in test set

window_length <- 24*7

model1 <- with(train, lprq(ws100, p, h = 1, tau=0.5))
model2 <- with(train, rq(p~bs(ws100,Boundary.knots = c(0,20),df=10)+
                           bs(dir100,Boundary.knots = pi/c(-2,2),df=18),
                         tau=0.5))

roll_results <- data.frame(Window=1:floor(nrow(test)/window_length),
                           Model1_mae=NA,
                           Model2_mae=NA)
model1_pred_test <- model2_pred_test <- rep(NA,nrow(test))

for(window in 1:floor(nrow(test)/window_length)){
  
  test_i <- (window-1)*floor(nrow(test)/window_length) + 1:window_length
  
  model1_pred_test[test_i] <- approx(model1$xx, model1$fv, test[test_i,]$ws100)$y
  roll_results[window,2] <- mean(abs(test$p[test_i]-model1_pred_test[test_i]),na.rm=T)
  
  model2_pred_test[test_i] <- predict(model2, newdata = test[test_i,])
  roll_results[window,3] <- mean(abs(test$p[test_i]-model2_pred_test[test_i]),na.rm=T)
  
}

plot(roll_results$Window,roll_results$Model1_mae,type="b",ylim = c(0,0.2),
     xlab="Window",ylab="Mean Absolute Error")
lines(roll_results$Window,roll_results$Model2_mae,type="b",col=2)




MAE <- NULL
MAE_roll <- NULL
for(window in 1:floor(nrow(test)/window_length)){
  for(i in 1:250) {
    
    # Individual Windows
    test_i <- (window-1)*floor(nrow(test)/window_length) + 1:window_length
    
    ind <- sample(test_i, replace = TRUE)
    testshort <- test[ind,]
    MAE <- rbind(MAE,  c(Window=window,
                         model1 = mean(abs(model1_pred_test[ind]-testshort$p),na.rm=T),
                         model2 = mean(abs(model2_pred_test[ind]-testshort$p),na.rm=T)))
    

    # Rolling Accumulation of Windows
    test_i <-  1:(window_length*window)
    ind <- sample(test_i, replace = TRUE)
    testshort <- test[ind,]
    MAE_roll <- rbind(MAE_roll,  c(Window=window,
                         model1 = mean(abs(model1_pred_test[ind]-testshort$p),na.rm=T),
                         model2 = mean(abs(model2_pred_test[ind]-testshort$p),na.rm=T)))
    
        
  }
}

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


