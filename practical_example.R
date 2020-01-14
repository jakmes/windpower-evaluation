###############################
## Suplementary material for Section 4.4 of "Evaluation of Wind Power Forecasts - An up-to-date view" by Jakob W. Messner, 
## Pierre Pinson, Jethro Browell, MathiasB Bjerregard, and Irene Schicker.
## Author: Jethro Browell (jethro.browell@strath.ac.uk)
## Encoding: UTF-8
## License: GPL-3
###############################

## Load required packages
library("quantreg")
library("qgam")
library("forecast")
library("rstudioapi") # only used for setting working directory

# Set working directory
setwd(dirname(getSourceEditorContext()$path))

# Set seed to reproduce results as they appear in the paper.
set.seed(1)

## Load and prepare data ####
# 
# A subset of data from GEFcom2014 is used. Originally published here: Tao Hong, Pierre Pinson, Shu Fan, 
# Hamidreza Zareipour, Alberto Troccoli, Rob J. Hyndman, "Probabilistic energy forecasting: Global Energy Forecasting 
# Competition 2014 and beyond", International Journal of Forecasting,Volume 32, Issue 3,2016,Pages 896-913,ISSN 0169-2070,
# https://doi.org/10.1016/j.ijforecast.2016.02.001

data <- read.csv("Task15_W_Zone1_10/Task15_W_Zone1.csv")
names(data) <- c("zone_id", "timestamp", "p", "u10", "v10", "u100", "v100")

# Add wind speed and direction features
data$ws10 <- with(data, sqrt(u10^2 + v10^2))
data$ws100 <- with(data, sqrt(u100^2 + v100^2))
data$dir100 <- with(data, atan(u100/v100))
data <- na.omit(data)

# Partition into training and test data
train_length <- 8760L
train <- data[1:train_length,]
test <- tail(data, -train_length)


## Off-line: k-fold Cross-validation ####

# Number of CV folds and stratification
n_fold <- 6L
stratification <-2L
validation_is <- rep(rep(1:n_fold,each=floor(train_length/(stratification*n_fold))),stratification)

# Initialise containers for results and models
cv_results <- data.frame(Fold=1:n_fold,
                         Benchmark=NA,
                         Proposed=NA)
Benchmark_pred <- Proposed_pred <- rep(NA,train_length)

# Perform k-fold cross-validation
for(fold in 1:n_fold){
  
  validation_i <- which(validation_is==fold)
  
  # Benchmark: LOESS with Wind Speed at 100m only
  Benchmark <- with(train[-validation_i,], lprq(ws100, p, h = 1, tau=0.5))
  Benchmark_pred[validation_i] <- approx(Benchmark$xx, Benchmark$fv, train[validation_i,]$ws100)$y
  cv_results[fold,2] <- mean(abs(train$p[validation_i]-Benchmark_pred[validation_i]),na.rm = T)
  
  # Proposed: Spline regression with Wind Speed and direction
  Proposed <- qgam(p~s(ws100,bs="cr",k=20)+s(dir100,bs="cc",k = 36),
                 qu=0.5,
                 data=train[-validation_i,])
  Proposed_pred[validation_i] <- predict(Proposed, newdata = train[validation_i,])
  cv_results[fold,3] <- mean(abs(train$p[validation_i]-Proposed_pred[validation_i]),na.rm=T)
  
}
cv_results[7,1] <- "All"
cv_results[7,2] <- mean(abs(train$p-Benchmark_pred),na.rm=T)
cv_results[7,3] <- mean(abs(train$p-Proposed_pred),na.rm=T)

# Bootstrap re-sampling of MAE
MAE <- NULL
for(i in 1:250) {
  ind <- sample(train_length, replace = TRUE)
  testshort <- train[ind,]
  MAE <- rbind(MAE,  c(Benchmark = mean(abs(Benchmark_pred[ind]-testshort$p),na.rm=T),
                       Proposed = mean(abs(Proposed_pred[ind]-testshort$p),na.rm=T)))
}
# Skill score
SS <- 1-MAE/matrix(rep(MAE[,1],2), nrow = nrow(MAE))

# Plot for paper
par(mfrow=c(1,2))
boxplot(MAE*100, ylab = "MAE")
boxplot(100*SS, ylab = "MAE Skill Score", ylim = c(0,5)); par(mfrow=c(1,1))

## Results for paper
require(xtable)
print(xtable(cv_results,digits = 3,display = rep("g",4)),
      include.rownames=FALSE)

## DM Test
e1 <- train$p-Benchmark_pred
e2 <- train$p-Proposed_pred
notNA <- which(!is.na(e1) & !is.na(e2))
dm.test(e1[notNA],e2[notNA], h = 30,power=1)$p.value



## On-line: Rolling evaluation in test set ####

# Set window length
window_length <- 24*7

# Fit models
Benchmark <- with(train, lprq(ws100, p, h = 1, tau=0.5))
Proposed <- qgam(p~s(ws100,k=20,bs="cr")+s(dir100,bs="cc",k = 36),
               qu=0.5,
               data=train)

# Initialise containers for results
roll_results <- data.frame(Window=1:floor(nrow(test)/window_length),
                           Benchmark_mae=NA,
                           Proposed_mae=NA)
Benchmark_pred_test <- Proposed_pred_test <- rep(NA,nrow(test))

# Run sliding window evaluation
for(window in 1:floor(nrow(test)/window_length)){
  
  test_i <- (window-1)*window_length + 1:window_length
  
  Benchmark_pred_test[test_i] <- approx(Benchmark$xx, Benchmark$fv, test[test_i,]$ws100)$y
  roll_results[window,2] <- mean(abs(test$p[test_i]-Benchmark_pred_test[test_i]),na.rm=T)
  
  Proposed_pred_test[test_i] <- predict(Proposed, newdata = test[test_i,])
  roll_results[window,3] <- mean(abs(test$p[test_i]-Proposed_pred_test[test_i]),na.rm=T)
  
}

## Results for Paper
colnames(roll_results) <- c("Week","Benchmark","Proposed")
roll_results[,2:3] <- round(roll_results[,2:3], digits = 3)
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
                         Benchmark = mean(abs(Benchmark_pred_test[ind]-testshort$p),na.rm=T),
                         Proposed = mean(abs(Proposed_pred_test[ind]-testshort$p),na.rm=T)))
    
    
    # Rolling Aggregate of Windows
    test_i <-  1:(window_length*window)
    ind <- sample(test_i, replace = TRUE)
    testshort <- test[ind,]
    MAE_roll <- rbind(MAE_roll,  c(Window=window,
                                   Benchmark = mean(abs(Benchmark_pred_test[ind]-testshort$p),na.rm=T),
                                   Proposed = mean(abs(Proposed_pred_test[ind]-testshort$p),na.rm=T)))
    
    
  }
}


## Skill scores

# Individual windows
SS <- cbind(Window=MAE[,1],1-MAE[,2:3]/matrix(rep(MAE[,2],2), nrow = nrow(MAE)))

# Rolling accumulation of results
SS_roll <- cbind(Window=MAE_roll[,1],1-MAE_roll[,2:3]/matrix(rep(MAE_roll[,2],2), nrow = nrow(MAE_roll)))

## Plot for paper
par(mfrow=c(1,2))
boxplot(Proposed~Window,SS[which(SS[,1]<=30),],
        main="Individual Week Skill Score",xlab="Week",ylab = "MAE Skill Score",
        ylim=c(-0.15,0.15),axes=F); axis(1); axis(2); box()
abline(h=0,lty=2,lwd=2)
boxplot(Proposed~Window,SS_roll[which(SS[,1]<=30),],
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
  
  test_i <- 1:(window_length*window)

  e1 <- Benchmark_pred_test[test_i]-test$p[test_i]
  e2 <- Proposed_pred_test[test_i]-test$p[test_i]
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
