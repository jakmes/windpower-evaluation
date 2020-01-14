###################################################
## R code to reproduce the results in "Evaluation of Wind Power Forecasts - An up-to-date view" by Jakob W. Messner, 
## Pierre Pinson, Jethro Browell, MathiasB Bjerregard, and Irene Schicker.
## Author: Jakob W. Messner (jakob.messner@posteo.net)
## Encoding: UTF-8
## License: GPL-3
###################################################



###################################################
### code chunk number 1: Figure 1
###################################################
par(mfrow = c(1,3))
a <- seq(-4,4,0.1)
plot(a, a^2, type = "l", xlab = expression(epsilon), ylab = expression(rho(epsilon)))
plot(a, abs(a), type = "l", xlab = expression(epsilon), ylab = expression(rho(epsilon)))
plot(a, a*(0.3-(a<0)), type = "l", xlab = expression(epsilon), ylab = expression(rho(epsilon)), ylim = c(0,4))


###################################################
### code chunk number 2: Load data, fit models and plot Figure 2
###################################################
library("xtable")
library("quantreg")
## The full GEFCom2014 data set is available as appendix to 
# Hong T, Pinson P, Fan S, Zareipour H, Troccoli A, Hyndman RJ. "Probabilistic Energy Forecasting: 
# Global Energy Forecasting Competition 2014 and Beyond". International Journal of Forecasting 2016; 32(3): 896-913. 
# doi:10.1016/j.ijforecast.2016.02.001 
# or at https://www.dropbox.com/s/pqenrr2mcvl0hk9/GEFCom2014.zip?dl=0
## The subset that we use in the following is also contained in this git repository
data <- read.csv("Task15_W_Zone1_10/Task15_W_Zone1.csv")
names(data) <- c("zone_id", "timestamp", "p", "u10", "v10", "u100", "v100")
data$ws10 <- with(data, sqrt(u10^2 + v10^2))
data$ws100 <- with(data, sqrt(u100^2 + v100^2))
data <- na.omit(data)

## split into training and test data
train <- data[1:10000,]
test <- tail(data, -10000)

## fit models
ols <- loess(p ~ ws100, train, span = 0.5, degree = 1)
median <- with(train, lprq(ws100, p, h = 1, tau=0.5))
quantile <- with(train, lprq(ws100, p, h = 1, tau=0.3))
  
## generate predictions
ols_pred <- predict(ols, newdata = test)
median_pred <- approx(median$xx, median$fv, test$ws100)$y
quantile_pred <- approx(quantile$xx, quantile$fv, test$ws100)$y
  
    

## Figure 2: plot observed and predicted power
plot(tail(test$p, 24), type = "l", xlab = "lead time [h]", ylab = "normalized power generation", lwd = 2, ylim = c(0,1))
lines(tail(ols_pred, 24), col = 2, lwd = 2)
lines(tail(median_pred, 24), col = 3, lwd = 2)
lines(tail(quantile_pred, 24), col = 4, lwd = 2)

legend("topleft", col = 1:4, lwd = 2, c("measurements", "quadratic", "absolute", "quantile"), bty = "n")


###################################################
### code chunk number 3: Table 1
###################################################

## define quantile score function
qs <- function(pred, obs = test$p, tau) {
    err <- pred-obs
    mean(err*((obs<=pred)-tau))
  }

## build data frame of error measures for different models
errors <- data.frame(MSE = c(
quadratic = mean((ols_pred-test$p)^2),
absolute  = mean((median_pred-test$p)^2),
quantile  = mean((quantile_pred-test$p)^2)),
MAE = c(
quadratic = mean(abs(ols_pred-test$p)),
absolute  = mean(abs(median_pred-test$p)),
quantile  = mean(abs(quantile_pred-test$p))),
QS = c(
quadratic = qs(ols_pred,tau = 0.3),
absolute  = qs(median_pred,tau = 0.3),
quantile  = qs(quantile_pred,tau = 0.3))
) 
errors <- round(errors, digits = 4)

## highlight best model
for(i in 1:3) errors[which.min(errors[,i]),i] <- paste0("\\textbf{", errors[which.min(errors[,i]),i], "}")

## output latex table
latextab <- xtable(errors, digits = 4, label = "tab:motivationscores", caption = "Different evaluation measures for the three local linear models with quadratic, absolute and quantile loss function. The best model for each score is highlighted in bold")
print(latextab, , sanitize.text.function=identity)


###################################################
### code chunk number 4: Table 2
###################################################
## reduced test data set
ind <- 1:200
testshort <- test[ind,]

## build data frame with error measures
errors <- data.frame(MSE = c(
quadratic = mean((ols_pred[ind]-testshort$p)^2),
absolute  = mean((median_pred[ind]-testshort$p)^2),
quantile  = mean((quantile_pred[ind]-testshort$p)^2)),
MAE = c(
quadratic = mean(abs(ols_pred[ind]-testshort$p)),
absolute  = mean(abs(median_pred[ind]-testshort$p)),
quantile  = mean(abs(quantile_pred[ind]-testshort$p))),
QS = c(
quadratic = qs(ols_pred[ind], testshort$p),
absolute  = qs(median_pred[ind], testshort$p),
quantile  = qs(quantile_pred[ind], testshort$p))
) 
errors <- round(errors, digits = 4)

## highlight best model
for(i in 1:3) errors[which.min(errors[,i]),i] <- paste0("\\textbf{", errors[which.min(errors[,i]),i], "}")

## output latex table
latextab <- xtable(errors, digits = 6, label = "tab:motivationscoresshort", caption = paste0("Same as Table~\\ref{tab:motivationscores} but only derived from the first ", nrow(testshort), " time steps of the test data set"))
print(latextab, , sanitize.text.function=identity)


###################################################
### code chunk number 5: Table 3
###################################################
test$fc <- test$ws100<3
test$p0 <- test$p<=0
conttab1 <- table(test$p0, test$fc)
conttab2 <- table(test$p0, median_pred<=0)



###################################################
### code chunk number 6: Figure 3
###################################################
set.seed(1)
## define  and plot time series
y <- c(0.2,0.3,0.5,0.9,1,1,1,0.4,0.5)
plot(y, ylab="wind power", xlab = "", xaxt="n", lwd = 2, type = "l", ylim = c(0,1))
title(xlab = "Time", line = 1)
abline(h=mean(y), col = 1, lwd = 2, lty = 2)

## bias free forecast
fc1 <- sample(y)+0.01
lines(fc1, col = 2, lwd = 2)
abline(h=mean(fc1), col = 2, lwd = 2, lty = 2)

## shifted forecast
fc2 <- y-0.2
lines(fc2, col = 4, lwd = 2)
abline(h=mean(fc2), col = 4, lwd = 2, lty = 2)

legend("bottom", lwd=2, col = c(1:2,4,1:2,4), lty = c(1,1,1,2,2,2), 
       c("observations", "forecast 1", "forecast 2", "average observations", "average forecast 1", "average forecast 2"), 
       bty = "n", ncol = 2)





###################################################
### code chunk number 7: Figure 4
###################################################
## fit and predict model
mod <- glm(p>0.4~ws100, train, family=binomial)
pred <- predict(mod, type = "response", newdata = test)

## define thresholds to derive hit and false alarm rate
th <- seq(0,1, 0.01)

## derive hit rate and false alarm rate for all th
a <- NULL
for(i in 1:length(th)) {
  hr <- sum(pred>th[i] & test$p>0.4)/sum(test$p>0.4)
  far <- sum(pred>th[i] & test$p<=0.4)/sum(test$p<0.4)
  a <- rbind(a, c(hr,far))
}

## plot
plot(a[,1]~a[,2], type = "l", xlab = "FAR", ylab = "HR", lwd = 2, xlim = c(0.0336,1-0.036), ylim = c(0.036,1-0.036))
polygon(c(a[,2], 0,1), c(a[,1], 0,0), col = rgb(1,0,0, alpha = 0.1))
abline(0,1, lwd = 0.6)




###################################################
### code chunk number 8: Figure 5
###################################################
library("SpecsVerification")
par(mar=c(5,4,2,4))
ReliabilityDiagram(pred, test$p>0.4, plot=TRUE)


###################################################
### code chunk number 9: Figure 6
###################################################
set.seed(1)
hist(runif(500,1,11), freq = FALSE, xlab = "rank", ylab = "relative frequency", breaks = seq(1, 11, 1), main = "")
abline(h=1/10, col = 2)


###################################################
### code chunk number 10: Figure 7
###################################################
a <- seq(0,1,0.05)
plot(a, pnorm(a, 0.3, 0.1), type = "l", ylab = "Cumulative distribution function", xlab = "normalized wind power")
polygon(c(a,1,0.35,0.35,0), c(pnorm(a, 0.3, 0.1),1,1,0,0), col = 2, density = 20)
lines(a, pnorm(a, 0.3, 0.1), lwd = 2)
lines(c(0,0.35,0.35,1), c(0,0,1,1), col = 2, lwd = 2)


###################################################
### code chunk number 11: Figure 8
###################################################
library("scoringRules")
a <- seq(0,1,0.05)
plot(a, crps(a, family="norm", mean = 0.3, sd = 0.1), type = "l", xlab = "normalized wind power", 
     ylab = "score contribution", ylim = c(-1.6,1.5), lwd = 2)
lines(a, logs(a, family="norm", mean = 0.3, sd = 0.1), col = 2, lwd = 2)
legend("bottomright", col = 1:2, lwd = 2, c("CRPS", "IS"), bty = "n")


###################################################
### code chunk number 12: Figure 10
###################################################
par(mfrow = c(1,2))
MSE <- NULL
for(i in 1:20) {
  ind <- ((i-1)*200+1):(i*200)
  testshort <- test[ind,]
  MSE <- rbind(MSE,  c(quadratic = mean((ols_pred[ind]-testshort$p)^2), absolute = mean((median_pred[ind]-testshort$p)^2), quantile = mean((quantile_pred[ind]-testshort$p)^2)))
}


boxplot(MSE, ylab = "mean squared error", ylim = c(0.01, 0.1))

MSE <- NULL
for(i in 1:10) {
  ind <- ((i-1)*400+1):(i*400)
  testshort <- test[ind,]
  MSE <- rbind(MSE,  c(quadratic = mean((ols_pred[ind]-testshort$p)^2), absolute = mean((median_pred[ind]-testshort$p)^2), quantile = mean((quantile_pred[ind]-testshort$p)^2)))
}


boxplot(MSE, ylab = "mean squared error", ylim = c(0.01, 0.1))



###################################################
### code chunk number 13:Figure 11
###################################################
## derive skill scores
SS <- 1-MSE/matrix(rep(MSE[,1],3), nrow = nrow(MSE))

## plot
par(mfrow=c(1,2))
boxplot(SS, ylab = "mean squared error skill score", ylim = c(-0.7, 0.1), las = 2)
abline(h=0, lwd=0.7)

## bootstrap averages
MSE <- NULL
for(i in 1:250) {
  ind <- sample(nrow(test), replace = TRUE)
  testshort <- test[ind,]
  MSE <- rbind(MSE,  c(quadratic = mean((ols_pred[ind]-testshort$p)^2), absolute = mean((median_pred[ind]-testshort$p)^2), quantile = mean((quantile_pred[ind]-testshort$p)^2)))
}
SS <- 1-MSE/matrix(rep(MSE[,1],3), nrow = nrow(MSE))
boxplot(SS, ylab = "mean squared error skill score", ylim = c(-0.7, 0.1))
abline(h=0, lwd=0.7)



###################################################
### code chunk number 14: Table 5
###################################################
library("forecast")
E <- data.frame(
quadratic = (ols_pred-test$p),
absolute  = (median_pred-test$p),
quantile  = (quantile_pred-test$p))

## Diebold Mariano test
tab <- data.frame(absolute = dm.test(E$quadratic, E$absolute, h = 30)$p.value, 
  quantile = dm.test(E$quadratic, E$quantile, h = 30)$p.value)
rownames(tab) <- "quadratic"
 
#tab <- round(tab, digits = 4)
xtable(tab, display = rep("g",3), label = "tab:dmtest", caption = "p-values from two sided Diebold-Mariano tests for equality of the squared error for the absolute and quantile loss models compared to the quadratic loss model. See Section~\ref{sec:context} for details on the test setup. Lower values signify more significant differences.")


