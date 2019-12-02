### R code from vignette source 'evaluation.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: evaluation.Rnw:136-141
###################################################
par(mfrow = c(1,3))
a <- seq(-4,4,0.1)
plot(a, a^2, type = "l", xlab = expression(epsilon), ylab = expression(rho(epsilon)))
plot(a, abs(a), type = "l", xlab = expression(epsilon), ylab = expression(rho(epsilon)))
plot(a, a*(0.3-(a<0)), type = "l", xlab = expression(epsilon), ylab = expression(rho(epsilon)), ylim = c(0,4))


###################################################
### code chunk number 2: evaluation.Rnw:147-187
###################################################
library("xtable")
library("quantreg")
data <- read.csv("Task15_W_Zone1_10/Task15_W_Zone1.csv")
names(data) <- c("zone_id", "timestamp", "p", "u10", "v10", "u100", "v100")
data$ws10 <- with(data, sqrt(u10^2 + v10^2))
data$ws100 <- with(data, sqrt(u100^2 + v100^2))
data <- na.omit(data)


train <- data[1:10000,]
test <- tail(data, -10000)

if(file.exists("motivationmodels.rda")) {
  load("motivationmodels.rda")
} else {
  ols <- loess(p ~ ws100, train, span = 0.5)
  median <- with(train, lprq(ws100, p, h = 1, tau=0.5))
  quantile <- with(train, lprq(ws100, p, h = 1, tau=0.3))
  
  ols_pred <- predict(ols, newdata = test)
  median_pred <- approx(median$xx, median$fv, test$ws100)$y
  quantile_pred <- approx(quantile$xx, quantile$fv, test$ws100)$y
  
    
  save(ols, median, quantile, ols_pred, median_pred, quantile_pred, file = "motivationmodels.rda")
}



#plot(p~ws100, data, col = gray(0.2, alpha = 0.3), xlab = "wind speed forecast", ylab = "power generation")
#lines(predict(ols, newdata = data.frame(ws100 = quantile$xx))~quantile$xx, col = 2, lwd = 2)
#lines(median$fv~median$xx, col = 3, lwd = 2)
#lines(quantile$fv~quantile$xx, col = 4, lwd = 2)
#legend("bottomright", lwd = 2, col = 2:4, c("quadratic", "absolute", "quantile"), bty = "n")
plot(tail(test$p, 24), type = "l", xlab = "lead time [h]", ylab = "normalized power generation", lwd = 2, ylim = c(0,1))
lines(tail(ols_pred, 24), col = 2, lwd = 2)
lines(tail(median_pred, 24), col = 3, lwd = 2)
lines(tail(quantile_pred, 24), col = 4, lwd = 2)

legend("topleft", col = 1:4, lwd = 2, c("measurements", "quadratic", "absolute", "quantile"), bty = "n")


###################################################
### code chunk number 3: evaluation.Rnw:195-219
###################################################


qs <- function(pred, obs = test$p, tau = 0.4) {
    err <- pred-obs
    mean(err*((obs<=pred)-tau))
  }

errors <- data.frame(MSE = c(
quadratic = mean((ols_pred-test$p)^2),
absolute  = mean((median_pred-test$p)^2),
quantile  = mean((quantile_pred-test$p)^2)),
MAE = c(
quadratic = mean(abs(ols_pred-test$p)),
absolute  = mean(abs(median_pred-test$p)),
quantile  = mean(abs(quantile_pred-test$p))),
QS = c(
quadratic = qs(ols_pred),
absolute  = qs(median_pred),
quantile  = qs(quantile_pred))
) 
errors <- round(errors, digits = 4)
for(i in 1:3) errors[which.min(errors[,i]),i] <- paste0("\\textbf{", errors[which.min(errors[,i]),i], "}")
latextab <- xtable(errors, digits = 4, label = "tab:motivationscores", caption = "Different evaluation measures for the three local linear models with quadratic, absolute and quantile loss function. The best model for each score is highlighted in bold")
print(latextab, , sanitize.text.function=identity)


###################################################
### code chunk number 4: evaluation.Rnw:233-253
###################################################
ind <- 1:200
testshort <- test[ind,]

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
for(i in 1:3) errors[which.min(errors[,i]),i] <- paste0("\\textbf{", errors[which.min(errors[,i]),i], "}")
latextab <- xtable(errors, digits = 6, label = "tab:motivationscoresshort", caption = paste0("Same as Table~\\ref{tab:motivationscores} but only derived from the first ", nrow(testshort), " time steps of the test data set"))
print(latextab, , sanitize.text.function=identity)


###################################################
### code chunk number 5: evaluation.Rnw:263-268
###################################################
test$fc <- test$ws100<3
test$p0 <- test$p<=0
conttab1 <- table(test$p0, test$fc)
conttab2 <- table(test$p0, median_pred<=0)



###################################################
### code chunk number 6: fig_bias
###################################################
set.seed(1)
y <- c(0.2,0.3,0.5,0.9,1,1,1,0.4,0.5)
plot(y, ylab="wind power", xlab = "", xaxt="n", lwd = 2, type = "l", ylim = c(0,1))
title(xlab = "Time", line = 1)
abline(h=mean(y), col = 1, lwd = 2, lty = 2)

fc1 <- sample(y)+0.01
lines(fc1, col = 2, lwd = 2)
abline(h=mean(fc1), col = 2, lwd = 2, lty = 2)

fc2 <- y-0.2
lines(fc2, col = 4, lwd = 2)
abline(h=mean(fc2), col = 4, lwd = 2, lty = 2)

legend("bottom", lwd=2, col = c(1:2,4,1:2,4), lty = c(1,1,1,2,2,2), c("observations", "forecast 1", "forecast 2", "average observations", "average forecast 1", "average forecast 2"), bty = "n", ncol = 2)





###################################################
### code chunk number 7: fig_roc
###################################################
mod <- glm(p>0.4~ws100, train, family=binomial)
pred <- predict(mod, type = "response", newdata = test)
th <- seq(0,1, 0.01)
a <- NULL
for(i in 1:length(th)) {
  hr <- sum(pred>th[i] & test$p>0.4)/sum(test$p>0.4)
  far <- sum(pred>th[i] & test$p<=0.4)/sum(test$p<0.4)
  a <- rbind(a, c(hr,far))
}
plot(a[,1]~a[,2], type = "l", xlab = "FAR", ylab = "HR", lwd = 2, xlim = c(0.0336,1-0.036), ylim = c(0.036,1-0.036))
polygon(c(a[,2], 0,1), c(a[,1], 0,0), col = rgb(1,0,0, alpha = 0.1))
abline(0,1, lwd = 0.6)




###################################################
### code chunk number 8: fig_reldiag
###################################################
library("SpecsVerification")
par(mar=c(5,4,2,4))
ReliabilityDiagram(pred, test$p>0.4, plot=TRUE)


###################################################
### code chunk number 9: fig_rankhist
###################################################
set.seed(1)
hist(runif(500,1,11), freq = FALSE, xlab = "rank", ylab = "relative frequency", breaks = seq(1, 11, 1), main = "")
abline(h=1/10, col = 2)


###################################################
### code chunk number 10: fig_crps
###################################################
a <- seq(0,1,0.05)
plot(a, pnorm(a, 0.3, 0.1), type = "l", ylab = "Cumulative distribution function", xlab = "normalized wind power")
polygon(c(a,1,0.35,0.35,0), c(pnorm(a, 0.3, 0.1),1,1,0,0), col = 2, density = 20)
lines(a, pnorm(a, 0.3, 0.1), lwd = 2)
lines(c(0,0.35,0.35,1), c(0,0,1,1), col = 2, lwd = 2)


###################################################
### code chunk number 11: fig_probloss
###################################################
library("scoringRules")
a <- seq(0,1,0.05)
plot(a, crps(a, family="norm", mean = 0.3, sd = 0.1), type = "l", xlab = "normalized wind power", ylab = "score contribution", ylim = c(-1.6,1.5), lwd = 2)
lines(a, logs(a, family="norm", mean = 0.3, sd = 0.1), col = 2, lwd = 2)
legend("bottomright", col = 1:2, lwd = 2, c("CRPS", "IS"), bty = "n")


###################################################
### code chunk number 12: evaluation.Rnw:856-877
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
### code chunk number 13: evaluation.Rnw:907-922
###################################################
par(mfrow=c(1,2))
SS <- 1-MSE/matrix(rep(MSE[,1],3), nrow = nrow(MSE))
boxplot(SS, ylab = "mean squared error skill score", ylim = c(-0.7, 0.1), las = 2)
abline(h=0, lwd=0.7)

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
### code chunk number 14: evaluation.Rnw:971-983
###################################################
library("forecast")
E <- data.frame(
quadratic = (ols_pred-test$p),
absolute  = (median_pred-test$p),
quantile  = (quantile_pred-test$p))

tab <- data.frame(absolute = dm.test(E$quadratic, E$absolute, h = 30)$p.value, 
  quantile = dm.test(E$quadratic, E$quantile, h = 30)$p.value)
rownames(tab) <- "quadratic"
 
#tab <- round(tab, digits = 4)
xtable(tab, display = rep("g",3), label = "tab:dmtest", caption = "p-values from two sided Diebold-Mariano tests for equality of the squared error for the absolute and quantile loss models compared to the quadratic loss model. See Section~\ref{sec:context} for details on the test setup. Lower values signify more significant differences.")


