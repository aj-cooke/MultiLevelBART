library(dbarts)
library(dplyr)
library(MASS)
library(tidyr)
library(ggplot2)
source('multilevel_sim.R')

r.squared <- function(ytrue, ypred){
  return(1-(sum((ytrue-ypred)^2)/sum((ytrue-mean(ytrue))^2)))
}

mae <- function(ytrue, ypred){
  return(mean(abs(ytrue-ypred)))
}

minmax <- function(x, shrink) {
  mm <- (x- min(x)) /(max(x)-min(x))
  mm <- mm * (1-shrink*2) + shrink
  return(mm)
}

N <- 100
d <- 10

dat <- generate_mvn(N, d)
betas <- rnorm(d, 0, 1.2)
p_raw <- as.matrix(dat)%*%betas
dat$p.score <- minmax(p_raw, 0.1)
hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
dat$z <- rbinom(nrow(dat),1,dat$p.score)

bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10)
bart_pred <- fitted(bart_fit)
hist(bart_pred)

logistic <- glm(z ~ . -p.score, data = dat, family = binomial())
pred_log <- fitted(logistic)
hist(pred_log)

r.squared(dat$p.score, bart_pred)
r.squared(dat$p.score, pred_log)
mae(rank(dat$p.score), rank(bart_pred)) / max(rank(dat$p.score))
mae(rank(dat$p.score), rank(pred_log)) / max(rank(dat$p.score))

### Parameters that optimize for rank vs r^2, accuracy on z





