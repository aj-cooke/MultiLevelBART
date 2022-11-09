library(dbarts)
library(dplyr)
library(MASS)
library(tidyr)
library(ggplot2)
source('multilevel_sim.R')

r.squared <- function(ytrue, ypred){
  return(1-(sum((ytrue-ypred)^2)/sum((ytrue-mean(ytrue))^2)))
}

minmax <- function(x, shrink) {
  mm <- (x- min(x)) /(max(x)-min(x))
  mm <- mm * (1-shrink*2) + shrink
  return(mm)
}

cv_bart <- function(dat, a, b, k, nt){
  params <- crossing(a = a, b = b, k = k, nt = nt)
  results <- data.frame(a= c(), b=c(), k=c(), nt = c(), calibration_r2 = c(),
                        minimum = c(), maximum = c())
  for(i in 1:nrow(params)){  
    bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10,
                              n.trees = params$nt[i],
                              k = params$k[i],
                              power = params$b[i],
                              base = params$a[i])
    pred <- fitted(bart_fit)
    calibration_r2 <- r.squared(dat$p.score, pred)
    minimum <- min(pred)
    maximum <- max(pred)
    results <- rbind(results, cbind(params[i,], calibration_r2, minimum, maximum))
    print(params[i,], calibration_r2, minimum, maximum)
  }
  return(results)
}

### CV all params, varying N
results <- data.frame()

N <- c(100, 500, 1000, 5000, 10000)
vars <- 10

a <- c(0.95, 0.75, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
b <- c(2, 4, 6, 8, 10) # crashes if b too low
k <- c(1, 2, 4, 5, 6)
nt <- c(200)

for(i in 1:N){
  dat <- generate_mvn(i, vars)
  betas <- rnorm(vars, .25, .6)
  p_raw <- as.matrix(dat)%*%betas 
  dat$p.score <- minmax(p_raw, 0.1)
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  r <- cv_bart(dat, a, b, k, nt)
  r$N <- as.integer(i)
  r$vars <- as.integer(vars)
  results <- rbind(results, r)
}


######## Searching K
cv_bart_k <- function(dat, k){
  results <- data.frame(k=c(), calibration_r2 = c(),
                        minimum = c(), maximum = c())
  for(i in 1:length(k)){  
    bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10,
                              k = k[i])
    pred <- fitted(bart_fit)
    calibration_r2 <- r.squared(dat$p.score, pred)
    minimum <- min(pred)
    maximum <- max(pred)
    k_val = k[i]
    results <- rbind(results, cbind(k_val, calibration_r2, minimum, maximum))
  }
  return(results)
}

###### K Fixed

results2 <- data.frame()
ITERS <- 1000 ### Change this to run for X amount of times
N <- c(100, 200, 300, 500, 1000)
rep_n <- 1:ITERS
nv <- crossing(N=N, rep_n = rep_n)
k <- seq(0.1,10,0.1) 

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.integer(nv[i,'N']), 10)
  betas <- rnorm(10, 0, 1.2)
  loggits <- as.matrix(dat)%*%betas
  dat$p.score <- minmax(loggits, 0.1)
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  # default bart pred
  bart_fit_default <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10)
  pred_bd <- fitted(bart_fit_default)
  c_r2_bd <- r.squared(dat$p.score, pred_bd)
  min_bd <- min(pred_bd)
  max_bd <- max(pred_bd)
  # logit pred
  logistic <- glm(z ~ . -p.score, data = dat, family = binomial())
  pred_log <- fitted(logistic)
  c_r2_log <- r.squared(dat$p.score, pred_log)
  min_log <- min(pred_log)
  max_log <- max(pred_log)
  
  r <- cv_bart_k(dat, k)
  r$N <- as.integer(nv[i,'N'])
  r$rep_n <- as.integer(nv[i,'rep_n'])
  
  results2 <- rbind(results2, cbind(r, c_r2_bd, min_bd, max_bd, c_r2_log, min_log, max_log))
}

###### K Chisq
cv_bart_k_chisq <- function(dat, kdf){
  results <- data.frame(kdf=c(), calibration_r2 = c(),
                        minimum = c(), maximum = c())
  for(i in 1:length(kdf)){  
    bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10,
                              k = chi(degreesOfFreedom = kdf[i], scale = Inf))
    pred <- fitted(bart_fit)
    calibration_r2 <- r.squared(dat$p.score, pred)
    minimum <- min(pred)
    maximum <- max(pred)
    kdf_val = kdf[i]
    results <- rbind(results, cbind(kdf_val, calibration_r2, minimum, maximum))
  }
  return(results)
}

results3 <- data.frame()
ITERS <- 1000 ### Change this to run for X amount of times
N <- c(100, 200, 300, 500, 1000)
rep_n <- 1:ITERS
nv <- crossing(N=N, rep_n = rep_n)
kdf <- seq(1,10,0.1)

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.integer(nv[i,'N']), 10)
  betas <- rnorm(10, 0, 1.2)
  p_raw <- as.matrix(dat)%*%betas
  dat$p.score <- minmax(p_raw, 0.1)
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  # default bart pred
  bart_fit_default <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10)
  pred_bd <- fitted(bart_fit_default)
  c_r2_bd <- r.squared(dat$p.score, pred_bd)
  min_bd <- min(pred_bd)
  max_bd <- max(pred_bd)
  # logit pred
  logistic <- glm(z ~ . -p.score, data = dat, family = binomial())
  pred_log <- fitted(logistic)
  c_r2_log <- r.squared(dat$p.score, pred_log)
  min_log <- min(pred_log)
  max_log <- max(pred_log)
  
  r <- cv_bart_k_chisq(dat, kdf)
  r$N <- as.integer(nv[i,'N'])
  r$rep_n <- as.integer(nv[i,'rep_n'])
  
  results3 <- rbind(results3, cbind(r, c_r2_bd, min_bd, max_bd, c_r2_log, min_log, max_log))
}


### INTERACTIONS K fixed

results4 <- data.frame()
ITERS <- 1000 ### Change this to run for X amount of times
N <- c(100, 200, 300, 500, 1000)
rep_n <- 1:ITERS
nv <- crossing(N=N, rep_n = rep_n)
k <- seq(0.1,10,0.1) 

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.integer(nv[i,'N']), 10)
  dat2 <- dat
  dat2$I1 <- dat2[,1] * dat2[,2]
  dat2$I2 <- dat2[,3] * dat2[,4]
  dat2$I3 <- dat2[,5] * dat2[,6]
  betas <- rnorm(13, 0, 1.2)
  p_raw <- as.matrix(dat2)%*%betas
  dat$p.score <- minmax(p_raw, 0.1)
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  # default bart pred
  bart_fit_default <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10)
  pred_bd <- fitted(bart_fit_default)
  c_r2_bd <- r.squared(dat$p.score, pred_bd)
  min_bd <- min(pred_bd)
  max_bd <- max(pred_bd)
  # logit pred
  logistic <- glm(z ~ . -p.score, data = dat, family = binomial())
  pred_log <- fitted(logistic)
  c_r2_log <- r.squared(dat$p.score, pred_log)
  min_log <- min(pred_log)
  max_log <- max(pred_log)
  # logit interaction pred
  logistic_int <- glm(z ~ . + X1:X2 + X3:X4 + X5:X6 -p.score, data = dat, family = binomial())
  pred_log_int <- fitted(logistic_int)
  c_r2_log_int <- r.squared(dat$p.score, pred_log_int)
  min_log_int <- min(pred_log_int)
  max_log_int <- max(pred_log_int)
  
  r <- cv_bart_k(dat, k)
  r$N <- as.integer(nv[i,'N'])
  r$rep_n <- as.integer(nv[i,'rep_n'])
  
  results4 <- rbind(results4, cbind(r, c_r2_bd, min_bd, max_bd, c_r2_log, min_log, max_log,
                                    c_r2_log_int, min_log_int, max_log_int))
}

### INTERACTIONS K Chi-Squared

results5 <- data.frame()
ITERS <- 1000 ### Change this to run for X amount of times
N <- c(100, 200, 300, 500, 1000)
rep_n <- 1:ITERS
nv <- crossing(N=N, rep_n = rep_n)
kdf <- seq(0.1,10,0.1) 

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.integer(nv[i,'N']), 10)
  dat2 <- dat
  dat2$I1 <- dat2[,1] * dat2[,2]
  dat2$I2 <- dat2[,3] * dat2[,4]
  dat2$I3 <- dat2[,5] * dat2[,6]
  betas <- rnorm(13, 0, 1.2)
  p_raw <- as.matrix(dat2)%*%betas
  dat$p.score <- minmax(p_raw, 0.1)
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  # default bart pred
  bart_fit_default <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10)
  pred_bd <- fitted(bart_fit_default)
  c_r2_bd <- r.squared(dat$p.score, pred_bd)
  min_bd <- min(pred_bd)
  max_bd <- max(pred_bd)
  # logit pred
  logistic <- glm(z ~ . -p.score, data = dat, family = binomial())
  pred_log <- fitted(logistic)
  c_r2_log <- r.squared(dat$p.score, pred_log)
  min_log <- min(pred_log)
  max_log <- max(pred_log)
  # logit interaction pred
  logistic_int <- glm(z ~ . + X1:X2 + X3:X4 + X5:X6 -p.score, data = dat, family = binomial())
  pred_log_int <- fitted(logistic_int)
  c_r2_log_int <- r.squared(dat$p.score, pred_log_int)
  min_log_int <- min(pred_log_int)
  max_log_int <- max(pred_log_int)
  
  r <- cv_bart_k_chisq(dat, kdf)
  r$N <- as.integer(nv[i,'N'])
  r$rep_n <- as.integer(nv[i,'rep_n'])
  
  results5 <- rbind(results5, cbind(r, c_r2_bd, min_bd, max_bd, c_r2_log, min_log, max_log,
                                    c_r2_log_int, min_log_int, max_log_int))
}




