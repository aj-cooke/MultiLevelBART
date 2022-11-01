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

### CV for no groups
results <- data.frame()

N <- c(100, 500, 1000, 5000, 10000)
vars <- c(1, 3, 5, 10, 20)
nv <- crossing(N=N, vars=vars)

a <- c(0.95, 0.75, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
b <- c(2, 4, 6, 8, 10) # crashes if b too low
k <- c(1, 2, 4, 5, 6)
nt <- c(200)

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.vector(nv[i,'N'])$N, as.vector(nv[i,'vars'])$vars)
  betas <- rnorm(as.vector(nv[i,'vars'])$vars, .25, .6)
  loggits <- as.matrix(dat)%*%betas 
  dat$p.score <- as.vector(exp(loggits)/(1 + exp(loggits))) # should min-max instead of logit
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  r <- cv_bart(dat, a, b, k, nt)
  r$N <- as.integer(as.vector(nv[i,'N'])$N)
  r$vars <- as.integer(as.vector(nv[i,'vars'])$vars)
  results <- rbind(results, r)
}
####
### Visualize CV for each NV

p <- ggplot(grid_df, aes(x = a, y = calibration_r2)) +
  geom_point() +
  facet_grid(N ~ vars) +
  ggtitle('A vs. Calibration')
p

p <- ggplot(grid_df, aes(x = b, y = calibration_r2)) +
  geom_point() +
  facet_grid(N ~ vars) +
  ggtitle('B vs. Calibration')
p

p <- ggplot(grid_df, aes(x = k, y = calibration_r2)) +
  geom_point() +
  facet_grid(N ~ vars) +
  ggtitle('K vs. Calibration')
p

### parameters depend on individual datasets?
results2 <- data.frame()

N <- c(1000, 1000, 1000, 1000, 1000, 1000)
vars <- c(10, 10 , 10, 10, 10, 10)
nv <- data.frame(N=N, vars=vars)

a <- c(0.95, 0.75, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
b <- c(2, 4, 6, 8, 10)
k <- c(1, 2, 4, 5, 6) # add bigger grid, default with chisq, vary rows only, 100, 200, 300, 500, 1000, compare to logistic
nt <- c(200)

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.vector(nv[i,'N']), as.vector(nv[i,'vars']))
  betas <- rnorm(as.vector(nv[i,'vars']), .25, .6)
  loggits <- as.matrix(dat)%*%betas
  dat$p.score <- as.vector(exp(loggits)/(1 + exp(loggits))) # minmax
  
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  r <- cv_bart(dat, a, b, k, nt)
  r$iter <- i
  results2 <- rbind(results2, r)
}

p <- ggplot(results2, aes(x = a, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ iter) +
  ggtitle('A vs. Calibration')
p

p <- ggplot(results2, aes(x = b, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ iter) +
  ggtitle('B vs. Calibration')
p

p <- ggplot(results2, aes(x = k, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ iter) +
  ggtitle('K vs. Calibration')
p




######## Dig into K
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

results3 <- data.frame()

N <- c(100, 200, 300, 500, 1000)
rep_n <- c(1,2,3)
nv <- crossing(N=N, rep_n = rep_n)
k <- seq(0.1,10,0.1) # add bigger grid, default with chisq, vary rows only, 100, 200, 300, 500, 1000, compare to logistic

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
  
  results3 <- rbind(results3, cbind(r, c_r2_bd, min_bd, max_bd, c_r2_log, min_log, max_log))
}


p <- ggplot(results3, aes(x = k_val, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Calibration')
p

results3$default_diff <- results3$calibration_r2 - results3$c_r2_bd
results3$log_diff <- results3$calibration_r2 - results3$c_r2_log
results3$default_log_diff <- results3$c_r2_bd - results3$c_r2_log

p <- ggplot(results3 %>% filter(k_val > 1), aes(x = k_val, y = default_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. BART Default diff')
p

p <- ggplot(results3 %>% filter(k_val > 1), aes(x = k_val, y = log_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Logistic Diff')
p

res_base <- results3 %>%
  group_by(N, rep_n)%>%
  summarise(diff = mean(default_log_diff))

p <- ggplot(res_base, aes(x = N, y = diff, col = as.factor(rep_n))) +
  geom_line() +
  ggtitle('Default BART - Logisitc Calibration by N')
p

r100 <- results3 %>% filter(N == 100)
r200 <- results3 %>% filter(N == 200)
r500 <- results3 %>% filter(N == 500)


# will more chains help hyperprior perform better on low N?
### repeat for interaction
results4 <- data.frame()
N <- c(100, 200, 300, 500, 1000)
rep_n <- c(1,2,3)
nv <- crossing(N=N, rep_n = rep_n)
k <- seq(0.1,10,0.1) # add bigger grid, default with chisq, vary rows only, 100, 200, 300, 500, 1000, compare to logistic

for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.integer(nv[i,'N']), 10)
  dat2 <- dat
  dat2$I1 <- dat2[,1] * dat2[,2]
  dat2$I2 <- dat2[,3] * dat2[,4]
  dat2$I3 <- dat2[,5] * dat2[,6]
  betas <- rnorm(13, 0, 1.2)
  loggits <- as.matrix(dat2)%*%betas
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
  # logit interaction pred
  #logistic_int <- glm(z ~ .*. -p.score, data = dat, family = binomial())
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


p <- ggplot(results4, aes(x = k_val, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Calibration')
p

results4$default_diff <- results4$calibration_r2 - results4$c_r2_bd
results4$log_diff <- results4$calibration_r2 - results4$c_r2_log
results4$log_int_diff <- results4$calibration_r2 - results4$c_r2_log_int
results4$default_log_diff <- results4$c_r2_bd - results4$c_r2_log
results4$default_log_int_diff <- results4$c_r2_bd - results4$c_r2_log_int

p <- ggplot(results4 %>% filter(k_val > 1), aes(x = k_val, y = default_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. BART Default diff')
p

p <- ggplot(results4 %>% filter(k_val > 1), aes(x = k_val, y = log_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Logistic Diff')
p

p <- ggplot(results4 %>% filter(k_val > 1), aes(x = k_val, y = log_int_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Logistic Int. Diff')
p

res_base4 <- results4 %>%
  group_by(N, rep_n)%>%
  summarise(diff_log = mean(default_log_diff), 
            diff_log_int = mean(default_log_int_diff))

p <- ggplot(res_base4, aes(x = N, y = diff_log, col = as.factor(rep_n))) +
  geom_line() +
  ggtitle('Default BART - Logisitc Calibration by N')
p

p <- ggplot(res_base4, aes(x = N, y = diff_log_int, col = as.factor(rep_n))) +
  geom_line() +
  ggtitle('Default BART - Logisitc Int. Calibration by N')
p









### linear assignment but add random indicator splits







####

cv_bart_multilevel <- function(dat, dat_filtered, a, b, k, nt){
  params <- crossing(a = a, b = b, k = k, nt = nt)
  results_ind <- data.frame(a= c(), b=c(), k=c(), nt = c(), calibration_r2 = c(),
                            minimum = c(), maximum = c())
  results_group <- data.frame(a= c(), b=c(), k=c(), nt = c(), calibration_r2 = c(),
                              minimum = c(), maximum = c())
  
  for(i in 1:nrow(params)){
    # bart 
    bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10,
                              n.trees = params$nt[i],
                              k = params$k[i],
                              power = params$b[i],
                              base = params$a[i])
    ind_pred <- fitted(bart_fit)
    calibration_r2 <- r.squared(dat$p.score, ind_pred)
    minimum <- min(ind_pred)
    maximum <- max(ind_pred)
    results_ind <- rbind(results_ind, cbind(params[i,], calibration_r2, minimum, maximum))
    
    bart_fit <- dbarts::bart2(z ~ . -p.score -j, data = dat_filtered, n.chains = 10,
                              n.trees = params$nt[i],
                              k = params$k[i],
                              power = params$b[i],
                              base = params$a[i])
    group_pred <- fitted(bart_fit)
    calibration_r2 <- r.squared(dat_filtered$p.score, group_pred)
    minimum <- min(group_pred)
    maximum <- max(group_pred)
    results_group <- rbind(results_group, cbind(params[i,], calibration_r2, minimum, maximum))
    print(params[i,], calibration_r2, minimum, maximum)
  }
  return(list(results_group, results_ind))
}

results_multi <- list()

ind <- c(1000, 1000, 1000, 1000, 10000, 10000, 10000, 10000)
grp <- c(40, 40, 250, 250, 80, 80, 500, 500)

a <- c(0.95, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
b <- c(2, 4, 6, 8, 10) # crashes if b too low
k <- c(1, 2, 4, 5, 6)
nt <- c(200)

for(i in 1:length(ind)){
  dat <- gen_multilevel(ind[i], grp[i], 7, 3)
  colnames(dat) <- c('j', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'z', 'p.score')
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  r <- cv_bart_multilevel(dat, dat_filtered, a, b, k, nt)
  results_multi <- append(results_multi, r)
}







