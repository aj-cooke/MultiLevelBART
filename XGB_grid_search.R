library(xgboost) # use gbm package
library(dplyr)
library(MASS)
library(tidyr)
source('multilevel_sim.R')

r.squared <- function(ytrue, ypred){
  return(1-(sum((ytrue-ypred)^2)/sum((ytrue-mean(ytrue))^2)))
}

cv_xgb <- function(dat, max_depth, n_trees, lr, cs, g){
  params <- crossing(max_depth = max_depth, n_trees = n_trees, lr = lr, col_sample = cs, gamma = g)
  results <- data.frame(max_depth= c(), n_trees=c(), lr=c(), 
                        col_sample = c(),
                        gamma = c(),
                        calibration_r2 = c(),
                        minimum = c(), maximum = c())
  for(i in 1:nrow(params)){  
    xgb_train <- xgb.DMatrix(data = data.matrix(dat %>% dplyr::select(-p.score, -z)), label = dat$z)
    fit <- xgb.train(data = xgb_train, 
                     nrounds = as.integer(as.vector(params[i,'n_trees'])), 
                     max_depth = as.vector(params[i,'max_depth']),
                     eta = as.vector(params[i,'lr']),
                     colsample_bytree = as.vector(params[i, 'col_sample']),
                     gamma = as.vector(params[i, 'gamma']))
    pred <- predict(fit, xgb_train, type = 'prob')
    calibration_r2 <- r.squared(dat$p.score, pred)
    minimum <- min(pred)
    maximum <- max(pred)
    results <- rbind(results, cbind(params[i,], calibration_r2, minimum, maximum))
    
  }
  return(results)
}

results_xgb <- data.frame()

N <- c(1000)
vars <- c(10)
nv <- data.frame(N=N, vars=vars)

md <- c(2, 4, 6)
nt <- c(100, 500, 1000)
lr <- c(0.01, 0.05, 0.1, 0.2, 0.4)
cs <- c(0.5, 0.75, 1)
g <- c(0, 5, 10)


for(i in 1:nrow(nv)){
  dat <- generate_mvn(as.vector(nv[i,'N']), as.vector(nv[i,'vars']))
  betas <- rnorm(as.vector(nv[i,'vars']), .25, .6)
  loggits <- as.matrix(dat)%*%betas
  dat$p.score <- as.vector(exp(loggits)/(1 + exp(loggits)))
  hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')
  dat$z <- rbinom(nrow(dat),1,dat$p.score)
  r <- cv_xgb(dat, md, nt, lr, cs, g)
  r$iter <- i
  results_xgb <- rbind(results_xgb, r)
}


p <- ggplot(results_xgb, aes(x = md, y = calibration_r2)) +
  geom_point() +
  ggtitle('Max Depth vs. Calibration')
p

p <- ggplot(results_xgb, aes(x = nt, y = calibration_r2)) +
  geom_point() +
  ggtitle('N Trees vs. Calibration')
p
p <- ggplot(results_xgb, aes(x = lr, y = calibration_r2)) +
  geom_point() +
  ggtitle('Learning Rate vs. Calibration')
p
p <- ggplot(results_xgb, aes(x = cs, y = calibration_r2)) +
  geom_point() +
  ggtitle('Col Subsample vs. Calibration')
p
p <- ggplot(results_xgb, aes(x = g, y = calibration_r2)) +
  geom_point() +
  ggtitle('Gamma vs. Calibration')
p


