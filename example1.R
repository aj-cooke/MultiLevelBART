library(dbarts)
library(dplyr)
library(tidyr)
library(ggplot2)

# example 1 altering signal with one variable
expit <- function(x) exp(x)/(1+exp(x))
set.seed(7)
x <- cbind(rnorm(100),rnorm(100))
u <- 1*x[,1]+0*x[,2]
Y <- rbinom(100,1,expit(u))

fit1 <- bart2(Y ~ x)
fit2 <- bart2(Y ~x[,1])
fit3 <- rpart::rpart(Y ~ x)
tibble(p.score = expit(u), 
       prediction1 = fitted(fit1), 
       prediction2 = fitted(fit2)) %>% 
  pivot_longer(1:3) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)

x <- cbind(rnorm(100),rnorm(100))
u <- .5*x[,1]+0*x[,2]
Y <- rbinom(100,1,expit(u))

fit1 <- bart2(Y ~ x)
fit2 <- bart2(Y ~x[,1])
fit3 <- rpart::rpart(Y ~ x)
tibble(p.score = expit(u), 
       prediction1 = fitted(fit1), 
       prediction2 = fitted(fit2)) %>% 
  pivot_longer(1:3) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)


x <- cbind(rnorm(100),rnorm(100))
u <- .25*x[,1]+0*x[,2]
Y <- rbinom(100,1,expit(u))

fit1 <- bart2(Y ~ x)
fit2 <- bart2(Y ~x[,1])
fit3 <- rpart::rpart(Y ~ x)
tibble(p.score = expit(u), 
       prediction1 = fitted(fit1), 
       prediction2 = fitted(fit2)) %>% 
  pivot_longer(1:3) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)

# example 2 many weak signals
R <- matrix(NA, nrow = 5, ncol = 5)
diag(R) <- 1
R[lower.tri(R)] <- 0
R[upper.tri(R)] <- t(R[lower.tri(R)])
X <- MASS::mvrnorm(100, rep(0, 5), R)
u <- expit(X %*% c(rep(.2, 5)))
z <- rbinom(100, 1, u)
fit <- bart2(z ~ X)
tibble(p.score = u, 
       prediction1 = fitted(fit)) %>% 
  pivot_longer(1:2) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)

# lets add corelation

R <- matrix(NA, nrow = 5, ncol = 5)
diag(R) <- 1
R[lower.tri(R)] <- runif(5, min = -.5, max = .5)
R[upper.tri(R)] <- t(R[lower.tri(R)])
X <- MASS::mvrnorm(100, rep(0, 5), R)
u <- expit(X %*% c(rep(.2, 5)))
z <- rbinom(100, 1, u)
fit <- bart2(z ~ X)
tibble(p.score = u, 
       prediction1 = fitted(fit)) %>% 
  pivot_longer(1:2) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)

# example 3 many moderate signals
R <- matrix(NA, nrow = 5, ncol = 5)
diag(R) <- 1
R[lower.tri(R)] <- 0
R[upper.tri(R)] <- t(R[lower.tri(R)])
X <- MASS::mvrnorm(100, rep(0, 5), R)
u <- expit(X %*% c(rep(.4, 5)))
z <- rbinom(100, 1, u)
fit <- bart2(z ~ X)
tibble(p.score = u, 
       prediction1 = fitted(fit)) %>% 
  pivot_longer(1:2) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)

# lets add corelation

R <- matrix(NA, nrow = 5, ncol = 5)
diag(R) <- 1
R[lower.tri(R)] <- runif(5, min = -.5, max = .5)
R[upper.tri(R)] <- t(R[lower.tri(R)])
X <- MASS::mvrnorm(100, rep(0, 5), R)
u <- expit(X %*% c(rep(.4, 5)))
z <- rbinom(100, 1, u)
fit <- bart2(z ~ X)
tibble(p.score = u, 
       prediction1 = fitted(fit)) %>% 
  pivot_longer(1:2) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)


# adding interactions 
R <- matrix(NA, nrow = 5, ncol = 5)
diag(R) <- 1
R[lower.tri(R)] <- 0
R[upper.tri(R)] <- t(R[lower.tri(R)])
X <- MASS::mvrnorm(100, rep(0, 5), R)
temp <- rnorm(100)
interat <- cbind.data.frame(temp, X)
X.mat <- lm(temp ~ .^2, interat, x = TRUE)
X.mat <- X.mat$x
X.mat[, 2:ncol(X.mat)]
u <- expit(X.mat %*% c(rep(.4, 5),.4, .4,.4, rep(0, 8)))
z <- rbinom(100, 1, u)
fit <- bart2(z ~ X)
tibble(p.score = u, 
       prediction1 = fitted(fit)) %>% 
  pivot_longer(1:2) %>%  
  ggplot(aes(value, fill = name)) + 
  geom_histogram(position = 'identity', col = 'black', alpha = .4) + 
  facet_wrap(~name, ncol = 1)


