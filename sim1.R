library(dbarts)
library(rstanarm)
library(MASS)
library(dplyr)

# basic
set.seed(2)
R <- matrix(NA, nrow = 10, ncol = 10)
diag(R) <- 1
R[lower.tri(R)] <- 0
R[upper.tri(R)] <- t(R)[upper.tri(R)]
X <- MASS::mvrnorm(1000, rep(0, 10), Sigma = R)
j <- sample(1:40, 1000, replace = TRUE)
dat <- data.frame(X, j)
dat <- dat %>% 
  group_by(j) %>% 
  mutate_at(vars(1:10), mean) 

betas <- rnorm(10, .25, .6)

loggits <- as.matrix(dat[, 1:10])%*%betas 
dat$p.score <- as.vector(exp(loggits)/(1 + exp(loggits)))
hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')

dat$j <- as.factor(dat$j)

dat <- dat %>% 
  dplyr::select(p.score, j) %>% 
  distinct() %>% 
  mutate(z = rbinom(1, 1, p.score)) %>% 
  full_join(dat) %>% 
  ungroup()
  

dat_filtered <- dat %>% distinct()


# loggistic regression 
fit <- glm(z ~ . -p.score -j, family = binomial, data = dat)
hist(fitted(fit))


fit <- stan_glmer(z ~ . -p.score - j + (1|j), 
                  family = binomial, 
                  data = dat, 
                  cores = 4)
hist(fitted(fit))

fit <- glm(z ~ . -p.score -j, family = binomial, data = dat_filtered)
hist(fitted(fit))

fit <- glm(z ~ . -p.score, family = binomial, data = dat_filtered)
hist(fitted(fit))

fit <- glm(dat$z ~ X, family = binomial)
hist(fitted(fit))


fit <- glm(dat$z ~ as.factor(j), family = binomial)
hist(fitted(fit))

fit <- glm(dat$z ~ X + as.factor(j), family = binomial)
hist(fitted(fit))

# bart 
bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat, n.chains = 10)
hist(fitted(bart_fit))

bart_fit <- dbarts::bart2(z ~ . -p.score -j, data = dat, n.chains = 10)
hist(fitted(bart_fit))

bart_fit <- dbarts::bart2(z ~ . -p.score -j, data = dat_filtered, n.chains = 10)
hist(fitted(bart_fit))

bart_fit <- dbarts::bart2(z ~ . -p.score, data = dat_filtered, n.chains = 10)
hist(fitted(bart_fit))

# bart with 5 number summary 
tapply(X[, 1],dat$j, summary)
five_number <- lapply(1:ncol(X), function(i){
  temp <- tapply(X[,i], dat$j, summary) %>% 
    bind_rows()
  names(temp) <- c(paste0('min_', i),
                   paste0('q25_', i), 
                   paste0('median_', i), 
                   paste0('mean_', i), 
                   paste0('q75_', i), 
                   paste0('max_', i))
  
  return(temp)
  })

five_number <- five_number %>% 
  bind_cols() %>% 
  mutate_all(as.double)
ordered_z <- dat %>% 
  group_by(j) %>% 
  select(z, j) %>% 
  arrange(j) %>% 
  distinct() %>% 
  ungroup() %>% 
  select(z) 

five_number$z <- ordered_z$z
rm(ordered_z)

bart_fit <- dbarts::bart2(z ~ . , data = five_number, n.chains = 10)
hist(fitted(bart_fit))
