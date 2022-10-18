library(dbarts)
library(rstanarm)
library(dplyr)
source('multilevel_sim.R')

dat <- gen_multilevel(1000, 40, 7, 3)
colnames(dat) <- c('j', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'z', 'p.score')
hist(dat$p.score, main = 'True propensity score', xlab = 'p.score')

dat_filtered <- dat %>% 
  group_by(j) %>% 
  summarise_at(vars(2:12), mean)

X <- as.matrix(dat[,2:11])
j <- dat$j

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
  mutate_all(as.double) %>% 
  dplyr::select(-contains('mean'))

ordered_z <- dat %>% 
  group_by(j) %>% 
  dplyr::select(z, j) %>% 
  arrange(j) %>% 
  distinct() %>% 
  ungroup() %>% 
  dplyr::select(z) 

five_number$z <- ordered_z$z
rm(ordered_z)

bart_fit <- dbarts::bart2(z ~ . , data = five_number, n.chains = 10)
hist(fitted(bart_fit))

