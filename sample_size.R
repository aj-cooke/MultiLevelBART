
# basic dpg requiors you have MASS installed
library(dplyr)
library(dbarts)
library(ggplot2)


results <- list()
groups <- c(50, 100, 200, 500)
cases <- c(500, 1000,2000, 5000)
run <- tidyr::crossing(cases, groups)

for (i in 1:nrow(run)) {
  n <- run[[i, 1]]
  g <- run[[i, 2]]
  set.seed(2)
  R <- matrix(NA, nrow = 10, ncol = 10)
  diag(R) <- 1
  R[lower.tri(R)] <- 0
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  X <- MASS::mvrnorm(n, rep(0, 10), Sigma = R)
  j <- sample(1:g, n, replace = TRUE)
  dat <- data.frame(X, j)
  dat <- dat %>% 
    group_by(j) %>% 
    mutate_at(vars(1:10), mean) 
  
  betas <- rnorm(10, .3, .2)
  
  loggits <- as.matrix(dat[, 1:10])%*%betas 
  dat$p.score <- as.vector(exp(loggits)/(1 + exp(loggits)))
  
  rand <- lapply(1:1000, function(j){
    set.seed(j)
    dat <- dat %>% 
      dplyr::select(p.score, j) %>% 
      distinct() %>% 
      mutate(z = rbinom(1, 1, p.score)) %>% 
      full_join(dat) %>% 
      ungroup()
    
    dat_filtered <- dat %>% distinct()
    fit <- bart2(z~ .-p.score -j, dat_filtered)
    data.frame(p.score = dat_filtered$p.score, 
               prdiction = fitted(fit), 
               group_size = g, 
               sample_size = n)
    
  })
  
  results[[i]] <- rand
  
}


df_results <- results %>% 
  bind_rows()
df_results %>% 
  mutate(dgp = paste0(group_size,' groups ', sample_size, ' cases')) %>% 
  group_by(p.score, dgp) %>% 
  mutate(prediction = mean(prdiction)) %>% 
  select(prediction, p.score, dgp) %>% 
  distinct() %>% 
  ggplot(aes(prediction, p.score, col = dgp)) + 
  geom_smooth(method = 'lm', se = F) + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) + 
  theme_bw() + facet_wrap(~dgp) + 
  theme(legend.position = 'none')


