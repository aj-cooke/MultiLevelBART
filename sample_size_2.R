# no grouping just sample size
sample_size <- c(50, 100, 200, 500, 1000, 5000)
results <- list()
for (i in 1:length(sample_size)) {
  set.seed(2)
  n <- sample_size[i]
  R <- matrix(NA, nrow = 10, ncol = 10)
  diag(R) <- 1
  R[lower.tri(R)] <- 0
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  X <- MASS::mvrnorm(n, rep(0, 10), Sigma = R)
  
  betas <- rnorm(10, .2, .1)
  
  p <- X%*%betas 
  p.score <- pnorm(p)
  
  rand <- lapply(1:1000, function(i){
    set.seed(i)
    z <- rbinom(n, 1, p.score)
    fit <- bart2(z~X)
    prediction <- fitted(fit)
    data.frame(p.score, prediction, n = length(z))
  })
  
  results[[i]] <- rand
}

results_df <- results %>% bind_rows()

results_df %>% 
  group_by(n, p.score) %>% 
  mutate(prediction = mean(prediction)) %>% 
  ggplot(aes(prediction, p.score, col = as.factor(n))) + 
  geom_smooth(method = 'lm', se = F) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw()

readr::write_rds(results, 'sample_size_study2.rds')
