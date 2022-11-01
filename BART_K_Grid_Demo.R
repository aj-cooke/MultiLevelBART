library(ggplot2)
library(dplyr)

bart_int_k_search <- read.csv('BART_int_k_search.csv')
bart_k_search <- read.csv('BART_k_search.csv')

### Linear

bart_k_search$default_diff <- bart_k_search$calibration_r2 - bart_k_search$c_r2_bd
bart_k_search$log_diff <- bart_k_search$calibration_r2 - bart_k_search$c_r2_log
bart_k_search$default_log_diff <- bart_k_search$c_r2_bd - bart_k_search$c_r2_log

p <- ggplot(bart_k_search, aes(x = k_val, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Calibration')
p

p <- ggplot(bart_k_search %>% filter(k_val > 1), aes(x = k_val, y = default_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. BART Default diff')
p

p <- ggplot(bart_k_search %>% filter(k_val > 1), aes(x = k_val, y = log_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Logistic Diff')
p

k_base <- bart_k_search %>%
  group_by(N, rep_n)%>%
  summarise(diff = mean(default_log_diff))

p <- ggplot(k_base, aes(x = N, y = diff, col = as.factor(rep_n))) +
  geom_line() +
  ggtitle('Default BART - Logisitc Calibration by N')
p

### Interaction

p <- ggplot(bart_int_k_search, aes(x = k_val, y = calibration_r2)) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Calibration')
p

bart_int_k_search$default_diff <- bart_int_k_search$calibration_r2 - bart_int_k_search$c_r2_bd
bart_int_k_search$log_diff <- bart_int_k_search$calibration_r2 - bart_int_k_search$c_r2_log
bart_int_k_search$log_int_diff <- bart_int_k_search$calibration_r2 - bart_int_k_search$c_r2_log_int
bart_int_k_search$default_log_diff <- bart_int_k_search$c_r2_bd - bart_int_k_search$c_r2_log
bart_int_k_search$default_log_int_diff <- bart_int_k_search$c_r2_bd - bart_int_k_search$c_r2_log_int

p <- ggplot(bart_int_k_search %>% filter(k_val > 1), aes(x = k_val, y = default_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. BART Default diff')
p

p <- ggplot(bart_int_k_search %>% filter(k_val > 1), aes(x = k_val, y = log_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Logistic Diff')
p

p <- ggplot(bart_int_k_search %>% filter(k_val > 1), aes(x = k_val, y = log_int_diff, col = as.factor(rep_n))) +
  geom_point() +
  facet_wrap(~ N) +
  ggtitle('K vs. Logistic Int. Diff')
p

res_int_k <- bart_int_k_search %>%
  group_by(N, rep_n)%>%
  summarise(diff_log = mean(default_log_diff), 
            diff_log_int = mean(default_log_int_diff))

p <- ggplot(res_int_k, aes(x = N, y = diff_log, col = as.factor(rep_n))) +
  geom_line() +
  ggtitle('Default BART - Logisitc Calibration by N')
p

p <- ggplot(res_int_k, aes(x = N, y = diff_log_int, col = as.factor(rep_n))) +
  geom_line() +
  ggtitle('Default BART - Logisitc Int. Calibration by N')
p