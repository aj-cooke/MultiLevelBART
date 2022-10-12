dat <- read.csv('~/Downloads/classroom.csv') # change path to fit your env
y <- dat$mathgain
dat <- dat %>% dplyr::select(-mathgain, -mathknow, -classid, -yearstea)
dat <- dat %>% 
  group_by(schoolid) %>% 
  mutate(avg_mathprep = mean(mathprep), 
         avg_mathkind = mean(mathkind))


# rank schools
variables <- c('avg_mathkind', 'avg_mathprep')
for (i in 1:length(variables)) {
  var <- variables[i]
  dat$temp <- dat[[var]]
  dat <-dat %>% 
    group_by(schoolid) %>%
    dplyr::select(temp, schoolid) %>% 
    distinct() %>% 
    ungroup() %>% 
    mutate(rank = dense_rank(temp)) %>% 
    right_join(dat)
  names(dat)[which(names(dat) == 'rank')] <- paste0(var, '_rank')
  dat <- dat %>% dplyr::select(-temp)
}

# now house pov (opposite directioin)
dat <- dat %>% 
  group_by(schoolid) %>%
  dplyr::select(housepov, schoolid) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(housepov_rank = dense_rank(desc(housepov))) %>% 
  right_join(dat)
  

dat %>% 
  count(housepov_rank < 40, 
        avg_mathprep_rank < 40, 
        avg_mathkind_rank < 40)

dat <- dat %>% 
  mutate(p.score = 
           case_when(
             avg_mathkind_rank < 40 & avg_mathprep_rank < 40 & housepov_rank < 40 ~ .8, 
             avg_mathkind_rank < 40 & avg_mathprep_rank < 40 ~ .7, 
             avg_mathkind_rank < 40 & housepov_rank < 40 ~ .65, 
             housepov_rank < 40 & avg_mathprep_rank < 40 ~ .6, 
             avg_mathkind_rank < 40 ~ .4, 
             avg_mathprep_rank < 40 ~ .3, 
             housepov_rank < 40 ~ .5, 
             TRUE ~ .2))


hist(dat$p.score)

set.seed(2)
dat <- dat %>% 
  group_by(schoolid) %>% 
  mutate(z = rbinom(1, 1, p.score)) %>% 
  ungroup()

# bart with "raw data"
fit <-dbarts::bart2(z~ housepov + mathkind + mathprep, data = dat)
hist(fitted(fit))

fit <-dbarts::bart2(z~ housepov + avg_mathkind + avg_mathprep, data = dat)
hist(fitted(fit))

# let give it just the means
dat_filtered <- dat %>% 
  select(z, housepov, avg_mathprep, avg_mathkind) %>% 
  distinct()
fit <-dbarts::bart2(z~ ., data = dat_filtered)
hist(fitted(fit))

# lets give it the ranks (big time cheating)
dat_rank <- dat %>% 
  select(z, housepov_rank, avg_mathprep_rank, avg_mathkind_rank) %>% 
  distinct()
fit <-dbarts::bart2(z~ ., data = dat_rank)
hist(fitted(fit))




