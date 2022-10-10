library(dplyr)
library(ggplot2)
library(MASS)
library(mvtnorm)

# correlation matrix is causing errors if tol is at default of 1e-6
# increasing tol to 1 changes correlations in resulting data from input

generate_mvn <- function(n, k, mu = "none", R = "none"){
  if(R == "none"){
    R <- matrix(nrow = k, ncol = k)
    diag(R) <- rep(1, k)
    off_diag <- c()
    
    for(i in rev(1:(k-1))){
      off_diag <- c(off_diag, runif(i, -0.99, 0.99))
    }
    
    R[lower.tri(R)] <- off_diag
    R[upper.tri(R)] <- t(R)[upper.tri(t(R))]
  }
  
  if(mu == "none"){
    mu <- rep(0, k)
  }
  return(MASS::mvrnorm(n, mu = mu, Sigma = R, tol = 1))
}

# best way to make group level variables correlated with individual level?

# ONE APPROACH:
# Make individual level dataset 
# Randomly assign groups 
# within groups, "de normalize" all variables with random mu and sigma
# re normalize columns across groups
# make each group level some function + error of an aggregated individual level variable

group_assign <- function(x, group_list){
  return(which(group_list == min(group_list[group_list >= x])))
}

# propensity is arbitrary number min-max scaled. Could change to expit
linear_treatment_group <- function(data){
  grouped <- data %>%
    group_by(group) %>%
    summarise(across(
      .cols = where(is.numeric), 
      .fns = list(Mean = mean), na.rm = TRUE, 
      .names = "{col}_{fn}"
    ))
  grouped <- data.frame(grouped)
  betas <- rnorm(ncol(grouped) - 1)
  grouped$propensity <- as.vector(as.matrix(grouped[,2:ncol(grouped)]) %*% betas)
  grouped$propensity <- (grouped$propensity-min(grouped$propensity))/(max(grouped$propensity)-min(grouped$propensity))
  grouped$z <- rbinom(nrow(grouped), 1, grouped$propensity)
  data <- merge(data, subset(grouped, select = c("group", "z")), by = "group", all.x = T, sort = F)
  return(data)
}

linear_treatment_ind <- function(data){
  betas <- rnorm(ncol(data)-1)
  data$propensity <- as.vector(as.matrix(data[,2:ncol(data)]) %*% betas)
  data$propensity <- (data$propensity-min(data$propensity))/(max(data$propensity)-min(data$propensity))
  data$z <- rbinom(nrow(data), 1, data$propensity)
  return(data)
}

# can enhance by recursively splitting.
# This approach may require making every variable share a positive relationship with treatment so splits can be made intuitive
discrete_treatment_group <- function(data){
  grouped <- data %>%
    group_by(group) %>%
    summarise(across(
      .cols = where(is.numeric), 
      .fns = list(Mean = mean), na.rm = TRUE, 
      .names = "{col}_{fn}"
    ))
  grouped <- data.frame(grouped)
  
  selection_col <- sample.int(ncol(grouped),1)
  cutoff <- quantile(grouped[,selection_col], probs = runif(1,0.25,0.75))
  grouped$z <- if_else(grouped[,selection_col] >= cutoff, 1, 0)
  data <- merge(data, subset(grouped, select = c("group", "z")), by = "group", all.x = T, sort = F)
  return(data)
}

# can add saving rules, sometimes errors if splits go haywire but usually ok
discrete_rec <- function(data, depth = 0, max_depth = 2, treat = 0){
  print(data)
  print(depth)
  if(depth < max_depth){
    selection_col <- sample.int(ncol(data),1)
    data[,"selection_col"] <- data[,selection_col]
    cutoff <- quantile(data$selection_col, probs = runif(1,0.25,0.75))
    right_df <- data %>% filter(selection_col < cutoff) %>% dplyr::select(-selection_col)
    left_df <- data %>% filter(selection_col >= cutoff) %>% dplyr::select(-selection_col)
    return(rbind(discrete_rec(right_df, depth + 1, max_depth), 
                 discrete_rec(left_df, depth + 1, max_depth, treat = 1)))
  }
  else{
    data$z <- treat
    return(data)
  }
}

discrete_rec <- function(data, depth = 0, max_depth = 2, treat = 0){
  if(depth < max_depth){
    selection_col <- sample.int(ncol(data),1)
    data[,"selection_col"] <- data[,selection_col]
    cutoff <- quantile(data$selection_col, probs = runif(1,0.25,0.75))
    right_df <- data %>% filter(selection_col < cutoff) %>% dplyr::select(-selection_col)
    left_df <- data %>% filter(selection_col >= cutoff) %>% dplyr::select(-selection_col)
    if(nrow(right_df) == 0){
      return(discrete_rec(left_df, depth + 1, max_depth, treat = 1))
    }
    else if(nrow(left_df) == 0){
      return(discrete_rec(right_df, depth + 1, max_depth))
    }
    else{
      return(rbind(discrete_rec(right_df, depth + 1, max_depth), 
                   discrete_rec(left_df, depth + 1, max_depth, treat = 1)))
    }
  }
  else{
    data$z <- treat
    return(data)
  }
}

discrete_rec_group <- function(data){
  grouped <- data %>%
    group_by(group) %>%
    summarise(across(
      .cols = where(is.numeric), 
      .fns = list(Mean = mean), na.rm = TRUE, 
      .names = "{col}_{fn}"
    ))
  grouped <- data.frame(grouped)
  grouped <- discrete_rec(grouped)
  data <- merge(data, subset(grouped, select = c("group", "z")), by = "group", all.x = T, sort = F)
  return(data)
}

gen_multilevel <- function(n, n_groups, k_ind, k_group, group_p = "none", assignment = "group"){
  data <- data.frame(generate_mvn(n, k_ind)) # first just ind level
  
  # randomly assign groups
  
  if(group_p == "none"){
    group_p <- runif(n_groups)
    group_p <- cumsum(group_p / sum(group_p))
  }
  data$index <- 1:nrow(data) / nrow(data)
  data$group <- lapply(data$index, group_assign, group_p)
  data <- subset(data, select = -index)
  
  # break out into group to de-standardize
  
  for(i in unique(data$group)){
    mus <- runif(k_ind)
    sigs <- runif(k_ind)
    data[data$group == i, 1:k_ind] = sweep(data[data$group == i, 1:k_ind], 2, sigs, "*")
    data[data$group == i, 1:k_ind] = sweep(data[data$group == i, 1:k_ind], 2, mus, "+")
  }
  
  # re-standardize
  data[,1:k_ind] = scale(data[,1:k_ind])
  
  # make group level vars a function of the mean of a random ind var
  
  for(i in 1:k_group){
    name <- paste0("GV", i)
    g_coef <- rnorm(1)
    var <- sample(x=1:k_ind, size=1)
    grouped <- data %>% group_by(group) %>% summarise(vals = mean(X1)) # i
    grouped$vals <- grouped$vals*g_coef + rnorm(n_groups, 0, 0.2)
    grouped[,name] = grouped[,'vals']
    grouped <- data.frame(grouped)
    grouped <- subset(grouped, select = c("group", name))
    data <- merge(data, grouped, by.x = "group", by.y = "group", all.x = T, sort = F)
  }
  
  # assign treatment
  if(assignment == "group"){data <- linear_treatment_group(data)}
  else if(assignment %in% c("ind", "individual")){data <- linear_treatment_ind(data)}
  else if(assignment %in% c("discrete")){data <- discrete_treatment_group(data)}
  else if(assignment %in% c("rec_ind")){data <- discrete_rec(data)}
  else if(assignment %in% c("rec_group")){data <- discrete_rec_group(data)}
  return(data)
}

# TESTING

data1 = gen_multilevel(10000, 6, 7,3)
data2 = gen_multilevel(10000, 6, 7,3, assignment = "ind")
data3 = gen_multilevel(10000, 6, 7,3, assignment = "discrete")
data4 = gen_multilevel(10000, 6, 7,3, assignment = "rec_ind")
data5 = gen_multilevel(10000, 6, 7,3, assignment = "rec_group")


num_ind_invalid <- 0
num_group_invalid <- 0

num_ind_broke <- 0
num_group_broke <- 0

for(i in 1:100){
  data4 = gen_multilevel(10000, 6, 7,3, assignment = "rec_ind")
  print('ind')
  data5 = gen_multilevel(10000, 6, 7,3, assignment = "rec_group")
  print('group')
  if(nrow(data4 %>% filter(z==1)) == 0){num_ind_invalid = num_ind_invalid+1}
  if(nrow(data4 %>% filter(z==0)) == 0){num_ind_invalid = num_ind_invalid+1}
  if(nrow(data5 %>% filter(z==1)) == 0){num_group_invalid = num_group_invalid+1}
  if(nrow(data5 %>% filter(z==0)) == 0){num_group_invalid = num_group_invalid+1}
  
  if(nrow(data4) != 10000){num_ind_broke = num_ind_broke + 1}
  if(nrow(data5) != 10000){num_group_broke = num_group_broke + 1}
  print(i)
}


