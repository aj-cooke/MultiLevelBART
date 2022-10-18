library(dplyr)
library(ggplot2)
library(MASS)
library(mvtnorm)


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
  return(data.frame(MASS::mvrnorm(n, mu = mu, Sigma = R, tol = 1)))
}

group_assign <- function(data, n_groups){
  group <- sample(1:n_groups, nrow(data), replace = TRUE)
  data <- cbind(data, group)
  return(data)
}

add_group_vars <- function(data, group_vars = 3){
  for(i in 1:group_vars){
    grouped <- data %>% group_by(group) %>% summarise_at(ncol(data)-i,mean)
    data <- merge(data[, !colnames(data) %in% c(paste0('X', ncol(data)-i))], grouped, by = 'group', all.x = T)
  }
  return(data)
}

expit <- function(x){return(exp(x)/(1+exp(x)))}

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
  grouped$propensity <- expit(grouped$propensity)
  grouped$z <- rbinom(nrow(grouped), 1, grouped$propensity)
  data <- merge(data, subset(grouped, select = c("group", "z", "propensity")), by = "group", all.x = T, sort = F)
  return(data)
}

linear_treatment_ind <- function(data){
  betas <- rnorm(ncol(data)-1)
  data$propensity <- as.vector(as.matrix(data[,2:ncol(data)]) %*% betas)
  data$propensity <- expit(data$propensity)
  data$z <- rbinom(nrow(data), 1, data$propensity)
  return(data)
}

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

gen_multilevel <- function(n, n_groups, k_ind, k_group, assignment = "group", rand_intercepts = F){
  data <- data.frame(generate_mvn(n, k_ind+k_group)) # first just ind level
  
  # randomly assign groups and group variables
  
  data <- group_assign(data, n_groups)
  
  # break out into group to de-standardize
  if(rand_intercepts == T){
    for(i in unique(data$group)){
      mus <- runif(k_ind)
      sigs <- runif(k_ind)
      data[data$group == i, 1:k_ind] = sweep(data[data$group == i, 1:k_ind], 2, sigs, "*")
      data[data$group == i, 1:k_ind] = sweep(data[data$group == i, 1:k_ind], 2, mus, "+")
    }
  }  
    # re-standardize
  data[,1:k_ind] = scale(data[,1:k_ind])
  # make group level vars
  
  data <- add_group_vars(data, k_group)
  
  # assign treatment
  if(assignment == "group"){data <- linear_treatment_group(data)}
  else if(assignment %in% c("ind", "individual")){data <- linear_treatment_ind(data)}
  else if(assignment %in% c("discrete")){data <- discrete_treatment_group(data)}
  else if(assignment %in% c("rec_ind")){data <- discrete_rec(data)}
  else if(assignment %in% c("rec_group")){data <- discrete_rec_group(data)}
  return(data)
}

# TESTING

# data1 = gen_multilevel(10000, 6, 7,3)
# data2 = gen_multilevel(10000, 6, 7,3, assignment = "ind")
# data3 = gen_multilevel(10000, 6, 7,3, assignment = "discrete")
# data4 = gen_multilevel(10000, 6, 7,3, assignment = "rec_ind")
# data5 = gen_multilevel(10000, 6, 7,3, assignment = "rec_group")

