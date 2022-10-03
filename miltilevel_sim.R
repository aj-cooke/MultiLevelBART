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

data <- generate_mvn(10000, 10)

# best way to make group level variables correlated with individual level?

# ONE APPROACH:
# Make individual level dataset X
# Randomly assign groups X
# within groups, "de normalize" all variables with random mu and sigma
# re normalize columns across groups
# make each group level some function + error of an aggregated individual level variable

group_assign <- function(x, group_list){
  return(which(group_list == min(group_list[group_list >= x])))
}

gen_multilevel <- function(n, n_groups, k_ind, k_group, group_p = "none"){
  data <- data.frame(generate_mvn(n, k_ind)) # first just ind level
  
  # randomly assign groups
  
  if(group_p == "none"){
    group_p <- runif(n_groups)
    group_p <- cumsum(group_p / sum(group_p))
    print(group_p)
  }
  data$index <- 1:nrow(data) / nrow(data)
  data$group <- lapply(data$index, group_assign, group_p)
  
  # break out into group to de-standardize
  
  for( i in unique(data$group)){
    mus <- runif(k_ind)
    sigs <- runif(k_ind)
    data[data$group == i, 1:k_ind] = sweep(data[data$group == i, 1:k_ind], 2, sigs, "*")
    data[data$group == i, 1:k_ind] = sweep(data[data$group == i, 1:k_ind], 2, mus, "+")
  }
  
  # re-standardize
  data[,1:k_ind] = scale(data[,1:k_ind])
  
  # make group level vars a function of the mean of a random ind var
  
   
  
  return(data)
}


data = gen_multilevel(100, 4, 3,3)





