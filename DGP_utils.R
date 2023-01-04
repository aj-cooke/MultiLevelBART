library(dplyr)

########################## P SCORE UTILITIES ##############################

shrink <- function(x, shrink){return((1-shrink*2)*x+shrink)}
expit <- function(x){return(exp(x)/(1+exp(x)))}
minmax <- function(x) {return((x- min(x)) /(max(x)-min(x)))}

############ TRANSFORMATIONS ON ENTIRE DGP ####################

exp_dgp <- function(data){
  betas <- rnorm(ncol(data), 0, 1)
  return(exp(as.matrix(data)%*%betas))
}

sqrt_dgp <- function(data){
  betas <- rnorm(ncol(data), 0, 1)
  res <- as.matrix(data) %*% betas
  signs <- if_else(res < 0 , -1, 1)
  return(sqrt(abs(res))*signs)
}

log_dgp <- function(data){
  betas <- rnorm(ncol(data), 0, 1)
  res <- as.matrix(data) %*% betas
  signs <- if_else(res < 0 , -1, 1)
  return(log(abs(res))*signs)
}

linear_dgp <- function(data){
  betas <- rnorm(ncol(data), 0, 1)
  return(as.matrix(data) %*% betas)
}

parse_dgp <- function(equation, data, transform="linear"){
  equation <- gsub(" ", "", equation, fixed = TRUE)
  lr <- unlist(strsplit(equation, "~", fixed = TRUE))
  rh <- unlist(strsplit(as.character(lr[2]), "+", fixed = T))
  out_col <- lr[1]
  used_cols <- c()
  for(i in rh){
    # polynomial
    if(substr(i, 1, 2) == "I("){
      pnom <- as.numeric(substr(i,nchar(i), nchar(i)))
      col <- substr(i, 3, nchar(i)-3)
      data[,paste0(col, "^", pnom)] <- data[,col]^pnom
      used_cols <- c(used_cols, paste0(col, "^", pnom))
    }
    # * interaction
    else if(grepl("*", i, fixed=TRUE)){
      loc <- unlist(strsplit(i, "*", fixed = T))
      c1 <- loc[1]
      c2 <- loc[2]
      data[,paste0(c1,"*",c2)] <- data[,c1]*data[,c2]
      used_cols <- c(used_cols, c1, c2, paste0(c1,"*",c2))
    }
    # : interaction
    else if(grepl(":", i, fixed=TRUE)){
      loc <- unlist(gregexpr(':', i))[1]
      c1 <- substr(i,1,loc-1)
      c2 <- substr(i, loc+1, nchar(i))
      data[,paste0(c1,"*",c2)] <- data[,c1]*data[,c2]
      used_cols <- c(used_cols, paste0(c1,"*",c2))
    }
    # sin
    else if((nchar(i) > 5) & (substr(i, 1, 4) == "sin(")){
      col <- substr(i, 5, nchar(i)-1)
      data[,paste0("sin_", col)] <- sin(data[,col])
      used_cols <- c(used_cols, paste0("sin_", col))
    }
    # sqrt
    else if((nchar(i) > 6) & (substr(i, 1, 5) == "sqrt(")){
      col <- substr(i, 6, nchar(i)-1)
      signs <- if_else(data[,col] < 0, -1, 1)
      data[,paste0("sqrt_", col)] <- sqrt(abs(data[,col]))* signs
      used_cols <- c(used_cols, paste0("sqrt_", col))
    }
    # log
    else if((nchar(i) > 5) & (substr(i, 1, 4) == "log(")){
      col <- substr(i, 5, nchar(i)-1)
      signs <- if_else(data[,col] < 0, -1, 1)
      data[,paste0("log_", col)] <- log(abs(data[,col]))* signs
      used_cols <- c(used_cols, paste0("log_", col))
    }
    # exp
    else if((nchar(i) > 5) & (substr(i, 1, 4) == "exp(")){
      col <- substr(i, 5, nchar(i)-1)
      data[,paste0("exp_", col)] <- exp(data[,col])
      used_cols <- c(used_cols, paste0("exp_", col))
    }
    else{used_cols <- c(used_cols, i)}
  }
  dat_dgp <- subset(data, select = used_cols)
  if(transform == "exp"){out <- exp_dgp(dat_dgp)}
  else if(transform == "sqrt"){out <- sqrt_dgp(dat_dgp)}
  else if(transform == "log"){out <- log_dgp(dat_dgp)}
  else if(transform == "linear"){out <- linear_dgp(dat_dgp)}
  else{out <- linear_dgp(dat_dgp)}
  data[,out_col] <- as.vector(out)
  return(data)
}

############ TRANSFORMATIONS ON COLUMNS #######################

add_ints <- function(data, num=3){
  for(i in 1:num){
    c1 <- sample(1:ncol(data), 1)
    c2 <- sample(c(1:c1-1, c1+1:ncol(data)), 1)
    data[,paste0("interaction_", i)] = data[,c1]*data[,c2]
  }
  return(data[,(ncol(data)-(num-1)):ncol(data)])
}

add_polynoms <- function(data, p2 = 2, p3 = 1){
  p2cols <- sample(1:ncol(data),p2)
  p3cols <- sample(p2cols, p3)
  for(i in 1:length(p2cols)){
    data[,paste0("poly2_", i)] <- data[,p2cols[i]]^2
    if(i<=length(p3cols)){
      data[,paste0("poly3_", i)] <- data[,p3cols[i]]^2
    }
  }
  return(data[,(ncol(data)-(p2+p3-1)):ncol(data)])
}

add_indicators <- function(data, num = 2){
  cols <- sample(1:ncol(data), num)
  for(i in 1:length(cols)){
    cut <- rnorm(1, 0.5, 0.166)
    cut <- if_else(cut <= 0.05, 0.05, if_else(cut>=0.95, 0.95, cut))
    cut <- as.numeric(quantile(data[,cols[i]],cut))
    data[,paste0("ind_", i)] <- if_else(data[,cols[i]] >= cut, 1, 0)
  }
  return(data[,(ncol(data)-(num-1)):ncol(data)])
}

add_trig <- function(data, num=2){
  cols <- sample(1:ncol(data), num)
  for(i in 1:length(cols)){
    beta1 <- rnorm(1, 0, 4)
    beta2 <- rnorm(1, 1, 0.3333)
    data[,paste0("sin_", i)] <- beta1*sin(beta2*data[,cols[i]])
  }
  return(data[,(ncol(data)-(num-1)):ncol(data)])  
}

################### Adding Bias

remove_cor <- function(data, ranks){
  cmr <- as.data.frame(rank(colMeans(abs(cor(data)))))
  colnames(cmr) <- c('rank')
  cmr <- cmr %>% filter(!rank %in% ranks)
  keep <- rownames(cmr)
  return(subset(data, select = keep))
}

remove_groups <- function(data, cols){
  for(i in cols){
    if(is.numeric(data[,i])){
      percs <- runif(2,1,100)
      min_p <- min(percs)
      max_p <- max(percs)
      data$ind <- 0
      data$ind[ntile(data[,i], 100) > min_p & ntile(data[,i], 100) < max_p] <- 1
      data <- data %>% filter(ind == 1)
    }
    else{
      cat <- sample(unique(data[,i]),1)
      data$ind <- 0
      data$ind[data[,i] == cat] <- 1
      data <- data %>% filter(ind == 1)
    }
  }
  data <- subset(data, select = -ind)
  return(data)
}


############### DEMO
# N <- 1000
# data <- data.frame("X1" = rnorm(N), "X2" = rnorm(N), "X3" = rnorm(N), "X4" = rnorm(N), "X5" = rnorm(N))
# 
# data <- cbind(data, add_ints(data[1:5]))
# data <- cbind(data, add_polynoms(data[1:5]))
# data <- cbind(data, add_indicators(data[1:5]))
# data <- cbind(data, add_trig(data[1:5]))
# 
# ################## can mix and match any of these:
# data$p.score1 <- minmax(sqrt_dgp(data))
# data$p.score2 <- expit(log_dgp(data))
# data$p.score3 <- minmax(exp_dgp(data))
# 
# #### Parse DGP
# data <- data.frame("X1" = rnorm(N), "X2" = rnorm(N), "X3" = rnorm(N), "X4" = rnorm(N), "X5" = rnorm(N))
# 
# data1 <- parse_dgp("Y0 ~ X1*X2 + log(X3) + X4 + I(X4)^2 + sqrt(X5) + sin(X1) + exp(X2)", data, "linear")
# data2 <- parse_dgp("Y0 ~ X1*X2 + log(X3) + X4 + I(X4)^2 + sqrt(X5) + sin(X1) + exp(X2)", data, "exp")
# data3 <- parse_dgp("Y0 ~ X1*X2 + log(X3) + X4 + I(X4)^2 + sqrt(X5) + sin(X1) + exp(X2)", data, "sqrt")
# data4 <- parse_dgp("Y0 ~ X1*X2 + log(X3) + X4 + I(X4)^2 + sqrt(X5) + sin(X1) + exp(X2)", data, "log")

# remove_cor(data, c(1,2))
# remove_groups(data, c("X1", "X2", "c1"))



