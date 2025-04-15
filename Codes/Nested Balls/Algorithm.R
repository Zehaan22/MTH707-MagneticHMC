#### defining the concentric L-1 Balls target

# Define the posterior distribution
target <- function(q) {
  mu_vals <- list(c(0, 0), c(8, 0), c(-8, 0), c(0, 8), c(0, -8))
  r_vals <- c(20, 4, 4, 4, 4)
  sigma <- 0.5
  
  sum_exp <- 0  # Initialize sum
  
  for (i in 1:5) {
    norm_l1 <- sum(abs(q - mu_vals[[i]]))  # Compute L1 norm
    exponent <- -((norm_l1 - r_vals[i])^2) / (2 * sigma^2)  # Compute exponent
    sum_exp <- sum_exp + exp(exponent)  # Accumulate exponentials
  }
  
  return(sum_exp)  # Return the unnormalized density
}

# Log of the posterior
log_target <- function(q) {
  mu_vals <- list(c(0, 0), c(8, 0), c(-8, 0), c(0, 8), c(0, -8))
  r_vals <- c(20, 4, 4, 4, 4)
  sigma <- 0.5
  
  sum_exp <- 0  # Initialize sum
  
  for (i in 1:5) {
    norm_l1 <- sum(abs(q - mu_vals[[i]]))  # Compute L1 norm
    exponent <- -((norm_l1 - r_vals[i])^2) / (2 * sigma^2)  # Compute exponent
    sum_exp <- sum_exp + exp(exponent)  # Accumulate exponentials
  }
  
  log_density <- log(sum_exp)  # Take log of sum
  return(log_density)
}

# Gradient of log posterior
grad_log <- function(q) {
  mu_vals <- list(c(0, 0), c(8, 0), c(-8, 0), c(0, 8), c(0, -8))
  r_vals <- c(20, 4, 4, 4, 4)
  sigma <- 0.5
  
  sum_exp <- 0  # Store sum of exponentials
  grad_sum <- c(0, 0)  # Initialize gradient vector
  
  for (i in 1:5) {
    norm_l1 <- sum(abs(q - mu_vals[[i]]))  # Compute L1 norm
    exponent <- -((norm_l1 - r_vals[i])^2) / (2 * sigma^2)  # Compute exponent
    term <- exp(exponent)  # Compute exponential term
    
    sum_exp <- sum_exp + term  # Accumulate exponentials
    
    # Compute gradient
    sign_q <- sign(q - mu_vals[[i]])  # Derivative of L1 norm
    grad_sum <- grad_sum + term * (2 * (norm_l1 - r_vals[i]) / (2 * sigma^2)) * sign_q
  }
  
  return(grad_sum / sum_exp)  # Normalize by posterior sum
}

########## Leapfrog 
conformalLF <- function(p,x,
                        L,eps,
                        G){
  ## Progressing the position and momentum L*eps time ahead
  for(j in 1:L){
    # Leap Frog Steps
    ## Step-1
    p <- p - (eps/2)*grad_log(x)
    
    ## Step-2
    x <- x + (1/G)*(exp(G*eps/2) - 1)*p
    p <- exp(G*eps/2)*p
    
    ## Step-1
    p <- p - (eps/2)*grad_log(x)
  }
  ## Time incremented State
  return(list(p,x))
}


################# Good Starting Values for parameters ####
find_start_e <- function(x, max_iter = 5){
  e <- 1
  p <- rnorm(2)
  g0 <- 0.1
  
  new_st <- conformalLF(p,x,1,e,g0)
  
  log_st_ratio <- (log_target(new_st[[1]]) + log_target(new_st[[2]])) -
    (log_target(p) + log_target(x))
  
  if (log_st_ratio > log(0.5)){
    a <- 1
  }else{
    a <- -1
  }
  
  count <- 1
  while((log_st_ratio) < log(0.5) & count < max_iter){
    e <- (2^a)*e
    new_st <- conformalLF(new_st[[1]],new_st[[2]],1,e,g0)
    log_st_ratio <- (log_target(new_st[[1]]) + log_target(new_st[[2]])) - 
      (log_target(p) + log_target(x))
    count <- count + 1
  }
  return(e)
}

find_start_G <- function(x, max_iter = 5){
  G <- 1
  p <- rnorm(2)
  e0 <- 0.1
  
  new_st <- conformalLF(p,x,1,e0,G)
  log_st_ratio <- (log_target(new_st[[1]]) + log_target(new_st[[2]])) - 
    (log_target(p) + log_target(x))
  if (log_st_ratio > log(0.5)){
    a <- 1
  }else{
    a <- -1
  }
  
  count <- 1
  while(log_st_ratio < log(0.5) & count < max_iter){
    G <- (2^a)*G
    new_st <- conformalLF(new_st[[1]],new_st[[2]],1,e0,G)
    st_ratio <- (log_target(new_st[[1]]) + log_target(new_st[[2]])) - 
      (log_target(p) + log_target(x))
    count <- count + 1
  }
  return(G)
}

################# The exploration routine #################

BuildTree <- function(x,p,e,g,
                      del_max = 1e-6, max_iter = 10){
  
  ## Initialise values
  x_pls <- x
  x_min <- x
  p_pls <- p
  p_min <- p
  
  s = 1
  j = 1
  
  while(s < max_iter){
    ## Sample a direction
    v <- 2*rbinom(1,1,0.5) - 1
    
    ## Negative motion
    if(v == -1){
      for(i in 1:j){
        new_st <- conformalLF(p_min,x_min,1,e,g)
        p_min <- new_st[[1]]
        x_min <- new_st[[2]]
        ## Check the condition
        if(sum((x_min - x_pls)*p_min) <= del_max){
          return(c((i + j),v))
        }
      }
    }else{ ## Positive motion
      for(i in 1:j){
        new_st <- conformalLF(p_pls,x_pls,1,e,g)
        p_pls <- new_st[[1]]
        x_pls <- new_st[[2]]
        ## Check the condition
        if(sum((x_pls - x_min)*p_pls) <= del_max){
          return(c((i + j),v))
        }
      }
    }
    j <- 2*j
    s <- s + 1
  }
  return(c(j,v))
}

################# The Bishop and Knight Algorithm #################

KnBSampler <- function(N = 1e5){ ## only argument is sample size
  
  ## Initialising Values
  samples <- matrix(0, nrow = N, ncol = 2)
  samples[1,] <- c(1,1) ## Starting point
  momentum <- matrix(rnorm(2*N), nrow = N, ncol = 2)
  accept <- 0
  
  ## Initialising the parameters
  e <- find_start_e(samples[1,])
  g <- find_start_G(samples[1,])
  
  ## Setting hyper-params
  mu1 <- log(10*e)
  mu2 <- log(10*g)
  e0 <- 1
  g0 <- 0.5
  H0 <- 0
  H <- H0
  gamma <- 0.05
  t0 <- 10
  kappa <- 0.75
  Madapt <- 1e3
  
  L_tot <- 0
  a <- 0
  
  ## Tuning Loop
  for(i in 2:Madapt){
    print(i)
    p <- momentum[i,]
    x <- samples[i,]
    ret <- BuildTree(x,p,e,g)
    
    ## Adaptation of e and G
    H <- (1 - 1/(i + t0))*H + (1/(i + t0))*(0.65 - min(1,exp(a)))
    
    ## Epsilon
    e <- exp(mu1 - sqrt(i)/gamma*H)
    e0 <- exp(i^(-kappa)*log(e) + (1 - i^(-kappa))*log(e0))
    
    ## Gamma
    sign <- 1
    if (g < 0){
      g <- -g 
      sign <- -1
    }
    g <- exp(mu2 - sqrt(i)/gamma*H)
    g0 <- exp(i^(-kappa)*log(e) + (1 - i^(-kappa))*log(g0))
    g <- sign*g
    
    ## L
    L <- ret[1]
    v <- ret[2]
    L_tot <- L_tot + L
    
    ## The leaping frog
    ret <- conformalLF(p,x,L,v*e,g)
    p_prop <- ret[[1]]
    x_prop <- ret[[2]]
    
    ## Energy Calc
    eng.curr <- -log_target(x) + sum(p^2)/2
    eng.prop <- -log_target(x_prop) + sum(p_prop^2)/2
    
    ## acceptance ratio 
    a <- (eng.curr - eng.prop)
    
    ## Acceptance and Flips
    if(log(runif(1)) <= a){
      samples[i,] <- x_prop
      accept <- accept + 1
      g <- -g
    }
    else{
      samples[i,] <- samples[i-1,]
    }
    
  }
  
  ## Setting the parameters
  e <- e0
  L <- L_tot/Madapt
  
  if(g < 0){
    g <- -g0
  }else{
    g <- g0
  }
  
  print(paste("Final Epsilon: ", e0))
  print(paste("Final Gamma: ", g0))
  print(paste("Average L: ", L_tot/N))
  
  ## The main loop
  for(i in 2:N){
    p <- momentum[i,]
    x <- samples[i-1,]
    ret <- BuildTree(x,p,e,g)
    
    ## The leaping frog
    ret <- conformalLF(p,x,L,v*e,g)
    p_prop <- ret[[1]]
    x_prop <- ret[[2]]
    
    ## Energy Calc
    eng.curr <- -log_target(x) + sum(p^2)/2
    eng.prop <- -log_target(x_prop) + sum(p_prop^2)/2
    
    ## acceptance ratio 
    a <- (eng.curr - eng.prop)
    
    ## Acceptance and Flips
    if(log(runif(1)) <= a){
      samples[i,] <- x_prop
      accept <- accept + 1
      g <- -g
    }
    else{
      samples[i,] <- samples[i-1,]
    }
  }
  
  ## Printing acceptance rate
  print(paste("Acceptance Rate: ",accept/N))
  
  return(samples)
}
