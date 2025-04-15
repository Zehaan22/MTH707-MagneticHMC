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


## Defining the sampler
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



################# The Bishop and Knight Algorithm #################

MHMCSampler <- function(L, eps, G,
                        N = 1e4){ ## only argument is sample size
  
  ## Initialising Values
  samples <- matrix(0, nrow = N, ncol = 2)
  samples[1,] <- c(0, 0) ## Starting point
  momentum <- matrix(rnorm(2*N), nrow = N, ncol = 2)
  accept <- 0
  
  ## Setting the parameters
  e <- eps
  g <- G
  
  ## The main loop
  for(i in 2:N){
    p <- momentum[i,]
    x <- samples[i-1,]
    
    ## The leaping frog
    ret <- conformalLF(p,x,L,e,g)
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
