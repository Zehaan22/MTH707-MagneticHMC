#### defining the concentric L-1 Balls target

# Define the posterior distribution
target <- function(q) {
  r_vals <- c(4, 8, 16)
  sigma <- 0.5
  q_norm1 <- sum(abs(q))  # L1 norm
  
  density <- sum(exp(- (q_norm1 - r_vals)^2 / (2 * sigma^2)))
  return(density)
}

# Log of the posterior
log_target <- function(q) {
  r_vals <- c(4, 8, 16)
  sigma <- 0.5
  q_norm1 <- sum(abs(q))  # L1 norm
  
  log_density <- log(sum(exp(- (q_norm1 - r_vals)^2 / (2 * sigma^2))))
  return(log_density)
}

# Gradient of log posterior
grad_log <- function(q) {
  r_vals <- c(4, 8, 16)
  sigma <- 0.5
  q_norm1 <- sum(abs(q))  # L1 norm
  
  exp_terms <- exp(- (q_norm1 - r_vals)^2 / (2 * sigma^2))
  weight_terms <- (-(q_norm1 - r_vals) / sigma^2) * exp_terms
  
  # Weighted sum for gradient
  grad_weight <- sum(weight_terms) / sum(exp_terms)
  
  # Compute gradient for each component
  grad_q <- grad_weight * sign(q)  # Sign function for L1 norm gradient
  return(grad_q)
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
