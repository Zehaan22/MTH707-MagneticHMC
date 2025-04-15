## Coding a toy example for the current proposal of Standard HMC

################### Target Specification N_2(0,I)
target <- function(x){
  # Standard Normal Target
  return((1/sqrt(2*pi))*exp(-sum(x^2)/2))
}

log_target <- function(x){
  # Standard Normal log Target
  return(-sum(x^2)/2)
}

grad_log_target <- function(x){
  # Standard Normal log Target
  return(-x)
}
###############################################


################### Magnetic HMC Algorithm
## Conformal Leap Frog
LeapFrog <- function(p,q,
                        L,eps
                        ){
  ## Progressing the position and momentum L*eps time ahead
  for(j in 1:L){
    # Leap Frog Steps
    ## Step-1
    p <- p + eps/2*grad_log_target(q)
    
    ## Step-2
    q <- q + eps*p
    
    ## Step-3
    p <- p + eps/2*grad_log_target(q)
  }
  ## Time incremented State
  return(list(p,q))
}

## Main Algorithm
Standard_HMC <- function(
    N = 1e4,
    k = 2,
    L = 10,
    eps = 0.1
){
  
  ## Initialize a starting point 
  q <- matrix(nrow = N, ncol = k)
  q[1,] <- rep(0,k)
  
  ## Generate momentums 
  p <- matrix(rnorm(N*k),
              nrow = N, ncol = k)
  
  ## Generate the chain
  accept <- 0
  for(i in 2:N){
    
    ## "re-sampling" momentum
    p_curr <- p[i,]
    q_curr <- q[i-1,]
    
    ## Conformal Leapfrog
    res <- LeapFrog(p_curr,q_curr,
                       L,eps)
    p_prop <- res[[1]]
    q_prop <- res[[2]]
    
    ## Energy Calc
    eng.curr <- -log_target(q_curr) + (sum(p_curr^2)/2)
    eng.prop <- -log_target(q_prop) + (sum(p_prop^2)/2)
    
    ## acceptance ratio 
    a <- (eng.curr - eng.prop)
    
    if(log(runif(1)) <= a){
      q[i,] <- q_prop
      accept <- accept + 1
    }
    else{
      q[i,] <- q[i-1,]
    }
  }
  
  ## Final Result
  print(paste("Acceptance Ratio =", accept/N))
  
  return(q)
}
###############################################




