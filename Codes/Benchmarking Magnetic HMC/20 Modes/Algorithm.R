#### defining the 20 modal target

## Means taken from Kou et. al. 2006
means <- matrix(c(
  2.18, 5.76, 
  8.67, 9.59, 
  4.24, 8.48, 
  8.41, 1.68, 
  3.93, 8.82, 
  3.25, 3.47, 
  1.70, 0.50, 
  4.59, 5.60, 
  6.91, 5.81, 
  6.87, 5.40, 
  5.41, 2.65, 
  2.70, 7.88, 
  4.98, 3.70, 
  1.14, 2.39, 
  8.33, 9.50, 
  4.93, 1.50, 
  1.83, 0.09, 
  2.26, 0.31, 
  5.54, 6.86, 
  1.69, 8.11
), ncol = 2, byrow = FALSE)

target <- function(x, var = 0.05){
  ## Defining the likelihood for the 20 modal target
  like <- 0
  for(i in 1:20){
    den <- exp(-0.5 * ((x[1] - means[i,1])^2 + (x[2] - means[i,2])^2) / var)
    like <- like + den
  }
  like <- (1/20) * like
  like <- (1/sqrt(2*pi*var)) * like
  return(like)
}

## Visulaising the location of the modes
par(mfrow = c(1,1))
plot(means, col = "red", pch = 19, cex = 2,
     main = "Visualising the 20 modes",
     xlab = "X", ylab = "Y",
     cex.main = 2)

## Definfing log target and grad log

log_target <- function(x, var = 0.05){
  ## Compute log of the unnormalized target
  return(log(target(x, var)))
}

grad_log <- function(x, var = 0.05){
  ## Defining the gradient of the log target
  
  # Compute squared distances and exponent terms
  diffs <- means - matrix(x, nrow = nrow(means), ncol = ncol(means), byrow = TRUE)
  squared_norms <- rowSums(diffs^2)
  exp_terms <- exp(-squared_norms / (2 * var))
  
  # Compute gradient
  weights <- exp_terms / sum(exp_terms)
  gradient <- colSums(weights * diffs) / var
  
  return(gradient)
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
  samples[1,] <- c(4.98, 3.70) ## Starting point
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
