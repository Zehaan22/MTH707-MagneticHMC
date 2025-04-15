#### KNBS for bi-modal gaussian

################# Target specification ##################
## Bi-modal Normal

target <- function(
    x, 
    mean1 = 0, 
    mean2 = 5){
  
  ## Density of Bi-modal Normal
  return((1/sqrt(2*pi))*(exp(-0.5*(x-mean1)^2) + exp(-0.5*(x-mean2)^2)))
}

log_target <- function(x){
  return(log(target(x)))
}

grad_log_target <- function(
    x,
    mean1 = 0,
    mean2 = 5){
  return((-1/target(x))*((x-mean1)*exp(-0.5*(x-mean1)^2) + (x-mean2)*exp(-0.5*(x-mean2)^2)))
}

################# Leapfrog ###############################
e0 <- 1
g0 <- 0.5

conformalLF <- function(p,x,
                        L,eps,
                        G){
  ## Progressing the position and momentum L*eps time ahead
  for(j in 1:L){
    # Leap Frog Steps
    ## Step-1
    p <- p - (eps/2)*grad_log_target(x)
    
    ## Step-2
    x <- x + (1/G)*(exp(G*eps/2) - 1)*p
    p <- exp(G*eps/2)*p
    
    ## Step-1
    p <- p - (eps/2)*grad_log_target(x)
  }
  ## Time incremented State
  return(list(p,x))
}

################# Good Starting Values for parameters ####
find_start_e <- function(x, max_iter = 10){
  e <- 1
  p <- rnorm(1)
  new_st <- conformalLF(p,x,1,e,g0)
  st_ratio <- (dnorm(new_st[[1]])*target(new_st[[2]]))/(dnorm(p)*target(x))
  if (st_ratio > 0.5){
    a <- 1
  }else{
    a <- -1
  }
  
  count <- 1
  while(st_ratio^a > 2^(-a) & count < max_iter){
    e <- (2^a)*e
    new_st <- conformalLF(new_st[[1]],new_st[[2]],1,e,g0)
    st_ratio <- (dnorm(new_st[[1]])*target(new_st[[2]]))/(dnorm(p)*target(x))
    count <- count + 1
  }
  return(e)
}

find_start_G <- function(x, max_iter = 10){
  G <- g0
  p <- rnorm(1)
  new_st <- conformalLF(p,x,1,e0,G)
  st_ratio <- (dnorm(new_st[[1]])*target(new_st[[2]]))/(dnorm(p)*target(x))
  if (st_ratio > 0.5){
    a <- 1
  }else{
    a <- -1
  }
  
  count <- 1
  while(st_ratio^a > 2^(-a) & count < max_iter){
    G <- (2^a)*G
    new_st <- conformalLF(new_st[[1]],new_st[[2]],1,e0,G)
    st_ratio <- (dnorm(new_st[[1]])*target(new_st[[2]]))/(dnorm(p)*target(x))
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
        if((x_min - x_pls)*p_min <= del_max){
          return(c((i + j),v))
        }
      }
    }else{ ## Positive motion
      for(i in 1:j){
        new_st <- conformalLF(p_pls,x_pls,1,e,g)
        p_pls <- new_st[[1]]
        x_pls <- new_st[[2]]
        ## Check the condition
        if((x_pls - x_min)*p_pls <= del_max){
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

BnKSampler <- function(N = 1e5){ ## only argument is sample size
  
  ## Initialising Values
  samples <- numeric(N)
  samples[1] <- 0 ## Starting point
  momentum <- rnorm(N)
  accept <- 0
  
  ## Initialising the parameters
  e <- find_start_e(0)
  g <- find_start_G(0)
  
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
  Madapt <- N/100
  
  L_tot <- 0
  a <- 0
  
  ## The main loop
  for(i in 2:N){
    p <- momentum[i]
    x <- samples[i-1]
    ret <- BuildTree(x,p,e,g)
    
    ## Adaptation of e and G
    if (i < Madapt){
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
      
      
    }else{
      e <- e0
      L <- L_tot/Madapt
      
      if(g < 0){
        g <- -g0
      }else{
        g <- g0
      }
    }
    
    ## The leaping frog
    ret <- conformalLF(p,x,L,v*e,g)
    p_prop <- ret[[1]]
    x_prop <- ret[[2]]
    
    ## Energy Calc
    eng.curr <- -log_target(x) + (p^2)/2
    eng.prop <- -log_target(x_prop) + (p_prop^2)/2
    
    ## acceptance ratio 
    a <- (eng.curr - eng.prop)
    
    ## Acceptance and Flips
    if(log(runif(1)) <= a){
      samples[i] <- x_prop
      accept <- accept + 1
      g <- -g
    }
    else{
      samples[i] <- samples[i-1]
    }
  }
  
  ## Printing acceptance rate
  print(paste("Acceptance Rate: ",accept/N))
  print(paste("Final Epsilon: ", e0))
  print(paste("Final Gamma: ", g0))
  print(paste("Average L: ", L_tot/N))
  
  return(samples[(Madapt+1):N])
}