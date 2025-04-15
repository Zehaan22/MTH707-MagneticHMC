## Loading the algorithm
source("Algorithm.R")
library(mcmcse)
## Running the algorithm
chain <- MHMCSampler(L = 100 , eps = 0.005, G = 0.5,
                      N = 5e4
                      )
## Plotting the results
plot.ts(chain)
plot(chain)
acf(chain)


## Plotting the monte carlo samples
posterior <- function(q) {
  mu_vals <- list(c(0, 0), c(8, 0), c(-8, 0), c(0, 8), c(0, -8))
  r_vals <- c(20, 4, 4, 4, 4)
  sigma <- 0.5
  
  sum_exp <- 0  # Initialize sum
  
  for (i in 1:5) {
    norm_l1 <- sum(abs(q - mu_vals[[i]]))  # Compute L1 norm
    exponent <- -((norm_l1 - r_vals[i])^2) / (2 * sigma^2)
    sum_exp <- sum_exp + exp(exponent)
  }
  
  return(sum_exp)  # Unnormalized posterior
}

generate_samples <- function(n_samples, dim = 2, proposal_sd = 10, max_iter = 1e6) {
  samples <- matrix(NA, nrow = n_samples, ncol = dim)
  count <- 0
  iter <- 0
  
  # Estimate a rough upper bound (you can make this tighter)
  max_posterior <- posterior(c(0, 0)) + 1  # conservative upper bound
  
  while (count < n_samples && iter < max_iter) {
    iter <- iter + 1
    
    # Proposal from a broad Gaussian
    q_prop <- rnorm(dim, mean = 0, sd = proposal_sd)
    u <- runif(1, 0, max_posterior)
    
    if (u < posterior(q_prop)) {
      count <- count + 1
      samples[count, ] <- q_prop
    }
  }
  
  if (iter == max_iter) {
    warning("Max iterations reached before collecting all samples.")
  }
  
  return(samples[1:count, ])
}

samples <- generate_samples(5000)

plot(samples, pch = 19, col = rgb(0, 0, 1, 0.3),
     main = "Nested L1 Balls",
     xlab = "X", ylab = "Y",
     cex.main = 2)
points(chain[1e4:5e4,], pch = 19, col = rgb(0, 1, 0, 0.2))
legend("topright", legend = c("MHMC samples", "Monte Carlo samples"), 
       col = c(rgb(0, 1, 0), rgb(0, 0, 1)), pch = c(19, 19))

ess(chain)
