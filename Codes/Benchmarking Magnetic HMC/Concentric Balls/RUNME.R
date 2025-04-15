## Loading the algorithm
source("Algorithm.R")

library(mcmcse)

chain <- MHMCSampler(N = 5e4,
                      L = 200,
                      eps = 0.006,
                      G = 0.3)
## Plotting the results

plot.ts(chain)
plot(chain)
acf(chain)

## Plotting the monte carlo samples
posterior <- function(q) {
  r_vals <- c(4, 8, 16)
  sigma <- 0.5
  q_norm1 <- sum(abs(q))  # L1 norm
  density <- sum(exp(- (q_norm1 - r_vals)^2 / (2 * sigma^2)))
  return(density)
}

# Rejection sampling function
generate_samples <- function(n_samples, dim = 2, proposal_sd = 10, max_iter = 1e6) {
  samples <- matrix(NA, nrow = n_samples, ncol = dim)
  count <- 0
  iter <- 0
  
  # Estimate a rough upper bound for the posterior (since it's unnormalized)
  max_posterior <- posterior(rep(8, dim))  # peak around L1 ~ 8
  
  while (count < n_samples && iter < max_iter) {
    iter <- iter + 1
    
    # Proposal from a uniform or normal distribution
    q_prop <- rnorm(dim, mean = 0, sd = proposal_sd)  # adjust sd as needed
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

# Visualize
plot(samples, pch = 19, col = rgb(0, 0, 1, 0.3), 
     main = "Concentric L1 Balls",
     xlab = "X", ylab = "Y",
     cex.main = 2)
points(chain[1e4:5e4,], pch = 19, col = rgb(0, 1, 0, 0.2))
legend("topright", legend = c("MHMC samples", "Monte Carlo samples"), 
       col = c(rgb(0, 1, 0), rgb(0, 0, 1)), pch = c(19, 19))

ess(chain)
