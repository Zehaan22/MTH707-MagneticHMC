## Loading the algorithm
source("Algorithm.R")

## libraries
library(mcmcse)

## Answer to the universe
set.seed(42)

## Running the sampler
chain <- BnKSampler(N = 1e6)

## Looking at the actual density
x <- seq(-5,5,0.01)

par(mfrow = c(1,1))
plot(density(chain),
     main = "Density Plot for KnighNBishop Algorithm",
     cex.main = 2,
     lwd = 2)
lines(x,dnorm(x),type = "l",col = "red",
      lwd = 2,
      lt = 2)
legend("topright",c("Algorithm","True Density"),
       col = c("black","red"),
       lwd = 2,
       lty = 1:2)

## Assessing the quality of samples
par(mfrow = c(2,1))
acf(chain,main = "ACF Plot for KnighNBishop Algorithm",
    cex.main = 3)

## Looking at the trace plots
plot.ts(chain,
        main = "Trace Plot for KnighNBishop Algorithm",
        cex.main = 2)

## ESS values
ess(chain)

## "Acceptance Rate:  0.650025"