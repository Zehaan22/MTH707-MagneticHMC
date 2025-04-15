## Loading the algorithm
source("Algorithm.R")

## Getting the chain
chain <- MHMCSampler(L = 200, eps = 0.002, G = 0.3,
                      N = 1e4)

plot(means, pch = 19, cex = 2,
     col = "red", xlab = "X1",
     ylab = "X2", main = "Kou's Benchmarking Target",
     cex.main = 2)
points(chain, pch = 19, cex = 0.1,
       col = "blue")
points(means, pch = 19, cex = 2,
       col = "red")


## Getting the chain
chain <- MHMCSampler(L = 300, eps = 0.0025, G = 0.3,
                     N = 5e4)

plot(means, pch = 19, cex = 2,
     col = "red", xlab = "X1",
     ylab = "X2", main = "Kou's Benchmarking Target",
     cex.main = 2)
points(chain, pch = 19, cex = 0.1,
       col = "blue")
points(means, pch = 19, cex = 2,
       col = "red")


## Getting the chain
chain <- MHMCSampler(L = 200, eps = 0.006, G = 0.3,
                     N = 1e4)

plot(means, pch = 19, cex = 2,
     col = "red", xlab = "X1",
     ylab = "X2", main = "Magnetic HMC",
     cex.main = 2)
points(chain, pch = 19, cex = 0.1,
       col = "blue")
points(means, pch = 19, cex = 2,
       col = "red")