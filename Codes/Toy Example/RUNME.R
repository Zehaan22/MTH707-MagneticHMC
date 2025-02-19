## Loading the algorithm
source("Algorithm.R")


## Making the chain
chain <- magnetic_HMC(N = 1e6)

# Trace Plots
plot.ts(chain)

## Density Plots
seq.x <- seq(-5,5, length = 1e4)
densities <- sapply(seq.x, dnorm)
par(mfrow = c(2,1))
for(i in 1:2){
  plot(
    density(chain[,i]), 
    col = "black", 
    type = "l",
    main = paste("Density Plot for Component",i),
    xlab = "",
    ylab = "Density",
    ylim = c(0,0.5))
  lines(seq.x, 
        densities,
        col = "red")
  print(i)
}

## ACF Plots
acf(chain)

