####### Loading the algorithm
source("Algorithm.R")


###### Getting the chain
chain <- KnBSampler(N = 1e4)

plot.ts(chain)
acf(chain)

plot(chain,
     xlab = "X1", ylab = "X2", main = "Comparing Samples"
     )

## Plotting the monte-carlo samples
points(rnorm(5e3),rnorm(5e3),
       col = "red")

### Comparing against theoretical density
library(ggplot2)


### Visualising the Trajectory in 2-d plane
par(mfrow = c(1,1))

density1 <- matrix(rnorm(2*5e4 , mean = 0), ncol = 2)

density_df <- as.data.frame(density1)
colnames(density_df) <- c("x", "y")

## True Target x Momentum
plot.gg <- ggplot() +
  stat_density_2d(data = density_df, aes(x = x, y = y, fill = after_stat(level)), geom = "polygon", alpha = 0.3) +
  geom_point(aes(x = chain[,1], y = chain[,2]), color = "coral2", size = 0.6, alpha = 0.6)+
  scale_fill_viridis_c(option = "magma") +
  theme_void() +  # Removes background grid and axes
  theme(legend.position = "none",  # Removes legend
        legend.title = element_blank())  # Removes legend title
plot.gg

## Plotting the densities

x_seq <- seq(-3, 3, length.out = 1e3)

par(mfrow = c(2,1))
plot(density(chain[,1]),
     xlab = "X", ylab = "Density",
     main = "Density Plot for Component 1",
     cex.main = 2, col = "coral2")
lines(x_seq, dnorm(x_seq), lwd = 2, lt = 2,
     xlab = "X", ylab = "Density",
     main = "Comparing the densities",
     cex.main = 2)
legend("topright", legend = c("Sampled Density", "True Density"),
       col = c("coral2", "black"), lty = c(1, 2), lwd = 2)

plot(density(chain[,2]),
     xlab = "X", ylab = "Density",
     main = "Density Plot for Component 2",
     cex.main = 2, col = "coral2")
lines(x_seq, dnorm(x_seq), lwd = 2, lt = 2,
      xlab = "X", ylab = "Density",
      main = "Comparing the densities",
      cex.main = 2)
legend("topright", legend = c("Sampled Density", "True Density"),
       col = c("coral2", "black"), lty = c(1, 2), lwd = 2)


