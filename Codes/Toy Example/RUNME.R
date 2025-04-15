## Loading the algorithm
source("Algorithm.R")


## Making the chain
chain <- magnetic_HMC(N = 1e6)

# Trace Plots
plot.ts(chain)

## Density Plots
seq.x <- seq(-3,3, length = 1e4)
densities <- sapply(seq.x, dnorm)
par(mfrow = c(2,1))
for(i in 1:2){
  plot(
    density(chain[,i]), 
    col = "gold2", 
    type = "l",
    main = paste("Density Plot for Component",i),
    xlab = "",
    ylab = "Density",
    ylim = c(0,0.5),
    xlim = c(-3,3),
    lwd = 2)
  lines(seq.x, 
        densities,
        col = "red",
        lwd = 1.5, 
        lt = 2)
  legend("topright",
         legend = c("Theoretical","Sampled"),
         col = c("red", "gold2"),lwd = c(1.5,2), lty = c(2,1))
}

## ACF Plots
acf(chain)

## True Target x Momentum
par(mfrow = c(1,1))

density1 <- matrix(rnorm(2*5e4 , mean = 0), ncol = 2)

density_df <- as.data.frame(density1)
colnames(density_df) <- c("x", "y")

plot.gg <- ggplot() +
  stat_density_2d(data = density_df, aes(x = x, y = y, fill = after_stat(level)), geom = "polygon", alpha = 0.3) +
  scale_fill_viridis_c(option = "magma") +
  theme_void() +  # Removes background grid and axes
  theme(legend.position = "none",  # Removes legend
        legend.title = element_blank())  # Removes legend title
plot.gg

## Visualising the samples
pts_1 <- sample(chain[,1], 5e3)
pts_2 <- sample(chain[,2], 5e3)
plot.gg + 
  geom_point(aes(x = pts_1, y = pts_2), color = "gold2", size = 0.6, alpha = 0.5)