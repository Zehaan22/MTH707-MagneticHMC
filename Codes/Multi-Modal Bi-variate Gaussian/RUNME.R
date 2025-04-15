## Loading the algorithm
source("Algorithm.R")

## Answer to the universe
set.seed(42)

## Making the chain
chain <- magnetic_HMC(
  N = 1e4,
  k = 2,
  G = 0.5,
  L = 10,
  eps = 0.1
)

# Trace Plots
plot.ts(chain,
        main = "Trace Plot",
        cex.main = 2,
        col = "coral")

## Density Plots
mean1 = c(0,0) 
mean2 = c(5,5)
seq.x <- seq(-5,10, length = 1e5)

densities <- 0.5*dnorm(seq.x, 0) + 0.5*dnorm(seq.x, 5)
par(mfrow = c(2,1))
for(i in 1:2){
  plot(
    density(chain[,i]), 
    col = "black", 
    type = "l",
    main = paste("Density Plot for Component",i),
    xlab = "",
    cex.main = 2,
    ylab = "Density",
    ylim = c(0,0.5))
  lines(seq.x, 
        densities,
        col = "red")
  print(i)
}

## ACF Plots
par(mfrow = c(2,1))
acf(chain[,1], main = "ACF Plots", cex.main = 2, ylab = "Component-1")
acf(chain[,2], main = "", ylab = "Component-2")

## Visualizing the actual Samples
density1 <- matrix(rnorm(2*5e4 , mean = 0), ncol = 2)
density2 <- matrix(rnorm(2*5e4 , mean = 5), ncol = 2)

density_df <- as.data.frame(rbind(density1, density2))
colnames(density_df) <- c("x", "y")

## True Target
plot.gg <- ggplot() +
  stat_density_2d(data = density_df, aes(x = x, y = y, fill = after_stat(level)), geom = "polygon", alpha = 0.3) +
  scale_fill_viridis_c(option = "magma") +
  theme_void() +  # Removes background grid and axes
  theme(legend.position = "none",
        legend.title = element_blank())  # Removes legend title
plot.gg

## Samples 
plot.gg + geom_point(data = as.data.frame(chain), aes(x = V1, y = V2), color = "coral", alpha = 0.05, size = 0.1) 

## Visualing the trajectory
## True Target
plot.gg <- ggplot() +
  stat_density_2d(data = density_df, aes(x = x, y = y, fill = after_stat(level)), geom = "polygon", alpha = 0.3) +
  scale_fill_viridis_c(option = "magma") +
  theme_void() +  # Removes background grid and axes
  theme(legend.position = "none",
        legend.title = element_blank())  # Removes legend title
plot.gg


## Generating Trajectory points
n = 1e4
k = 2 
L = 30
eps = 0.1
G = 0.5

# Current State with a new momentum sample
p <- rnorm(2)
qt <- matrix(0, nrow = L+1, ncol = k)
qt[1,] <- c(1, 1)
count <- 2

## Adding the point to the plot
plot.gg <- plot.gg + geom_point(aes(x = qt[1,1], y = qt[1,2]), size = 3, shape = 21, fill = "lightgreen")
plot.gg ## Starting point

## Conformal LF
## Progressing the position and momentum L*eps time ahead
q <- qt[1,]
for(j in 1:L){
  # Leap Frog Steps
  ## Step-1
  p <- p - (eps/2)*grad_log_target(q)
  
  ## Step-2
  q <- q + (1/G)*(exp(G*eps/2) - 1)*p
  p <- exp(G*eps/2)*p
  
  ## Step-1
  p <- p - (eps/2)*grad_log_target(q)
  qt[count,] <- q
  count <- count + 1
}

# Making the final plot
plot.gg <- plot.gg +
  geom_point(aes(x = qt[2:L+1,1], y = qt[2:L+1,2]), size = 3, shape = 21, fill = "white")+
  geom_path(aes(x = qt[1:L+1,1], y = qt[1:L+1,2]), linewidth = 2,  color = "coral")

plot.gg
