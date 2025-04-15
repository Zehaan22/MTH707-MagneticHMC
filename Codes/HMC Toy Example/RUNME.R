################### Loading the Algorithm
source("Algorithm.R")

library(ggplot2)

################### Solution to the universe
set.seed(42)

################### Samples
chain <- Standard_HMC(N = 1e5,
                       k = 2,
                       L = 10,
                       eps = 0.1)

################### Visualising Samples

## Trace plot
par(mfrow = c(2,1))
plot(chain[,1],type = "l",col = "blue",
     xlab = "Iteration",ylab = "",
     main = "Trace Plot of Component 1")
plot(chain[,2],type = "l",col = "blue",
     xlab = "Iteration",ylab = "",
     main = "Trace Plot of Component 2")
## Density plot

## Determining the theoretical density
x <- seq(-3,3,length = 1e4)
y <- dnorm(x)

par(mfrow = c(2,1))
plot(density(chain[,1]),col = "blue", alpha = 0.5,
     main = "Density Plot of Component 1", lwd = 2,
     xlim = c(-3,3))
lines(x,y,col = "red", lt = 2, lwd = 1.5)
legend("topright",
       legend = c("Theoretical","Sampled"),
       col = c("red", "blue"),lwd = c(1,2), lty = c(2,1))

plot(density(chain[,2]),col = "blue", alpha = 0.5,
     main = "Density Plot of Component 2", lwd = 2,
     xlim = c(-3,3))
lines(x,y,col = "red", lt = 2, lwd = 1.5)
legend("topright",
       legend = c("Theoretical","Sampled"),
       col = c("red", "blue"),lwd = c(1,2), lty = c(2,1))

### Visualising the Trajectory in 2-d plane
par(mfrow = c(1,1))

density1 <- matrix(rnorm(2*5e4 , mean = 0), ncol = 2)

density_df <- as.data.frame(density1)
colnames(density_df) <- c("x", "y")

## True Target x Momentum
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
  geom_point(aes(x = pts_1, y = pts_2), color = "blue", size = 0.6, alpha = 0.5)

## Starting Point
x_0 <- 0
p_0 <- 1

plot.gg <- plot.gg +
  geom_point(aes(x = x_0, y = p_0), color = "blue", size = 3)
plot.gg

## Leap Frog Steps
eps = 0.1
L = 40
p = p_0
q = x_0

moms <- numeric(L)
pos <- numeric(L)

count <- 1


for(j in 1:L){
  # Leap Frog Steps
  ## Step-1
  p <- p + eps/2*grad_log_target(q)
  
  ## Step-2
  q <- q + eps*p
  
  ## Step-3
  p <- p + eps/2*grad_log_target(q)
  
  # Storing the values
  moms[count] <- p
  pos[count] <- q
  
  count <- count + 1
}

plot.gg +
  geom_point(aes(x = pos, y = moms), color = "red", size = 3) 
