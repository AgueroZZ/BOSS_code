library(BayesGP)
set.seed(123)
lower = 0
upper = 10
x_grid <- seq(lower, upper, length.out = 1000)

func_generator <- function(f1, f2) {
  func <- function(x, a) {
    sapply(x, function(xi) {
      if (xi <= a) {
        f1(xi)
      } else {
        # print( - f2(a) + f1(a))
        f2(xi) - f2(a) + f1(a)
      }
    })
  }
  func
}
my_func <- func_generator(f1 = function(x) x*log(x^2 + 1), f2 = function(x) 3.3*x)
# Simulate observations from a regression model with a piecewise linear function
simulate_observation <- function(a, func, x_grid, measurement = 3) {
  # Generate n random x values between 0 and 1
  x <- rep(x_grid, each = measurement)
  n <- length(x)
  
  # Initialize y
  y <- numeric(n)
  
  # Loop through each x to compute the corresponding y value based on the piecewise function
  fx <- func(x = x, a = a)
  
  # Add random noise e from a standard normal distribution
  e <- rnorm(n, mean = 0, sd = 0.3)
  y <- fx + e
  
  return(data.frame(x, y))
}

# Knot value
a <- 6.5
# Simulate the data
plot(my_func(x = seq(0.1,10,length.out = 100), a = a)~seq(0.1,10,length.out = 100), type = "l")
data <- simulate_observation(a = a, func = my_func, x_grid = x_grid, measurement = 1)
plot(y~x, data)
save(data, file = "data.rda")







