library(BayesGP)
library(npreg)
library(INLA)

set.seed(123)
source_location <- "function/"
source(file = paste0(source_location, "bo_functions.R"))

load("data.rda")
plot(y~x, data)

lower = 0.1
upper = 9.9
a <- 6.5
eval_once <- function(alpha){
  a_fit <- alpha
  data$x1 <- ifelse(data$x <= a_fit, data$x, a_fit)
  data$x2 <- ifelse(data$x > a_fit, (data$x - a_fit), 0)
  mod <- model_fit(formula = y ~ f(x1, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0) + f(x2, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0), 
                   data = data, method = "aghq", family = "Gaussian", aghq_k = 3
  )
  (mod$mod$normalized_posterior$lognormconst)
}
x_vals <- seq(lower, upper, length.out = 1000)

# Initialize the progress bar
begin_time <- Sys.time()
total <- length(x_vals)
pb <- txtProgressBar(min = 0, max = total, style = 3)
# Initialize exact_vals if needed
exact_vals <- c()
# Loop with progress bar update
for (i in 1:total) {
  xi <- x_vals[i]
  # Your existing code
  exact_vals <- c(exact_vals, eval_once(xi))
  # Update the progress bar
  setTxtProgressBar(pb, i)
}
# Close the progress bar
close(pb)
exact_grid_result <- data.frame(x = x_vals, exact_vals = exact_vals)
exact_grid_result$exact_vals <- exact_grid_result$exact_vals - max(exact_grid_result$exact_vals)
exact_grid_result$fx <- exp(exact_grid_result$exact_vals)
end_time <- Sys.time()
end_time - begin_time
# Time difference of 1.27298 hours
save(exact_grid_result, file = "exact_grid_result.rda")

# Calculate the differences between adjacent x values
dx <- diff(exact_grid_result$x)
# Compute the trapezoidal areas and sum them up
integral_approx <- sum(0.5 * (exact_grid_result$fx[-1] + exact_grid_result$fx[-length(exact_grid_result$fx)]) * dx)
exact_grid_result$pos <- exact_grid_result$fx / integral_approx
plot(exact_grid_result$x, exact_grid_result$pos, type = "l", col = "red", xlab = "x (0-10)", ylab = "density", main = "Posterior")
abline(v = a, col = "purple")
grid()
save(exact_grid_result, file = "exact_grid_result.rda")


# Plot smoothed version:
# surrogate <- function(xvalue, data_to_smooth){
#   predict(ss(x = data_to_smooth$x, y = data_to_smooth$exact_vals, m = 3, lambda = (1e-10)/nrow(data_to_smooth)), x = xvalue)$y
# }
surrogate <- function(xvalue, data_to_smooth){
  predict(ss(x = data_to_smooth$x, y = data_to_smooth$exact_vals, df = length(unique(data_to_smooth$x)), m = 2, all.knots = TRUE), x = xvalue)$y
}

# Convert to the internal scale:
exact_grid_result_internal <- data.frame(x = (exact_grid_result$x - lower)/(upper-lower),
                                         exact_vals = exact_grid_result$exact_vals + log(upper - lower))

exact_grid_result_internal_smooth <- data.frame(x = exact_grid_result_internal$x)
exact_grid_result_internal_smooth$exact_vals <- surrogate(xvalue = exact_grid_result_internal$x, data_to_smooth = exact_grid_result_internal)

# Convert back:
exact_grid_result_smooth <- data.frame(x = (exact_grid_result_internal_smooth$x)*(upper - lower) + lower, exact_vals = exact_grid_result_internal_smooth$exact_vals - log(upper - lower))
exact_grid_result_smooth$exact_vals <- exact_grid_result_smooth$exact_vals - max(exact_grid_result_smooth$exact_vals)
exact_grid_result_smooth$fx <- exp(exact_grid_result_smooth$exact_vals)
dx <- diff(exact_grid_result_smooth$x)
integral_approx <- sum(0.5 * (exact_grid_result_smooth$fx[-1] + exact_grid_result_smooth$fx[-length(exact_grid_result_smooth$fx)]) * dx)
exact_grid_result_smooth$pos <- exact_grid_result_smooth$fx / integral_approx
save(exact_grid_result_smooth, file = "exact_grid_result_smooth.rda")

plot(x = exact_grid_result$x, y = exact_grid_result$exact_vals, col = "black", cex = 0.5, type = "p", ylim = c(-50,0))
lines(exact_grid_result_smooth$x, (exact_grid_result_smooth$exact_vals), type = "l", col = "red")
abline(v = exact_grid_result_smooth$x[which.max(exact_grid_result_smooth$exact_vals)], col = "green", lty = "dashed")
abline(v = exact_grid_result$x[which.max(exact_grid_result$exact_vals)], col = "blue", lty = "dashed")
abline(v = a, col = "purple", lty = "dashed")
grid()

plot(exact_grid_result_smooth$x, exp(exact_grid_result_smooth$exact_vals), type = "l", col = "red", xlab = "x", ylab = "y", main = "Plot of ff$fn")
abline(v = a, col = "purple")
grid()

