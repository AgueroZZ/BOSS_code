library(BayesGP)
library(npreg)
library(ggplot2)
library(INLA)
source_location <- "function/"
source(file = paste0(source_location, "bo_functions.R"))
load("results/data.rda")
load(file = "results/exact_grid_result.rda")
lower = 0.1; upper = 9.9
set.seed(123)
obtain_aghq <- function(f, k = 100, startingvalue = 0, optresult = NULL){
  if(!is.null(optresult)){
    return(aghq::aghq(ff = ff, k = k, startingvalue = startingvalue, optresults = optresult))
  }
  else{
    ff <- list(fn = f, gr = function(x) numDeriv::grad(f, x), he = function(x) numDeriv::hessian(f, x))
    return(aghq::aghq(ff = ff, k = k, startingvalue = startingvalue))
  }
}

### 1. AGHQ on exact grid
lf_data <- data.frame(x = exact_grid_result$x/10,
                     lfx = exact_grid_result$exact_vals)

## Convert to the real line:
lg_data <- data.frame(y = qnorm(lf_data$x),
                      lgy = exact_grid_result$exact_vals + dnorm(qnorm(lf_data$x), log = T))
surrogate <- function(xvalue, data_to_smooth){
  predict(ss(x = data_to_smooth$y, y = data_to_smooth$lgy, df = length(unique(data_to_smooth$y)), m = 2, all.knots = TRUE), x = xvalue)$y
}
fn <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lg_data))

grid_opt_list = list(ff = list(fn = fn, gr = function(x) numDeriv::grad(fn, x), he = function(x) numDeriv::hessian(fn, x)), 
                     mode = lg_data$y[which.max(lg_data$lgy)])

grid_opt_list$hessian = -grid_opt_list$ff$he(grid_opt_list$mode)

# y_vec = seq(-2,2, by = 0.01)
# plot(fn(y_vec) ~ y_vec, type = "l")

## Compute the runtime:
start_time <- Sys.time()
aghq_result_grid <- obtain_aghq(f = fn, k = 10, optresult = grid_opt_list)
end_time <- Sys.time()
end_time - start_time

### Compute the quadrature from BOSS
load(file = "results/BO_data_to_smooth.rda")
data_to_smooth <- BO_result_original_list[[4]]
lf_data_BO <- data.frame(x = data_to_smooth$x_original/10,
                      lfx = data_to_smooth$y)
## Convert to the real line:
lg_data_BO <- data.frame(y = qnorm(lf_data_BO$x),
                      lgy = data_to_smooth$y + dnorm(qnorm(lf_data_BO$x), log = T))

fn_BO <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lg_data_BO))
grid_opt_list_BO = list(ff = list(fn = fn_BO, gr = function(x) numDeriv::grad(fn_BO, x), he = function(x) numDeriv::hessian(fn_BO, x)), 
                     mode = lg_data_BO$y[which.max(lg_data_BO$lgy)])

grid_opt_list_BO$hessian = -grid_opt_list_BO$ff$he(grid_opt_list_BO$mode)
## Compute the runtime:
start_time <- Sys.time()
aghq_result_BOSS <- obtain_aghq(f = fn_BO, k = 10, optresult = grid_opt_list_BO)
end_time <- Sys.time()
end_time - start_time


### Compare the quadrature
quad_exact <- aghq_result_grid$normalized_posterior$nodesandweights
quad_BO <- aghq_result_BOSS$normalized_posterior$nodesandweights

plot(weights ~ theta1, type = "o", col = "red", data = quad_exact, ylim = c(0, 1), xlim = c(0,1))
points(weights ~ theta1, type = "o", col = "blue", data = aghq_result_BOSS$normalized_posterior$nodesandweights)

save(quad_exact, file = "quad_exact.rda")
save(quad_BO, file = "quad_BO.rda")





