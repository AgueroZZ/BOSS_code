library(BayesGP)
library(npreg)
library(ggplot2)
library(INLA)

set.seed(123)
noise_var = 1e-8
source_location <- "function/"
source(file = paste0(source_location, "bo_functions.R"))

load("data.rda")
plot(y~x, data)
eval_once <- function(alpha){
  a_fit <- alpha
  data$x1 <- ifelse(data$x <= a_fit, data$x, a_fit)
  data$x2 <- ifelse(data$x > a_fit, (data$x - a_fit), 0)
  mod <- model_fit(formula = y ~ f(x1, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0) + f(x2, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0), 
                   data = data, method = "aghq", family = "Gaussian", aghq_k = 3
  )
  (mod$mod$normalized_posterior$lognormconst)
}
surrogate <- function(xvalue, data_to_smooth){
  predict(ss(x = data_to_smooth$x, y = data_to_smooth$y, df = length(unique(data_to_smooth$x)), m = 2, all.knots = TRUE), x = xvalue)$y
}
lower = 0.1
upper = 9.9
objective_func <- eval_once


#### Runtime analysis:
eval_num <- seq(from = 15, to = 35, by = 5)
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()

for (i in 1:length(eval_num)) {
  eval_number <- eval_num[i]
  begin_time <- Sys.time()
  result_ad <- BO_adap(func = objective_func,
                       # update_step = round(eval_number/2),
                       update_step = 10,
                       number_eval = (eval_number - 5),
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 1000,
                       noise_var = noise_var,
                       num_initial = 5,
                       grid_points = 1000)
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth
  
  # surrogate2 <- function(xvalue, data_to_smooth){
  #   pp <- predict_gp(data = data_to_smooth, x_pred = xvalue, noise_var = noise_var, choice_cov = square_exp_cov_generator(length_scale = result_ad$length_scale, signal_var = result_ad$signal_var))
  #   pp$mean
  # }

  ff <- list()
  # ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
  ff$fn <- function(x) as.numeric(surrogate(x, data_to_smooth = data_to_smooth))
  # ff$fn <- function(x) as.numeric(surrogate2(x, data_to_smooth = data_to_smooth))
  
  x_vals <- (seq(from = lower, to = upper, length.out = 1000) - lower)/(upper - lower)
  # y_vals <- qnorm(x_vals)
  # fn_vals <- sapply(y_vals, ff$fn)
  fn_vals <- sapply(x_vals, ff$fn)
  # obj <- function(y) {exp(ff$fn(y))}
  obj <- function(x) {exp(ff$fn(x))}
  lognormal_const <- log(integrate(obj, lower = 0, upper = 1)$value)
  post_x <- data.frame(y = x_vals, pos = exp(fn_vals - lognormal_const))
  BO_result_list[[i]] <- data.frame(x = (lower + x_vals*(upper - lower)), pos = post_x$pos /(upper - lower))
  end_time <- Sys.time()
  rel_runtime[i] <- as.numeric((end_time - begin_time), units = "mins")/1.961008
}
save(BO_result_original_list, file = "BO_data_to_smooth.rda")
save(BO_result_list, file = "BO_result_list.rda")

#### Plot the Rel-Runtime:
plot(rel_runtime ~ eval_num, type = "o")


#### Comparison:
load(file = "exact_grid_result.rda")
load(file = "exact_grid_result_smooth.rda")
load(file = "mcmc_samps.rda")
burnin <- 1000
thinning <- 3
mcmc_samps_selected <- mcmc_samps[-c(1:burnin)][seq(1, length(mcmc_samps[-c(1:burnin)]), by=thinning)]

for (i in 1:length(eval_num)) {
  ggplot() +
    geom_histogram(aes(x = mcmc_samps_selected, y = ..density..), bins = 300, alpha = 0.8, fill = "skyblue") +
    geom_line(data = BO_result_list[[i]], aes(x = x, y = pos), color = "red", size = 1) +
    geom_line(data = exact_grid_result_smooth, aes(x = x, y = pos), color = "black", size = 0.5, linetype = "dashed") +
    ggtitle(paste0("Comparison Posterior Density: B = ", eval_num[i])) +
    xlab("Value") +
    ylab("Density") +
    lims(x = c(5,8)) +
    lims(y = c(0,10)) +
    theme(text=element_text(size=10)) +
    theme_minimal() 
    ggsave(filename = paste0("figures/change_point_BO_", eval_num[i], ".pdf"), width = 5, height = 5)
}




#### Compute the KL distance:
Compute_KL <- function(x, qx, px){
  to_kept <- which(px > 0)
  x <- x[to_kept]
  qx <- qx[to_kept]
  px <- px[to_kept]
  # px <- px + .Machine$double.eps
  # qx <- qx + .Machine$double.eps
  dx <- diff(x)
  left <- c(0,dx)
  right <- c(dx,0)
  0.5 * sum(left * log(px/qx) * px) + 0.5 * sum(right * log(px/qx) * px)
}
KL_vec <- c()
for (i in 1:length(eval_num)) {
  KL_vec[i] <- Compute_KL(x = exact_grid_result_smooth$x, px = exact_grid_result_smooth$pos, qx = BO_result_list[[i]]$pos)
}
plot(log(KL_vec) ~ eval_num, type = "o")



#### Compute the KS distance:
Compute_KS <- function(x, qx, px){
  dx <- c(diff(x),0)
  max(abs(cumsum(qx * dx) - cumsum(px * dx)))
}
KS_vec <- c()
for (i in 1:length(eval_num)) {
  KS_vec[i] <- Compute_KS(x = exact_grid_result_smooth$x, px = exact_grid_result_smooth$pos, qx = BO_result_list[[i]]$pos)
}
png(filename = "figures/ks_compare.png", height = 500, width = 500)
plot((KS_vec) ~ eval_num, type = "o", ylab = "KS", xlab = "eval number: B")
dev.off()
