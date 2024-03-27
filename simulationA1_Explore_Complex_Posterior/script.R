library(npreg)
library(ggplot2)
library(aghq)

set.seed(123)
noise_var = 1e-8
source_location <- "function/"
source(file = paste0(source_location, "bo_functions.R"))
num_initial <- 3
surrogate <- function(xvalue, data_to_smooth){
  predict(ss(x = data_to_smooth$x, y = data_to_smooth$y, df = length(unique(data_to_smooth$x)), m = 2, all.knots = TRUE), x = xvalue)$y
}
lower = 0
upper = 10

############################################################
##### Function to integrate with aghq: #####################
##### here f is a function mapping from the entire real line:
integrate_aghq <- function(f, k = 100, startingvalue = 0){
  ff <- list(fn = f, gr = function(x) numDeriv::grad(f, x), he = function(x) numDeriv::hessian(f, x))
  aghq(ff = ff, k = k, startingvalue = startingvalue)$normalized_posterior$lognormconst
}


#########################################
#########################################
#########################################
### Easy Example:
log_prior <- function(x){
  # dnorm(x = x, mean = 1, log = T, sd = 2)
  1
}
log_likelihood <- function(x){
  # dnorm(x = x, mean = 3, log = T, sd = 2)
  x*sin(x)
}
eval_once <- function(x){
  log_prior(x) + log_likelihood(x)
}
eval_once_mapped <- function(y){
  eval_once(pnorm(y) * (upper - lower) + lower) + dnorm(y, log = T) + log(upper - lower)
}
x <- seq(0.01,9.99, by = 0.01)
y <- qnorm((x - lower)/(upper - lower))
true_log_norm_constant <- integrate_aghq(f = function(y) eval_once_mapped(y))
true_log_post_mapped <- function(y) {eval_once_mapped(y) - true_log_norm_constant}
plot((true_log_post_mapped(y)) ~ y, type = "l", cex.lab = 2.5, cex.axis = 2.5)
true_log_post <- function(x) {true_log_post_mapped(qnorm((x - lower)/(upper - lower))) - dnorm(qnorm((x - lower)/(upper - lower)), log = T) - log(upper - lower)}
integrate(function(x) exp(true_log_post(x)), lower = 0, upper = 10)

### Apply BO:
objective_func <- eval_once
eval_num <- seq(5, 100, by = 5)
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()
for (i in 1:length(eval_num)) {
  eval_number <- eval_num[i]
  result_ad <- BO_adap(func = objective_func,
                       update_step = 10,
                       number_eval = eval_number - num_initial,
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 100,
                       noise_var = noise_var,
                       initial = NULL,
                       num_initial = num_initial,
                       grid_points = 1000)
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth

  ff <- list()
  ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
  fn_vals <- sapply(y, ff$fn)

  lognormal_const <- integrate_aghq(f = ff$fn)
  post_y <- data.frame(y = y, pos = exp(fn_vals - lognormal_const))
  post_x <- data.frame(x = pnorm(post_y$y) * (upper - lower) + lower, post = (post_y$pos / dnorm(post_y$y))/(upper - lower) )
  
  BO_result_list[[i]] <- post_x
  
  to_plot_data <- BO_result_list[[i]]
  to_plot_data$logpos <- log(to_plot_data$pos)
  to_plot_data_obs <- BO_result_original_list[[i]]
  y_min <- -15
  y_max <- 0
  png(filename = paste0("figures/easy/", "BOSS_approxi_B_", eval_number, "_easy.png"), width = 800, height = 800)
  mar.default <- c(5,4,4,2)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(to_plot_data$logpos ~ to_plot_data$x, 
       type = "l", lty = "dashed", col = "red", lwd = 3,
       ylab = "Log Post", xlab = expression(alpha),
       ylim = c(y_min, y_max), cex.lab = 2.5, cex.axis = 2.5)
  lines((true_log_post(x)) ~ x, lwd = 1)
  y_offset <- 0.03 * (y_max - y_min) # adjust offset as needed
  for(x_val in to_plot_data_obs$x_original) {
    segments(x_val, y_min, x_val, y_min - y_offset, col = "red")
  }
  dev.off()
  
  png(filename = paste0("figures/easy/", "pos_BOSS_approxi_B_", eval_number, "_easy.png"), width = 800, height = 800)
  mar.default <- c(5,4,4,2)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(exp(to_plot_data$logpos) ~ to_plot_data$x, 
       type = "l", lty = "dashed", col = "red", lwd = 2,
       ylab = "Post", xlab = expression(alpha), ylim = c(0,2), cex.lab = 2.5, cex.axis = 2.5)
  lines(exp(true_log_post(x)) ~ x, lwd = 1)
  y_min <- -0.03
  y_offset <- 0.01 * (2) # adjust offset as needed
  for(x_val in to_plot_data_obs$x_original) {
    segments(x_val, y_min, x_val, y_min - y_offset, col = "red")
  }
  dev.off()
}

#### Compute the KL distance:
Compute_KL <- function(x, logpx, logqx){
  dx <- diff(x)
  left <- c(0,dx)
  right <- c(dx,0)
  0.5 * sum(left * (logpx - logqx) * exp(logpx)) + 0.5 * sum(right * (logpx - logqx) * exp(logpx))
}
KL_vec <- c()
for (i in 1:length(eval_num)) {
  KL_vec[i] <- Compute_KL(x = x, logpx = true_log_post(x), logqx = log(BO_result_list[[i]]$pos))
}
png(filename = "figures/easy/kl_compare_easy.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot((KL_vec) ~ eval_num, type = "o", ylab = "KL", xlab = "eval number: B", cex.lab = 2.5, cex.axis = 2.5)
dev.off()

#### Compute the KS distance:
Compute_KS <- function(x, qx, px){
  dx <- c(diff(x),0)
  max(abs(cumsum(qx * dx) - cumsum(px * dx)))
}
KS_vec <- c()
for (i in 1:length(eval_num)) {
  KS_vec[i] <- Compute_KS(x = x, px = exp(true_log_post(x)), qx = BO_result_list[[i]]$pos)
}
png(filename = "figures/easy/ks_compare_easy.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot((KS_vec) ~ eval_num, type = "o", ylab = "KS", xlab = "eval number: B", ylim = c(0,1), cex.lab = 2.5, cex.axis = 2.5)
dev.off()

png(filename = "figures/easy/log_ks_compare_easy.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(KS_vec) ~ eval_num, type = "o", ylab = "Log KS", xlab = "eval number: B", ylim = c(-10,0), cex.lab = 2.5, cex.axis = 2.5)
dev.off()




#########################################
#########################################
#########################################
### Medium Example:
log_prior <- function(x){
  # dnorm(x = x, mean = 1, log = T, sd = 2)
  1
}
log_likelihood <- function(x){
  # dnorm(x = x, mean = 3, log = T, sd = 2)
  log(x+1)*sin(x*2) - x*cos(x*2)
}
eval_once <- function(x){
  log_prior(x) + log_likelihood(x)
}
eval_once_mapped <- function(y){
  eval_once(pnorm(y) * (upper - lower) + lower) + dnorm(y, log = T) + log(upper - lower)
}
x <- seq(0.01,9.99, by = 0.01)
y <- qnorm((x - lower)/(upper - lower))
true_log_norm_constant <- integrate_aghq(f = function(y) eval_once_mapped(y), k = 100)
true_log_post_mapped <- function(y) {eval_once_mapped(y) - true_log_norm_constant}
plot((true_log_post_mapped(y)) ~ y, type = "l")
true_log_post <- function(x) {true_log_post_mapped(qnorm((x - lower)/(upper - lower))) - dnorm(qnorm((x - lower)/(upper - lower)), log = T) - log(upper - lower)}
integrate(function(x) exp(true_log_post(x)), lower = 0, upper = 10)

### Apply BO:
objective_func <- eval_once
eval_num <- seq(5, 100, by = 5)
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()
for (i in 1:length(eval_num)) {
  eval_number <- eval_num[i]
  result_ad <- BO_adap(func = objective_func,
                       update_step = 10,
                       number_eval = eval_number - num_initial,
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 100,
                       noise_var = noise_var,
                       initial = NULL,
                       num_initial = num_initial,
                       grid_points = 1000)
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth
  
  ff <- list()
  ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
  fn_vals <- sapply(y, ff$fn)
  
  lognormal_const <- integrate_aghq(f = ff$fn)
  post_y <- data.frame(y = y, pos = exp(fn_vals - lognormal_const))
  post_x <- data.frame(x = pnorm(post_y$y) * (upper - lower) + lower, post = (post_y$pos / dnorm(post_y$y))/(upper - lower) )
  
  BO_result_list[[i]] <- post_x
  
  to_plot_data <- BO_result_list[[i]]
  to_plot_data$logpos <- log(to_plot_data$pos)
  to_plot_data_obs <- BO_result_original_list[[i]]
  y_min <- -20
  y_max <- 5
  png(filename = paste0("figures/medium/", "BOSS_approxi_B_", eval_number, "_med.png"), width = 800, height = 800)
  mar.default <- c(5,4,4,2)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(to_plot_data$logpos ~ to_plot_data$x, 
       type = "l", lty = "dashed", col = "red", lwd = 3,
       ylab = "Log Post", xlab = expression(alpha),
       ylim = c(y_min, y_max), cex.lab = 2.5, cex.axis = 2.5)
  lines((true_log_post(x)) ~ x, lwd = 1)
  y_offset <- 0.03 * (y_max - y_min) # adjust offset as needed
  for(x_val in to_plot_data_obs$x_original) {
    segments(x_val, y_min, x_val, y_min - y_offset, col = "red")
  }
  dev.off()
  
  y_min <- 0
  y_max <- 2.5
  png(filename = paste0("figures/medium/", "pos_BOSS_approxi_B_", eval_number, "_med.png"), width = 800, height = 800)
  mar.default <- c(5,4,4,2)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(exp(to_plot_data$logpos) ~ to_plot_data$x, 
       type = "l", lty = "dashed", col = "red", lwd = 2,
       ylab = "Post", xlab = expression(alpha), ylim = c(y_min,y_max), cex.lab = 2.5, cex.axis = 2.5)
  lines(exp(true_log_post(x)) ~ x, lwd = 1)
  y_offset <- 0.01 * (y_max - y_min) # adjust offset as needed
  for(x_val in to_plot_data_obs$x_original) {
    segments(x_val, y_min, x_val, y_min - y_offset, col = "red")
  }
  dev.off()
}

#### Compute the KL distance:
Compute_KL <- function(x, logpx, logqx){
  dx <- diff(x)
  left <- c(0,dx)
  right <- c(dx,0)
  0.5 * sum(left * (logpx - logqx) * exp(logpx)) + 0.5 * sum(right * (logpx - logqx) * exp(logpx))
}
KL_vec <- c()
for (i in 1:length(eval_num)) {
  KL_vec[i] <- Compute_KL(x = x, logpx = true_log_post(x), logqx = log(BO_result_list[[i]]$pos))
}
png(filename = "figures/medium/kl_compare_med.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot((KL_vec) ~ eval_num, type = "o", ylab = "KL", xlab = "eval number: B", cex.lab = 2.5, cex.axis = 2.5)
dev.off()

#### Compute the KS distance:
Compute_KS <- function(x, qx, px){
  dx <- c(diff(x),0)
  max(abs(cumsum(qx * dx) - cumsum(px * dx)))
}
KS_vec <- c()
for (i in 1:length(eval_num)) {
  KS_vec[i] <- Compute_KS(x = x, px = exp(true_log_post(x)), qx = BO_result_list[[i]]$pos)
}
png(filename = "figures/medium/ks_compare_med.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot((KS_vec) ~ eval_num, type = "o", ylab = "KS", xlab = "eval number: B", ylim = c(0,1), cex.lab = 2.5, cex.axis = 2.5)
dev.off()

png(filename = "figures/medium/log_ks_compare_med.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(KS_vec) ~ eval_num, type = "o", ylab = "Log KS", xlab = "eval number: B", ylim = c(-10,0), cex.lab = 2.5, cex.axis = 2.5)
dev.off()



#########################################
#########################################
#########################################
### Very Hard Example:
log_prior <- function(x){
  1
}
log_likelihood <- function(x){
  log(x + 1) * (sin(x * 4) + cos(x * 2))
}
eval_once <- function(x){
  log_prior(x) + log_likelihood(x)
}
eval_once_mapped <- function(y){
  eval_once(pnorm(y) * (upper - lower) + lower) + dnorm(y, log = T) + log(upper - lower)
}
x <- seq(0.01,9.99, by = 0.01)
y <- qnorm((x - lower)/(upper - lower))
# true_log_norm_constant <- integrate_aghq(f = function(y) eval_once_mapped(y), k = 100)
true_log_norm_constant <- log(integrate(f = function(y) exp(eval_once_mapped(y)), lower = -Inf, upper = Inf)$value)
true_log_post_mapped <- function(y) {eval_once_mapped(y) - true_log_norm_constant}
plot((true_log_post_mapped(y)) ~ y, type = "l")
plot(exp(true_log_post_mapped(y)) ~ y, type = "l")
true_log_post <- function(x) {true_log_post_mapped(qnorm((x - lower)/(upper - lower))) - dnorm(qnorm((x - lower)/(upper - lower)), log = T) - log(upper - lower)}
integrate(function(x) exp(true_log_post(x)), lower = 0, upper = 10)

### Apply BO:
objective_func <- eval_once
eval_num <- seq(5, 100, by = 5)
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()
for (i in 1:length(eval_num)) {
  eval_number <- eval_num[i]
  result_ad <- BO_adap(func = objective_func,
                       update_step = 10,
                       number_eval = eval_number - num_initial,
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 100,
                       noise_var = noise_var,
                       initial = NULL,
                       num_initial = num_initial,
                       grid_points = 1000)
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth
  
  ff <- list()
  ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
  fn_vals <- sapply(y, ff$fn)
  # lognormal_const <- integrate_aghq(f = ff$fn, k = 100)
  lognormal_const <- log(integrate(f = function(y) exp(ff$fn(y)), lower = -Inf, upper = Inf)$value)
  
  post_y <- data.frame(y = y, pos = exp(fn_vals - lognormal_const))
  post_x <- data.frame(x = pnorm(post_y$y) * (upper - lower) + lower, post = (post_y$pos / dnorm(post_y$y))/(upper - lower) )
  BO_result_list[[i]] <- post_x
  to_plot_data <- BO_result_list[[i]]
  to_plot_data$logpos <- log(to_plot_data$pos)
  to_plot_data_obs <- BO_result_original_list[[i]]
  y_min <- -10
  y_max <- 5
  png(filename = paste0("figures/very_hard/", "BOSS_approxi_B_", eval_number, "_very_hard.png"), width = 800, height = 800)
  mar.default <- c(5,4,4,2)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(to_plot_data$logpos ~ to_plot_data$x, 
       type = "l", lty = "dashed", col = "red", lwd = 3,
       ylab = "Log Post", xlab = expression(alpha),
       ylim = c(y_min, y_max), cex.lab = 2.5, cex.axis = 2.5)
  lines((true_log_post(x)) ~ x, lwd = 1)
  y_offset <- 0.03 * (y_max - y_min) # adjust offset as needed
  for(x_val in to_plot_data_obs$x_original) {
    segments(x_val, y_min, x_val, y_min - y_offset, col = "red")
  }
  dev.off()
  
  y_min <- 0
  y_max <- 3
  png(filename = paste0("figures/very_hard/", "pos_BOSS_approxi_B_", eval_number, "_very_hard.png"), width = 800, height = 800)
  mar.default <- c(5,4,4,2)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(exp(to_plot_data$logpos) ~ to_plot_data$x, 
       type = "l", lty = "dashed", col = "red", lwd = 2,
       ylab = "Post", xlab = expression(alpha), ylim = c(y_min,y_max), cex.lab = 2.5, cex.axis = 2.5)
  lines(exp(true_log_post(x)) ~ x, lwd = 1)
  y_offset <- 0.01 * (y_max - y_min) # adjust offset as needed
  for(x_val in to_plot_data_obs$x_original) {
    segments(x_val, y_min, x_val, y_min - y_offset, col = "red")
  }
  dev.off()
}

#### Compute the KL distance:
Compute_KL <- function(x, logpx, logqx){
  dx <- diff(x)
  left <- c(0,dx)
  right <- c(dx,0)
  0.5 * sum(left * (logpx - logqx) * exp(logpx)) + 0.5 * sum(right * (logpx - logqx) * exp(logpx))
}
KL_vec <- c()
for (i in 1:length(eval_num)) {
  KL_vec[i] <- Compute_KL(x = x, logpx = true_log_post(x), logqx = log(BO_result_list[[i]]$pos))
}
png(filename = "figures/very_hard/kl_compare_very_hard.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot((KL_vec) ~ eval_num, type = "o", ylab = "KL", xlab = "eval number: B", cex.lab = 2.5, cex.axis = 2.5)
dev.off()

#### Compute the KS distance:
Compute_KS <- function(x, qx, px){
  dx <- c(diff(x),0)
  max(abs(cumsum(qx * dx) - cumsum(px * dx)))
}
KS_vec <- c()
for (i in 1:length(eval_num)) {
  KS_vec[i] <- Compute_KS(x = x, px = exp(true_log_post(x)), qx = BO_result_list[[i]]$pos)
}
png(filename = "figures/very_hard/ks_compare_very_hard.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot((KS_vec) ~ eval_num, type = "o", ylab = "KS", xlab = "eval number: B", ylim = c(0,1), cex.lab = 2.5, cex.axis = 2.5)
dev.off()

png(filename = "figures/very_hard/log_ks_compare_very_hard.png", height = 800, width = 800)
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(KS_vec) ~ eval_num, type = "o", ylab = "Log KS", xlab = "eval number: B", ylim = c(-10,0), cex.lab = 2.5, cex.axis = 2.5)
dev.off()
