library(BayesGP)
library(npreg)
library(ggplot2)

set.seed(123)
noise_var = 1e-8
num_initial <- 3

source_location <- "function/"
source(file = paste0(source_location, "bo_functions.R"))

load("data.rda")
plot(y~x, data)

log_prior <- function(alpha){
  # dexp(x = alpha, rate = 5, log = T)
  dnorm(x = alpha, mean = 3, log = T, sd = 0.5)
}

eval_once <- function(alpha){
  a_fit <- (2*pi)/alpha
  x <- data$x
  data$cosx <- cos(a_fit * x)
  data$sinx <- sin(a_fit * x)
  data$cos2x <- cos(2*a_fit * x)
  data$sin2x <- sin(2*a_fit * x)
  mod <- model_fit(formula = y ~ cosx + sinx + cos2x + sin2x + f(x = indx, model = "IID", 
                                                                 sd.prior = list(param = 1)),
                   data = data, method = "aghq", family = "Poisson", aghq_k = 4
  )
  (mod$mod$normalized_posterior$lognormconst) + log_prior(alpha)
}


surrogate <- function(xvalue, data_to_smooth){
  predict(ss(x = data_to_smooth$x, y = data_to_smooth$y, df = length(unique(data_to_smooth$x)), m = 2, all.knots = TRUE), x = xvalue)$y
}

lower = 0.5
upper = 4.5
a <- 1.5
objective_func <- eval_once


#### Runtime analysis:
eval_num <- seq(from = 10, to = 80, by = 5)
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()

load(file = "exact_grid_result.rda")

for (i in 1:length(eval_num)) {
  n_grid <- nrow(exact_grid_result)
  eval_number <- eval_num[i]
  begin_time <- Sys.time()
  result_ad <- BO_adap(func = objective_func,
                       # update_step = round(eval_number/2),
                       update_step = 10,
                       number_eval = eval_number - num_initial,
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 1000,
                       noise_var = noise_var,
                       initial = NULL,
                       num_initial = num_initial,
                       grid_points = n_grid)
  end_time <- Sys.time()
  rel_runtime[i] <- as.numeric((end_time - begin_time), units = "mins")/1.344585
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth

  ff <- list()
  ff$fn <- function(x) as.numeric(surrogate(x, data_to_smooth = data_to_smooth))
  x_vals <- (seq(from = lower, to = upper, length.out = n_grid) - lower)/(upper - lower)
  fn_vals <- sapply(x_vals, ff$fn)
  obj <- function(x) {exp(ff$fn(x))}
  lognormal_const <- log(integrate(obj, lower = 0, upper = 1)$value)
  post_x <- data.frame(y = x_vals, pos = exp(fn_vals - lognormal_const))
  BO_result_list[[i]] <- data.frame(x = (lower + x_vals*(upper - lower)), pos = post_x$pos /(upper - lower))
}
save(BO_result_list, file = "BO_result_list.rda")

#### Plot the Rel-Runtime:
png(filename = "figures/runtime_compare.png", height = 800, width = 800)
plot(rel_runtime ~ eval_num, type = "o", ylab = "rel-runtime", xlab = "eval number: B", cex.lab = 2.0, cex.axis = 2.0)
dev.off()

#### Comparison:
load(file = "exact_grid_result.rda")
load(file = "exact_grid_result_smooth.rda")
# load(file = "mcmc_samps.rda")
# burnin <- 1000
# thinning <- 3
# mcmc_samps_selected <- mcmc_samps[-c(1:burnin)][seq(1, length(mcmc_samps[-c(1:burnin)]), by=thinning)]

for (i in 1:length(eval_num)) {
  ggplot() +
    # geom_histogram(aes(x = mcmc_samps_selected, y = ..density..), bins = 100, alpha = 0.5, fill = "blue") +
    geom_line(data = BO_result_list[[i]], aes(x = x, y = pos), color = "red", size = 1) +
    geom_line(data = exact_grid_result_smooth, aes(x = x, y = pos), color = "black", size = 1, linetype = "dashed") + 
    ggtitle(paste0("Comparison Posterior Density: B = ", eval_num[i])) +
    xlab("Value") +
    ylab("Density") +
    theme_minimal() +
    theme(text = element_text(size = 20), axis.text = element_text(size = 25)) + # only change the lab and axis text size
    # lims(y = c(0,15))
    lims(y = range(exact_grid_result_smooth$pos))
  ggsave(filename = (paste0("figures/Comparison Posterior Density: B = ", eval_num[i], " .pdf")),
         width = 8, height = 8)
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
png(filename = "figures/kl_compare.png", height = 800, width = 800)
plot((KL_vec) ~ eval_num, type = "o", ylab = "KL", xlab = "eval number: B", cex.lab = 2.0, cex.axis = 2.0)
dev.off()

Compute_KL(x = exact_grid_result_smooth$x, px = exact_grid_result_smooth$pos, qx = BO_result_list[[i]]$pos)


#### Compute the KS distance:
Compute_KS <- function(x, qx, px){
  dx <- c(diff(x),0)
  max(abs(cumsum(qx * dx) - cumsum(px * dx)))
}
KS_vec <- c()
for (i in 1:length(eval_num)) {
  KS_vec[i] <- Compute_KS(x = exact_grid_result_smooth$x, px = exact_grid_result_smooth$pos, qx = BO_result_list[[i]]$pos)
}
png(filename = "figures/ks_compare.png", height = 800, width = 800)
plot((KS_vec) ~ eval_num, type = "o", ylab = "KS", xlab = "eval number: B", cex.lab = 2.0, cex.axis = 2.0)
dev.off()
