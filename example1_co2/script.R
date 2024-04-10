.libPaths( c( .libPaths(), "~/./lib") )

#### PATH:
working_path <- getwd()
source_path <- paste0(working_path,"/source/")
figure_path <- paste0(working_path,"/figures/")
result_path <- paste0(working_path,"/results/")

source(file = paste0(source_path, "00_load.R"))
source(file = paste0(source_path, "bo_functions.R"))

### Read in the full data:
cUrl = paste0("http://scrippsco2.ucsd.edu/assets/data/atmospheric/",
              "stations/flask_co2/daily/daily_flask_co2_mlo.csv")
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header = FALSE, sep = ",",
                  skip = 69, stringsAsFactors = FALSE, col.names = c("day",
                                                                     "time", "junk1", "junk2", "Nflasks", "quality",
                                                                     "co2"))
co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M",
                     tz = "UTC")
co2s$day = as.Date(co2s$date)
timeOrigin = as.Date("1960/03/30")
co2s$timeYears = round(as.numeric(co2s$day - timeOrigin)/365.25,
                         3)
co2s$dayInt = as.integer(co2s$day)
allDays = seq(from = min(co2s$day), to = max(co2s$day),
              by = "7 day")
observed_dataset <- co2s %>% filter(!is.na(co2s$co2)) %>% dplyr::select(c("co2", "timeYears"))
observed_dataset$quality <- ifelse(co2s$quality > 0, 1, 0)
observed_dataset <- observed_dataset %>% filter(quality == 0)
### Create the covariate for trend and the 1-year seasonality 
observed_dataset$t1 <- observed_dataset$timeYears
observed_dataset$t2 <- observed_dataset$timeYears
observed_dataset$t3 <- observed_dataset$timeYears
observed_dataset$t4 <- observed_dataset$timeYears


lower <- 2
upper <- 6
eval_once <- function(a){
  ### reduce the number of harmonic of 2*pi/a to just one. (previously two)
  mod_once <- model_fit(formula = co2 ~ f(x = t1, model = "IWP", order = 2, sd.prior = list(param = 30, h = 10), k = 30, initial_location = median(observed_dataset$t1)) +
                          f(x = t2, model = "sGP", sd.prior = list(param = 1, h = 10), a = (2*pi), m = 2, k = 30) + 
                          f(x = t3, model = "sGP", sd.prior = list(param = 1, h = 10), a = (2*pi/a), m = 1, k = 30),
                        data = observed_dataset,
                        control.family = list(sd.prior = list(param = 1)),
                        family = "Gaussian", aghq_k = 3)
  mod_once$mod$normalized_posterior$lognormconst
}
begin_time <- Sys.time()
BO_result <- BO_adap(eval_once, update_step = 10, number_eval = 30, 
        lower = lower, upper = upper, noise_var = 1e-6, 
        length_scale = 1, signal_var = 1, initial = NULL,
        num_initial = 5,
        delta = 0.01, grid_points = 100)
end_time <- Sys.time()
end_time - begin_time
save(file = "BO_result.rda", BO_result)

data_to_smooth <- BO_result$result[order(BO_result$result$x), ]
plot(exp(y) ~ x_original, data_to_smooth, type = "o")

# Construct the surrogate:
# surrogate <- function(xvalue, data_to_smooth){
#   predict(ss(x = data_to_smooth$x, y = data_to_smooth$y, df = length(unique(data_to_smooth$x)), m = 2, all.knots = TRUE), x = xvalue)$y
# }
surrogate <- function(xvalue, data_to_smooth){
  pp <- predict_gp(data = data_to_smooth, x_pred = xvalue, noise_var = 1e-6, choice_cov = square_exp_cov_generator(length_scale = BO_result$length_scale, signal_var = BO_result$signal_var))
  pp$mean
}
integrate_aghq <- function(f, k = 100, startingvalue = 0){
  ff <- list(fn = f, gr = function(x) numDeriv::grad(f, x), he = function(x) numDeriv::hessian(f, x))
  aghq::aghq(ff = ff, k = k, startingvalue = startingvalue)$normalized_posterior$lognormconst
}

x <- seq(2.01, 5.99, by = 0.01)
y <- qnorm((x - lower)/(upper - lower))
ff <- list()
ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
fn_vals <- sapply(y, ff$fn)
plot(exp(fn_vals) ~ y, type = "l")

lognormal_const <- integrate_aghq(f = ff$fn, k = 10)
post_y <- data.frame(y = y, pos = exp(fn_vals - lognormal_const))
post_x <- data.frame(x = pnorm(post_y$y) * (upper - lower) + lower, post = (post_y$pos / dnorm(post_y$y))/(upper - lower) )

### Plot the posterior of alpha:
tikzDevice::tikz(file = paste0("co2_pi_alpha.tex"), width = 5, height = 5, standAlone = TRUE)
plot(post_x$post ~ post_x$x, type = "l", ylab = "Post", 
     xlab = '$\\alpha$' , cex.lab = 1.5, cex.axis = 1.5)
for(x_val in data_to_smooth$x_original) {
  segments(x_val, -0.02, x_val, -0.05, col = "red")
}
dev.off()
system(paste0("pdflatex ", "co2_pi_alpha.tex"))
file.remove(paste0("co2_pi_alpha.aux"))
file.remove(paste0("co2_pi_alpha.log"))
file.remove(paste0("co2_pi_alpha.tex"))

### Compute the HPD:
post_x$x <- round(post_x$x, digits = 2)
my_cdf <- cumsum(post_x$post * c(diff(post_x$x), 0))
my_cdf[which(post_x$x == 3.47)] - my_cdf[which(post_x$x == 3.11)]
plot(post_x$post ~ post_x$x, type = "l", ylab = "Post", 
     xlab = expression(alpha), cex.lab = 2.0, cex.axis = 2.0)
abline(v = 3.47, col = "purple")
abline(v = 3.11, col = "purple")

my_cdf[which(post_x$x == 3.82)] - my_cdf[which(post_x$x == 2.82)]
plot(post_x$post ~ post_x$x, type = "l", ylab = "Post", 
     xlab = expression(alpha), cex.lab = 2.0, cex.axis = 2.0)
abline(v = 3.82, col = "purple")
abline(v = 2.82, col = "purple")

### Approximate Quadrature rules from BOSS:
aghq_boss <- function(f, k = 100, data_to_smooth){
  ff <- list(fn = f, gr = function(x) numDeriv::grad(f, x), he = function(x) numDeriv::hessian(f, x))
  opt_result <- list(ff = ff, mode = qnorm(data_to_smooth$x[which.max(data_to_smooth$y)]))
  opt_result$hessian = -matrix(ff$he(opt_result$mode))
  aghq::aghq(ff = ff, k = k, optresults = opt_result, startingvalue = opt_result$mode)
}
aghq_boss_obj <- aghq_boss(f = ff$fn, k = 10, data_to_smooth)

nodes <- aghq_boss_obj[[1]]$grid$nodes
nodes_converted <- as.numeric(pnorm(nodes)*(upper - lower) + lower)
L <- as.numeric(sqrt(aghq_boss_obj$normalized_posterior$grid$features$C))
weights <- as.numeric(aghq_boss_obj[[1]]$nodesandweights$weights)
lambda <- weights * exp(aghq_boss_obj[[1]]$nodesandweights$logpost_normalized)

### Fit new models at these nodes:
fit_once <- function(a){
  ### reduce the number of harmonic of 2*pi/a to just one. (previously two)
  mod_once <- model_fit(formula = co2 ~ f(x = t1, model = "IWP", order = 2, sd.prior = list(param = 30, h = 10), k = 30, initial_location = median(observed_dataset$t1)) +
                          f(x = t2, model = "sGP", sd.prior = list(param = 1, h = 10), a = (2*pi), m = 2, k = 30) + 
                          f(x = t3, model = "sGP", sd.prior = list(param = 1, h = 10), a = (2*pi/a), m = 1, k = 30),
                        data = observed_dataset,
                        control.family = list(sd.prior = list(param = 1)),
                        family = "Gaussian", aghq_k = 3)
  mod_once
}
set.seed(123)
for (i in 1:length(nodes_converted)){
  mod <- fit_once(nodes_converted[i]) 
  save(mod, file = paste0("models/", "model", "_", i,".rda"))
}


### Inferring the latent field:
num_samples <- round(3000 * lambda, 0)

t1_samples <- data.frame(t = seq(0,62, by = 0.01))
t2_samples <- data.frame(t = seq(0,62, by = 0.01))
t3_samples <- data.frame(t = seq(0,62, by = 0.01))

for (i in 1:length(num_samples)) {
  load(file = paste0("models/model_", i, ".rda"))
  pred <- predict(mod, variable = "t1", only.samples = T, newdata = data.frame(t1 = seq(0,62, by = 0.01)))[,-1][,1:num_samples[i]]
  t1_samples <- cbind(t1_samples, pred)
  pred <- predict(mod, variable = "t2", only.samples = T, newdata = data.frame(t2 = seq(0,62, by = 0.01)))[,-1][,1:num_samples[i]]
  t2_samples <- cbind(t2_samples, pred)
  pred <- predict(mod, variable = "t3", only.samples = T, newdata = data.frame(t3 = seq(0,62, by = 0.01)))[,-1][,1:num_samples[i]]
  t3_samples <- cbind(t3_samples, pred)
}

t1_summary <- extract_mean_interval_given_samps(t1_samples)
t2_summary <- extract_mean_interval_given_samps(t2_samples)
t3_summary <- extract_mean_interval_given_samps(t3_samples)

t_all_samples <- t1_samples + t2_samples + t3_samples
t_all_samples[,1] <- t1_samples[,1]
t_all_summary <- extract_mean_interval_given_samps(t_all_samples)
t_all_summary$time <- (t_all_summary$x * 365.25) + timeOrigin

## Plot the results to pdf
tikzDevice::tikz(file = paste0("co2_overall.tex"), width = 5, height = 5, standAlone = TRUE)
plot(t_all_summary$q0.5 ~ t_all_summary$time, type = "l", 
     lty = "solid", ylab = "CO2", xlab = "year", cex.lab = 1.5, cex.axis = 1.5)
lines(t_all_summary$q0.975 ~ t_all_summary$time, col = "red", lty = "dashed")
lines(t_all_summary$q0.025 ~ t_all_summary$time, col = "red", lty = "dashed")
dev.off()
system(paste0("pdflatex ", "co2_overall.tex"))
file.remove(paste0("co2_overall.aux"))
file.remove(paste0("co2_overall.log"))
file.remove(paste0("co2_overall.tex"))

t1_summary$time <- (t1_summary$x * 365.25) + timeOrigin

tikzDevice::tikz(file = paste0("co2_trend.tex"), width = 5, height = 5, standAlone = TRUE)
plot(t1_summary$q0.5 ~ t1_summary$time, type = "l", lty = "solid", ylab = "CO2", 
     xlab = "year", cex.lab = 1.5, cex.axis = 1.5)
lines(t1_summary$q0.975 ~ t1_summary$time, col = "red", lty = "dashed")
lines(t1_summary$q0.025 ~ t1_summary$time, col = "red", lty = "dashed")
dev.off()
system(paste0("pdflatex ", "co2_trend.tex"))
file.remove(paste0("co2_trend.aux"))
file.remove(paste0("co2_trend.log"))
file.remove(paste0("co2_trend.tex"))


ts_samples <- t2_samples + t3_samples
ts_samples[,1] <- t2_samples[,1]
ts_summary <- extract_mean_interval_given_samps(ts_samples)
ts_summary$time <- (ts_summary$x * 365.25) + timeOrigin

tikzDevice::tikz(file = paste0("co2_seasonal.tex"), width = 5, height = 5, standAlone = TRUE)
plot(ts_summary$q0.5 ~ ts_summary$time, type = "l", lty = "solid", ylab = "CO2", xlab = "year",
     xlim = as.Date(c("1985-01-01","2000-01-01")), ylim = c(725, 740), cex.lab = 1.5, cex.axis = 1.5)
lines(ts_summary$q0.975 ~ ts_summary$time, col = "red", lty = "dashed")
lines(ts_summary$q0.025 ~ ts_summary$time, col = "red", lty = "dashed")
dev.off()
system(paste0("pdflatex ", "co2_seasonal.tex"))
file.remove(paste0("co2_seasonal.aux"))
file.remove(paste0("co2_seasonal.log"))
file.remove(paste0("co2_seasonal.tex"))






