library(tidyverse)
library(BayesGP)
library(npreg)

set.seed(123)
noise_var = 1e-7
source_location <- "function/"
figure_location <- "figures/"
source(file = paste0(source_location, "bo_functions.R"))
surrogate <- function(xvalue, data_to_smooth){
  predict(ss(x = data_to_smooth$x, y = data_to_smooth$y, df = length(unique(data_to_smooth$x)), m = 2, all.knots = TRUE), x = xvalue)$y
}

cUrl = paste0("https://raw.githubusercontent.com/akarlinsky/world_mortality/main/world_mortality.csv")
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
world_death = read.table(cFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)

### West EU: NL
NL_death <- world_death %>% filter(country_name == "Netherlands")
NL_death$date <- make_date(year = NL_death$year) + weeks(NL_death$time)
NL_death$x <- as.numeric(NL_death$date)/365.25; 
ref_val <- min(NL_death$x)
NL_death$x <- NL_death$x - ref_val
plot(NL_death$deaths ~ NL_death$date)

fit_once <- function(alpha, data){
  a_fit <- alpha
  data$x1 <- ifelse(data$x <= a_fit, data$x, a_fit); data$xx1 <- cos(data$x1*2*pi); data$xxx1 <- sin(data$x1*2*pi)
  data$x2 <- ifelse(data$x > a_fit, (data$x - a_fit), 0) ; data$xx2 <- cos(data$x2*2*pi); data$xxx2 <- sin(data$x1*2*pi)
  
  mod <- model_fit(formula = deaths ~ xx1 + xx2 + xxx1 + xxx2 +
                     f(x1, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = "left") + f(x2, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = "left"),
                   data = data, method = "aghq", family = "Poisson", aghq_k = 3
  )
  mod
}
eval_once <- function(alpha, data = NL_death){
  mod <- fit_once(alpha = alpha, data = data)
  (mod$mod$normalized_posterior$lognormconst)
} 
lower = 0.5
upper = 7.2
objective_func <- eval_once

#### Runtime analysis:
eval_num <- 20
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()
output_dir <- "figures/"
for (i in 1:length(eval_num)) {
  eval_number <- eval_num[i]
  begin_time <- Sys.time()
  result_ad <- BO_adap(func = objective_func,
                       # update_step = round(eval_number/2),
                       update_step = 10,
                       number_eval = eval_number,
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 1000,
                       noise_var = noise_var,
                       initial = 5,
                       grid_points = 1000)
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth
  
  ff <- list()
  # ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
  ff$fn <- function(x) as.numeric(surrogate(x, data_to_smooth = data_to_smooth))
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
for (i in 1:length(eval_num)) {
  base_name <- paste0("NL_change_point: B = ", eval_num[i])
  tex_file_name <- paste0(output_dir, base_name, ".tex")
  pdf_file_name <- paste0(output_dir, base_name, ".pdf")
  BO_result_list[[i]]$year <- as.Date((BO_result_list[[i]]$x + ref_val)*365.25)
  tikzDevice::tikz(tex_file_name, width = 5, height = 5, standAlone = TRUE)
  ggplot() +
    # geom_histogram(aes(x = mcmc_samps_selected, y = ..density..), bins = 100, alpha = 0.5, fill = "blue") +
    geom_line(data = BO_result_list[[i]], aes(x = year, y = pos), color = "red", size = 1) +
    # geom_line(data = exact_grid_result_smooth, aes(x = x, y = pos), color = "black", size = 1, linetype = "dashed") +
    # ggtitle(paste0("Change-point in NL: B = ", eval_num[i])) +
    xlab("") +
    ylab("Density") +
    scale_x_date(
      limits = as.Date(c("2019-12-01", "2021-01-01")),
      date_labels = "%b %y", 
      date_breaks = "3 month"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 15), axis.text = element_text(size = 15))  # only change the lab and axis text size
    # lims(y = range(exact_grid_result_smooth$pos))
  dev.off()
  # Run pdflatex to generate the PDF
  if (file.exists(tex_file_name)) {
    system(sprintf('pdflatex -output-directory=%s "%s"', output_dir, tex_file_name))
    
    # Check if the PDF was created
    if (file.exists(pdf_file_name)) {
      # Delete the .tex, .aux, and .log files to clean up
      file.remove(tex_file_name)
      file.remove(paste0(output_dir, base_name, ".aux"))
      file.remove(paste0(output_dir, base_name, ".log"))
      file.remove(paste0(output_dir, base_name, ".tex"))
    } else {
      warning("PDF file was not created: ", pdf_file_name)
    }
  } else {
    warning("TeX file was not created: ", tex_file_name)
  }
  
  # ggsave(filename = (paste0("figures/NL_change_point: B = ", eval_num[i], " .pdf")))
}
### Which day is most likely: "2020-02-06"
as.Date((BO_result_list[[i]]$x[which.max(BO_result_list[[i]]$pos)] + ref_val)*365.25)

### Take a look at the fit:
my_alpha_NL <- (BO_result_list[[i]]$x[which.max(BO_result_list[[i]]$pos)])
mod_NL <- fit_once(alpha = my_alpha_NL, data = NL_death)

### Predict before COVID:
fixed_result_1 <- do.call(cbind, mod_NL$design_mat_fixed[which(names(mod_NL$fixed_samp_indexes) %in% c("xx1", "xxx1"))]) %*% t(sample_fixed_effect(mod_NL, variables = c("xx1", "xxx1")))
fixed_result_2 <- do.call(cbind, mod_NL$design_mat_fixed[which(names(mod_NL$fixed_samp_indexes) %in% c("xx2", "xxx2"))]) %*% t(sample_fixed_effect(mod_NL, variables = c("xx2", "xxx2")))
fixed_result <- fixed_result_1 + fixed_result_2
smooth_result <- predict(mod_NL, variable = "x1", newdata = mod_NL$instances[[1]]@data, only.samples = T)
pre_result <- fixed_result + smooth_result[,-1]
pre_sum <- data.frame(mean = pre_result %>% apply(1, mean), upper = pre_result %>% apply(1, quantile, p = 0.975), lower = pre_result %>% apply(1, quantile, p = 0.025))
pre_sum$x <- smooth_result$x
pre_sum <- distinct(pre_sum, x, .keep_all = TRUE)
pre_sum$year <- as.Date((pre_sum$x + ref_val)*365.25)
matplot(x = pre_sum$year, y = exp(pre_sum[,1:3]), col = c(1,2,2), lty = c(1,2,2), type = "l", xlim = as.Date(c("2015-01-01", "2023-01-01")))

### Predict after COVID:
smooth_result_2 <- predict(mod_NL, variable = "x2", newdata = mod_NL$instances[[1]]@data, only.samples = T, include.intercept = FALSE)
smooth_result_2$x <- smooth_result_2$x + my_alpha_NL
post_result <- rbind(fixed_result,fixed_result) + rbind(smooth_result[,-1], smooth_result[,-1]) + rbind(smooth_result_2[,-1], smooth_result_2[,-1])
post_result$x <- c(smooth_result$x, smooth_result_2$x)
post_result <- distinct(post_result, by = x, .keep_all = TRUE)

post_sum <- data.frame(mean = post_result %>% apply(1, mean), upper = post_result %>% apply(1, quantile, p = 0.975), lower = post_result %>% apply(1, quantile, p = 0.025))
post_sum$x <- post_result$x 
post_sum$year <- as.Date((post_sum$x + ref_val)*365.25)
pdf(file = "figures/NL_fit.pdf", height = 5, width = 5)
matplot(x = post_sum$year, y = exp(post_sum[,1:3]),
        col = c(1,2,2), lty = c(1,2,2), type = "l", 
        xlim = as.Date(c("2015-01-01", "2023-01-01")),
        xlab = "", ylab = "Weekly Death",
        main = "", cex.lab = 1.5, cex.axis = 1.5)
points(NL_death$deaths ~ NL_death$date, col = "black", cex = 0.2)
abline(v = as.Date("2020-02-06"), col = "purple")
dev.off()

tikzDevice::tikz("figures/NL_curve.tex", width = 5, height = 5, standAlone = TRUE)
matplot(x = post_sum$year, y = exp(post_sum[,1:3]),
        col = c(1,2,2), lty = c(1,2,2), type = "l", 
        xlim = as.Date(c("2015-01-01", "2023-01-01")),
        xlab = "", ylab = "Weekly Death",
        main = "", cex.lab = 1.5, cex.axis = 1.5, lwd = 2)
points(NL_death$deaths ~ NL_death$date, col = "black", cex = 0.2)
abline(v = as.Date("2020-02-06"), col = "purple")
dev.off()
system("pdflatex -output-directory=./figures ./figures/NL_curve.tex")
file.remove("./figures/NL_curve.tex")
file.remove("./figures/NL_curve.aux")
file.remove("./figures/NL_curve.log")



### East EU: BG
BG_death <- world_death %>% filter(country_name == "Bulgaria")
BG_death$date <- make_date(year = BG_death$year) + weeks(BG_death$time)
BG_death$x <- as.numeric(BG_death$date)/365.25; 
ref_val <- min(BG_death$x)
BG_death$x <- BG_death$x - ref_val
plot(BG_death$deaths ~ BG_death$date)

eval_once <- function(alpha, data = BG_death){
  mod <- fit_once(alpha = alpha, data = data)
  (mod$mod$normalized_posterior$lognormconst)
} 
lower = 0.5
upper = 7.2
objective_func <- eval_once

#### Runtime analysis:
rel_runtime <- c()
BO_result_list <- list()
BO_result_original_list <- list()
for (i in 1:length(eval_num)) {
  eval_number <- eval_num[i]
  begin_time <- Sys.time()
  result_ad <- BO_adap(func = objective_func,
                       # update_step = round(eval_number/2),
                       update_step = 10,
                       number_eval = eval_number,
                       delta = 0.01,
                       lower = lower,
                       upper = upper,
                       length_scale = 0.1,
                       signal_var = 1000,
                       noise_var = noise_var,
                       initial = 5,
                       grid_points = 1000)
  
  data_to_smooth <- result_ad$result
  BO_result_original_list[[i]] <- data_to_smooth
  
  ff <- list()
  # ff$fn <- function(y) as.numeric(surrogate(pnorm(y), data_to_smooth = data_to_smooth) + dnorm(y, log = TRUE))
  ff$fn <- function(x) as.numeric(surrogate(x, data_to_smooth = data_to_smooth))
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
for (i in 1:length(eval_num)) {
  BO_result_list[[i]]$year <- as.Date((BO_result_list[[i]]$x + ref_val)*365.25)
  tikzDevice::tikz(paste0("figures/BG_change_point: B = ", eval_num[i], ".tex"), width = 5, height = 5, standAlone = TRUE)
  ggplot() +
    # geom_histogram(aes(x = mcmc_samps_selected, y = ..density..), bins = 100, alpha = 0.5, fill = "blue") +
    geom_line(data = BO_result_list[[i]], aes(x = year, y = pos), color = "red", size = 1) +
    # geom_line(data = exact_grid_result_smooth, aes(x = x, y = pos), color = "black", size = 1, linetype = "dashed") +
    # ggtitle(paste0("Change-point in BG: B = ", eval_num[i])) +
    xlab("") +
    ylab("Density") +
    # xlim(as.Date(c("2019-12-01", "2021-01-01"))) +
    scale_x_date(
      limits = as.Date(c("2019-12-01", "2021-01-01")),
      date_labels = "%b %y", 
      date_breaks = "3 month"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 15), axis.text = element_text(size = 15))  # only change the lab and axis text size
  # lims(y = range(exact_grid_result_smooth$pos))
  dev.off()
  # Run pdflatex to generate the PDF
  if (file.exists(paste0("figures/BG_change_point: B = ", eval_num[i], ".tex"))) {
    system(sprintf('pdflatex -output-directory=figures "figures/BG_change_point: B = %s.tex"', eval_num[i]))
    
    # Check if the PDF was created
    if (file.exists(paste0("figures/BG_change_point: B = ", eval_num[i], ".pdf"))) {
      # Delete the .tex, .aux, and .log files to clean up
      file.remove(paste0("figures/BG_change_point: B = ", eval_num[i], ".tex"))
      file.remove(paste0("figures/BG_change_point: B = ", eval_num[i], ".aux"))
      file.remove(paste0("figures/BG_change_point: B = ", eval_num[i], ".log"))
    } else {
      warning("PDF file was not created: ", paste0("figures/BG_change_point: B = ", eval_num[i], ".pdf"))
    }
  } else {
    warning("TeX file was not created: ", paste0("figures/BG_change_point: B = ", eval_num[i], ".tex"))
  }
  
  # ggsave(filename = (paste0("figures/BG_change_point: B = ", eval_num[i], " .pdf")))
}
### Which day is most likely: "2020-08-28"
as.Date((BO_result_list[[i]]$x[which.max(BO_result_list[[i]]$pos)] + ref_val)*365.25)

### Take a look at the fit:
my_alpha_BG <- (BO_result_list[[i]]$x[which.max(BO_result_list[[i]]$pos)])
mod_BG <- fit_once(alpha = my_alpha_BG, data = BG_death)

### Predict before COVID:
fixed_result_1 <- do.call(cbind, mod_BG$design_mat_fixed[which(names(mod_BG$fixed_samp_indexes) %in% c("xx1", "xxx1"))]) %*% t(sample_fixed_effect(mod_BG, variables = c("xx1", "xxx1")))
fixed_result_2 <- do.call(cbind, mod_BG$design_mat_fixed[which(names(mod_BG$fixed_samp_indexes) %in% c("xx2", "xxx2"))]) %*% t(sample_fixed_effect(mod_BG, variables = c("xx2", "xxx2")))
fixed_result <- fixed_result_1 + fixed_result_2
smooth_result <- predict(mod_BG, variable = "x1", newdata = mod_BG$instances[[1]]@data, only.samples = T)
pre_result <- fixed_result + smooth_result[,-1]
pre_sum <- data.frame(mean = pre_result %>% apply(1, mean), upper = pre_result %>% apply(1, quantile, p = 0.975), lower = pre_result %>% apply(1, quantile, p = 0.025))
pre_sum$x <- smooth_result$x
pre_sum <- distinct(pre_sum, x, .keep_all = TRUE)
pre_sum$year <- as.Date((pre_sum$x + ref_val)*365.25)
matplot(x = pre_sum$year, y = exp(pre_sum[,1:3]), col = c(1,2,2), lty = c(1,2,2), type = "l", xlim = as.Date(c("2015-01-01", "2023-01-01")))

### Predict after COVID:
smooth_result_2 <- predict(mod_BG, variable = "x2", newdata = mod_BG$instances[[1]]@data, only.samples = T, include.intercept = FALSE)
smooth_result_2$x <- smooth_result_2$x + my_alpha_BG
post_result <- rbind(fixed_result,fixed_result) + rbind(smooth_result[,-1], smooth_result[,-1]) + rbind(smooth_result_2[,-1], smooth_result_2[,-1])
post_result$x <- c(smooth_result$x, smooth_result_2$x)
post_result <- distinct(post_result, by = x, .keep_all = TRUE)

post_sum <- data.frame(mean = post_result %>% apply(1, mean), upper = post_result %>% apply(1, quantile, p = 0.975), lower = post_result %>% apply(1, quantile, p = 0.025))
post_sum$x <- post_result$x 
post_sum$year <- as.Date((post_sum$x + ref_val)*365.25)
pdf(file = "figures/BG_fit.pdf", height = 5, width = 5)
matplot(x = post_sum$year, y = exp(post_sum[,1:3]),
        col = c(1,1,1), lty = c(1,2,2), type = "l", 
        xlim = as.Date(c("2015-01-01", "2023-01-01")),
        xlab = "year", ylab = "All-Cause Death",
        main = "BG")
points(BG_death$deaths ~ BG_death$date, col = "red", cex = 0.2)
abline(v = as.Date("2020-08-28"), col = "purple")
dev.off()

tikzDevice::tikz("figures/BG_curve.tex", width = 5, height = 5, standAlone = TRUE)
matplot(x = post_sum$year, y = exp(post_sum[,1:3]),
        col = c(1,2,2), lty = c(1,2,2), type = "l", 
        xlim = as.Date(c("2015-01-01", "2023-01-01")),
        xlab = "", ylab = "Weekly Death",
        main = "", cex.lab = 1.5, cex.axis = 1.5, lwd = 2)
points(BG_death$deaths ~ BG_death$date, col = "black", cex = 0.2)
abline(v = as.Date("2020-08-28"), col = "purple")
dev.off()
system("pdflatex -output-directory=./figures ./figures/BG_curve.tex")
file.remove("./figures/BG_curve.tex")
file.remove("./figures/BG_curve.aux")
file.remove("./figures/BG_curve.log")

