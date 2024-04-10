library(BayesGP)
library(npreg)
library(ggplot2)
library(INLA)
source_location <- "function/"
source(file = paste0(source_location, "bo_functions.R"))
load("results/data.rda")
load("results/quad_exact.rda")
load("results/quad_BO.rda")

g1 <- function(x){x*log((x^2) + 1)}
g2 <- function(x){3.3*x + 3.035}


### Function that simulate from Quad, then make inference of the function
fit_once <- function(alpha, data){
  a_fit <- alpha
  data$x1 <- ifelse(data$x <= a_fit, data$x, a_fit)
  data$x2 <- ifelse(data$x > a_fit, (data$x - a_fit), 0)
  mod <- model_fit(formula = y ~ f(x1, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0) + f(x2, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0), 
                   data = data, method = "aghq", family = "Gaussian", aghq_k = 3
  )
  mod
}
sim_quad <- function(n, quad){
  prob_vec = quad$weights * exp(quad$logpost_normalized)
  freq <- as.numeric(rmultinom(n = 1, size = n, prob = prob_vec))
  samp_theta <- rep(quad$theta1, times = freq)
  samp_alpha <- pnorm(samp_theta)*10
  samp_alpha
}
fit_all_mod <- function(alpha_samps){
  alpha_samps_table <- table(alpha_samps)
  result_list <- list(alpha = as.numeric(names(alpha_samps_table)), mod = list(), M = as.numeric(alpha_samps_table))
  for (i in 1:length(alpha_samps_table)) {
    M <- as.numeric(alpha_samps_table[i])
    alpha <- as.numeric(names(alpha_samps_table)[i])
    mod <- fit_once(alpha = alpha, data = data)
    result_list$mod[[i]] <- mod
  }
  result_list
}
infer_g1 <- function(xvec, all_mod){
  all_samps <- matrix(nrow = length(xvec), ncol = 0)
  for (i in 1:length(all_mod$M)) {
    mod <- all_mod$mod[[i]]
    samples_g1 <- predict(mod, variable = "x1",
                          newdata = data.frame(x1 = xvec), 
                          only.sample = T)[,-1]
    all_samps <- cbind(all_samps, samples_g1[,(1:all_mod$M[i])])
  }
  all_samps
}
infer_g2 <- function(xvec, all_mod){
  all_samps <- matrix(nrow = length(xvec), ncol = 0)
  for (i in 1:length(all_mod$M)) {
    mod <- all_mod$mod[[i]]
    samples_g2 <- predict(mod, variable = "x2",
                          newdata = data.frame(x2 = xvec), 
                          only.sample = T)[,-1]
    all_samps <- cbind(all_samps, samples_g2[,(1:all_mod$M[i])])
  }
  all_samps
}

### Using exact grid:
set.seed(123)
exact_all_mod <- fit_all_mod(sim_quad(n = 3000, quad = quad_exact))
exact_samples_g1 <- infer_g1(xvec = seq(0, 10, by = 0.01), all_mod = exact_all_mod)
exact_samples_g2 <- infer_g2(xvec = seq(0, 10, by = 0.01), all_mod = exact_all_mod)
exact_samples_g1_sum <- BayesGP::extract_mean_interval_given_samps(samps = cbind(seq(0, 10, by = 0.01), exact_samples_g1))
exact_samples_g2_sum <- BayesGP::extract_mean_interval_given_samps(samps = cbind(seq(0, 10, by = 0.01), exact_samples_g2))

## Plot g1:
png(filename = "figures/grid_g1.png", height = 500, width = 500)
plot(q0.5~x, data = exact_samples_g1_sum, type = "l", col = "blue", ylab = "g1")
lines(q0.975~x, data = exact_samples_g1_sum, lty = "dashed", col = "red")
lines(q0.025~x, data = exact_samples_g1_sum, lty = "dashed", col = "red")
lines(g1(x)~x, exact_samples_g1_sum, col = "black")
dev.off()

## Plot g2: (initialized)
png(filename = "figures/grid_g2.png", height = 500, width = 500)
plot(q0.5~x, data = exact_samples_g2_sum, type = "l", col = "blue", ylab = "g2")
lines(q0.975~x, data = exact_samples_g2_sum, lty = "dashed", col = "red")
lines(q0.025~x, data = exact_samples_g2_sum, lty = "dashed", col = "red")
lines(I(g2(x)-g2(0))~x, exact_samples_g1_sum, col = "black")
dev.off()

### Using BOSS approximation:
set.seed(123)
BO_all_mod <- fit_all_mod(sim_quad(n = 3000, quad = quad_BO))
BO_samples_g1 <- infer_g1(xvec = seq(0, 10, by = 0.01), all_mod = BO_all_mod)
BO_samples_g2 <- infer_g2(xvec = seq(0, 10, by = 0.01), all_mod = BO_all_mod)
BO_samples_g1_sum <- BayesGP::extract_mean_interval_given_samps(samps = cbind(seq(0, 10, by = 0.01), BO_samples_g1))
BO_samples_g2_sum <- BayesGP::extract_mean_interval_given_samps(samps = cbind(seq(0, 10, by = 0.01), BO_samples_g2))

## Plot g1:
png(filename = "figures/BO_g1.png", height = 500, width = 500)
plot(q0.5~x, data = BO_samples_g1_sum, type = "l", col = "blue", ylab = "g1")
lines(q0.975~x, data = BO_samples_g1_sum, lty = "dashed", col = "red")
lines(q0.025~x, data = BO_samples_g1_sum, lty = "dashed", col = "red")
lines(g1(x)~x, BO_samples_g1_sum, col = "black")
dev.off()


## Plot g2: (initialized)
png(filename = "figures/BO_g2.png", height = 500, width = 500)
plot(q0.5~x, data = BO_samples_g2_sum, type = "l", col = "blue", ylab = "g2")
lines(q0.975~x, data = BO_samples_g2_sum, lty = "dashed", col = "red")
lines(q0.025~x, data = BO_samples_g2_sum, lty = "dashed", col = "red")
lines(I(g2(x)-g2(0))~x, BO_samples_g1_sum, col = "black")
dev.off()


### Comparison in one figure:
exact_samples_color <- rgb(1, 0, 0, alpha = 0.2) # Red with transparency
BO_samples_color <- rgb(0, 0, 1, alpha = 0.2) # Blue with transparency
mar.default <- c(5,4,4,2)
par(mar = mar.default + c(0, 1, 0, 0))
png(filename = "figures/compare_g1.png", height = 500, width = 500)
plot(q0.5~x, data = exact_samples_g1_sum, type = "l", col = "red", ylab = "", lty = "dashed", cex.lab = 2.0, cex.axis = 2.0)
polygon(c(exact_samples_g1_sum$x, rev(exact_samples_g1_sum$x)),
        c(exact_samples_g1_sum$q0.025, rev(exact_samples_g1_sum$q0.975)),
        col = exact_samples_color, border = NA)
lines(g1(x)~x, exact_samples_g1_sum, col = "black")
polygon(c(BO_samples_g1_sum$x, rev(BO_samples_g1_sum$x)),
        c(BO_samples_g1_sum$q0.025, rev(BO_samples_g1_sum$q0.975)),
        col = BO_samples_color, border = NA)
lines(q0.5~x, data = BO_samples_g1_sum, col = "blue", lty = "dashed")
dev.off()

# pdf output through tikzDevice
tikzDevice::tikz(file = paste0("figures/compare_g1.tex"), 
                 width = 5, height = 5, standAlone = TRUE)
plot(q0.5~x, data = exact_samples_g1_sum, type = "l", col = "red", ylab = "", lty = "dashed", cex.lab = 2.0, cex.axis = 2.0)
polygon(c(exact_samples_g1_sum$x, rev(exact_samples_g1_sum$x)),
        c(exact_samples_g1_sum$q0.025, rev(exact_samples_g1_sum$q0.975)),
        col = exact_samples_color, border = NA)
lines(g1(x)~x, exact_samples_g1_sum, col = "black")
polygon(c(BO_samples_g1_sum$x, rev(BO_samples_g1_sum$x)),
        c(BO_samples_g1_sum$q0.025, rev(BO_samples_g1_sum$q0.975)),
        col = BO_samples_color, border = NA)
lines(q0.5~x, data = BO_samples_g1_sum, col = "blue", lty = "dashed")
dev.off()
system("pdflatex -output-directory=./figures ./figures/compare_g1.tex")
file.remove("./figures/compare_g1.tex")
file.remove("./figures/compare_g1.aux")
file.remove("./figures/compare_g1.log")



png(filename = "figures/compare_g2.png", height = 500, width = 500)
plot(q0.5~x, data = exact_samples_g2_sum, type = "l", col = "red", ylab = "", lty = "dashed", cex.lab = 2.0, cex.axis = 2.0)
polygon(c(exact_samples_g2_sum$x, rev(exact_samples_g2_sum$x)),
        c(exact_samples_g2_sum$q0.025, rev(exact_samples_g2_sum$q0.975)),
        col = exact_samples_color, border = NA)
lines(I(g2(x)-g2(0))~x, exact_samples_g1_sum, col = "black")
polygon(c(BO_samples_g2_sum$x, rev(BO_samples_g2_sum$x)),
        c(BO_samples_g2_sum$q0.025, rev(BO_samples_g2_sum$q0.975)),
        col = BO_samples_color, border = NA)
lines(q0.5~x, data = BO_samples_g2_sum, col = "blue", lty = "dashed")
dev.off()

tikzDevice::tikz(file = paste0("figures/compare_g2.tex"), 
                 width = 5, height = 5, standAlone = TRUE)
plot(q0.5~x, data = exact_samples_g2_sum, type = "l", col = "red", ylab = "", lty = "dashed", cex.lab = 2.0, cex.axis = 2.0)
polygon(c(exact_samples_g2_sum$x, rev(exact_samples_g2_sum$x)),
        c(exact_samples_g2_sum$q0.025, rev(exact_samples_g2_sum$q0.975)),
        col = exact_samples_color, border = NA)
lines(I(g2(x)-g2(0))~x, exact_samples_g1_sum, col = "black")
polygon(c(BO_samples_g2_sum$x, rev(BO_samples_g2_sum$x)),
        c(BO_samples_g2_sum$q0.025, rev(BO_samples_g2_sum$q0.975)),
        col = BO_samples_color, border = NA)
lines(q0.5~x, data = BO_samples_g2_sum, col = "blue", lty = "dashed")
dev.off()
system("pdflatex -output-directory=./figures ./figures/compare_g2.tex")
file.remove("./figures/compare_g2.tex")
file.remove("./figures/compare_g2.aux")
file.remove("./figures/compare_g2.log")

