library(tidyverse)
library(tikzDevice)
library(rstan)
library(INLA)
library(inlabru)

# load BOSS algorithm functions
source('source/bo_functions_nd.R')
# simulate the data
r <- seq(0, 20, by = 0.1)
beta <- 10
a <- 2
b <- 2
c <- -2.5

set.seed(1234)
Ir <- beta*(1 + (r/a)^b)^c
lr <- log(Ir) + rnorm(length(r), 0, 0.5)

plot(r, lr)

data <- data.frame(r, lr)

# fit the non-linear regression model through inlabru
a_fun <- function(u){
  qunif(pnorm(u), 0.1, 5)
}

b_fun <- function(u){
  qunif(pnorm(u), 0.1, 4)
}

cmp <- ~ a(1, model="linear", mean.linear=0, prec.linear=1) +
  b(1, model="linear", mean.linear=0, prec.linear=1) + 
  c(1) + Intercept(1)

form <- lr ~ Intercept + c*log(1 + (r/a_fun(a))^b_fun(b))

system.time(
fit <- bru(cmp, formula = form, data = data, family = 'gaussian')
)

# save results
saveRDS(fit, 'inlabru_nlreg.RDS')

# get joint posterior of (R, beta) from inlabru
joint_samp <- inla.posterior.sample(10000, fit, selection = list(a = 1, b = 1), seed = 12345)
joint_samp <- do.call('rbind', lapply(joint_samp, function(x) matrix(x$latent, ncol = 2)))

inla.joint.samps <- data.frame(a = a_fun(joint_samp[,1]), b = b_fun(joint_samp[,2]))

# plot joint posterior of (R, beta) from inlabru
tikz(file = "joint_post_R_beta_inlabru.tex", standAlone=T, width = 5, height = 4)
ggplot(inla.joint.samps, aes(a, b)) + stat_density_2d(
  geom = "raster",
  aes(fill = after_stat(density)), n = 500,
  contour = FALSE) +
  geom_point(data = data.frame(a = a_fun(fit$summary.fixed$mode[1]), b = b_fun(fit$summary.fixed$mode[2])), color = 'red', shape = 1, size =0.5) + 
  geom_point(data = data.frame(a = 2, b = 2), color = 'orange', size =0.5) +
  coord_fixed() + scale_fill_viridis_c(name = 'Density') + theme_minimal() + xlab('$R$') + ylab('$\\beta$') + xlim(c(0.1, 5)) + ylim(c(0.1, 4))
dev.off()
system('pdflatex joint_post_R_beta_inlabru.tex')

# run MCMC to fit the non-linear regression model
set.seed(1234)
system.time(
MCMC_fit <- stan(
  file = "nlreg.stan",  # Stan program
  data = list(x = r, y = lr, N = length(r)),    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 100000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  algorithm = 'NUTS'             # no progress shown
))

# thin the samples
MCMC_samp <- as.data.frame(MCMC_fit)
MCMC_samp_thin <- MCMC_samp[seq(1, 396000, by = 8),]

# specify the objective function for BOSS: unnormalized log posterior of (R, beta)
eval_func <- function(par, x = r, y = lr){
  a <- par[1]
  b <- par[2]
  sim_data <- data.frame(r = log(1 + (r/a)^b), Ir = y)
  fit <- inla(formula = Ir ~ r, data = sim_data, family = 'gaussian')
  
  unname(fit$mlik[2,])
}

# run BOSS
set.seed(12345)
system.time(
res_opt <- BO_adap_optim(eval_func, update_step = 5, number_eval = 100, D = 2, 
                   lower = c(0.1, 0.1), upper = c(5, 4), length_scale = 0.2, signal_var = 2000, num_initial = 10, delta = 0.0001)
)

# save BOSS results. 
saveRDS(res_opt, 'BOSS_nlreg.RDS')
# NOTE: Since INLA executes parallel processing automatically, setting seed cannot guarantee exactly the same results.
# To obtain the results and figures in the paper, please directly read in the following saved result from BOSS.
res_opt <- readRDS('BOSS_nlreg.RDS')

# construct covariance function for GP regression
square_exp_cov <- square_exp_cov_generator_nd(length_scale = res_opt$length_scale, signal_var = res_opt$signal_var)

# get BOSS design points to fit GP for surrogate log posterior density 
data_to_smooth <- list()
unique_data <- unique(res_opt$result[,c('x.1','x.2','y')])
data_to_smooth$x <- as.matrix(unique_data[,c('x.1','x.2')])
data_to_smooth$y <- unique_data$y - mean(unique_data$y)

surrogate <- function(xvalue, data_to_smooth, cov){
  predict_gp(data_to_smooth, x_pred = xvalue, choice_cov = cov, noise_var = 5e-10)$mean
}

# grid method to normalize posterior density for plotting joint posterior of (R, beta)
ff <- list()
ff$fn <- function(x) as.numeric(surrogate(x, data_to_smooth, square_exp_cov))
x.1 <- (seq(from = 0.1, to = 5, length.out = 500) - 0.1)/4.9
x.2 <- (seq(from = 0.1, to = 4, length.out = 500) - 0.1)/3.9
x_vals <- expand.grid(x.1, x.2)
names(x_vals) <- c('x.1','x.2')
x_original <- t(t(x_vals)*(c(5, 4) - c(0.1, 0.1)) + c(0.1, 0.1)) 

fn_vals <- apply(x_vals, 1, function(x) ff$fn(matrix(x, ncol = 2))) + mean(res_opt$result$y)

# normalize
lognormal_const <- log(sum(exp(fn_vals))*0.0098*0.0078)
post_x <- data.frame(x_original, pos = exp(fn_vals - lognormal_const))

# plot joint posterior of (R, beta) from BOSS
tikz(file = "joint_post_R_beta.tex", standAlone=T, width = 5, height = 4)
ggplot(post_x, aes(x.1,x.2)) + geom_raster(aes(fill = pos)) + 
  geom_point(data = data.frame(x.1 = post_x$x.1[which.max(post_x$pos)], x.2 = post_x$x.2[which.max(post_x$pos)]), color = 'red', shape = 1, size =0.5) +
  geom_point(data = data.frame(x.1 = 2, x.2 = 2), color = 'orange', size =0.5) + coord_fixed() + scale_fill_viridis_c(name = 'Density') + theme_minimal() + xlab('$R$') + ylab('$\\beta$')
dev.off()
system('pdflatex joint_post_R_beta.tex')

# get posterior samples of (R, beta) from BOSS
dx <- unique(post_x$x.1)[2] - unique(post_x$x.1)[1]
dy <- unique(post_x$x.2)[2] - unique(post_x$x.2)[1]
set.seed(123456)
sample_idx <- rmultinom(1:250000, size = 49500, prob = post_x$pos)
sample_x <- data.frame(post_x, n = sample_idx)

samples <- data.frame(do.call(rbind, apply(sample_x, 1, function(x) cbind(runif(x[4], x[1], x[1]+dx), runif(x[4], x[2], x[2] + dy)))))

# plot joint posterior of (R, beta) from MCMC
tikz(file = "joint_post_R_beta_MCMC.tex", standAlone=T, width = 5, height = 4)
ggplot(MCMC_samp_thin, aes(a, b)) + stat_density_2d(
  geom = "raster",
  aes(fill = after_stat(density)), n = 500,
  contour = FALSE) + 
  geom_point(data = data.frame(a = post_x$x.1[which.max(post_x$pos)], b = post_x$x.2[which.max(post_x$pos)]), color = 'red', shape = 1, size =0.5) +
  geom_point(data = data.frame(a = 2, b = 2), color = 'orange', size =0.5) + coord_fixed() + scale_fill_viridis_c(name = 'Density') + theme_minimal() + xlab('$R$') + ylab('$\\beta$') + xlim(c(0.1, 5)) + ylim(c(0.1, 4))
dev.off()
system('pdflatex joint_post_R_beta_MCMC.tex')

# get posterior samples of (R, beta) from inlabru
set.seed(1234)
inla.samples.a <- a_fun(inla.rmarginal(49500, fit$marginals.fixed$a))
inla.samples.b <- b_fun(inla.rmarginal(49500, fit$marginals.fixed$b))

# get marginal of R 
R_marginal <- data.frame(R = c(inla.samples.a, samples[,1], MCMC_samp_thin$a), 
                         method = rep(c('inlabru', 'BOSS (ours)', 'MCMC'),  c(length(inla.samples.a), length(samples[,1]), length(MCMC_samp_thin$a))))

# plot marginal of R
tikz(file = "R_marginal.tex", standAlone=T, width = 4, height = 2)
ggplot(R_marginal, aes(R)) + geom_density(aes(color = method)) + theme_minimal() + xlab('$R$')
dev.off()
system('pdflatex R_marginal.tex')

# get marginal of beta
beta_marginal <- data.frame(beta = c(inla.samples.b, samples[,2], MCMC_samp_thin$b), 
                         method = rep(c('inlabru', 'BOSS (ours)', 'MCMC'),  c(length(inla.samples.b), length(samples[,2]), length(MCMC_samp_thin$b))))

# plot marginal of R
tikz(file = "beta_marginal.tex", standAlone=T, width = 4, height = 2)
ggplot(beta_marginal, aes(beta)) + geom_density(aes(color = method)) + theme_minimal() + xlab('$\\beta$')
dev.off()
system('pdflatex beta_marginal.tex')

