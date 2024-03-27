library(BayesGP)
set.seed(123)
source("01_simulate.R")


eval_once <- function(alpha){
  a_fit <- alpha
  data$x1 <- ifelse(data$x <= a_fit, data$x, a_fit)
  data$x2 <- ifelse(data$x > a_fit, (data$x - a_fit), 0)
  mod <- model_fit(formula = y ~ f(x1, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0) + f(x2, model = "IWP", order = 2, sd.prior = list(param = 1, h = 1), initial_location = 0), 
                   data = data, method = "aghq", family = "Gaussian", aghq_k = 3
  )
  (mod$mod$normalized_posterior$lognormconst)
}

### MCMC with symmetric (RW) proposal:
propose <- function(ti, step_size = 0.1){
  ti + rnorm(n = 1, sd = step_size)
}

## prior: uniform prior over [0,10]
prior <- function(t){
  ifelse(t <= 0 | t >= 10, 0, 0.1)
}

### Acceptance rate:
acceptance_Rate <- function(ti1, ti2){
  ## To make the algorithm numerically stable.
  if(ti2 >= 9.9 | ti2 <= 0.1){
    0
  }
  else{
    min(1, exp(log(prior(ti2)+.Machine$double.eps) + eval_once(ti2) - log(prior(ti1)+.Machine$double.eps) - eval_once(ti1)))
  }
}

### Run MCMC:
run_mcmc <- function(ti0 = 2, M = 10, size = 0.3){
  samps <- numeric(M)
  ti1 <- ti0
  for (i in 1:M) {
    ti2 <- propose(ti = ti1, step_size = size)
    ui <- runif(1)
    Rate <- acceptance_Rate(ti1, ti2)
    if(ui <= acceptance_Rate(ti1, ti2)){
      cat(paste0(ti2, " is accepted by ", ti1, " at iteration ", i, "\n"))
      ti1 <- ti2
    }
    else{
      cat(paste0(ti2, " is rejected by ", ti1, " at iteration ", i, "\n"))
    }
    samps[i] <- ti1
  }
  samps
}
begin_runtime <- Sys.time()
mcmc_samps <- run_mcmc(ti0 = 2, M = 10000, size = 0.5)
end_runtime <- Sys.time()
end_runtime - begin_runtime

save(mcmc_samps, file = "mcmc_samps.rda")

plot(mcmc_samps)
burnin <- 500
hist(mcmc_samps[-c(1:burnin)], breaks = 30)
thinning <- 3
mcmc_samps_selected <- mcmc_samps[-c(1:burnin)][seq(1, length(mcmc_samps[-c(1:burnin)]), by=thinning)]
hist(mcmc_samps_selected, breaks = 30)

dens <- density(mcmc_samps_selected, from = 0, to = 10)
plot(dens, xlim = c(0, 10), main = "Kernel Density Estimation")

### For run-time:
begin_runtime <- Sys.time()
mcmc_samps <- run_mcmc(ti0 = 2, M = 100, size = 0.5)
end_runtime <- Sys.time()
end_runtime - begin_runtime
