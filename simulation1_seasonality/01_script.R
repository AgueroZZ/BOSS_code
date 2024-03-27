library(BayesGP)
library(tidyverse)

lower = 0.5
upper = 4.5

### Simulate data:
set.seed(123)
n <- 100
x <- runif(n = n, min = 0, max = 10)
true_func <- function(x, true_alpha = 1.5){1 + 
    0.5 * cos((2*pi*x)/true_alpha) - 1.3 * sin((2*pi*x)/true_alpha) +
    1.1 * cos((4*pi*x)/true_alpha) + 0.3 * sin((4*pi*x)/true_alpha)}
log_mu <- true_func(x) + rnorm(n, sd = 2)
y <- rpois(n = n, lambda = exp(log_mu))
data <- data.frame(y = y, x = x, indx = 1:n, log_mu = log_mu)
plot(log_mu ~ x, type = "o", data = arrange(data, x))
lines(true_func(x) ~ x, data = arrange(data, x), lty = "dashed", col = "red")
plot(y ~ x, type = "p", data = arrange(data, x))
save(data, file = "data.rda")

log_prior <- function(alpha){
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


### MCMC with symmetric (RW) proposal:
propose <- function(ti, step_size = 0.1){
  ti + rnorm(n = 1, sd = step_size)
}


### Acceptance rate:
acceptance_Rate <- function(ti1, ti2){
  ## To make the algorithm numerically stable.
  if(ti2 >= 5 | ti2 <= 0.1){
    0
  }
  else{
    min(1, exp(log_prior(ti2) + eval_once(ti2) - log_prior(ti1) - eval_once(ti1)))
  }
}

### Run MCMC:
run_mcmc <- function(ti0 = 2, M = 10, size = 0.1){
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
# begin_runtime <- Sys.time()
# mcmc_samps <- run_mcmc(ti0 = 2, M = 10000, size = 0.5)
# end_runtime <- Sys.time()
# end_runtime - begin_runtime
# save(mcmc_samps, file = "mcmc_samps.rda")
# 
# plot(mcmc_samps)
# burnin <- 500
# hist(mcmc_samps[-c(1:burnin)], breaks = 30)
# thinning <- 3
# mcmc_samps_selected <- mcmc_samps[-c(1:burnin)][seq(1, length(mcmc_samps[-c(1:burnin)]), by=thinning)]
# hist(mcmc_samps_selected, breaks = 30)
# 
# dens <- density(mcmc_samps_selected, from = 0.1, to = 5)
# plot(dens, xlim = c(0, 10), main = "Kernel Density Estimation")



### For run-time:
begin_runtime <- Sys.time()
mcmc_samps <- run_mcmc(ti0 = 2, M = 100, size = 0.1)
end_runtime <- Sys.time()
end_runtime - begin_runtime
# 41.87028 secs
# 0.697838 mins



