# Square Exponential covariance function
square_exp_cov_generator <- function(length_scale = 1, signal_var = 1) {
  square_exp_cov <- function(x, x_prime) {
    # Ensure x and x_prime are column vectors
    x <- matrix(x, ncol = 1)
    x_prime <- matrix(x_prime, ncol = 1)
    
    # Compute the pairwise squared distances
    sq_dist <- as.matrix(dist(rbind(x, x_prime)))^2
    
    # Extract the blocks corresponding to x and x_prime
    n <- nrow(x)
    m <- nrow(x_prime)
    sq_dist <- sq_dist[(n + 1):(n + m), 1:n]
    
    # Compute the covariance
    cov_matrix <- signal_var * exp(-sq_dist / (2 * length_scale^2))
    dimnames(cov_matrix) <- NULL
    if (n == 1 || m == 1) {
      dim(cov_matrix) <- c(m, n)
    }
    
    return(cov_matrix)
  }
  return(square_exp_cov)
}

# Square Exponential covariance function's derivative
square_exp_cov_deriv_generator <- function(length_scale = 1, signal_var = 1) {
  square_exp_cov_deriv <- function(x, x_prime) {
    n <- length(x)
    m <- length(x_prime)
    # Compute the pairwise squared distances
    dist <- outer(x, x_prime, FUN = "-") # x - x_prime
    sq_dist <- dist^2
    
    # Compute the covariance
    cov_matrix <- signal_var * exp(-sq_dist / (2 * length_scale^2)) * (dist/(length_scale^2))
    dimnames(cov_matrix) <- NULL
    if (n == 1 || m == 1) {
      dim(cov_matrix) <- c(m, n)
    }
    
    return(cov_matrix)
  }
  return(square_exp_cov_deriv)
}

# Square Exponential covariance function's hessian
square_exp_cov_hess_generator <- function(length_scale = 1, signal_var = 1) {
  square_exp_cov_hess <- function(x, x_prime) {
    n <- length(x)
    m <- length(x_prime)
    # Compute the pairwise squared distances
    dist <- outer(x, x_prime, FUN = "-") # x - x_prime
    sq_dist <- dist^2
    
    # Compute the covariance
    cov_matrix <- signal_var * exp(-sq_dist / (2 * length_scale^2)) * ((dist/(length_scale^2))^2) +
      signal_var * exp(-sq_dist / (2 * length_scale^2)) * (-1/(length_scale^2))
    dimnames(cov_matrix) <- NULL
    if (n == 1 || m == 1) {
      dim(cov_matrix) <- c(m, n)
    }
    
    return(cov_matrix)
  }
  return(square_exp_cov_hess)
}

# Function to compute conditional mean and variance of a GP
predict_gp <- function(data, x_pred, noise_var = 1e-10, choice_cov) {
  # Extract x and y from the data
  x_obs <- data$x
  y_obs <- data$y
  
  # Compute covariance matrices
  K_obs_obs <- choice_cov(x_obs, x_obs)
  K_obs_pred <- choice_cov(x_obs, x_pred)
  K_pred_pred <- choice_cov(x_pred, x_pred)
  
  # Add noise to the diagonal of observed covariance matrix
  # noise_var <- max(abs(min(eigen(K_obs_obs, only.values = T)$values)) + .Machine$double.eps, noise_var)
  K_obs_obs <- K_obs_obs + noise_var * diag(length(x_obs))
  
  # Compute conditional mean and variance
  cond_mean <- K_obs_pred %*% solve(K_obs_obs, y_obs)
  cond_var <- K_pred_pred - K_obs_pred %*% solve(K_obs_obs, t(K_obs_pred))
  
  return(data.frame(x = x_pred, mean = as.vector(cond_mean), var = diag(cond_var)))
}

# Function to compute the derivative (first or second) of the GP's posterior mean function
predict_gp_deriv <- function(data, x_pred, noise_var = 1e-10, choice_cov_deriv, choice_cov) {
  # Extract x and y from the data
  x_obs <- data$x
  y_obs <- data$y
  
  # Compute covariance matrices
  K_obs_obs <- choice_cov(x_obs, x_obs)
  K_obs_pred <- choice_cov_deriv(x_obs, x_pred)

  # Add noise to the diagonal of observed covariance matrix
  K_obs_obs <- K_obs_obs + noise_var * diag(length(x_obs))
  
  # Compute conditional mean and variance
  cond_mean <- K_obs_pred %*% solve(K_obs_obs, y_obs)

  return(data.frame(x = x_pred, mean = as.vector(cond_mean)))
}


# Compute likelihood using Square Exponential:
compute_like <- function(x, y, length_scale, signal_var, noise_var){
  square_exp_cov <- square_exp_cov_generator(length_scale = length_scale, signal_var = signal_var)
  C <- square_exp_cov(x = x, x_prime = x)
  # noise_var <- max(abs(min(eigen(C, only.values = T)$values)) + .Machine$double.eps, noise_var)
  Q <- solve(C + noise_var * diag(length(x)))
  # return(log(det(Q)))
  # return(-log(det(C)))
  like <- as.numeric((-tcrossprod(tcrossprod(x = y, Q), y)/2) + ((log(det(Q)))/2) - (log(2*pi)*length(x)/2))
  # like <- as.numeric((-tcrossprod(tcrossprod(x = y, Q), y)/2) - ((log(det(C)))/2) - (log(2*pi)*length(x)/2))
  like
}

#### Making BO adaptive:
BO_adap <- function(func, update_step = 10, number_eval = 10, lower = 0, upper = 1, noise_var = 1e-6, 
                    length_scale = 1, signal_var = 1, initial = NULL,
                    num_initial = 5,
                    delta = 0.01, grid_points = 100){
  
  # Define the grid
  grid <- seq(from = lower, to = upper, length.out = grid_points)
  grid_transformed <- (grid-lower)/(upper - lower)
  
  # Initialize xvec and yvec to store evaluations
  xvec <- c()
  xvec_trans <- c()
  yvec <- c()
  
  # Set initial point if not provided
  if(is.null(initial)){
    initial = median(grid)
  }
  
  # Consider multiple initial val
  if(num_initial > 1){
    initial <- as.numeric(quantile(grid, p = seq(0,1, length.out = num_initial)))
  }
  
  for (initial_val in initial) {
    # Add the initial point to xvec and evaluate the function
    xvec_trans <- c(xvec_trans, (initial_val - lower)/(upper - lower))
    xvec <- c(xvec, initial_val)
    yvec <- c(yvec, func(initial_val))
  }
  # Assign the reference value
  rel <- max(yvec)
  yvec <- yvec - rel
  
  # Remove the initial point from the grid
  grid <- grid[!grid %in% xvec]
  grid_transformed <- (grid-lower)/(upper - lower)
  
  choice_cov <- square_exp_cov_generator(length_scale = length_scale, signal_var = signal_var)
  
  # Perform Bayesian Optimization
  for (i in 1:number_eval) {
    print(paste("Iteration:", i))
    newdata <- data.frame(x = xvec_trans, x_original = xvec, y = yvec)
    
    if(i %% update_step == 0){
      print(paste("Time to update the parameters!"))
      signal_var = var(newdata$y)
      length_scale_vec <- seq(from = 0.01, to = 1, by = 0.01)
      like <- c()
      for (j in 1:length(length_scale_vec)) {
        like[j] <- compute_like(y = newdata$y, x = newdata$x, signal_var = signal_var, length_scale = length_scale_vec[j], noise_var = noise_var)
      }
      length_scale <- length_scale_vec[which.max(like)]
      print(paste("The new length.scale:", length_scale))
      print(paste("The new signal_var:", signal_var))
      choice_cov <- square_exp_cov_generator(length_scale = length_scale, signal_var = signal_var)
    }
    
    # Update the grid by removing already evaluated points
    grid <- grid[!grid %in% newdata$x_original]
    grid_transformed <- (grid-lower)/(upper - lower)
    
    # Update the GP
    fnew <- predict_gp(newdata, grid_transformed, choice_cov = choice_cov, noise_var = noise_var)
    
    # Compute the UCB acquisition function
    beta <- 2*log(((i + num_initial)^2)*(pi^2)/(6*delta))
    UCB <- fnew$mean + sqrt(beta) * sqrt(fnew$var)
    
    # Select the next point to evaluate
    next_point <- grid[which.max(UCB)]
    xvec <- c(xvec, next_point)
    xvec_trans <- (xvec - lower)/(upper - lower)
    yvec <- c(yvec, (func(next_point) - rel))
    
    # Print some information for debugging
    print(paste("Next point:", next_point))
    print(paste("Function value:", yvec[i + 1]))
  }
  
  # Return the result
  return(list(result = data.frame(x = xvec_trans, x_original = xvec, y = yvec),
              length_scale = length_scale, signal_var = signal_var))
}


