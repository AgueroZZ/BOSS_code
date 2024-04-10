# This function computes the pairwise squared distances for matrices
compute_sq_dist <- function(X, Y) {
  if(is.vector(X)) X <- matrix(X, ncol = 1)
  if(is.vector(Y)) Y <- matrix(Y, ncol = 1)
  XX <- rowSums(X^2)
  YY <- rowSums(Y^2)
  XY <- tcrossprod(X, Y)
  sq_dist <- matrix(XX, ncol=nrow(Y), nrow=nrow(X)) + 
    t(matrix(YY, ncol=nrow(X), nrow=nrow(Y))) - 
    2 * XY
  return(sq_dist)
}

# Square Exponential covariance function
square_exp_cov_generator_nd <- function(length_scale = 1, signal_var = 1) {
  square_exp_cov <- function(x, x_prime) {
    
    # Check if x and x_prime are the same, if so only compute half the matrix
    same_input <- identical(x, x_prime)
    
    # Compute the pairwise squared distances for multivariate data
    sq_dist <- compute_sq_dist(x, x_prime)
    
    # If inputs are the same, make the matrix symmetric
    if (same_input) {
      upper_tri <- upper.tri(sq_dist)
      sq_dist[upper_tri] <- t(sq_dist)[upper_tri]
    }
    
    # Compute the covariance
    cov_matrix <- signal_var * exp(-sq_dist / (2 * length_scale^2))
    
    return(t(cov_matrix))
  }
  
  return(square_exp_cov)
}

# Square Exponential covariance function's derivative
square_exp_cov_deriv_generator_nd <- function(length_scale = 1, signal_var = 1) {
  square_exp_cov_deriv_nd <- function(x, x_prime) {
    if(is.vector(x)) x <- matrix(x, ncol = 1)
    if(is.vector(x_prime)) x_prime <- matrix(x_prime, ncol = 1)

    n <- nrow(x)
    m <- nrow(x_prime)
    # Compute the pairwise squared distances
    dist <- vector('list', ncol(x)) # x - x_prime
    for(i in 1:ncol(x)){
      dist[[i]] <- outer(x[,i], x_prime[,i], '-')
    }
    names(dist) <- apply(rbind(1:ncol(x)), 1, function(n) paste0('dim', n))
    sq_dist <- t(compute_sq_dist(x, x_prime))

    # Compute the covariance
    cov_matrix <- lapply(dist, function(z) signal_var * exp(-sq_dist / (2 * length_scale^2)) * (z/(length_scale^2)))

    return(cov_matrix)
  }
  return(square_exp_cov_deriv_nd)
}
# 


# Function to compute conditional mean and variance of a GP
predict_gp <- function(data, x_pred, noise_var = 1e-10, choice_cov) {
  # Extract x and y from the data
  x_obs <- data$x
  y_obs <- data$y
  if(is.vector(x_obs)){
    N <- length(x_obs)
    }
  else{
    N <- nrow(x_obs)
  }
  
  # Compute covariance matrices
  K_obs_obs <- choice_cov(x_obs, x_obs)
  K_obs_pred <- choice_cov(x_obs, x_pred)
  K_pred_pred <- choice_cov(x_pred, x_pred)
  
  # Add noise to the diagonal of observed covariance matrix
  # noise_var <- max(abs(min(eigen(K_obs_obs, only.values = T)$values)) + .Machine$double.eps, noise_var)
  K_obs_obs <- K_obs_obs + noise_var * diag(N)
  
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
  square_exp_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
  C <- square_exp_cov(x = x, x_prime = x)
  # noise_var <- max(abs(min(eigen(C, only.values = T)$values)) + .Machine$double.eps, noise_var)
  Q <- solve(C + noise_var * diag(nrow(x)))
  # return(log(det(Q)))
  # return(-log(det(C)))
  like <- as.numeric((-tcrossprod(tcrossprod(x = y, Q), y)/2) + ((log(det(Q)))/2) - (log(2*pi)*length(x)/2))
  # like <- as.numeric((-tcrossprod(tcrossprod(x = y, Q), y)/2) - ((log(det(C)))/2) - (log(2*pi)*length(x)/2))
  like
}

#### Making BO adaptive:
BO_adap <- function(func, update_step = 10, number_eval = 10, D = 1, lower = rep(0, D), upper = rep(1, D), noise_var = 1e-6, 
                    length_scale = 1, signal_var = 1, initial = NULL, 
                    num_initial = 5, delta = 0.01, grid_points = 100){
  
  # Define the grid
  
  # Check if dimensions of lower_bounds and upper_bounds match
  if(length(lower) != D | length(upper) != D) {
    stop("lower_bounds and upper_bounds must have the same length as function input dimension")
  }
  
  # Generate points
  grid <- matrix(runif(grid_points * D, lower, upper),
                   nrow = grid_points, ncol = D, byrow = T)
  grid_transformed <- t((t(grid) - lower)/(upper - lower))
  
  # Initialize xvec and yvec to store evaluations
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()
  
  # Set initial point if not provided
  if(is.null(initial)){
    initial = colMeans(grid)
  }
  
  # Consider multiple initial val
  if(num_initial > 1){
    initial <- matrix(unname(apply(grid, 2, function(x) quantile(x, p = seq(0.01,0.99, length.out = 5)))), ncol = D)
  }
  
  for (i in 1:nrow(initial)) {
    # Add the initial point to xvec and evaluate the function
    xmat_trans <- rbind(xmat_trans, (initial[i,] - lower)/(upper - lower))
    xmat <- rbind(xmat, initial[i,])
    yvec <- c(yvec, func(initial[i,]))
  }
  # Assign the reference value
  rel <- max(yvec)
  yvec <- yvec - rel
  
  # Remove the initial point from the grid
  grid <- matrix(grid[is.na(match(do.call(paste, data.frame(grid)), do.call(paste, data.frame(xmat)))),], ncol = D)
  grid_transformed <- t((t(grid)-lower))/(upper - lower)
  
  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
  
  # Perform Bayesian Optimization
  for (i in 1:number_eval) {
    print(paste("Iteration:", i))
    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
    
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
      choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
    }
    
    # Update the grid by removing already evaluated points
    grid <- matrix(grid[is.na(match(do.call(paste, data.frame(grid)), do.call(paste, data.frame(newdata$x_original)))),], ncol = D)
    grid_transformed <- (grid-lower)/(upper - lower)
    
    # Update the GP
    fnew <- predict_gp(newdata, grid_transformed, choice_cov = choice_cov, noise_var = noise_var)
    
    # Compute the UCB acquisition function
    beta <- 2*log(((i + num_initial)^2)*(pi^2)/(6*delta))
    UCB <- fnew$mean + sqrt(beta) * sqrt(fnew$var)
    
    # Select the next point to evaluate
    next_point <- grid[which.max(UCB),]
    xmat <- rbind(xmat, next_point)
    xmat_trans <- t((t(xmat) - lower)/(upper - lower))
    yvec <- c(yvec, (func(next_point) - rel))
    
    # Print some information for debugging
    print(paste("Next point:", next_point))
    print(paste("Function value:", yvec[i + 1]))
  }
  
  # Return the result
  return(list(result = data.frame(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var))
}


UCB <- function(x, data, cov, nv, D, d){
  fnew <- predict_gp(data, x, choice_cov = cov, noise_var = nv)
  
  # Compute the UCB acquisition function
  beta <- 2*log((D^2)*(pi^2)/(6*d))
  return(as.numeric(-fnew$mean - sqrt(beta) * sqrt(fnew$var)))
}

#### Making BO adaptive:
BO_adap_optim <- function(func, update_step = 10, number_eval = 10, D = 1, lower = rep(0, D), upper = rep(1, D), noise_var = 1e-6, 
                    length_scale = 1, signal_var = 1, initial = NULL, 
                    num_initial = 5, delta = 0.01){
  
  # Check if dimensions of lower_bounds and upper_bounds match
  if(length(lower) != D | length(upper) != D) {
    stop("lower_bounds and upper_bounds must have the same length as function input dimension")
  }
  
  # Initialize xvec and yvec to store evaluations
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()
  
  initial <- matrix(runif(num_initial * D, lower, upper),
                    nrow = num_initial, ncol = D, byrow = T)

  for (i in 1:nrow(initial)) {
    # Add the initial point to xvec and evaluate the function
    xmat_trans <- rbind(xmat_trans, (initial[i,] - lower)/(upper - lower))
    xmat <- rbind(xmat, initial[i,])
    yvec <- c(yvec, func(initial[i,]))
  }
  # Assign the reference value
  rel <- max(yvec)
  yvec <- yvec - rel
  
  # Remove the initial point from the grid
  # grid <- matrix(grid[is.na(match(do.call(paste, data.frame(grid)), do.call(paste, data.frame(xmat)))),], ncol = D)
  # grid_transformed <- t((t(grid)-lower))/(upper - lower)
  
  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
  
  # Perform Bayesian Optimization
  for (i in 1:number_eval) {
    print(paste("Iteration:", i))
    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)
    
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
      choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
    }
    
    # Update the grid by removing already evaluated points
    # grid <- matrix(grid[is.na(match(do.call(paste, data.frame(grid)), do.call(paste, data.frame(newdata$x_original)))),], ncol = D)
    # grid_transformed <- (grid-lower)/(upper - lower)
    
    # Update the GP
    print('Maximize Acquisition Function')
    initialize_UCB <- matrix(runif(D, lower, upper),
                             nrow = 1, ncol = D, byrow = T)
    initialize_UCB <- t((t(initialize_UCB) - lower)/(upper - lower))
    next_point <- optim(initialize_UCB, function(x) UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata, cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta), 
                        control = list(maxit = 100), lower = rep(0, D), upper = rep(1, D), method = 'L-BFGS-B')$par
    
    # Select the next point to evaluate
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point)
    yvec <- c(yvec, (func(next_point) - rel))
    
    # Print some information for debugging
    print(paste("Next point:", next_point))
    print(paste("Function value:", yvec[i + 1]))
  }
  
  # Return the result
  return(list(result = data.frame(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var))
}


