# Symmetrize a matrix

symmetrize <- function(X) {
  return((X + t(X)) / 2)
}

# Generate a 2x2 rotation matrix

rotation_matrix <- function(theta) {
  matrix(
    data = c(cos(theta), sin(theta), -sin(theta), cos(theta)),
    nrow = 2,
    ncol = 2
  )
}

# Construct X = rotation_matrix(theta) . diag(lambda1, lambda2) . rotation_matrix(theta)^T

construct_X <- function(theta, lambda1, lambda2) {
  R_theta <- rotation_matrix(theta)
  diag_lambda <- diag(c(lambda1, lambda2))
  X <- R_theta %*% diag_lambda %*% t(R_theta)
  symmetrize(X)
}

#' Generate vector autoregression
#'
#' @param n sample size
#' @param M matrix of autoregressive coefficients
#' @param Sigma covariance matrix
#' @param K degrees of freedom
generate_VAR1_collection <- function(n, M, Sigma, K = 1L) {
  # M <- M_list[[i]]
  # Sigma <- Sigma_list[[j]]
  d <- ncol(M)
  K <- as.integer(K)
  stopifnot(K >= 0)
  stopifnot(nrow(M) == d, ncol(Sigma) == d, nrow(Sigma) == d)
  # Initialize list to store X[t], each element is a K x d matrix
  X_list <- vector("list", n)

  # Generate the first time step X[1] as a K x d matrix
  X_t <- matrix(0, nrow = K, ncol = d)
  for (k in 1:K) {
    X_t[k, ] <- mvtnorm::rmvnorm(
      n = 1,
      mean = rep(0, d),
      sigma = Sigma
    )
  }
  X_list[[1]] <- X_t

  # Generate subsequent time steps using the VAR(1) structure
  for (t in 2:n) {
    X_t <- matrix(0, nrow = K, ncol = d)
    for (k in 1:K) {
      X_t[k, ] <- c(
        M %*%
          X_list[[t - 1]][k, ]
      ) +
        c(mvtnorm::rmvnorm(1, mean = rep(0, d), sigma = Sigma))
    }
    X_list[[t]] <- X_t
  }

  return(X_list)
}

# Function to generate the WAR(1) process as a sum of outer products of VAR(1) processes
# [LEO] rather than i,j, pass mean and covariance matrix
# extract dimension d from it.
generate_WAR1 <- function(i, j, K, n) {
  # Initialize list to store Y[t], each element is a d x d matrix
  Y_list <- vector("list", n)

  # Generate the VAR(1) collection as a list of matrices
  X_list <- generate_VAR1_collection(i, j, K, n)
  d <- ncol(X_list[[1]])
  # Compute Y[t] as the sum of outer products of X[k, t, ] over k from 1 to K
  for (t in 1:n) {
    Y_t <- matrix(0, nrow = d, ncol = d)
    for (k in 1:K) {
      X_kt <- X_list[[t]][k, ]
      Y_t <- Y_t + X_kt %*% t(X_kt)
    }
    Y_list[[t]] <- Y_t
  }

  return(Y_list)
}

# Function to solve the Riccati equation Sigma_inf = M Sigma_inf M^T + Sigma
solve_riccati <- function(
  M,
  Sigma,
  tol = 1e-8,
  max_iter = 1000L
) {
  Sigma_inf <- Sigma # Initial guess for Sigma_inf
  for (i in 1:max_iter) {
    Sigma_new <- M %*% Sigma_inf %*% t(M) + Sigma
    if (norm(Sigma_new - Sigma_inf, type = "F") < tol) break
    Sigma_inf <- Sigma_new
  }
  # Symmetrization step to ensure the result is symmetric
  symmetrize(Sigma_inf)
}

# LaplacesDemon::dwishart
dwishart <- function(Omega, nu, S, log = FALSE) {
  if (!is.matrix(Omega)) Omega <- matrix(Omega)
  if (!is.positive.definite(Omega))
    stop("Matrix Omega is not positive-definite.")
  if (!is.matrix(S)) S <- matrix(S)
  if (!is.positive.semidefinite(S))
    stop("Matrix S is not positive-semidefinite.")
  if (!identical(dim(Omega), dim(S)))
    stop("The dimensions of Omega and S differ.")
  if (nu < nrow(S)) stop("The nu parameter is less than the dimension of S.")
  k <- nrow(Omega)
  gamsum <- 0
  for (i in 1:k) {
    gamsum <- gamsum + lgamma((nu + 1 - i) / 2)
  }
  dens <- -((nu * k) / 2) *
    log(2) -
    ((k * (k - 1)) / 4) * log(pi) -
    gamsum -
    (nu / 2) * log(det(S)) +
    ((nu - k - 1) / 2) * logdet(Omega) -
    (tr(as.inverse(S) %*% Omega) / 2)
  if (log == FALSE) dens <- exp(dens)
  return(dens)
}

# Centered SMN density

G <- function(Y, b) {
  # b > 0, Y > 0 is a d x d symmetric matrix
  d <- nrow(Y)
  return(
    exp(-sum(diag(Y %*% Y)) / (2 * b)) /
      ((2 * pi * b)^(d * (d + 1) / 4)) *
      2^(d * (d - 1) / 4)
  )
}


LG <- function(X, b, S) {
  # Recenter logm(S)
  Y <- logm(S) - logm(X)

  # Calculate the G term
  G_term <- G(Y, b) / det(S)

  # Compute the eigenvalues of S
  eigenvalues <- eigen(S, only.values = TRUE)$values
  d <- length(eigenvalues)

  # Compute the product factor
  product_factor <- 1
  for (i in 1:d) {
    if (i < d) {
      for (j in (i + 1):d) {
        if (abs(eigenvalues[i] - eigenvalues[j]) > .Machine$double.eps) {
          # > 0
          numerator <- log(eigenvalues[i]) - log(eigenvalues[j])
          denominator <- eigenvalues[i] - eigenvalues[j]
          product_factor <- product_factor * abs(numerator / denominator)
        } else {
          product_factor <- product_factor * (1 / eigenvalues[i])
        }
      }
    }
  }

  # Return the product of G_term and product_factor
  return(G_term * product_factor)
}


hat_f <- function(XX, S, b, method = "WK") {
  estimator <- switch(
    method,
    "WK" = mean(sapply(
      XX,
      function(X) LaplacesDemon::dwishart(X, 1 / b + d + 1, b * S)
    )),
    "LG" = mean(sapply(XX, function(X) LG(X, b, S))),
    warning("Invalid method. Should be 'WK' or 'LG'.")
  )
  return(estimator)
}


LSCV <- function(XX, b, i, j, K, method, tolerance = tol1) {
  integrand <- function(vars) {
    # Construct SPD matrix S
    S <- construct_X(vars[1], vars[2], vars[3])

    # Compute the squared difference between the estimator and the target density
    diff_squared <- (hat_f(XX, S, b, method) - f(i, j, K, S))^2
    jacobian_value <- abs(vars[2] - vars[3]) / 4

    # Return the product of diff_squared and the Jacobian factor
    return(diff_squared * jacobian_value)
  }

  # Define the limits of integration
  lower_limit <- c(0, delta, delta)
  upper_limit <- c(2 * pi, 1 / delta, 1 / delta)

  # Perform the numerical integration using the cubature package
  result <- adaptIntegrate(
    integrand,
    lowerLimit = lower_limit,
    upperLimit = upper_limit,
    tol = tolerance
  )

  # Return the result of the integral
  return(result$integral)
}


# Least Squares Cross Validation (LSCV) criterion (Monte Carlo version)

LSCV_MC <- function(XX, b, ii, jj, K, method) {
  n <- length(XX)
  d <- nrow(XX[[1]])
  r_d <- d * (d + 1) / 2

  sum1 <- 0
  sum2 <- 0

  if (method == "WK") {
    second_term <- 0 # initialization of second term using h-block cross-validation

    for (i in 1:n) {
      n_h <- 0 # number of indices j for a given i in h-block cross-validation
      sum2 <- 0 # for h-block cross-validation

      for (j in 1:n) {
        Xi <- XX[[i]]
        Xj <- XX[[j]]
        Sij <- Xi + Xj

        # Avoid numerical errors due to Stirling's formula
        term2 <- exp(
          (-d / b - r_d) *
            log(2) -
            (r_d / 2) * log(b) +
            lmvgamma(1 / b + (d + 1) / 2, d) -
            2 * lmvgamma(1 / (2 * b) + (d + 1) / 2, d)
        )

        # Combined determinants to avoid numerical errors
        term3 <- det(solve(Sij %*% Sij, Xi %*% Xj))^(1 / (2 * b))
        term4 <- det(Sij)^(-(d + 1) / 2)

        sum1 <- sum1 + term2 * term3 * term4

        if (abs(j - i) > h(n)) {
          # h-block cross-validation
          sum2 <- sum2 +
            LaplacesDemon::dwishart(Xi, nu = 1 / b + d + 1, S = b * Xj)
          n_h <- n_h + 1 # update the number of indices j for a given i in h-block cross-validation
        }
      }
      second_term <- second_term + 2 * sum2 / (n * n_h) # second term update for h-block cross-validation
    }
    first_term <- (2^(d / b) * b^(-r_d / 2)) * sum1 / (n^2)
  } else if (method == "LG") {
    second_term <- 0 # initialization of second term using h-block cross-validation

    for (i in 1:n) {
      n_h <- 0 # number of indices j for a given i in h-block cross-validation
      sum2 <- 0 # for h-block cross-validation

      for (j in 1:n) {
        Xi <- XX[[i]]
        Xj <- XX[[j]]
        Xi_log <- logm(Xi)
        Xj_log <- logm(Xj)

        etr_term <- exp(tr(
          (-Xi_log^2 - Xj_log^2 + (Xi_log + Xj_log)^2 / 2) / (2 * b)
        ))
        sum1 <- sum1 + etr_term

        if (abs(j - i) > h(n)) {
          # h-block cross-validation
          sum2 <- sum2 + G(Xj_log - Xi_log, b)
          n_h <- n_h + 1 # update the number of indices j for a given i in h-block cross-validation
        }
      }
      second_term <- second_term + 2 * sum2 / (n * n_h) # second term update for h-block cross-validation
    }
    first_term <- sum1 / ((2 * pi * b)^(r_d / 2) * 2^(d / 2) * n^2)
  } else {
    stop("Invalid method. Should be 'WK' or 'LG'.")
  }

  return(first_term - second_term)
}


LSCV_MC_integral <- function(XX, b, j, method) {
  # Total number of Monte Carlo samples
  n_MC_total <- ceiling(1000 / cores_per_node) * cores_per_node
  n_MC_per_core <- n_MC_total / cores_per_node
  # guaranteed to be an integer

  # Set up a parallel cluster to distribute the computation
  cl <- setup_parallel_cluster()

  # Perform parallel MC estimation with local sampling
  estimates <- parLapply(cl, 1:cores_per_node, function(core_id) {
    set.seed(123 + core_id) # Optional: different seed per core
    theta_vals <- runif(
      n_MC_per_core,
      min = 0,
      max = 2 * pi
    )
    lambda1_vals <- runif(
      n_MC_per_core,
      min = delta,
      max = 1 / delta
    )
    lambda2_vals <- runif(
      n_MC_per_core,
      min = delta,
      max = 1 / delta
    )

    sapply(1:n_MC_per_core, function(i) {
      S <- construct_X(
        theta = theta_vals[i],
        lambda1 = lambda1_vals[i],
        lambda2 = lambda2_vals[i]
      )
      diff_squared <- (hat_f(XX, S, b, method) - f(j, S))^2
      jacobian_value <- abs(lambda1_vals[i] - lambda2_vals[i]) / 4
      return(diff_squared * jacobian_value)
    })
  })

  stopCluster(cl)

  # Aggregate results
  estimates_vec <- unlist(estimates)
  volume <- 2 * pi * (1 / delta - delta)^2
  result <- mean(estimates_vec) * volume

  return(result)
}


# Density function for matrix-type-II Beta distribution

dmatrixbeta_typeII <- function(X, a, b) {
  # X is an SPD matrix of size d x d, and a,b > (d - 1)/2
  d <- nrow(X)
  numerator <- mvgamma(a + b, d)
  denominator <- mvgamma(a, d) * mvgamma(b, d)
  det_X <- det(X)
  det_I_plus_X <- det(diag(d) + X)

  density <- (numerator / denominator) *
    det_X^(a - (d + 1) / 2) *
    det_I_plus_X^(-(a + b))
  return(density)
}

# Likelihood Cross Validation (LCV) criterion

LCV <- function(XX, b, method) {
  n <- length(XX)
  lcv_values <- numeric(n)

  for (i in 1:n) {
    XX_minus_i <- XX[-i]
    density_estimate <- hat_f(XX_minus_i, XX[[i]], b, method)
    lcv_values[i] <- log(density_estimate)
  }

  return(mean(lcv_values))
}
