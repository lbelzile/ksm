#' Random generation from first-order vector autoregressive model
#'
#' Given a matrix of coefficients \code{M} and a covariance matrix \code{Sigma},
#' simulate \code{K} vectors from a first-order autoregressive process.
#'
#' @param n sample size
#' @param M matrix of autoregressive coefficients
#' @param Sigma covariance matrix
#' @param K integer, degrees of freedom
#' @param burnin number of iterations discarded
#' @param order order of autoregressive process, only \code{1} is supported at current.
#' @useDynLib ksm, .registration=TRUE
#' @return a list of length \code{n} containing matrices of size \code{K} by \code{d}
#' @keywords internal
#' @export
#' @examples
#' M <- matrix(c(0.3, -0.3, -0.3, 0.3), nrow = 2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' rVAR(n = 100, M = M, Sigma = Sigma, K = 10)
rVAR <- function(
  n,
  M,
  Sigma,
  K = 1L,
  order = 1L,
  burnin = 25L
) {
  d <- ncol(M)
  K <- as.integer(K)
  n <- as.integer(n)
  order <- as.integer(order)
  if (order != 1L) {
    stop("Only first order autoregressive model is supported at current.")
  }
  burnin = as.integer(burnin)
  stopifnot(n > 1, burnin >= 0)
  if (isTRUE(any(abs(eigen(M, only.values = TRUE)$values) >= 1))) {
    warning("VAR(1) model is not stationary.")
  }
  stopifnot(K >= 0)
  stopifnot(nrow(M) == d, ncol(Sigma) == d, nrow(Sigma) == d)
  # Initialize list to store X[t], each element is a K x d matrix
  X <- vector("list", n)
  # Generate the first time step X[1] as a K x d matrix
  X[[1]] <- rmnorm(K, mean = rep(0, d), vcov = Sigma)
  # Generate subsequent time steps using the VAR(1) structure
  for (t in 2:(n + burnin)) {
    X[[t]] <- X[[t - 1]] %*% t(M) + rmnorm(K, mean = rep(0, d), vcov = Sigma)
  }
  return(tail(X, n))
}

#' Random matrix generation from first-order autoregressive Wishart process
#'
#' Given a matrix of coefficients \code{M} and a covariance matrix \code{Sigma},
#' simulate \code{n} random matrices from a first-order autoregressive Wishart process
#' by simulating from cross-products of vector autoregressions
#'
#' @param n sample size
#' @param M matrix of autoregressive coefficients
#' @param Sigma covariance matrix
#' @param K integer, degrees of freedom
#' @param burnin number of iterations discarded
#' @param order order of autoregressive process, only \code{1} is supported at current.
#' @export
#' @return an array of size \code{d} by \code{d} by \code{n} containing the samples
#' @references C. Gourieroux, J. Jasiak, and R. Sufana (2009). The Wishart Autoregressive process of multivariate stochastic volatility, \emph{Journal of Econometrics}, 150(\bold{2}), 167-181, <doi:10.1016/j.jeconom.2008.12.016>.
#' @examples
#' M <- matrix(c(0.3, -0.3, -0.3, 0.3), nrow = 2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' rWAR(n = 10, M = M, Sigma = Sigma, K = 5)
rWAR <- function(n, M, Sigma, K = 1L, order = 1L, burnin = 25L) {
  d <- ncol(Sigma)
  order <- as.integer(order)
  if (order != 1L) {
    stop("Only first order autoregressive model is supported at current.")
  }
  # Generate the VAR(1) collection as a list of matrices
  X_list <- rVAR(n = n, M = M, Sigma = Sigma, K = K, burnin = burnin)
  array(
    data = c(sapply(X_list, function(x) {
      crossprod(x)
    })),
    dim = c(d, d, n)
  )
}


#' Bandwidth optimization for symmetric matrix kernels
#'
#' Given a sample of positive definite matrices,
#' perform numerical maximization of the \code{h}-block least square (\code{lscv}) or leave-one-out likelihood (\code{lcv}) cross-validation criteria using a root search.
#' @importFrom stats optimize
#' @param x sample of symmetric matrix observations from which to build the kernel density kernel
#' @param criterion optimization criterion, one of \code{lscv} for least square cross-validation at lag \code{h} or \code{lcv} for leave-one-out cross-validation.
#' @param kernel string, one of \code{Wishart}, \code{smlnorm} (log-Gaussian) or \code{smnorm} (Gaussian).
#' @param tol double, tolerance of optimization (root search)
#' @param h lag step for consideration of observations, for the case \code{criterion=lscv}
#' @export
#' @return double, the optimal bandwidth up to \code{tol}
bandwidth_optim <- function(
  x,
  criterion = c("lscv", "lcv"),
  kernel = c("Wishart", "smlnorm", "smnorm"),
  tol = 1e-4,
  h = 1L
) {
  criterion <- match.arg(criterion)
  kernel <- match.arg(kernel)
  if (criterion == "lscv" & kernel == "smnorm") {
    stop(
      "Least square cross validation not currently implemented for matrix normal kernel estimator."
    )
  }
  if (criterion == "lscv") {
    if (kernel == "Wishart") {
      optfun <- function(band) {
        fn <- lscv_kern_Wishart(x = x, b = band, h = h)
        -exp(fn[1]) + exp(fn[2])
      }
    } else if (kernel == "smlnorm") {
      optfun <- function(band) {
        fn <- lscv_kern_smlnorm(x = x, b = band, h = h)
        -exp(fn[1]) + exp(fn[2])
      }
    }
  } else {
    # LCV
    optfun <- function(band) {
      c(lcv_kdens_symmat(x = x, b = band, kernel = kernel)$lcv)
    }
  }

  optimize(
    f = optfun,
    lower = 1e-4,
    upper = 10,
    maximum = TRUE,
    tol = tol
  )$maximum
}
