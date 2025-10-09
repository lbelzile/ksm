# Tests for Wishart package
library(tinytest)
# Check that error if x or S is not positive definite
source(
  "~/Documents/Dropbox/Publications/Ongoing/BGOR_2024_Wishart_asym_kernels/Code/functions.R"
)
library(Wishart)

# Check dimension of output returned by rWishart
d <- 10L
S <- Wishart::symmetrize(diag(d) + matrix(0.5, d, d))
df <- d + 1
b <- 2
n <- 100L
set.seed(6543)
samp <- rWishart(
  n = n,
  df = df,
  S = S
)
for (i in 1:n) {
  samp[,, i] <- symmetrize(cov2cor(samp[,, i]))
}
sampfull <- samp
expect_equal(dim(samp), c(d, d, n))
expect_equal(
  current = c(dWishart(x = samp[,, 1, drop = FALSE], df = df, S = S)),
  LaplacesDemon::dwishart(Omega = symmetrize(samp[,, 1]), nu = df, S = S)
)

expect_equal(
  current = c(dinvWishart(
    x = samp[,, 1, drop = FALSE],
    df = df,
    S = S,
    log = TRUE
  )),
  LaplacesDemon::dinvwishart(
    log = TRUE,
    Sigma = symmetrize(samp[,, 1]),
    nu = df,
    S = S
  )
)


expect_equal(
  LaplacesDemon::dwishart(
    symmetrize(samp[,, 1]),
    nu = 1 / b + d + 1,
    S = symmetrize(b * samp[,, 1]),
    log = TRUE
  ),
  c(dWishart(
    samp[,, 1, drop = FALSE],
    df = 1 / b + d + 1,
    S = b * samp[,, 1],
    log = TRUE
  ))
)

isTRUE(all(
  apply(samp, 3, function(x) {
    tail(eigen(x, only.values = TRUE)$values, 1)
  }) >
    0
))

# Check matrix-variate normal
expect_equal(
  c(dsmnorm(x = samp, b = b, M = S, log = TRUE)),
  log(apply(samp, 3, function(x) {
    G(x - S, b = b)
  }))
)

expect_equal(
  log(LG(X = S, b = b, S = samp[,, 1])),
  c(dsmlnorm(x = samp[,, 1, drop = FALSE], b = b, M = S, log = TRUE))
)

xg <- rexp(n = 10, rate = 0.1)
expect_equal(
  c(CholWishart::lmvgamma(xg, p = 5)),
  c(mgamma(xg, p = 5, log = TRUE))
)
expect_equal(
  Wishart:::lmgamma(3.2, 4),
  c(CholWishart::lmvgamma(3.2, 4))
)


expect_equal(
  c(dmbeta2(samp, shape1 = 15, shape2 = 10)),
  log(apply(samp, 3, function(x) {
    dmatrixbeta_typeII(x, 15, 10)
  }))
)

M2 <- matrix(c(0.3, -0.3, -0.3, 0.3), nrow = 2)
# Covariance matrices
Sigma2 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

expect_equal(
  Riccati(S = Sigma2, M = M2)$solution,
  solve_riccati(M = M2, Sigma = Sigma2)
)

## Kernels and Leave-one-out cross validation
## This does not work anymore for some reason
# expect_equal(
#   mean(sapply(1:10, function(i) {
#     kdens_Wishart(x = samp[,, i, drop = FALSE], xs = samp[,, -i], b = b)
#   })),
#   lcv_kern_Wishart(x = samp[,, 1:10], b = b)
# )
#
# expect_equal(
#   mean(sapply(1:10, function(i) {
#     kdens_smlnorm(x = samp[,, i, drop = FALSE], xs = samp[,, -i], b = b)
#   })),
#   lcv_kern_smlnorm(x = samp[,, 1:10], b = b)
# )

u <- runif(100)
expect_equal(
  exp(meanlog(log(u))),
  mean(u)
)


# Reduce sample size since R functions are slow

samp <- samp[,, 1:10]
XX <- list()
for (i in 1:dim(samp)[3]) {
  XX[[i]] <- symmetrize(samp[,, i])
}


expect_equal(
  log(hat_f(XX[-1], S = XX[[1]], b = b, method = "WK")),
  kdens_Wishart(
    x = samp[,, 1, drop = FALSE],
    xs = samp[,, -1],
    b = b
  )
)

expect_equal(
  log(hat_f(XX[-1], S = XX[[1]], b = b, method = "LG")),
  kdens_smlnorm(x = samp[,, 1, drop = FALSE], xs = samp[,, -1], b = b)
)

expect_equal(
  LCV(XX = XX, b = b, method = "LG"),
  lcv_kern_smlnorm(x = samp[,, 1:10], b = b)
)
expect_equal(
  LCV(XX = XX, b = b, method = "WK"),
  lcv_kern_Wishart(x = samp[,, 1:10], b = b)
)


expect_equal(
  c(lscv_kern_Wishart(x = samp, b = b, h = 1)),
  LSCV_MC(XX = XX, b = b, method = "WK", h = 1)
)

expect_equal(
  c(lscv_kern_Wishart(x = samp, b = b, h = 2)),
  LSCV_MC(XX = XX, b = b, method = "WK", h = 2)
)
expect_equal(
  c(lscv_kern_smlnorm(x = samp, b = b, h = 1)),
  LSCV_MC(XX = XX, b = b, method = "LG", h = 1)
)

expect_equal(
  c(lscv_kern_smlnorm(x = samp, b = b, h = 2)),
  LSCV_MC(XX = XX, b = b, method = "LG", h = 2)
)

microbenchmark::microbenchmark(
  c(lcv_kern_Wishart(x = sampfull, b = b)),
  c(lcv_kern_smlnorm(x = sampfull, b = b)),
  c(lcv_kern_smnorm(x = sampfull, b = b))
)

microbenchmark::microbenchmark(
  c(lscv_kern_Wishart(x = sampfull, b = b, h = 1)),
  c(lscv_kern_smlnorm(x = sampfull, b = b, h = 1))
)

microbenchmark::microbenchmark(
  kdens_Wishart(x = sampfull, xs = sampfull, b = b),
  kdens_smlnorm(x = sampfull, xs = sampfull, b = b),
  kdens_smnorm(x = sampfull, xs = sampfull, b = b)
)


# Define the integrand for cubature
integrand <- function(vars, shape) {
  # Construct SPD matrix X from the parameters theta, lambda1, lambda2
  X <- array(
    Wishart::rotation_scaling(vars[1], c(vars[2], vars[3])),
    dim = c(2, 2, 1)
  )
  # Compute the density using the dmatrixbeta_typeII function
  density_value <- dmbeta2(X, shape[1], shape[2], log = FALSE)

  # Return the density multiplied by the Jacobian adjustment for polar coordinates
  jacobian_value <- abs(vars[2] - vars[3]) / 4
  return(density_value * jacobian_value)
}

# Perform the numerical integration using adaptIntegrate from the cubature package
result_test <- cubature::adaptIntegrate(
  integrand,
  lowerLimit = c(0, 0, 0),
  upperLimit = c(2 * pi, Inf, Inf),
  tol = 1e-4,
  shape = 2 + rexp(2)
)
expect_equal(1, result_test$integral, tolerance = 1e-4)
