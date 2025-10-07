# Tests for Wishart package
library(tinytest)
# Check that error if x or S is not positive definite
source("~/Documents/Dropbox/Rpackage/Wishart/R/functions.R")
Rcpp::sourceCpp("~/Documents/Dropbox/Rpackage/Wishart/src/functions.cpp")

# Check dimension of output returned by rWishart
d <- 10L
S <- Wishart:::symmetrize(diag(d) + matrix(0.5, d, d))
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
  lmgamma(3.2, 4),
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
  solve.Riccati(S = Sigma2, M = M2)$solution,
  solve_riccati(M = M2, Sigma = Sigma2)
)

## Kernels and Leave-one-out cross validation
## This does not work anymore for some reason
# expect_equal(
#   mean(sapply(1:10, function(i) {
#     kdens_Wishart(x = samp[,, i, drop = FALSE], pts = samp[,, -i], b = b)
#   })),
#   lcv_kern_Wishart(x = samp[,, 1:10], b = b)
# )
#
# expect_equal(
#   mean(sapply(1:10, function(i) {
#     kdens_smlnorm(x = samp[,, i, drop = FALSE], pts = samp[,, -i], b = b)
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
    pts = samp[,, -1],
    b = b
  )
)

expect_equal(
  log(hat_f(XX[-1], S = XX[[1]], b = b, method = "LG")),
  kdens_smlnorm(x = samp[,, 1, drop = FALSE], pts = samp[,, -1], b = b)
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
  c(lcv_kern_smlnorm(x = sampfull, b = b))
)

microbenchmark::microbenchmark(
  c(lscv_kern_Wishart(x = sampfull, b = b, h = 1)),
  c(lscv_kern_smlnorm(x = sampfull, b = b, h = 1))
)

microbenchmark::microbenchmark(
  kdens_Wishart(x = sampfull, pts = sampfull, b = b),
  kdens_smlnorm(x = sampfull, pts = sampfull, b = b)
)
