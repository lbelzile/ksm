## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5, 
  fig.height = 3,
  fig.align='center',
  global.par = TRUE
)

## -----------------------------------------------------------------------------
set.seed(2025)
library(ksm)
d <- 4L
# Equicorrelation matrix
S <- diag(rep(0.5, d)) + matrix(0.5, d, d)
# Generate simulated data from an inverse Wishart (d by d by n array).
samp <- rinvWishart(n = 100, df = 10, S)
dim(samp)
# Optimize the bandwidth for kernel using leave-one-out cross validation
b <- bandwidth_optim(
  x = samp, 
  kernel = "smlnorm", 
  criterion = "lcv")
# Evaluate kernel with a new psd matrix data point
newsamp <- rinvWishart(n = 2, df = 10, S = S)
# Kernel density (default to log scale)
kdens_symmat(
  x = newsamp, 
  xs = samp, 
  kernel = "smlnorm", 
  b = b)

## -----------------------------------------------------------------------------
data(realvar, package = "ksm")
dim(realvar)
# Select suitable lag for least square cross validation
h <- ceiling(dim(realvar)[3]^0.25)
# Calculate optimal bandwidth
bopt <- bandwidth_optim(
  x = realvar, 
  kernel = "Wishart", 
  criterion = "lscv",
  h = h)
# Evaluate at multiple values of bandwidth
bseq <- seq(0.5 * bopt, 2 * bopt, length.out = 101)
lscv <- ksm::lscv_kdens_symmat(
  x = realvar, 
  b = bseq, 
  h = h, 
  kernel = "Wishart")
with(lscv, 
     plot(x = b, 
          y = lscv,
          type = "l", 
          panel.first = {
            abline(v = bopt, lty = 2, col = "grey")
            },
          xlab = "bandwidth",
          ylab = "least square cross validation"
          )
     )

## -----------------------------------------------------------------------------
# Check that Wishart density integrates to unity
(int <- integrate_spd(
  dim = 2L,
  neval = 1e5L,
  f = function(x, S){
   dWishart(x, df = 10, S = S, log = FALSE)},
  S = diag(2),lb = 0, method = "hcubature"))
isTRUE(all.equal(int$integral, 1, tol = 2*int$error))

## -----------------------------------------------------------------------------
# Integrated squared error
ISE <- function(
  S,
  x,
  bandwidth,
  model = 1:6,
  kernel = c("Wishart", "smlnorm", "smnorm"),
  ...
) {
  model <- match.arg(as.character(model), choices = 1:6)
  kernel <- match.arg(kernel)
  dim <- ncol(x)
  # Compute the squared difference between the kernel and the target density
  (
    ksm::kdens_symmat(
      x = S,
      xs = x,
      b = bandwidth,
      kernel = kernel,
      log = FALSE
    ) -
      ksm::simu_fdens(S, model = model, d = dim)
  )^2
}

# Calculate numerical integral
d <- 2L
xs <- simu_rdens(n = 100, model = 2, d = d)
b <- bandwidth_optim(
  x = xs, 
  criterion = "lcv", 
  kernel = "smnorm")
# Computing ISE
ISE(
  S = xs[,,1, drop = FALSE], 
  x = xs, 
  bandwidth = b, 
  model = 2, 
  kernel = "smnorm")
# Next integrate over the space of psd matrices
ISE_int <- try(
            ksm::integrate_spd(
              f = function(S) {
                ISE(
                  S,
                  x = xs,
                  kernel = "smnorm",  # (only for lcv)
                  bandwidth = b,
                  model = 2
                )
              },
              dim = d,
              method = "hcubature",
              lb = 0,
              ub = Inf,
              neval = 1e6L
            ),
            silent = TRUE
          )
ISE_int

