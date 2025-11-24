# Tests for Wishart package
library(tinytest)
library(ksm)
# Check dimension of output returned by rWishart
d <- 4L
S <- ksm::symmetrize(diag(d) + matrix(0.5, d, d))
# Check symmetrization
expect_true(
  isSymmetric(
    ksm::symmetrize(
      cbind(c(2, 1), c(0.4, 3))
    )
  )
)

# Check simulations

M <- matrix(c(0.3, -0.3, -0.3, 0.3), nrow = 2)
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
samp <- rWAR(n = 10, M = M, Sigma = Sigma)
# Check that default degrees of freedom yield a
# valid set of WAR processes
expect_true(isTRUE(all(apply(samp, 3, function(x) {
  min(eigen(x, only.values = TRUE)$values) > 0
}))))


# Check that all models return numeric solution
# h = 2; kernel = "Wishart"; criterion = "lscv"

x <- simu_rdens(n = 20, model = 1, d = 3)
for (h in 1:2) {
  for (criterion in c("lscv", "lcv")) {
    for (kernel in c("Wishart", "smlnorm", "smnorm")) {
      b <- try(
        bandwidth_optim(
          x = x,
          criterion = criterion,
          kernel = kernel,
          h = h
        ),
        silent = TRUE
      )
      if (criterion == "lscv" & kernel == "smnorm") {
        print(expect_true(inherits(b, "try-error")))
      } else {
        print(expect_true(is.numeric(b)))
        # Compare the method for optim with the lcv wrapper
        band <- sort(c(b, seq(0.5 * b, 2 * b, length.out = 11)))
        if (criterion == "lcv") {
          crit <- ksm::lcv_kdens_symmat(
            x = x,
            b = band,
            h = h,
            kernel = kernel
          )
          # plot(c(crit$b), crit$lcv)
          criter <- crit$lcv[crit$b %in% crit$bandwidth]
        } else if (criterion == "lscv") {
          crit <- ksm::lscv_kdens_symmat(
            x = x,
            b = band,
            h = h,
            kernel = kernel
          )
          criter <- crit$lscv[crit$b %in% crit$bandwidth]
        }
        print(expect_true(crit$bandwidth == b))
        criter2 <- switch(
          kernel,
          Wishart = switch(
            criterion,
            lcv = lcv_kern_Wishart(x = x, b = b, h = h),
            lscv = lscv_kern_Wishart(x = x, b = b, h = h)
          ),
          smlnorm = switch(
            criterion,
            lcv = lcv_kern_smlnorm(x = x, b = b, h = h),
            lscv = lscv_kern_smlnorm(x = x, b = b, h = h)
          ),
          smnorm = lcv_kern_smnorm(x = x, b = b, h = h)
        )
        print(expect_equal(
          criter,
          criter2
        ))
      }
    }
  }
}

# Check that this returns a solution
expect_true(Riccati(M, Sigma, tol = 1e-08)$convergence)

# Should return vector
expect_true(is.array(simu_rdens(100, 2, 3)))
expect_true(is.vector(simu_fdens(x, 2, 3)))

# Check kernel densities (log, wrapper for symmat, positive densities)

dens <- ksm::kdens_smlnorm(x = x, xs = x, b = 2, FALSE)
expect_true(is.vector(dens))
expect_equal(
  exp(ksm::kdens_smlnorm(x = x, xs = x, b = 2, TRUE)),
  dens
)
expect_equal(
  ksm::kdens_symmat(x = x, xs = x, b = 2, kernel = "smlnorm", log = FALSE),
  dens
)
expect_true(isTRUE(all(dens >= 0)))

dens <- ksm::kdens_smnorm(x = x, xs = x, b = 2, FALSE)
expect_true(is.vector(dens))
expect_equal(
  exp(ksm::kdens_smnorm(x = x, xs = x, b = 2, TRUE)),
  dens
)
expect_equal(
  ksm::kdens_symmat(x = x, xs = x, b = 2, kernel = "smnorm", log = FALSE),
  dens
)
expect_true(isTRUE(all(dens >= 0)))

dens <- ksm::kdens_Wishart(x = x, xs = x, b = 2, FALSE)
expect_true(is.vector(dens))
expect_equal(
  exp(ksm::kdens_Wishart(x = x, xs = x, b = 2, TRUE)),
  dens
)
expect_equal(
  ksm::kdens_symmat(x = x, xs = x, b = 2, kernel = "Wishart", log = FALSE),
  dens
)
expect_true(isTRUE(all(dens >= 0)))
