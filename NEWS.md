# mig 1.1  (Release date 2026-06-06)

## New:

- Function `lscv_kdens_symmat`, mimicking `lcv_kdens_symmat`, to evaluate the least square cross validation criterion at multiple bandwidth
- Symmetric matrix normal least square cross validation criterion now implemented (`lscv_kern_smnorm`)
- New unit test codes for optimization and kernels
- Vignette now included

## Bug fixes:

- Function `lcv_kdens_symmat` now returns the bandwidth maximizing the criterion

## Changes:

- Leave-one-out cross validation methods now allow to lag `h` (default to 1) for selection criterion.
- Degrees of freedom for rWAR defaults to the dimension of the scale matrix provided, to ensure positive definiteness of the output.
- Simulation model 6 default degrees of freedom increased
- Change wrappers to return by default Rcpp vectors (since `Rcpp::wrap` returns column vectors as matrices).


# ksm 1.0 (Release date: 2025-10-23)

Initial release
