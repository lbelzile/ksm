// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Symmetrize matrix
//'
//'  Given an input matrix, symmetrize by taking average of lower and upper triangular components as \eqn{A + A^\top}.
//'
//' @param A square matrix
//' @return symmetrized version of \code{A}
//' @export
// [[Rcpp::export]]
arma::mat symmetrize(const arma::mat &A){
 arma::mat B = (A.t() + A) / 2;
 return B;
}

//' Log of sum with precision
//' @param x vector of log components to add
//' @return double with the log sum of elements
//' @keywords internal
//' @export
// [[Rcpp::export]]
double sumlog(const arma::vec &x){
  double x_max = arma::max(x);
  double logsum_x = x_max + std::log(arma::accu(arma::exp(x - x_max)));
  return logsum_x;
}

//' Log of mean with precision
//'
//' @keywords internal
//' @param x vector of log components to add
//' @return double with the log mean of elements
// [[Rcpp::export]]
double meanlog(const arma::vec &x){
  arma::uword m = x.n_elem;
  double logmean_x = sumlog(x) - std::log(m);
  return logmean_x;
}

//' Signed sum with precision
//'
//' Given a vector of signs and log arguments,
//' return their sum with precision avoiding
//' numerical underflow assuming that the sum is strictly positive.
//' @keywords internal
//' @param x vector of log components
//' @param sgn sign of the operator
//' @return log sum of elements
//' @export
// [[Rcpp::export]]
double sumsignedlog(const arma::vec &x, Rcpp::LogicalVector sgn){
  int n = sgn.length();
  if((int) x.n_elem != n){
    Rcpp::stop("Invalid arguments: \"x\" and \"sgn\" must match.");
  }
  arma::vec b = arma::abs(x);
  double bmax = b.max();
  double neg = 0;
  double pos = 0;
  for(int i = 0; i < n; i++){
    if(sgn[i] == 1){
     pos += std::exp(b(i) - bmax);
    } else if(sgn[i] == 0){
      neg += std::exp(b(i) - bmax);
    }
  }
    double out = bmax + std::log(pos - neg);
    return out;
}


//' Rotation matrix for 2 dimensional problems
//'
//' @param par double for angle in \eqn{[0,2\pi).}
//' @return a 2 by 2 rotation matrix
//' @export
//' @keywords internal
// [[Rcpp::export]]
arma::mat rotation2d(double par){
 arma::mat R(2,2);
 R(0,0) = std::cos(par);
 R(1,1) = R(0,0);
 R(1,0) = std::sin(par);
 R(0,1) = -R(1,0);
 return R;
}

//' Rotation matrix for 3 dimensional problems
//'
//' @param par vector of length 3 containing elements \eqn{\phi \in [0, 2*\pi)}, \eqn{\theta \in [0, \pi]} and \eqn{\psi \in [0, 2\pi)}, in that order.
//' @return a 3 by 3 rotation matrix
//' @export
//' @keywords internal
// [[Rcpp::export]]
 arma::mat rotation3d(arma::vec par){
   if(par.n_elem != 3){
     Rcpp::stop("Invalid input.");
   }
   // Order:
   // 1. \phi \in [0, 2*\pi),
   // 2. \theta \in [0, \pi],
   // 3, \psi \in [0, 2\pi)
   arma::vec cos_t = arma::cos(par);
   arma::vec sin_t = arma::sin(par);
   arma::mat R(3,3);
   R(0,0) = cos_t(0) * cos_t(1) * cos_t(2) - sin_t(0) * sin_t(2);
   R(1,0) = sin_t(0) * cos_t(1) * cos_t(2) + cos_t(0) * sin_t(2);
   R(2,0) = - sin_t(1)*cos_t(2);
   R(0,1) = - cos_t(0) * cos_t(1) * sin_t(2) - sin_t(0) * cos_t(2);
   R(1,1) = -sin_t(0)*cos_t(1) * sin_t(2) + cos_t(0) * cos_t(2);
   R(2,1) = sin_t(1) * sin_t(2);
   R(0,2) = cos_t(0) * sin_t(1);
   R(1,2) = sin_t(0) * sin_t(1);
   R(2,2) = cos_t(1);
   return R;
 }


//' Rotation matrix with scaling for Monte Carlo integration
//'
//' @param ang vector of length 1 containing  \eqn{\theta \in [0, 2\pi]} for
//' \eqn{d=2}, or of 3 containing elements \eqn{\phi \in [0, 2*\pi)}, \eqn{\theta \in [0, \pi]} and \eqn{\psi \in [0, 2\pi)} for \eqn{d=3}.
//' @param scale vector of length 2 (\eqn{d=2}) or 3 (\eqn{d=3}), strictly positive
//' @return a 2 by 2, or 3 by 3 scaling matrix, depending on inputs
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::mat rotation_scaling(arma::vec ang, arma::vec scale){
  if(ang.n_elem == 1){
    arma::mat R = rotation2d(ang(0));
    if((scale.n_elem != 2) | (arma::any(scale < 0))){
      Rcpp::stop("Invalid argument \"scale\".");
    }
    return symmetrize(R * arma::diagmat(scale) * R.t());
  } else if(ang.n_elem == 3){
    arma::mat R = rotation3d(ang);
    if((scale.n_elem != 3) | (arma::any(scale < 0))){
      Rcpp::stop("Invalid argument \"scale\".");
    }
    return symmetrize(R * arma::diagmat(scale) * R.t());
  } else{
    Rcpp::stop("Invalid arguments.");
  }
}


//' Density of Wishart random matrix
//'
//' @param x array of dimension \code{d} by \code{d} by \code{n}
//' @param S symmetric positive definite matrix of dimension \code{d} by \code{d}
//' @param df degrees of freedom
//' @param log logical; if \code{TRUE}, returns the log density
//' @return a vector of length \code{n} containing the log-density of the Wishart.
//' @export
// [[Rcpp::export('dWishart')]]
arma::vec dWishart(const arma::cube &x, double df, const arma::mat &S, bool log = false){
  arma::uword n = x.n_slices;
  arma::uword d = x.n_rows;
    if((x.n_cols != d) | (d != S.n_rows) | (d != S.n_cols)){
        Rcpp::stop("Non conformal size for \"S\" and sample \"x\": the latter should be an d by d by n cube, and S a square d by d matrix.");
  }
  if(!S.is_sympd()){
     Rcpp::stop("\"S\" is not positive definite.");
  }
  if(df < (double) d){
        Rcpp::stop("Invalid degrees of freedom \"df\".");
  }
  arma::vec logdens(n);
  double cst = 0;
  double logdetS = arma::log_det_sympd(S);
  arma::mat Sinv;
  for(arma::uword i=0; i < d; i++){
   cst -= std::lgamma(0.5 * (df - i));
  }
  cst -= 0.5 * df * d * std::log(2.0);
  cst -= 0.25 * d * (d - 1.0) * std::log(arma::datum::pi);
  cst -= 0.5 * df * logdetS;
  for(arma::uword i=0; i<n; i++){
   logdens(i) = cst + 0.5 * (df - d - 1) * arma::log_det_sympd(x.slice(i)) -
   	0.5 * arma::trace(arma::solve(S, x.slice(i), arma::solve_opts::likely_sympd));
  }
  if(log){
    return logdens;
  } else{
    return arma::exp(logdens);
  }
}

// [[Rcpp::export('dWishart_mat')]]
double dWishart_mat(const arma::mat &x, double df, const arma::mat &S, bool log = false){
  arma::uword d = x.n_rows;
  if(df < (double) d){
    Rcpp::stop("Invalid degrees of freedom \"df\".");
  }
  double cst = 0;
  arma::mat Sinv;
  for(arma::uword i = 0; i < d; i++){
    cst -= std::lgamma(0.5 * (df - i));
  }
  cst -= 0.5 * df * d * std::log(2.0);
  cst -= 0.25 * d * (d - 1.0) * std::log(arma::datum::pi);
  cst -= 0.5 * df *  arma::log_det_sympd(S);
  double logdens = cst + 0.5 * (df - d - 1) * arma::log_det_sympd(x) -
      0.5 * arma::trace(arma::solve(S, x, arma::solve_opts::likely_sympd));
if(log){
    return logdens;
  } else{
    return std::exp(logdens);
  }
}

//' Random matrix generation from Wishart distribution
//'
//' @param n [integer] sample size
//' @param df [double] degrees of freedom, positive
//' @param S [matrix] a \code{d} by \code{d} positive definite scale matrix
//' @return an array of dimension \code{d} by \code{d} by \code{n} containing the samples
//' @export
// [[Rcpp::export('rWishart')]]
arma::cube rWishart(int n, double df, const arma::mat &S){
  arma::uword d = S.n_cols;
      if(d != S.n_rows){
          Rcpp::stop("\"S\" should be a square d by d matrix.");
  }
      if(!S.is_sympd()){
     Rcpp::stop("\"S\" is not positive definite.");
  }
  if(df < (double) d){
        Rcpp::stop("Invalid degrees of freedom \"df\".");
  }
  arma::cube x(d, d, n);
  if(n < 1) {
      Rcpp::stop("Sample size must be positive.");
  }
  arma::mat R = chol(S);
  for(int i = 0; i < n; i++){
    x.slice(i) = wishrnd(S, df, R);
  }
   return x;
}

//' Symmetric matrix-variate normal density
//'
//' @param x [cube] array of dimension \code{d} by \code{d} by \code{n}
//' @param M [matrix] location matrix, positive definite
//' @param b [numeric] scale parameter, strictly positive
//' @param log [logical] if \code{TRUE} (default), returns the log density
//' @return a vector of length \code{n}
//' @export
// [[Rcpp::export]]
arma::vec dsmnorm(const arma::cube &x, double b, const arma::mat &M, bool log = true){
  arma::uword n = x.n_slices;
  arma::uword d = x.n_rows;
  if((x.n_cols != d) | (d != M.n_rows) | (d != M.n_cols)){
    Rcpp::stop("Non conformal size for \"M\" and sample \"x\": the latter should be an d by d by n cube, and M a square d by d matrix.");
  }
  // M = symmatu(M);
  // if(!M.is_sympd()){
  //   Rcpp::stop("\"M\" is not positive definite.");
  // }

  if(b <= 0){
    Rcpp::stop("\"b\" must be a positive scale.");
  }
  double logcst;
  logcst = 0.25 * d * (d + 1.0) * std::log(2.0 * arma::datum::pi * b) - 0.25 * d * (d - 1.0) * std::log(2);
  arma::vec logdens(n);
  for(arma::uword i = 0; i < n; i++){
   logdens(i) = - arma::accu(arma::pow(x.slice(i) - M, 2.0)) / (2.0 * b) - logcst;
  }
  if(log){
    return logdens;
  } else{
    return arma::exp(logdens);
  }
}

// [[Rcpp::export]]
double dsmnorm_mat(const arma::mat &x, double b, const arma::mat &M, bool log = true){
  int d = x.n_rows;
  if(b <= 0){
    Rcpp::stop("\"b\" must be a positive scale.");
  }
  double logcst;
  logcst = 0.25 * d * (d + 1.0) * std::log(2.0 * arma::datum::pi * b) - 0.25 * d * (d - 1.0) * std::log(2);
    double logdens = - arma::accu(arma::pow(x - M, 2.0)) / (2.0 * b) - logcst;
  if(log){
    return logdens;
  } else{
    return std::exp(logdens);
  }
}

//' Symmetric matrix-variate lognormal density
//'
//' Density of the lognormal matrix-variate density, defined through the matrix logarithm, with the Jacobian resulting from the transformation
//' @param x [cube] array of dimension \code{d} by \code{d} by \code{n}
//' @param M [matrix] location matrix, positive definite
//' @param b [numeric] scale parameter, strictly positive
//' @param log [logical] if \code{TRUE} (default), returns the log density
//' @return a vector of length \code{n}
//' @export
// [[Rcpp::export]]
arma::vec dsmlnorm(
    const arma::cube &x ,
    double b,
    const arma::mat &M,
    bool log = true){
  arma::uword n = x.n_slices;
  arma::uword d = x.n_rows;
  double logjac;
  arma::vec eigvals(d);
  arma::vec log_eigvals(d);
  arma::vec logdens(n);
  logdens.zeros();
  for(arma::uword k = 0; k < n; k++){
    logjac = 0;
    // Eigenvalues are stored in ASCENDING order
    eigvals = arma::eig_sym(x.slice(k));
    log_eigvals = eigvals;
    log_eigvals = arma::log(log_eigvals);
    // PSD matrix have positive eigenvalues
    for(arma::uword i = 0; i < d - 1; i++){
      for(arma::uword j = i + 1; j < d; j++){
        if(log_eigvals(i) != log_eigvals(j)){
          logjac += std::log(log_eigvals(j) - log_eigvals(i))-
            std::log(eigvals(j) - eigvals(i));
        } else{
          logjac += log_eigvals(j);
        }
      }
    }
   logdens(k) += -log_det_sympd(x.slice(k)) + logjac;
  }
  // Copy to apply function to each element
  arma::cube matlog_x = x; //(d,d,n);
  for(arma::uword k = 0; k < n; k++){
  // WARNING: applying logmat_sympd changes x.slice
    matlog_x.slice(k) = arma::logmat_sympd(matlog_x.slice(k));
  }
   // take matrix logarithm of both matrices
 logdens += dsmnorm_mat(
     matlog_x,
     b,
     arma::logmat_sympd(M),
     log = true);

   if(log){
     return logdens;
   } else{
     return arma::exp(logdens);
   }
 }

//' Multivariate gamma function
//'
//' @param x [vector] of points at which to evaluate the function
//' @param p [int] dimension of the multivariate gamma function, strictly positive.
//' @param log [logical] if \code{TRUE}, returns the log multivariate gamma function.
//' The function is defined as
//' \deqn{\gamma_p(x) = \pi^{p(p-1)/4}\prod_{i=1}^p \Gamma\{x + (1-i)/2\}.}
//' @export
// [[Rcpp::export]]
arma::vec mgamma(const arma::vec &x, int p, bool log = false){
 if(p < 1){
   Rcpp::stop("Invalid argument\"p\": must be a positive integer.");
 }
 arma::vec out(x.n_elem, arma::fill::value(0.25 * p * (p - 1) * std::log(arma::datum::pi)));
 for(int i = 1; i <= p; i++){
  out += arma::lgamma(x + 0.5 * (1 - i));
 }
 if(log){
   return out;
 } else{
   return arma::exp(out);
 }
}

// [[Rcpp::export]]
double lmgamma(double x, arma::uword p){
  if(p < 1){
    Rcpp::stop("Invalid argument\"p\": must be a positive integer.");
  }
  double out = 0.25 * p * (p - 1) * std::log(arma::datum::pi);
  for(arma::uword i = 0; i < p; i++){
    out += std::lgamma(x - 0.5 * i);
  }
  return out;
}

//' Matrix beta type II density function
//'
//' Given a random matrix \code{x}, compute the density
//' for arguments \code{shape1} and \code{shape2}
//' @param x cube of dimension \code{d} by \code{d} by \code{n} containing the random matrix samples
//' @param shape1 positive shape parameter, strictly larger than \eqn{(d-1)/2}.
//' @param shape2 positive shape parameter, strictly larger than \eqn{(d-1)/2}.
//' @param log [logical] if \code{TRUE} (default), returns the log density.
//' @return a vector of length \code{n}
//' @export
// [[Rcpp::export]]
arma::vec dmbeta2(const arma::cube &x, double shape1, double shape2, bool log = true) {
  arma::uword d = x.n_rows;
  arma::uword n = x.n_slices;
  if((d != x.n_cols) | (shape1 < (0.5 * (d - 1.0))) | (shape2 < (0.5 * (d - 1.0)))){
    Rcpp::stop("Invalid arguments: shape parameters must be larger than (d-1)/2");
  }
  double logcst = lmgamma(shape1 + shape2, d) - lmgamma(shape1, d) - lmgamma(shape2, d);
  arma::vec logdens(n);
  for(arma::uword k = 0; k < n; k++){
  logdens(k) = (shape1 - (0.5 * (d + 1.0))) * arma::log_det_sympd(x.slice(k)) -
    (shape1 + shape2) * arma::log_det_sympd(arma::eye(d, d) + x.slice(k)) + logcst;
  }
  if(log){
    return logdens;
  } else{
    return arma::exp(logdens);
  }
}
//' Solver for Riccati equation
//'
//' Given two matrices \code{M} and \code{S}, solve Riccati equation by iterative updating to find the solution \eqn{\mathbf{R}}, where the latter satisfies
//' \deqn{\mathbf{R}=\mathbf{M}\mathbf{R}\mathbf{M}^\top + \mathbf{S}}
//' until convergence (i.e., when the Frobenius norm is less than \code{tol}, or the maximum number of iterations \code{maxiter} is reached.
//' @param M matrix
//' @param S matrix
//' @param tol double for tolerance
//' @param maxiter integer, the maximum number of iterations
//' @export
//' @return a list containing
//' \itemize{
//' \item \code{solution} matrix solution to Riccati's equation
//' \item \code{error} numerical error
//' \item \code{niter} number of iteration
//' \item \code{convergence} bool indicating convergence (\code{TRUE}) if \code{niter < maxiter}
//' }
// [[Rcpp::export(Riccati)]]
Rcpp::List solvericcati(
    arma::mat M,
    arma::mat S,
    double tol = 1e-8,
    int maxiter = 1e4){
  arma::mat sol = S;
  double err;
  arma::mat update_sol(S.n_rows, S.n_cols);
  int niter = 0;
  for(int i = 0; i < maxiter; i++){
    niter ++;
    update_sol = M * sol * M.t() + S;
    err = arma::norm(update_sol - sol, "fro");
    if(err < tol){
      break;
    }
    sol = update_sol;
    }
  return Rcpp::List::create(
    Rcpp::Named("solution") = sol,
    Rcpp::Named("error") = err,
    Rcpp::Named("niter") = niter,
    Rcpp::Named("convergence") = niter < (maxiter - 1)
    );
}

// [[Rcpp::export]]
double dsmlnorm_mat(
    const arma::mat &x,
    const arma::mat &matlog_x,
    double b,
    const arma::mat &M,
    const arma::mat &matlog_M,
    bool log = true){
  arma::uword d = x.n_rows;
  double logjac = 0;
  arma::vec eigvals(d);
  arma::vec log_eigvals(d);
  double logdens = 0;
  // Eigenvalues are stored in ASCENDING order
  eigvals = arma::eig_sym(x);
  log_eigvals = eigvals;
  log_eigvals = arma::log(log_eigvals);
  // PSD matrix have positive eigenvalues
  for(arma::uword i = 0; i < d - 1; i++){
    for(arma::uword j = i + 1; j < d; j++){
      if(log_eigvals(i) != log_eigvals(j)){
        logjac += std::log(log_eigvals(j) - log_eigvals(i))-
          std::log(eigvals(j) - eigvals(i));
      } else{
        logjac += log_eigvals(j);
      }
    }
  }
  logjac -= arma::accu(log_eigvals); //log_det_sympd(x)
  logdens +=  logjac + dsmnorm_mat(
    matlog_x,
    b,
    matlog_M,
    log = true);
  if(log){
    return logdens;
  } else{
    return std::exp(logdens);
  }
}



//' Likelihood cross validation criterion for symmetric matrix lognormal kernel
//'
//' Given a cube \code{x} and a bandwidth \code{b}, compute
//' the leave-one-out cross validation criterion by taking out a slice
//' and evaluating the kernel at the holdout value.
//' @export
//' @inheritParams dsmlnorm
//' @return the value of the log objective function
//' @keywords internal
// [[Rcpp::export(lcv_kern_smlnorm)]]
double lcvkernsmlnorm(const arma::cube &x, const double &b){
  arma::uword n = x.n_slices;
  arma::uword d = x.n_cols;
  arma::cube logmat_x(d, d, n);
  for(arma::uword i = 0; i < n; i++){
   logmat_x.slice(i) = arma::logmat_sympd(x.slice(i));
  }
  arma::vec loo(n - 1);
  double logcrit = 0;
  int it;


  double logjac;
  arma::vec eigvals(d);
  arma::vec log_eigvals(d);
  double logcst;
  logcst = 0.25 * d * (d + 1.0) * std::log(2.0 * arma::datum::pi * b) - 0.25 * d * (d - 1.0) * std::log(2);

  for(arma::uword k = 0; k < n; k++){
    it = 0;
    logjac = 0;
    // Eigenvalues are stored in ASCENDING order
    eigvals = arma::eig_sym(x.slice(k));
    log_eigvals = eigvals;
    log_eigvals = arma::log(log_eigvals);
    // PSD matrix have positive eigenvalues
    for(arma::uword i = 0; i < d - 1; i++){
      for(arma::uword j = i + 1; j < d; j++){
        if(log_eigvals(i) != log_eigvals(j)){
          logjac += std::log(log_eigvals(j) - log_eigvals(i))-
            std::log(eigvals(j) - eigvals(i));
        } else{
          logjac += log_eigvals(j);
        }
      }
    }
    logjac -= arma::accu(log_eigvals);
    for(arma::uword j = 0; j < n; j++){
      if(k != j){
      loo(it) = - arma::accu(arma::pow(logmat_x.slice(k) - logmat_x.slice(j), 2.0)) * 0.5  / b;
        it++;
    }
    }
    logcrit += meanlog(loo) + logjac;
  }
  logcrit = logcrit / n - logcst;
  return logcrit;
}


//' Likelihood cross validation criterion for Wishart kernel
//'
//' Given a cube \code{x} and a bandwidth \code{b}, compute
//' the leave-one-out cross validation criterion by taking out a slice
//' and evaluating the kernel at the holdout value.
//'
//' @inheritParams dsmlnorm
//' @keywords internal
//' @export
//' @return the value of the log objective function
// [[Rcpp::export(lcv_kern_Wishart)]]
double lcvkernWishart(const arma::cube &x, const double &b){
   arma::uword n = x.n_slices;
   arma::uword d = x.n_cols;
   double bandwidth = 1.0 / b + d + 1.0;
   // This ensures we get a finite density if nu > d+1
   arma::vec loo(n - 1);
   double logcrit = 0;
   arma::vec logdetx(n);
   for(arma::uword i = 0; i < n; i++){
    logdetx(i) = arma::log_det_sympd(x.slice(i));
   }
   arma::mat Sinv(d,d);
   double logcst = 0;
   for(arma::uword i = 0; i < d; i++){
     logcst -= std::lgamma(0.5 * (bandwidth - i));
   }
   logcst -= 0.5 * bandwidth * d * (std::log(2.0) + std::log(b));
   logcst -= 0.25 * d * (d - 1.0) * std::log(arma::datum::pi);
   int it;
   for(arma::uword k = 0; k < n; k++){
     it = 0;
     loo.zeros();
     Sinv = arma::inv_sympd(x.slice(k));
     for(arma::uword j = 0; j < n; j++){
       if(k != j){
         // TODO: finish this to avoid recomputing unnecessarily everything
         // Determinants can be recycled
       loo(it) =  0.5 / b * (logdetx(j)  - arma::accu(Sinv % x.slice(j))) + logcst - bandwidth * 0.5 * logdetx(k);
      it++;
      }
     }
     logcrit += meanlog(loo);
   }
   logcrit = logcrit / n;
   return logcrit;
 }

//' Likelihood cross-validation for symmetric positive definite matrix kernels
//'
//' Given a cube of sample observations (consisting of random symmetric positive definite matrices), and a vector of candidate bandwidth parameters \code{b},
//' compute the leave-one-out likelihood cross-validation criterion and
//' return the bandwidth among the choices that minimizes the criterion.
//' @param x array of dimension \code{d} by \code{d} by \code{n}
//' @param b vector of candidate bandwidth, strictly positive
//' @param kernel string indicating the kernel, one of \code{Wishart} or \code{smlnorm}.
//' @importFrom utils tail
//' @import Rcpp
//' @export
//' @return a list with arguments
//' \itemize{
//' \item \code{lcv} vector of likelihood cross validation criterion
//' \item \code{b} vector of candidate bandwidth
//' \item \code{bandwidth} optimal bandwidth among candidates
//' \item \code{kernel} string indicating the choice of kernel function
//'}
// [[Rcpp::export(lcv_symmat)]]
Rcpp::List lcvsymmat(
    const arma::cube &x,
    const arma::vec b,
    std::string kernel = "Wishart"){
 int nl = (int) b.n_elem;
 Rcpp::NumericVector crit(nl);
 if(kernel == "Wishart"){
   for(int i = 0; i < nl; i++){
     crit[i] = lcvkernWishart(x, b(i));
   }
 } else if(kernel == "smlnorm"){
   for(int i = 0; i < nl; i++){
     crit[i] = lcvkernsmlnorm(x, b(i));
   }
 } else{
  Rcpp::stop("Invalid kernel choice.");
 }
 int bmin = Rcpp::which_min(crit);
 return Rcpp::List::create(
   Rcpp::Named("lcv") = crit,
   Rcpp::Named("b") = b,
   Rcpp::Named("bandwidth") = b(bmin),
   Rcpp::Named("kernel") = kernel
 );
}

// [[Rcpp::export(kdens_Wishart)]]
Rcpp::NumericVector kdensWishart(
    const arma::cube &x,
    const arma::cube &pts,
    const double &b,
    bool log = true){
  arma::uword d = x.n_rows;// [[Rcpp::export]]
  arma::uword n = x.n_slices;
  arma::uword m = pts.n_slices;
  arma::vec logdet_x(n);
  arma::vec logdet_pts(m);
  // pts = b * pts;
  if((x.n_cols != d) | (pts.n_cols != d) | (pts.n_rows != d)){
   Rcpp::stop("Invalid dimensions for \"x\" or \"pts\".");
  }
  if(b <= 0){
   Rcpp::stop("Invalid bandwidth parameter: must be strictly positive.");
  }
  double bandwidth = 1.0 / b + d + 1.0;
  Rcpp::NumericVector logdens(n);
  arma::vec logdens_k(m);

  for(arma::uword i = 0; i < n; i++){
    logdet_x(i) = arma::log_det_sympd(x.slice(i));
  }
  for(arma::uword k = 0; k < m; k++){
    logdet_pts(k) = arma::log_det_sympd(pts.slice(k));
  }
  arma::mat Sinv(d,d);
  double logcst = 0;
  for(arma::uword i = 0; i < d; i++){
    logcst -= std::lgamma(0.5 * (bandwidth - i));
  }
  logcst -= 0.5 * bandwidth * d * (std::log(2.0) + std::log(b));
  logcst -= 0.25 * d * (d - 1.0) * std::log(arma::datum::pi);
  for(arma::uword k = 0; k < n; k++){
    logdens_k.zeros();
    Sinv = arma::inv_sympd(x.slice(k));
    for(arma::uword j = 0; j < m; j++){
        // TODO: finish this to avoid recomputing unnecessarily everything
        // Determinants can be recycled
        logdens_k(j) =  0.5 / b * (logdet_pts(j)  - arma::accu(Sinv % pts.slice(j))) + logcst - bandwidth * 0.5 * logdet_x(k);
    }
  logdens[k] = meanlog(logdens_k);
  }
  if(log){
    return logdens;
  } else{
    return Rcpp::exp(logdens);
  }
}


// NumericVector kdensWishart(
//     const arma::cube &x,
//     const arma::cube &pts,
//     const double &b,
//     bool log = true){
//   arma::uword d = x.n_rows;// [[Rcpp::export]]
//   arma::uword n = x.n_slices;
//   arma::uword m = pts.n_slices;
//   // pts = b * pts;
//   if((x.n_cols != d) | (pts.n_cols != d) | (pts.n_rows != d)){
//     Rcpp::stop("Invalid dimensions for \"x\" or \"pts\".");
//   }
//   if(b <= 0){
//     Rcpp::stop("Invalid bandwidth parameter: must be strictly positive.");
//   }
//   double band = 1.0 / b + d + 1.0;
//   Rcpp::NumericVector logdens(n);
//   arma::vec logdens_k(m);
//   for(arma::uword i = 0; i < n; i++){
//     logdens_k.zeros();
//     for(arma::uword j = 0; j < m; j++){
//       logdens_k(j) = dWishart_mat(
//         pts.slice(j),
//         band,
//         b * x.slice(i),
//         true);
//     }
//     logdens[i] = meanlog(logdens_k);
//   }
//   if(log){
//     return logdens;
//   } else{
//     return Rcpp::exp(logdens);
//   }
// }


// [[Rcpp::export(kdens_smlnorm)]]
Rcpp::NumericVector kdenssmlnorm(
    const arma::cube &x,
    const arma::cube &pts,
    double b,
    bool log = true){
  arma::uword d = x.n_rows;
  arma::uword n = x.n_slices;
  arma::uword m = pts.n_slices;
  if((x.n_cols != d) | (pts.n_cols != d) | (pts.n_rows != d)){
    Rcpp::stop("Invalid dimensions for \"x\" or \"pts\".");
  }
  Rcpp::NumericVector logdens(n);

  double logjac;
  arma::vec eigvals(d);
  arma::vec log_eigvals(d);

  double logcst;
  logcst = 0.25 * d * (d + 1.0) * std::log(2.0 * arma::datum::pi * b) - 0.25 * d * (d - 1.0) * std::log(2);

   // The Jacobian of the transformation is the same for all points
  for(arma::uword k = 0; k < n; k++){
  logjac = 0;
  // Eigenvalues are stored in ASCENDING order
  eigvals = arma::eig_sym(x.slice(k));
  log_eigvals = eigvals;
  log_eigvals = arma::log(log_eigvals);
  // PSD matrix have positive eigenvalues
  for(arma::uword i = 0; i < d - 1; i++){
    for(arma::uword j = i + 1; j < d; j++){
      if(log_eigvals(i) != log_eigvals(j)){
        logjac += std::log(log_eigvals(j) - log_eigvals(i))-
          std::log(eigvals(j) - eigvals(i));
      } else{
        logjac += log_eigvals(j);
      }
    }
  }
  logjac -= arma::accu(log_eigvals); //log_det_sympd(x)
  logdens[k] = logjac;
  }
  arma::vec logdens_k(m);
  arma::cube logmat_pts = pts;
  for(arma::uword j = 0; j < m; j++){
    logmat_pts.slice(j) = arma::logmat_sympd(pts.slice(j));
  }
  arma::mat logmat_x_i(d,d);
  for(arma::uword i = 0; i < n; i++){
    logdens_k.zeros();
    logmat_x_i = arma::logmat_sympd(x.slice(i));
    for(arma::uword j = 0; j < m; j++){
      logdens_k(j) = - accu(arma::pow(logmat_x_i - logmat_pts.slice(j), 2.0)) / (2.0 * b);
    }
    logdens[i] = logdens[i] + meanlog(logdens_k) - logcst;
  }
  if(log){
    return logdens;
  } else{
    return Rcpp::exp(logdens);
  }
}
// LSCV
// Numerical integration

//
//
// NumericVector kdenssmlnormback(
//     const arma::cube &x,
//     const arma::cube &pts,
//     double b,
//     bool log = true){
//   arma::uword d = x.n_rows;
//   arma::uword n = x.n_slices;
//   arma::uword m = pts.n_slices;
//   if((x.n_cols != d) | (pts.n_cols != d) | (pts.n_rows != d)){
//     Rcpp::stop("Invalid dimensions for \"x\" or \"pts\".");
//   }
//   Rcpp::NumericVector logdens(n);
//   arma::vec logdens_k(m);
//   arma::cube logmat_pts = pts;
//   for(arma::uword j = 0; j < m; j++){
//     logmat_pts.slice(j) = arma::logmat_sympd(pts.slice(j));
//   }
//   arma::mat logmat_x_i(d,d);
//   for(arma::uword i = 0; i < n; i++){
//     logdens_k.zeros();
//     logmat_x_i = arma::logmat_sympd(x.slice(i));
//     for(arma::uword j = 0; j < m; j++){
//       logdens_k(j) = dsmlnorm_mat(
//         x.slice(i),
//         logmat_x_i,
//         b,
//         pts.slice(j),
//         logmat_pts.slice(j),
//         true);
//     }
//     logdens[i] = meanlog(logdens_k);
//   }
//   if(log){
//     return logdens;
//   } else{
//     return Rcpp::exp(logdens);
//   }
// }

//' Least square cross validation criterion for Wishart kernel
//'
//' Finite sample h-block leave-one-out approximation to the least
//' square criterion, omitting constant term.
//' @inheritParams lcv_kern_Wishart
//' @export
//' @param h separation vector; only pairs that are \eqn{|i-j| \leq h} apart are considered
//' @return a vector of length two containing the log of the summands
// [[Rcpp::export(lscv_kern_Wishart)]]
arma::vec lscvkernwishart(
    const arma::cube &x,
    const double &b,
    const int &h = 1){
 arma::uword n = x.n_slices;
 arma::uword d = x.n_rows;
 double bandwidth = 1.0 / b + d + 1.0;
 double rd = 0.5 * (1.0 + d) * d;
 if(b <= 0){
   Rcpp::stop("Invalid bandwidth parameter: must be strictly positive.");
 }
 arma::vec obj(2);
 double log2 = std::log(2);
 int nh;
 int nt = 0;
 arma::vec logkern(n);
 arma::vec logdets(n);
 arma::vec sumlogdets(n);
 arma::vec sumlogkern(n);
 arma::vec logdet_x(n);

 arma::mat Sinv(d,d);
 double logcst2 = 0;
 for(arma::uword i = 0; i < d; i++){
   logcst2 -= std::lgamma(0.5 * (bandwidth - i));
 }
 logcst2 -= 0.5 * bandwidth * d * (std::log(2.0) + std::log(b));
 logcst2 -= 0.25 * d * (d - 1.0) * std::log(arma::datum::pi);


 // Temporary solution (perhaps error-prone)
 // We average every ith round, then multiply by the number of summands
 // Keep everything on the log scale
 // Some discrepancy with the formula
 // TODO avoid computing every determinant twice needlessly
  double bcst = 1.0 / b + 0.5 * (1.0 + d);
 double logcst =  - 2.0 * std::log(n) +
   - rd * (std::log(b) + log2) +
   lmgamma(1.0 / b + 0.5 * (d + 1.0), d) -
   2.0 * lmgamma(0.5 * bandwidth, d);
 for(arma::uword i = 0; i < n; i++){
  logdet_x(i) =  log_det_sympd(x.slice(i));
 }
 for(arma::uword i = 0; i < n; i++){
   nh = 0;
   logkern.zeros();
   logdets.zeros();
   logdets(i) = logdet_x(i) / b -
     bcst * (log2 * d  + logdet_x(i));
   if(i > 0){
    for(arma::uword j = 0; j < i; j++){
      // This thing is symmetric! no need to compute off-diagonal entries twice
      // Likewise, det(AB) = det(A)det(B)
      // and A + B is positive definite if both A and B are pd
      logdets(j) = log2 + 0.5 / b * (logdet_x(i) + logdet_x(j)) -
        bcst * log_det_sympd(x.slice(i) + x.slice(j));
    }
   }
   Sinv = arma::inv_sympd(x.slice(i));
  for(arma::uword j = 0; j < n; j++){
    if(std::abs((int) (j-i)) >= (int) h){
      // The log kernel contribution isn't symmetric, however
     nt++;

      // TODO speed this up
     logkern(nh) =  0.5 / b * (logdet_x(j)  - arma::accu(Sinv % x.slice(j))) - bandwidth * 0.5 * logdet_x(i);
     nh++;
    }
  }
 // Compute partial sum with precision, but using the log trick
    sumlogdets(i) = sumlog(logdets.rows(0, i));
    sumlogkern(i) = sumlog(logkern.rows(0, nh - 1)) + logcst2;
  }
  //Rcpp::LogicalVector sgn = {1,0};
  obj(0) = sumlog(sumlogdets) + logcst;
  // cout << nt << "and difference" << n * (n - h) << std::endl;
  obj(1) = sumlog(sumlogkern) + log2 - std::log(nt); // nt = n * (n-h)
  //double out = sumsignedlog(obj, sgn);
  return obj;
 }



//' Least square cross validation criterion for log symmetric matrix normal kernel
//'
//' Finite sample h-block leave-one-out approximation to the least
//' square criterion, omitting constant term.
//' @inheritParams lcv_kern_Wishart
//' @export
//' @param h separation vector; only pairs that are \eqn{|i-j| \leq h}
//'  apart are considered
// [[Rcpp::export(lscv_kern_smlnorm)]]
 arma::vec lscvkernsmlnorm(
     const arma::cube &x,
     const double &b,
     const int &h = 1){
   if(b <= 0){
     Rcpp::stop("Invalid bandwidth parameter: must be strictly positive.");
   }
   arma::uword n = x.n_slices;
   arma::uword d = x.n_rows;
   double rd = 0.5 * (1.0 + d) * d;
   arma::cube logmat_x(d,d,n);
   arma::vec obj(2);
   double log2 = std::log(2);
   int nh;
   int nt = 0;
   arma::vec logkern(n);
   arma::vec logetr(n);
   arma::vec sumlogetr(n);
   arma::vec sumlogkern(n);
   arma::vec trmatlog(n);
   arma::mat logmatsum(d,d);
   // Save on calculations by saving some elements that are recycled
   for(arma::uword i = 0; i < n; i++){
     logmat_x.slice(i) = arma::logmat_sympd(x.slice(i));
     trmatlog(i) = arma::accu(arma::pow(logmat_x.slice(i), 2.0));
   }
   double logcst =
     0.5 * rd * (log2 + std::log(arma::datum::pi) + std::log(b)) +
     0.5 * d * log2 + 2 * std::log(n);

   double logcst2;
   logcst2 = 0.25 * d * (d + 1.0) * std::log(2.0 * arma::datum::pi * b) - 0.25 * d * (d - 1.0) * std::log(2);
     for(arma::uword i = 0; i < n; i++){
     nh = 0;
     logkern.zeros();
     logetr.zeros();
      // Diagonal terms (i, i) are exactly zero
     if(i > 0){
       for(arma::uword j = 0; j < i; j++){
         // This thing is symmetric! no need to compute off-diagonal entries twice
         logmatsum = logmat_x.slice(i) + logmat_x.slice(j);
         logetr(j) = log2 + 0.5 / b * (0.5 * arma::accu(arma::pow(logmatsum, 2.0)) - trmatlog(i) - trmatlog(j)); // multiplied by 2 because symmetric
       }
     }
     for(arma::uword j = 0; j < n; j++){
       if(std::abs((int) (j-i)) >= (int) h){
         nt++;
         logkern(nh) = - accu(arma::pow(logmat_x.slice(j) - logmat_x.slice(i), 2.0)) / (2.0 * b) - logcst2;
         nh++;
       }
     }
     // Compute partial sum with precision, but using the log trick
     sumlogetr(i) = sumlog(logetr.rows(0, i));
     sumlogkern(i) = sumlog(logkern.rows(0, nh - 1));
   }
   //Rcpp::LogicalVector sgn = {1,0};
   obj(0) = sumlog(sumlogetr) - logcst;
   obj(1) = sumlog(sumlogkern) + log2 - std::log(nt);
   //double out = sumsignedlog(obj, sgn);
   return obj;
 }


//' Random vector generation from the multivariate normal distribution
//'
//' Sampler derived using the eigendecomposition of the covariance
//' matrix \code{vcov}.
//'
//' @param n sample size
//' @param mean mean vector of length \code{d}
//' @param vcov a square positive definite covariance matrix, of the same dimension as \code{mean}.
//' @export
//' @return an \code{n} by \code{d} matrix of samples
//' @examples
//' rmnorm(n = 10, mean = c(0, 2), vcov = diag(2))
// [[Rcpp::export]]
Rcpp::NumericMatrix rmnorm(
    int n,
    const arma::vec &mean,
    const arma::mat &vcov){
   if ((vcov.n_rows != vcov.n_cols) | (mean.n_elem != vcov.n_cols)){
     Rcpp::stop("Incompatible arguments");
   }
   int d = (int) vcov.n_rows;
   arma::mat samp(n, d, arma::fill::randn);
   arma::vec eigval;
   arma::mat eigvec;
   //Covariance matrix must be symmetric - otherwise eig_sym throws error
   arma::eig_sym(eigval, eigvec, vcov);
   samp = samp*arma::diagmat(arma::sqrt(eigval))*trans(eigvec);
   for(int i = 0; i < d; i++){
     samp.col(i) += mean(i);
   }
   return Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(samp));
 }
