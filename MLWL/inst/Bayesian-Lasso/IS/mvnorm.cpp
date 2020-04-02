#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericVector dmvnorm_cholesky_inverse(const NumericMatrix & x, const NumericVector & mean, const Eigen::MatrixXd & cholesky_inverse){
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::MatrixXd xcentered(x_);
  double halflogdeterminant = cholesky_inverse.diagonal().array().log().sum();;
  double cst = - (-halflogdeterminant) - (x.cols() * 0.9189385);
  for(int j = 0; j < x.cols(); j++){
    for(int i = 0; i < x.rows(); i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = -0.5 * (cholesky_inverse.transpose() * xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = results(i) + cst;
  }
  return wrap(results);
}
