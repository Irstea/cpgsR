#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>
#include "cpgsEquality.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace RcppEigen;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::FullPivLU;



//' Complex Polytope Gibbs Sampling
//' This function draw uniform samples in a convex polytope with both equality and inequality constraints
//'
//' @param N the number of samples to generate
//' @param A a matrix of coefficients of inequality constants A.x<=b
//' @param b a vector of length equals to nrow(A)
//' @param C a matrix of coefficients of inequality constants C.x=v
//' @param v a vector of length equals to nrow(C)
//' @param x0 a vector of length equals to ncol(A) that should be in the polytope, for example returned by \code{\link{chebycenter}}
//'
//' @section Details:
//' This function is based on an initial matlab code developped called CPRND
//' (https://ch.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope)
//' It generates samples within the complex polytope defined by \eqn{A \cdot x \leqslant   b}
//'
//' @return a matrix with one row per sample and one column per parameter
//' @examples
//' n <- 20
//' A1 <- -diag(n)
//' b1 <- as.matrix(rep(0,n))
//' A2 <- diag(n)
//' b2 <- as.matrix(rep(1,n))
//' A <- rbind(A1,A2)
//' b <- rbind(b1,b2)
//' C <- rbind(c(1,1,rep(0,n-2)),c(0,0,1,1,rep(0,n-4)))
//' v <- matrix(rep(0.2,2),2)
//' X0 <- rep(0.1,n)
//' x <- cpgsEquality(1000,A,b,C,v,X0)
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd cpgsEquality(const int N,const Eigen::MatrixXd &A ,const Eigen::VectorXd &b,const Eigen::MatrixXd &C ,const Eigen::VectorXd &v,const Eigen::VectorXd &x0){
  int p=A.cols();
  int m=A.rows();
  int p2=C.cols();
  int m2=C.rows();
  int discard;

  VectorXd X(N);

  // Check input arguments
  if (m < (p+1) || b.size()!=m || x0.size()!=p){
    throw std::range_error("dimensions mismatch");
  }
  if (v.size()!=m2 || x0.size()!=p2){
    throw std::range_error("dimensions mismatch");
  }
  // Initialisation
  FullPivLU<MatrixXd> lu(C);
  MatrixXd Nt = lu.kernel().transpose();

  MatrixXd Abis=A*N;
  VectorXd bbis=b-A*x0;

  VectorXd x0bis=VectorXd::Zero(Nt.cols());
  MatrixXd x=cpgs(N,Abis,bbis,x0bis);
  for(int i=0;i<N;++i) X.row(i)=Nt*x.row(i)+x0;
  return X;
}
