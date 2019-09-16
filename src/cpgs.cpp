#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>


using namespace Rcpp;
using namespace RcppEigen;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

//' Complex Polytope Gibbs Sampling
//' This function is based on an initial matlab code developped by Tim Bernham - University of Queensland
//' (https://ch.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope
//' It generates samples within the complex polytope defined by \eqn{A \cdot x \leqslant   b}
//' @param N the number of samples to generate
//' @param A a matrix
//' @param b a vector of length equals to nrow(A)
//' @param x0 a vector of length equals to nrcol(A) that should be in the polytope, for example returned by \code{\link{chebycenter}}
//' @param verbose flag to display progress bar
//'
//' @return a matrix with one row per sample and one column per parameter
//' @example
//' n <- 20
//' A1 <- -diag(n)
//' b1 <- as.matrix(rep(0,n))
//' A2 <- diag(n)
//' b2 <- as.matrix(rep(1,n))
//' A <- rbind(A1,A2)
//' b <- rbind(b1,b2)
//' X0 <- chebycenter(A,b)
//' x <- cpgs(1000,A,b,X0)
//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpgs(const int N,const Eigen::MatrixXd &A ,const Eigen::VectorXd &b,const Eigen::VectorXd &x0, const bool verbose=false) {
  int p=A.cols();
  int m=A.rows();
  int discard;

  // Check input arguments
  if (m < (p+1) || b.size()!=m || x0.size()!=p){
    throw std::range_error("dimensions mismatch");
  }
  // Initialisation
  Eigen::MatrixXd X(N,p);
  int n=0;
  Eigen::MatrixXd x(p,1);
  Eigen::MatrixXd y(p,1);
  x=x0;

  // Initialize variables for keeping track of sample mean, covariance
  // and isotropic transform.
  Eigen::MatrixXd M(p,1);
  M.setZero();
  Eigen::MatrixXd S2(p,p);
  S2.setZero();

  // outer products.
  Eigen::MatrixXd S(p,p);
  Eigen::MatrixXd S0(p,p);
  S.setIdentity();

  IntegerVector index=Rcpp::seq(0,p-1);

  Eigen::MatrixXd T1(p,p);
  Eigen::MatrixXd T2(p,p);
  T1.setIdentity();
  Eigen::MatrixXd W(m,p);

  W = A;
  bool adapt=true;
  Eigen::MatrixXd d(m,1);
  Eigen::MatrixXd delta0(p,1);
  Eigen::MatrixXd delta1(p,1);
  Eigen::MatrixXd z(m,1);
  Eigen::MatrixXd L(p,p);
  Eigen::MatrixXd D(p,p);
  Eigen::VectorXd Dtmp(p);
  Eigen::VectorXd Dzero=VectorXd::Constant(p,1.0e-16);
  Eigen::LDLT<MatrixXd> ldltOfS(S.cols());
  int runup=0;

  while ((adapt==true) || (n < (N+discard))){               //sampling loop
    //std::random_shuffle(index.begin(), index.end()); //we change the order to limit the influence of initial ordering
    index=sample(index,p,false);
    y=x;
    NumericVector alea2=runif(p);
    // compute approximate stochastic transformation
    if ((!adapt) && n == 0){ //first true sample, we now make the isotropic transformation
      ldltOfS.compute(S.transpose());
      D=ldltOfS.vectorD().cwiseMax(Dzero).asDiagonal();
      L=ldltOfS.matrixL();
      T1=ldltOfS.transpositionsP().transpose()*L*D.sqrt();
      T2=T1.inverse();
      W = A*T1;
    }
    if (!adapt) y=T2*y; //otherwise y=I^-1 * y=y
    if (n>0 && n%50==0 && verbose) Rcpp::Rcout<<"sampling iter "<<n<<"/"<<discard+N<<std::endl;

    // choose p new components
    Eigen::VectorXd e(p);
    for (int ip=0;ip<p;++ip){
      int i=index[ip];
      //Find points where the line with the (p-1) components x_i
      //fixed intersects the bounding polytope.
      e.setZero();
      e(i)= 1;
      z = (W*e); //prevent any divisions by 0
      d=(b - W*y);
      d=d.cwiseQuotient(z);
      double tmin=-9e9;
      double tmax=9e9;
      for (int j=0;j<m;++j){
        if (z(j)<0 && tmin<d(j)) tmin=d(j);
        if (z(j)>0 && tmax>d(j)) tmax=d(j);
      }
      y(i)+=(tmin+(tmax-tmin)*alea2(i));
    }
    x=T1*y;
    if (!adapt){
      if (n>=discard){
        X.row(n-discard)= x.col(0);
      }
      ++n;
    }
    if (adapt && runup>1){  //part in the adapation phase to tune the isotropic transformation
      ++runup;
      // Incremental mean and covariance updates
      delta0 = x - M; // delta new point wrt old mean
      M+= delta0/(double)runup;     // sample mean
      delta1= x - M;      // delta new point wrt new mean

      S2 +=(runup-1)/(double)(runup*runup)*(delta0*delta0.transpose())+(delta1*delta1.transpose());
      S0 = S;
      S = S2/(double)(runup-1);           // sample covariance
      if ((S0-S).norm()/S0.norm()<0.05 && runup>=p) {
        adapt=false; //the covariance matrix is stable, adaptation stage is ok
        discard=runup;
        if (verbose) Rcpp::Rcout<<"end of adaptation phase after "<<runup<<" iterations"<<std::endl;
      } else{
        if (verbose) Rcpp::Rcout<<"adaptation phase iter "<<runup<<" - score "<<(S0-S).norm()/S0.norm()<<std::endl;
      }
    } else if (adapt) {
      ++runup;
    }
  }
  return X;
}

