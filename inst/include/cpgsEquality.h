#ifndef CPGSEQUALITY_H
#define CPGSEQUALITY_H
#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::MatrixXd cpgsEquality(const int N,const Eigen::MatrixXd &A ,const Eigen::VectorXd &b,const Eigen::MatrixXd &C ,const Eigen::VectorXd &v,const Eigen::VectorXd &x0);


#endif
