#' chebycenter
#' This code is a translation of the matlab code developped by Tim Bernham - University of Queensland
#' (https://ch.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope
#' It computes the centroid of the complex polytope defined by \eqn{A \cdot x \leqslant   b}
#'
#' @param A a matrix
#' @param b a vector of length equals to nrow(A)
#'
#' @return a vector corresponding to the centroid of the polytope
#' @importFrom linprog solveLP
#' @examples
#' n <- 20
#' A1 <- -diag(n)
#' b1 <- as.matrix(rep(0,n))
#' A2 <- diag(n)
#' b2 <- as.matrix(rep(1,n))
#' A <- rbind(A1,A2)
#' b <- rbind(b1,b2)
#' X0 <- chebycenter(A,b)
#' @export

chebycenter <- function(A, b) {
  n <- dim(A)[1]
  p <- dim(A)[2]

  an <- rowSums(A ^ 2) ^ 0.5
  A1 <- matrix(data = 0,
               nrow = n,
               ncol = p + 1)
  A1[, 1:p] <- A

  A1[, p + 1] <- an

  f <- matrix(data = 0,
              nrow = p + 1,
              ncol = 1)
  f[p + 1] <- -1



  d <-
    linprog::solveLP(
      cvec = f,
      bvec = as.numeric(b),
      Amat = A1,
      const.dir = rep("<=", n)
    )
  x <- d$solution
  return(x[-p - 1])
}
