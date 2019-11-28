#' @title Empirical copula
#'
#' @description This function computes the  empirical bivariate copula at a series of points.
#'
#' @param u         (nx2) data matrix of points.
#' @param U         (nx2) data matrix of pseudo-observations.
#'
#' @return \item{cdf}{Empirical copula values at u.}
#'
#' @author Bouchra R. Nasri, August 14, 2019
#'
#'
#'
#' @examples param <- c(0.8, 2.5, 0.7) ;
#' U <- SimNCSCop('Clayton', 250, param)
#' u = matrix(c(0.2,0.6,0.3,0.5,0.7,0.9),ncol=2,byrow=TRUE);
#' cdf=copulaEmp(u,U);
#'
#' @export

copulaEmp <- function(u,U)
{



  m = dim(u)[1];
  cdf = 0*c(1:m);


  U1 = U[,1];
  U2 = U[,2];

  for (j in 1:m){
    cdf[j] = mean ((U1<u[j,1]) * (U2<u[j,2]));
  }

  return(cdf)
}





