#' @title Initial values for estimation
#' @description This function computes initial values of non-centrality parameters and Kendall's tau
#'at selected points for the estimation non-central squared copula parameters. The results are not satisfactory. Do not use.
#'
#' @param U         (nx2) data matrix of pseudo-observations.
#' @param family    'Gaussian' , 't' , 'Clayton' , 'Frank' , 'Gumbel'.
#'
#' @return \item{paraml}{Initial  values for the non-centrality parameters and Kendall's tau to be
#'included in the EstNCSCop function.}
#'
#' @author Bouchra R. Nasri, August 14, 2019
#'
#'
#'
#' @examples param <- c(0.8, 2.5, 0.7) ;
#'U <- SimNCSCop('Clayton', 250, param)
#' param = initialValues(U, 'Clayton');
#'
#'
#' @export



 initialValues <- function(U, family='Clayton'){


n = dim(U)[1];


l=0;

u = matrix(c(0.05, 0.15, 0.15, 0.25, 0.25, 0.35, 0.35, 0.45, 0.45 ,0.55,0.55, 0.65, 0.65, 0.75, 0.75, 0.85, 0.85, 0.95),ncol=2,byrow=TRUE);

val = matrix(0,ncol=3,nrow=300);
cdf = matrix(0,ncol=300,nrow=9);
for(i in 1:10){
    for( j in 1:10){
        for(k in 1:3){
        l=l+1;
        a1 = i/4;
        a2 = j/4;
        tau = k/4;
        val[l,] = c(a1, a2, tau);
        cdf[,l] = NCSCopCdf(u, family, val[l,]);
        }
    }
}

cdf0 =  copulaEmp(u,U);

x = abs(cdf0-cdf);

bias = apply(x,2,max);
l=which.min(bias)

param = val[l,];
return(param)
}
