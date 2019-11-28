#' @title Gives the parameters of the copula family
#'
#' @description This function computes the parameter of the copula according to CRAN copula package where corresponding to the unconstrainted parameters alpha.
#'
#' @param family "Gaussian" , "t" , "Clayton" , "Frank" , "Gumbel"
#' @param alpha  unconstrainted parameters of the copula family
#' @author Bouchra R. Nasri, August 14, 2019
#'@return \item{theta}{Bivariate vector of constrained copula family parameters}
#'
#' @examples ParamCop('Clayton',0)
#'
#'@export
#'
ParamCop<- function(family,alpha){

  theta = c(0,1);
  theta[2] = NaN;
  switch(family,
         "Gaussian" = {
           theta[1] =  2 /(1+exp(-alpha[1]) )-1;
          },

         "t" = {
           theta[1] = 2 /(1+exp(-alpha[1])) -1;
           theta[2]   = 25/(1+exp(-alpha[2]));

         },

         "Clayton" = {
           theta[1] = exp(alpha[1]);

         },

         "Frank" = {
           theta[1] = alpha[1];

         },

         "Gumbel" = {
           theta[1] = exp(alpha[1])+1;

         }
  )

  return(theta)
}





