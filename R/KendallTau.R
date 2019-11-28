#' @title Kendall's tau of a copula
#'
#' @description This function computes the Kendall's tau of a copula family for a given a unconstrainted parameter alpha.
#'
#' @param family "Gaussian" , "t" , "Clayton" , "Frank" , "Gumbel"
#' @param alpha  unconstrainted parameters of the copula family
#'
#' @return \item{tau}{estimated Kendall's tau}
#' @return \item{theta}{estimated copula parameter (constrained)}
#'
#' @author Bouchra R. Nasri, August 14, 2019
#'
#' @examples KendallTau('Clayton',0)
#'
#' @export
#'
KendallTau<- function(family,alpha){

  theta = ParamCop(family,alpha);

  switch(family,
         "Gaussian" = {

         tau = copula::tau(normalCopula(theta[1]));
         },

         "t" = {
           tau = copula::tau(tCopula(theta[1]));
         },

         "Clayton" = {

           tau = copula::tau(claytonCopula(theta[1]));
         },

         "Frank" = {


           tau = copula::tau(frankCopula(theta[1]));
         },

         "Gumbel" = {

           tau = copula::tau(gumbelCopula(theta[1]));
         }
  )

  out=list(tau=tau,theta=theta)
}




