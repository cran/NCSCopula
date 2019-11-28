#' @title Unconstrained parameters
#'
#' @description This function computes the  unconstrainted parameter alpha for a given Kendall's tau
#'
#' @param family    'Gaussian' , 't' , 'Clayton' , 'Frank' , 'Gumbel'
#' @param tau        Kendall's tau of the copula family
#'
#' @author Bouchra R. Nasri, August 14, 2019
#' @return \item{alpha}{Unconstrainted parameter}
#'
#' @examples ParamTau('Clayton',0.5)
#'
#' @export
ParamTau<- function(family,tau){

  switch(family,
         "Gaussian" = {
           theta =  copula::iTau(copula::normalCopula(),  tau) ;
           alpha = log( (1+theta) / (1-theta) );
         },

         "t" = {
           theta =  copula::iTau(copula::normalCopula(), tau) ;
           alpha = log( (1+theta) / (1-theta) );
         },

         "Clayton" = {
           theta =  copula::iTau(copula::claytonCopula(),  tau);
           alpha = log(theta);
         },

         "Frank" = {
           theta =  copula::iTau(copula::frankCopula(),  tau);
           alpha=theta;
         },

         "Gumbel" = {
           theta =  copula::iTau(copula::gumbelCopula(),  tau );
           alpha = log(theta-1);
         }
  )

  return(alpha)
}
