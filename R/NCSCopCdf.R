#' @title Distribution function of a non-central squared copula
#'
#' @description This function computes the distribution function a non-central squared copula
#'
#' @param u         (nx2) data matrix of pseudo-observations.
#' @param family    'Gaussian' , 't' , 'Clayton' , 'Frank' , 'Gumbel'.
#' @param param   c(a1,a2,tau)  where  a1,a2  are the non-negative non-centrality
#                     parameters, tau is Kendall's tau of the copula family.
#' @param dof     degrees of freedom of the Student copula (if needed).
#' @return \item{cdf}{Non-central squared copula evaluated at points u.}
#' @author Bouchra R. Nasri, August 14, 2019

#'
#'
#'
#' @examples param = c(0.8,2.5,0.7);
#' u = matrix(c(0.2,0.6,0.3,0.5,0.7,0.9),ncol=2,byrow=TRUE);
#' cdf=NCSCopCdf(u,'Clayton',param);
#'
#' @export

NCSCopCdf <- function(u,family, param,dof=NULL)
{

  a1  = param[1];
  a2  = param[2];
  KendallTau= param[3];

  switch(family,
         "Gaussian" = {
           theta =  copula::iTau(copula::normalCopula(),  tau = KendallTau)
         },

         "t" = {
           theta =  copula::iTau(copula::normalCopula(),  tau = KendallTau)
         },

         "Clayton" = {
           theta =  copula::iTau(copula::claytonCopula(),  tau = KendallTau)
         },

         "Frank" = {
           theta =  copula::iTau(copula::frankCopula(),  tau = KendallTau)
         },

         "Gumbel" = {
           theta =  copula::iTau(copula::gumbelCopula(),  tau = KendallTau)
         }
  )


  d1 = a1^2;
  d2 = a2^2;

  U1 = u[,1];
  U2 = u[,2];


  G1 = sqrt(qchisq(U1,1,d1));  #Ga inverse
  G2 = sqrt(qchisq(U2,1,d2));

  ha = matrix(c((G1-a1),(-G1-a1) , (G2-a2), (-G2-a2)), ncol=4, byrow=FALSE);

  htilde = pnorm(ha)+ 1E-10;
  w0     = dnorm(ha);

  U11 = matrix(c(htilde[,1], htilde[,3]), ncol=2,byrow=FALSE);
  U12 = matrix(c(htilde[,1], htilde[,4]), ncol=2,byrow=FALSE);
  U21 = matrix(c(htilde[,2], htilde[,3]), ncol=2,byrow=FALSE);
  U22 = matrix(c(htilde[,2], htilde[,4]), ncol=2,byrow=FALSE);


  switch(family,
         "Gaussian" =
         {
           F11 = copula::pCopula(U11,copula::normalCopula(theta, dim = 2))+ 1E-20;
           F12 = copula::pCopula(U12,copula::normalCopula(theta, dim = 2))+ 1E-20;
           F21 = copula::pCopula(U21,copula::normalCopula(theta, dim = 2))+ 1E-20;
           F22 = copula::pCopula(U22,copula::normalCopula(theta, dim = 2))+ 1E-20;

         },

         "t" = {
           F11 = copula::pCopula(U11, copula::tCopula(theta, dim = 2,df = dof))+ 1E-20;
           F12 = copula::pCopula(U12, copula::tCopula(theta, dim = 2,df = dof))+ 1E-20;
           F21 = copula::pCopula(U21, copula::tCopula(theta, dim = 2,df = dof))+ 1E-20;
           F22 = copula::pCopula(U22, copula::tCopula(theta, dim = 2,df = dof))+ 1E-20;
         },

         "Clayton" = {
           F11 = copula::pCopula(U11,copula::claytonCopula(theta, dim = 2))+ 1E-20;
           F12 = copula::pCopula(U12,copula::claytonCopula(theta, dim = 2))+ 1E-20;
           F21 = copula::pCopula(U21,copula::claytonCopula(theta, dim = 2))+ 1E-20;
           F22 = copula::pCopula(U22,copula::claytonCopula(theta, dim = 2))+ 1E-20;
         },

         "Frank" = {
           F11 = copula::pCopula(U11,copula::frankCopula(theta, dim = 2))+ 1E-20;
           F12 = copula::pCopula(U12,copula::frankCopula(theta, dim = 2))+ 1E-20;
           F21 = copula::pCopula(U21,copula::frankCopula(theta, dim = 2))+ 1E-20;
           F22 = copula::pCopula(U22,copula::frankCopula(theta, dim = 2))+ 1E-20;
         },

         "Gumbel" = {
           F11 = copula::pCopula(U11,copula::gumbelCopula(theta, dim = 2))+ 1E-20;
           F12 = copula::pCopula(U12,copula::gumbelCopula(theta, dim = 2))+ 1E-20;
           F21 = copula::pCopula(U21,copula::gumbelCopula(theta, dim = 2))+ 1E-20;
           F22 = copula::pCopula(U22,copula::gumbelCopula(theta, dim = 2))+ 1E-20;
         }
  )


  cdf = F11-F12-F21+F22;


  return(cdf)
}
