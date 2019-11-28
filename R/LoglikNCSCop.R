#' @title Log-likelihood of a non-central squared copula
#'
#' @description This function computes the log-likelihood vector of a non-central squared copula
#'
#' @param alpha     unconstrained non-centrality parameters a1, a2, and unconstrained copula parameters.
#' @param U         (nx2) data matrix of pseudo-observations.
#' @param family    'Gaussian' , 't' , 'Clayton' , 'Frank' , 'Gumbel'.
#' @param p         number of different non-centrality parameters (0,1,2 default).
#'
#' @author Bouchra R. Nasri, August 14, 2019
#' @return \item{LL}{Vector of log-likelihoods}
#'
#'
#'
#' @examples alpha = c(log(0.2),log(5),log(2),log(12));
#'param = c(0.5,2.5,0.5);
#'data = SimNCSCop('Clayton', 250, param);
#'LL = LoglikNCSCop(alpha,data,'Clayton')
#'
#'
#'@export

LoglikNCSCop <- function(alpha,U,family, p=2)
{

  theta = ParamCop(family,alpha[(p+1):(p+2)]);
  if(p==1){
    a1 = 3/(1+exp(-alpha[1]));
    a2 = a1;
    }   else if (p==2){
    a1 = 3/(1+exp(-alpha[1]));
    a2 = 3/(1+exp(-alpha[2]));
    }

if(p>0){
  d1 = a1^2;
  d2 = a2^2;

  U1 = U[,1];
  U2 = U[,2];


  G1 = sqrt(qchisq(U1,1,d1));  #Ga inverse
  G2 = sqrt(qchisq(U2,1,d2));

  ha = matrix(c((G1-a1),(-G1-a1) , (G2-a2), (-G2-a2)), ncol=4, byrow=FALSE);

  htilde = pnorm(ha)+ 1E-10;
  w0     = dnorm(ha);

  U11 = matrix(c(htilde[,1], htilde[,3]), ncol=2,byrow=FALSE);
  U12 = matrix(c(htilde[,1], htilde[,4]), ncol=2,byrow=FALSE);
  U21 = matrix(c(htilde[,2], htilde[,3]), ncol=2,byrow=FALSE);
  U22 = matrix(c(htilde[,2], htilde[,4]), ncol=2,byrow=FALSE);
  w0011  = w0[,1] /(w0[,1]+w0[,2]);
  w0012  = 1-w0011;
  w0021 = w0[,3] /(w0[,3]+w0[,4]);
  w0022  = 1-w0021;

  w11 = w0011 * w0021;
  w12 = w0011 * w0022;
  w21 = w0012 * w0021;
  w22 = w0012 * w0022;

  switch(family,
         "Gaussian" =

         {
           f11 = copula::dCopula(U11,copula::normalCopula(theta[1], dim = 2))+ 1E-20;
           f12 = copula::dCopula(U12,copula::normalCopula(theta[1], dim = 2))+ 1E-20;
           f21 = copula::dCopula(U21,copula::normalCopula(theta[1], dim = 2))+ 1E-20;
           f22 = copula::dCopula(U22,copula::normalCopula(theta[1], dim = 2))+ 1E-20;

         },

         "t" = {
           f11 = copula::dCopula(U11, copula::tCopula(theta[1], dim = 2,df = theta[2]))+ 1E-20;
           f12 = copula::dCopula(U12, copula::tCopula(theta[1], dim = 2,df = theta[2]))+ 1E-20;
           f21 = copula::dCopula(U21, copula::tCopula(theta[1], dim = 2,df = theta[2]))+ 1E-20;
           f22 = copula::dCopula(U22, copula::tCopula(theta[1], dim = 2,df = theta[2]))+ 1E-20;
         },

         "Clayton" = {
           f11 = copula::dCopula(U11,copula::claytonCopula(theta[1], dim = 2))+ 1E-20;
           f12 = copula::dCopula(U12,copula::claytonCopula(theta[1], dim = 2))+ 1E-20;
           f21 = copula::dCopula(U21,copula::claytonCopula(theta[1], dim = 2))+ 1E-20;
           f22 = copula::dCopula(U22,copula::claytonCopula(theta[1], dim = 2))+ 1E-20;
         },

         "Frank" = {
           f11 = copula::dCopula(U11,copula::frankCopula(theta[1], dim = 2))+ 1E-20;
           f12 = copula::dCopula(U12,copula::frankCopula(theta[1], dim = 2))+ 1E-20;
           f21 = copula::dCopula(U21,copula::frankCopula(theta[1], dim = 2))+ 1E-20;
           f22 = copula::dCopula(U22,copula::frankCopula(theta[1], dim = 2))+ 1E-20;
         },

         "Gumbel" = {
           f11 = copula::dCopula(U11,copula::gumbelCopula(theta[1], dim = 2))+ 1E-20;
           f12 = copula::dCopula(U12,copula::gumbelCopula(theta[1], dim = 2))+ 1E-20;
           f21 = copula::dCopula(U21,copula::gumbelCopula(theta[1], dim = 2))+ 1E-20;
           f22 = copula::dCopula(U22,copula::gumbelCopula(theta[1], dim = 2))+ 1E-20;
         }
  )


  LL = log(w11 * f11 + w12 * f12 + w21 * f21 + w22 * f22);
}
else{

 switch(family,
       "Gaussian" =

       {
         LL = log(copula::dCopula(U,copula::normalCopula(theta[1], dim = 2))+ 1E-20);

       },

       "t" = {
         LL = log(copula::dCopula(U, copula::tCopula(theta[1], dim = 2,df = theta[2]))+ 1E-20);
       },

       "Clayton" = {
         LL = log(copula::dCopula(U,copula::claytonCopula(theta[1], dim = 2))+ 1E-20);
       },

       "Frank" = {
         LL = log(copula::dCopula(U,copula::frankCopula(theta[1], dim = 2))+ 1E-20);
       },

       "Gumbel" = {
         LL = log(copula::dCopula(U,copula::gumbelCopula(theta[1], dim = 2))+ 1E-20);
       }
    )
}
  return(LL)
}
