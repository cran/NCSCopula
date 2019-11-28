#' @title Simulation of a bivariate non-central squared copula
#'
#' @description This function simulates observations a bivariate non-central squared copula model.
#'
#' @param family    'Gaussian' , 't' , 'Clayton' , 'Frank' , 'Gumbel'.
#' @param n    number of simulated vectors.
#' @param param   c(a1,a2,tau)  where  a1,a2  are the non-negative non-centrality
#                     parameters, tau is Kendall's tau of the copula family.
#' @param DoF     degrees of freedom of the Student copula (if needed).
#' @return \item{U}{Simulated Data}
#'
#' @author Bouchra R. Nasri, August 14, 2019
#'
#' @examples param <- c(0.8, 2.5, 0.7) ;
#'U <- SimNCSCop('Clayton', 250, param)
#'
#' @export
SimNCSCop<-function(family, n, param, DoF = NULL){
  a1  = param[1];
  a2  = param[2];
  KendallTau= param[3];

  switch(family,
           "Gaussian" = {
             alpha =  copula::iTau(copula::normalCopula(),  tau = KendallTau)
           },

           "t" = {
             alpha =  copula::iTau(copula::normalCopula(),  tau = KendallTau)
           },

           "Clayton" = {
             alpha =  copula::iTau(copula::claytonCopula(),  tau = KendallTau)
           },

           "Frank" = {
             alpha =  copula::iTau(copula::frankCopula(),  tau = KendallTau)
           },

           "Gumbel" = {
             alpha =  copula::iTau(copula::gumbelCopula(),  tau = KendallTau)
           }
    )



  switch(family,
         "Gaussian" =

          {
            V = copula::rCopula(n,copula::normalCopula(alpha, dim = 2))

          },

         "t" = {
           V = copula::rCopula(n, copula::tCopula(alpha, dim = 2,df = DoF))
                      },

         "Clayton" = {
          V = copula::rCopula(n,copula::claytonCopula(alpha, dim = 2))
                   },

         "Frank" = {
          V = copula::rCopula(n,copula::frankCopula(alpha, dim = 2))
                },

         "Gumbel" = {
          V = copula::rCopula(n,copula::gumbelCopula(alpha, dim = 2))
         }
  )

Z = qnorm(V);

X1 = abs(Z[,1]+a1);
X2 = abs(Z[,2]+a2);

U1 = pnorm(X1-a1)-pnorm(-X1-a1);
U2 = pnorm(X2-a2)-pnorm(-X2-a2);


U = matrix(c(U1,U2),byrow = FALSE, ncol=2);


  return(U);
}

