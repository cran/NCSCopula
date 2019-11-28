#' @title Estimation of  a non-central squared copula model
#'
#' @description This function estimates the copula parameter and the non-centrality parameters of a
#' non-central squared copula
#'
#' @param y   (nx2) data matrix (observations or residuals) that will be transformed to pseudo-observations
#' @param family    'Gaussian' , 't' , 'Clayton' , 'Frank' , 'Gumbel'
#' @param p    number of non-centrality parameters to be estimated (p = 0,1,2)
#' @param InitialValues initial values c(a1,a2,tau) to start the estimation; otherwise pre-selected values will be used
#'
#' @author Bouchra R. Nasri, August 14, 2019
#' @return \item{theta}{Estimated parameter of the copula according to CRAN copula package}
#' @return \item{dof}{Estimated degrees of freedom, only for the Student copula}
#' @return \item{tau}{Estimated theoretical Kendall tau for the copula family}
#'
#' @references Section 5.1 of Nasri, Rémillard & Bouezmarni (2019). Semi-parametric copula-based models under non-stationarity,
#'Journal of Multivariate Analysis, 173, pages 347-365.
#'
#' @examples
#' \donttest{
#' param <- c(0.8, 2.5, 0.7) ;
#'U <- SimNCSCop('Clayton', 250, param)
#'estimation <- EstNCSCop(U,'Clayton')
#'}
#'
#'
#' @importFrom copula iTau pcopula dcopula tau normalCopula tCopula frankCopula gumbelCopula claytonCopula
#' @importFrom  stats qnorm pnorm dnorm qchisq
#' @export
#'

EstNCSCop<-function(y,family,p=2,InitialValues=NULL){

  n = dim(y)[1]; d = dim(y)[2]
  U = floor( apply(y,2,rank) )/ (n + 1)


#InitialValues

  alpha0 = ParamTau(family,0.5);  #unconstrained parameters

  if (is.null(InitialValues))
  {
    if (p==2){
      initialGuess = c(finv( 1.0),finv(1.5), alpha0, 0);    #a1 = 1.0, a2 = 1.5, tau = 0.5, nu = 12.5
    } else if (p==1){
      initialGuess =c(finv(1.5), alpha0, 0);    #a1 = 1.0, a2 = 1.5, tau = 0.5, nu = 12.5
    }
  else if (p==0){
    initialGuess = c(alpha0, 0);      # tau = 0.5,  nu = 12.5
  }
  }
  else{

  a10  = InitialValues[1];
  a20  = InitialValues[2];
  tau0 = InitialValues[3];
  nu0  = 12.5;

  initialGuess = c(finv(a10), finv(a20),ParamTau(family,tau0),-log(25/nu0 -1));
  }

  # Optimisation


  alpha =  stats::optim(initialGuess,function(x){ -sum(LoglikNCSCop(x,U,family, p))}, method = "Nelder-Mead",control=list(maxit=20000))$par

  #  alpha_new =  stats::optim(par= alpha,fun, method = "L-BFGS-B")$par



  dof = NaN;
  #theta = ParamCop(family,alpha[(p+1):(p+2)]);
 out = KendallTau(family,alpha[(p+1):(p+2)]);
  tau = out$tau;
  theta0 = out$theta;

  theta = theta0[1];
  dof = theta0[2];
  if(p>0){  a=f(alpha[1:p])}
  else{a=NaN}

  LL =  sum(LoglikNCSCop(alpha,U,family, p));
  out = list(tau=tau, a=a, theta=theta, dof=dof,LL=LL,alpha=alpha);
  return(out)
}

finv = function(x){ return( log(x/(3-x)) )}
f    = function(x){ return(3./(1+exp(-x)))}
