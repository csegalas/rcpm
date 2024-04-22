# test of the existence of a random changepoint in a mixed model
#' testRCPMM
#'
#' Realizes the supremum score test on longdata according to formu.
#'
#' @param longdata A longitudinal dataset containing all variables used in the formula \code{formu}
#' @param formu A formula object describing which variables are to be used. The formula has to be of the following form \code{markervar ~ scorevar | groupvar} for the function to work.
#' @param gamma A smoothing parameter for the transition on the changepoint date. 0.1 by default.
#' @param nbnodes Number of pseudo-adaptive Gaussian quadrature nodes used to compute the numeric integrals. 5 by default.
#' @param nbpert Number of perturbations used to compute the empirical p-value. 500 by default.
#' @param covariate An eventual covariate dependence of all the parameters in the model. Not implemented yet.
#'
#' @return The function returns a list with the computed empirical p-value and the observed test statistic.
#' @export
#'
#' @examples
testRCPMM <- function(longdata, formu, covariate = NULL, gamma = 0.1, nbnodes = 5, nbpert = 500){
  
  if(is.null(longdata)){stop("I can not do my job if a dataset is not provided in the longdata argument of the function.")}
  if(is.null(formu)){stop("Please indicate me in the formu argument which variables are to be used.")}
  if(!is.null(covariate)){stop("I can not deal with covariate at the moment. I am sorry for the inconvenience.")}
  
  scorevar = longdata[,all.vars(formu)[1]]
  timevar = longdata[,all.vars(formu)[2]]
  groupvar = longdata[,all.vars(formu)[3]]
  ngroupvar = rep(seq(length(unique(groupvar))), lapply(split(longdata, groupvar), function(x) return(dim(x)[[1]])))
  longdata <- cbind(longdata, ngroupvar)
  # adjustvar = longdata[,all.vars(covariate)]
  
  # null model estimation
  nullmodel <- nlme::lme(fixed = scorevar ~ 1 + timevar, # + adjustvar,
               random = list(ngroupvar = nlme::pdSymm(~ 1 + timevar)),
               na.action = na.omit,
               method = "ML")

  # estimated parameters from nlme
  Beta0 <- nullmodel$coefficients$fixed[[1]]
  Beta1 <- nullmodel$coefficients$fixed[[2]]
  sigma <- nullmodel$sigma
  sigmab0 <- as.numeric(nlme::VarCorr(nullmodel)[1,2])
  sigmab1 <- as.numeric(nlme::VarCorr(nullmodel)[2,2])
  covb01 <- as.numeric(nlme::VarCorr(nullmodel)[2,3]) * sigmab0 * sigmab1
  estis <- c(Beta0, Beta1, sigma)
  estiVarEA <- matrix(c(sigmab0**2, covb01, covb01, sigmab1**2), nrow = 2)
  remove(Beta0, Beta1, sigma, sigmab0, sigmab1, covb01)

  # estimating individual random effects
  init <- c(1,1)
  hatRE <- lapply(split(longdata, ngroupvar), function(x) return(hatRE(x, init, estis, estiVarEA, formu)))
  estiranef <- lapply(hatRE,function(x){return(c(x[[1]],0))})
  varbi <- lapply(hatRE,function(x){return(rbind(cbind(x[[2]],c(0,0)),c(0,0,1/2)))})
  cholvarbi <- lapply(varbi,chol)


  # computation of adaptGH nodes and weights
  nRE = 3
  nodes <- gauherm(nbnodes,nRE)$x
  weights <- gauherm(nbnodes,nRE)$w
  newnodes <- mapply(function(a,b){return(list(t(apply(nodes,1,function(x){return(sqrt(2)*b%*%x+a)}))))},estiranef,cholvarbi)
  newweights <- lapply(cholvarbi,function(a){return(weights*det(sqrt(2)*a))*apply(nodes,1,function(x){return(exp(t(x)%*%x))})})
  
  # calcul des qtes necessaires au calcul du score
  MATmus <- MATmu(estis,longdata,newnodes,nbnodes,3,ngroupvar,formu)

  # calcul de la stat de test observee sur l'echantillon simule
  tau <- c(round(mean(timevar, na.rm=TRUE),1),1)
  pert <- 1
  optobs <- optim(par=tau, fn=Score, estis=estis,estiVarEA=estiVarEA,longdata=longdata,newnodes=newnodes,newweights=newweights,nbnodes=nbnodes,nd=3,groupvar=ngroupvar, MATmus = MATmus, formu = formu, pert=pert, method=c("L-BFGS-B"), control=list(maxit=250, factr = 0.001, pgtol=0.001)) # optim : reltol001

  # perturbations
  res <- rep(NA, nbpert)
  n <- max(ngroupvar, na.rm = TRUE)
  for (i in 1:nbpert){
    pert <- rnorm(n,0,1)
    opt <- optim(par=tau, fn=Score, estis=estis,estiVarEA=estiVarEA,longdata=longdata,newnodes=newnodes,newweights=newweights,nbnodes=nbnodes,nd=3,groupvar=ngroupvar, MATmus = MATmus, formu = formu, pert=pert, method=c("L-BFGS-B"), control=list(maxit=250, factr = 0.001, pgtol=0.001)) # optim : reltol001
    res[i] <- opt$value
  }
  
  return(list("empirical p-value" = mean(-res > -optobs$value), "obs. stat" = optobs$value))
}
