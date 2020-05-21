# functions for individual random effects estimation
loglikmode <- function(param, estis, estiVarEA, data, scorevar, timevar){
  
  score = scorevar[!is.na(scorevar)]
  date = timevar[!is.na(scorevar)]
  
  b0 <- param[1]
  b1 <- param[2]
  
  Beta0 <- estis[1]
  Beta1 <- estis[2]
  sigma <- estis[3]
  
  MATmu <- Beta0 + b0 + (Beta1 + b1) * date
  
  return(-sum(dnorm(score, mean = MATmu, sigma, log = TRUE)) - mvtnorm::dmvnorm(c(b0,b1), c(0,0), estiVarEA, log = TRUE))
}
estisranef <- function(data, init, estis, estiVarEA, formu){
  
  opt <- optim(par = init, loglikmode, estis = estis, estiVarEA = estiVarEA, data = data, scorevar = data[,all.vars(formu)[1]], timevar = data[,all.vars(formu)[2]], hessian = TRUE, method=c("L-BFGS-B"))
  
  ranef <- opt$par
  ranefvar <- solve(opt$hessian)
  
  return(list(ranef,matrix(ranefvar, byrow=TRUE, nrow=2)))
}

# functions for gaussian quadrature
gauher <-  function(n) {
  m <- trunc((n + 1)/2)
  x <- w <- rep(-1, n)
  for (i in seq_len(m)) {
    z <- if (i == 1) {
      sqrt(2*n + 1) - 1.85575 * (2*n + 1)^(-0.16667)
    } else if (i == 2) {
      z - 1.14 * n^0.426 / z
    } else if (i == 3) {
      1.86 * z - 0.86 * x[1]
    } else if (i == 4) {
      1.91 * z - 0.91 * x[2]
    } else {
      2*z - x[i - 2]
    }
    for (its in seq_len(10)) {
      p1 <- 0.751125544464943
      p2 <- 0
      for (j in seq_len(n)) {
        p3 <- p2
        p2 <- p1
        p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
      }
      pp <- sqrt(2*n) * p2
      z1 <- z
      z <- z1 - p1/pp
      if (abs(z - z1) <= 3e-14) 
        break
    }
    x[i] <- z
    x[n + 1 - i] <- -z
    w[i] <- 2 / (pp * pp)
    w[n + 1 - i] <- w[i]
  }
  list(x = x, w = w)
}
gauherm <- function(n,m){
  grid <- gauher(n)
  idx <- as.matrix(expand.grid(rep(list(1:n),m)))
  points <- matrix(grid$x[idx],nrow(idx),m)
  weights <- apply(matrix(grid$w[idx],nrow(idx),m),1,prod)
  
  return(list(x = points, w = weights))
}

# functions for individual quantities
MATmuInd <- function(data, Beta0, Beta1, sigma, newnodes, nbnodes, nd=3, formu){
  
  scorevar = data[!is.na(data[,all.vars(formu)[1]]),all.vars(formu)[1]]
  timevar = data[!is.na(data[,all.vars(formu)[1]]),all.vars(formu)[2]]
  lgt = length(scorevar)
  
  if (lgt==0){
    return(0)
  } else {
    idind <- data[,"ngroupvar"][1]
    MATmu <- (Beta0 + newnodes[[idind]][,1]) %*% t(rep(1,lgt)) + (Beta1+newnodes[[idind]][,2]) %*% t(timevar)
    diffY_MATmu <- (rep(1,nbnodes**nd) %*% t(scorevar)) - MATmu
    expNorm <- exp(-0.5*(diffY_MATmu/sigma)**2)
    if (lgt == 1){
      # prodIndDiff <- <- rep(1,nq**nd) modifie pour avoir type matrix
      prodIndDiff <- matrix(rep(1,nbnodes**nd),nrow = nbnodes**nd, ncol = 1)
    } else if (lgt == 2){
      prodIndDiff <- sapply(1:lgt,function(x){return(expNorm[,-x])})
    } else { 
      prodIndDiff <- sapply(1:lgt,function(x){return(apply(expNorm[,-x],1,prod))})
    }
  } 
  return(list(MATmu, diffY_MATmu, expNorm, prodIndDiff))
}
MATmu <- function(estis, longdata, newnodes, nbnodes, nd=3, groupvar, formu){
  
  Beta0 <- estis[1]
  Beta1 <- estis[2]
  sigma <- estis[3]
  
  return(by(longdata, groupvar, MATmuInd, Beta0, Beta1, sigma, newnodes, nbnodes, nd, formu))
}

# functions for the test statistic
tst <- function(x, pert){
  # return(sum(x)**2/(sum(x**2)-length(x)*mean(x)**2))
  # return(sum(x)**2/(sum(x**2)))
  # return(sum(x * pert)**2/(sum(x**2)-length(x)*mean(x)**2))
  return(sum(x * pert)**2/(sum(x**2)))
}

Score <- function(tau,estis,estiVarEA,longdata,newnodes,newweights,nbnodes,nd,groupvar, MATmus, formu, pert){
  mutau <- tau[1]
  sigmatau <- tau[2]
  
  Beta0 <- estis[1]
  Beta1 <- estis[2]
  sigma <- estis[3]
  
  # U0 <- estis[4]
  # U01 <- estis[5]
  # U1 <- estis[6]
  # U <- matrix(c(U0,U01,0,0,U1,0,0,0,sigmatau), byrow = TRUE,nrow = 3)
  # 
  # B <- t(U) %*% U
  
  B <- rbind(cbind(estiVarEA,c(0,0)),c(0,0,sigmatau**2))
  
  return(-tst(by(longdata, groupvar, ScoreInd, Beta0, Beta1, sigma, mutau, sigmatau, B, MATmus, nbnodes, nd, newnodes, newweights, all.vars(formu)[1],all.vars(formu)[2],"ngroupvar"), pert)) # rcpp
  # return(-tst(by(data, grp, ScoreInd, sigma, mutau, sigmatau, B, MATmus, nq, nd, newnodes, newweights)))
}

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
  ngroupvar = rep(seq(length(unique(groupvar))), by(longdata,groupvar,function(x){return(dim(x)[[1]])}))
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
  hatRE <- by(longdata, ngroupvar, estisranef, init, estis, estiVarEA, formu)
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
