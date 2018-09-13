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

CPmodel <- function(a,b,c,d,t,gamma){
  return(a + b*t + c*sqrt((t-d)**2+gamma))
}

#' Random Change Point Mixed Model
#'
#' @param longdata A dataframe containing the variables used in the formula  \code{formu}
#' @param formu A formula object describing which variables are to be used. The formula has to be of the following form \code{markervar ~ scorevar | grouvpar} for the function to work.
#' @param covariate An optional string indicating a binary covariate to add on the fixed effects, i.e. intercept, mean slope, difference of slopes and changepoint date. The parameter \code{REadjust} indicates how this covariate influences the random effects variance structure. Default to NULL, i.e. no covariates.
#' @param REadjust An optional string value indicating how the random effects variance structure depends on \code{covariate}. "no" means that the structure doesn't depend upon \code{covariate}. "prop" indicates that the random effects variance structure is proportionnal according to \code{covariate} value. "yes" indicates that there is two different random effects variance structures, i.e. one for each level of \code{covariate}. Default to "no".
#' @param gamma A numeric parameter indicating how smooth the trajectory is on the changepoint date. Default to 0.1.
#' @param nbnodes A numeric parameter indicating how many nodes are to be used for the gaussian quadrature for numerical integration. Default to 10.
#'
#' @return The output contains several objects : \code{loglik} is the value of the log-likelihood at the optimum; \code{fixed} contains all fixed parameters estimates, standard errors, CIs, wald test statistic and corresponding pvalue when possible; \code{sdres} the estimated residual error; \code{VarEA} a 4x4 matrix or a list of 4x4 matrices - if there is some covariate for example - containing the estimated random effects covariance matrix; \code{optpar} the optimal parameters maximizing the log-likelihood; \code{covariate} the covariate declared in the function call; \code{REadjust} the string indicating how random effects structure is handled as declared in the function call, \code{invhessian} the covariance matrix containing all the standard errors and correlations of the parameter estimates;
#' @export
#'
#' @examples
rcpme <- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, param = NULL){
  
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (REadjust == "prop") stop("It has not been implemented yet. Sorry for the inconvenience...")
  if (!is.null(param)){
    if ((covariate == "NULL") & (length(param) != 12)) stop("Initial parameters vector must be a vector of size 12.")
    if ((covariate != "NULL") & (REadjust == "no") & (length(param) != 16)) stop("Initial parameters vector must be a vector of size 16.")
    if ((covariate != "NULL") & (REadjust == "yes") & (length(param) != 23)) stop("Initial parameters vector must be a vector of size 23.")
    if ((covariate != "NULL") & (REadjust == "prop") & (length(param) != 17)) stop("Initial parameters vector must be a vector of size 18.")
  }
  
  # =============================================
  
  scorevar = longdata[,all.vars(formu)[1]]
  timevar = longdata[,all.vars(formu)[2]]
  groupvar = longdata[,all.vars(formu)[3]]
  if (covariate != "NULL") adjustvar = longdata[,covariate]
  ngroupvar = rep(seq(length(unique(groupvar))), by(longdata,groupvar,function(x){return(dim(x)[[1]])}))
  longdata <- cbind(longdata, ngroupvar)

  ghcoeff <- gauher(nbnodes)
  nodes <- sqrt(2) * ghcoeff$x
  weights <- ghcoeff$w / sqrt(pi)
  
  if (is.null(param)){
    if (covariate != "NULL"){
      lmm <- nlme::lme(fixed = scorevar ~ 1 + timevar * adjustvar,
                       random = list(ngroupvar = nlme::pdSymm(~ 1 + timevar)),
                       na.action = na.omit,
                       method = "ML")
      
      param <- c(lmm$coefficients$fixed[1:2], -0.5, median(timevar), nlme::VarCorr(lmm)[c(3,1),2],0,0,1,0,1,1,lmm$coefficients$fixed[3:4],0,0)
      param <- as.numeric(param)
      
      if (REadjust == "yes"){
        param <- c(lmm$coefficients$fixed[1:2], -0.5, median(timevar), nlme::VarCorr(lmm)[c(3,1),2],0,0,1,0,1,1,lmm$coefficients$fixed[3:4],0,0,0,0,0,0,0,0,0)
        param <- as.numeric(param)
      }
      
      if (REadjust == "prop"){ # fonctionne pas pour le moment
        param <- c(lmm$coefficients$fixed[1:2], -0.5, median(timevar), nlme::VarCorr(lmm)[c(3,1),2],0,0,1,0,1,1,lmm$coefficients$fixed[3:4],0,0,0.1)
        param <- as.numeric(param)
      }
    }
    if (covariate == "NULL") {
      lmm <- nlme::lme(fixed = scorevar ~ 1 + timevar,
                       random = list(ngroupvar = nlme::pdSymm(~ 1 + timevar)),
                       na.action = na.omit,
                       method = "ML")  
      
      param <- c(lmm$coefficients$fixed[1:2], -0.5, median(timevar), nlme::VarCorr(lmm)[c(3,1),2],0,0,1,0,1,1)
      param <- as.numeric(param)
    }
  }
  

  # varsbeta0 <- if (all.vars(beta0)[1]=="1") all.vars(beta0)[-1] else all.vars(beta0)
  # varsbeta1 <- if (all.vars(beta1)[1]=="1") all.vars(beta1)[-1] else all.vars(beta1)
  # varsbeta2 <- if (all.vars(beta2)[1]=="1") all.vars(beta2)[-1] else all.vars(beta2)
  # varsbeta3 <- if (all.vars(mutau)[1]=="1") all.vars(mutau)[-1] else all.vars(mutau)
  
  opt <- optim(param,lvsblNCgen,data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}),nq=nbnodes,grp=ngroupvar,weights=weights, nodes=nodes, scorevar = all.vars(formu)[1], timevar = all.vars(formu)[2], covariate = covariate, REadjust = REadjust, hessian = TRUE, method="BFGS")

  # OUT : fixed parameters ===========================================================================================
  
  if (covariate != "NULL"){
    invhessian <- solve(opt$hessian)
    tab <- cbind(opt$par[c(1,13,2,14,3,15,4,16)],sqrt(diag(invhessian))[c(1,13,2,14,3,15,4,16)])
    tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2], tab[,1]/tab[,2],1-pchisq(tab[,1]**2/tab[,2]**2,df=1)); tab[c(5,7),c(5,6)] = NA;
    rownames(tab) <- c("beta0", paste("beta0*",covariate),"beta1", paste("beta1*",covariate),"beta2", paste("beta2*",covariate),"mutau", paste("mutau*",covariate))
    colnames(tab) <- c("par", "se(par)", "ICinf", "ICsup", "Wald stat.", "pvalue")
  }
  
  else {
    invhessian <- solve(opt$hessian)
    tab <- cbind(opt$par[c(1,2,3,4)],sqrt(diag(invhessian))[c(1,2,3,4)])
    tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2], tab[,1]/tab[,2], 1-pchisq(tab[,1]**2/tab[,2]**2,df=1))
    tab[c(3,4),c(5,6)] = NA
    rownames(tab) <- c("beta0", "beta1", "beta2", "mutau")
    colnames(tab) <- c("par", "se(par)", "ICinf", "ICsup", "Wald stat.", "pvalue")
  }
  
  # OUT : variance parameters ===========================================================================================
  
  if (REadjust == "no" | covariate == "NULL"){
    seps <- abs(opt$par[5]); U0 <- opt$par[6]; U01 <- opt$par[7]; U02 <- opt$par[8]; U1 <- opt$par[9]; U12 <- opt$par[10]; U2 <- opt$par[11]; stau <- abs(opt$par[12]);
    U <- matrix(c(U0,U01,U02,0,U1,U12,0,0,U2),nrow=3); 
    B <- t(U) %*% U;
    vars <- data.frame(c(seps, sqrt(diag(B)), stau, B[1,2:3], B[2,3]))
    rownames(vars) <- c("sdres", "sd0", "sd1", "sd2", "sdtau", "cov01", "cov02", "cov12")
    colnames(vars) <- c("par")
    VarEA <- matrix(0, nrow=4,ncol=4); VarEA[1:3,1:3] <- B; VarEA[4,4] <- stau**2;
  }
  
  if (REadjust == "yes"){
    seps <- abs(opt$par[5]); U0 <- opt$par[6]; U01 <- opt$par[7]; U02 <- opt$par[8]; U1 <- opt$par[9]; U12 <- opt$par[10]; U2 <- opt$par[11]; utau <- abs(opt$par[12]);
    U <- matrix(c(U0,U01,U02,0,U1,U12,0,0,U2),nrow=3); 
    B <- t(U) %*% U;
    V0 <- opt$par[6]+opt$par[17]; V01 <- opt$par[7]+opt$par[18]; V02 <- opt$par[8]+opt$par[19]; V1 <- opt$par[9]+opt$par[20]; V12 <- opt$par[10]+opt$par[21]; V2 <- opt$par[11]+opt$par[22]; vtau <- abs(opt$par[12]+opt$par[23]);
    V <- matrix(c(V0,V01,V02,0,V1,V12,0,0,V2),nrow=3); 
    BV <- t(V) %*% V;
    varsU <- c(seps, sqrt(diag(B)), utau, B[1,2:3], B[2,3])
    varsV <- c(seps, sqrt(diag(BV)), vtau, BV[1,2:3], BV[2,3])
    vars <- cbind(varsU,varsV);
    rownames(vars) <- c("sdres", "sd0", "sd1", "sd2", "sdtau", "cov01", "cov02", "cov12")
    colnames(vars) <- paste(covariate,levels(as.factor(adjustvar)))
    VarEA0 <- matrix(0, nrow=4,ncol=4); VarEA1 <- matrix(0, nrow=4,ncol=4);
    VarEA0[1:3,1:3] <- B; VarEA0[4,4] <- utau**2; VarEA1[1:3,1:3] <- BV; VarEA1[4,4] <- vtau**2;
    VarEA <- list(VarEA0, VarEA1); names(VarEA) <- paste(covariate,levels(as.factor(adjustvar)))
  }
  
  if (REadjust == "prop"){ # fonctionne pas pour le moment
    seps <- abs(opt$par[5]); U0 <- opt$par[6]; U01 <- opt$par[7]; U02 <- opt$par[8]; U1 <- opt$par[9]; U12 <- opt$par[10]; U2 <- opt$par[11]; utau <- abs(opt$par[12]);
    U <- matrix(c(U0,U01,U02,0,U1,U12,0,0,U2),nrow=3); 
    B <- t(U) %*% U;
    coeff <- abs(1+param[17])
    vars <- data.frame(c(seps, sqrt(diag(B)), utau, B[1,2:3], B[2,3], coeff))
    rownames(vars) <- c("sdres", "sd0", "sd1", "sd2", "sdtau", "cov01", "cov02", "cov12", "prop.")
    colnames(vars) <- c("par")
  }
  
  return(list("call" = as.list(match.call()), "Loglik" = opt$value, "formula" = formu, "fixed" = round(tab,3), "sdres"=vars[1,1], "VarEA" = VarEA, optpar= opt$par, "covariate" = covariate, "REadjust" = REadjust, "invhessian" = invhessian, "conv" = opt$convergence, "init" = param))
}

IndRePostDis <- function(re, data, rcpmeObj, scorevar, timevar){
  
  score = data[scorevar][!is.na(data[scorevar])]
  time = data[timevar][!is.na(data[scorevar])]
  lgt = length(score)
  
  if (rcpmeObj$covariate == "NULL"){
    betas <- rcpmeObj$fixed[,1]
    mu = CPmodel(a = betas[1]+re[1], b = betas[2]+re[2], c = betas[3]+re[3], d = betas[4]+re[4], t = time, gamma = 0.1)
  }
  
  if (rcpmeObj$covariate != "NULL"){
    betas <- rcpmeObj$fixed[,1]
    covInd = data[,rcpmeObj$covariate][1]
    mu = CPmodel(a = betas[1]+betas[2]*covInd+re[1], b = betas[3]+betas[4]*covInd+re[2], c = betas[5]+betas[6]*covInd+re[3], d = betas[7]+betas[8]*covInd+re[4], t = time, gamma = 0.1)
  }
  
  sdres = rcpmeObj$sdres; REadjust = rcpmeObj$REadjust;
  
  if (REadjust == "no" | rcpmeObj$covariate == "NULL") {estiVarEA <- rcpmeObj$VarEA;}
  if (REadjust == "yes") {estiVarEA <- rcpmeObj$VarEA[[covInd+1]];}
  
  return(-sum(dnorm(score, mean = mu, sd = sdres, log = TRUE)) - mvtnorm::dmvnorm(re, c(0,0,0,0), estiVarEA, log = TRUE))
} # coder en Rcpp pour aller plus vite ?!

#' Random Effects estimate
#'
#' @param rcpmeObj An object which contains the output of a \code{rcpme} function call.
#' @param var A boolean indicating if the covariance matrices of random effects estimation should be returned. Default to \code{FALSE}.
#'
#' @return A list with each item containing the subject specific random effects estimation and its covariance matrix if \code{var = TRUE}.
#' @export
#'
#' @examples
REestimate <- function(rcpmeObj, longdata, var = FALSE){
  
  formu <- rcpmeObj$formula
  scorevar = all.vars(formu)[1]
  timevar = all.vars(formu)[2]
  groupvar = all.vars(formu)[3]
  
  re <- rep(0,4)
  if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){opt<-optim(par = re, fn = IndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar);return("par"=opt$par)}))}
  if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){opt<-optim(par = re, fn = IndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, hessian = TRUE);return(list("par"=opt$par,"var"=solve(opt$hessian)))}))}

}

# by(paquidCEPDEM,paquidCEPDEM$ID,function(x){return(IndRePostDis(re=rep(0,4), x, rcpmeobj, "isa15", "delai"))})
# by(paquidCEPDEM,paquidCEPDEM$ID,function(x){opt<-optim(par = re, fn = IndRePostDis, rcpmeObj = rcpmeobj, data = x, scorevar = "isa15", timevar = "delai", hessian = TRUE);return(list("par"=opt$par,"hess"=opt$hessian))})

#' Bivariate Random Change Point Mixed Model
#'
#' @param longdata A dataframe containing the variables used in the formula  \code{formu}
#' @param formu A formula object describing which variables are to be used. The formula has to be of the following form \code{markervar1 + markervar2 ~ scorevar | grouvpar} for the function to work.
#' @param covariate An optional string indicating a binary covariate to add on the fixed effects, i.e. intercept, mean slope, difference of slopes and changepoint date. The parameter \code{REadjust} indicates how this covariate influences the random effects variance structure. Default to NULL, i.e. no covariates.
#' @param REadjust An optional string value indicating how the random effects variance structure depends on \code{covariate}. "no" means that the structure doesn't depend upon \code{covariate}. "prop" indicates that the random effects variance structure is proportionnal according to \code{covariate} value. "yes" indicates that there is two different random effects variance structures, i.e. one for each level of \code{covariate}. Default to "no".
#' @param gamma A numeric parameter indicating how smooth the trajectory is on the changepoint date. Default to 0.1.
#' @param nbnodes A numeric parameter indicating how many nodes are to be used for the gaussian quadrature for numerical integration. Default to 10.
#' @param adapt A boolean indicating whether adaptive gaussian quadrature should be used for numerical integration. Default to FALSE, the adaptive quadrature is still being implemented and TRUE will return an error.
#' @param param 
#'
#' @return The output contains several objects : \code{loglik} is the value of the log-likelihood at the optimum; \code{fixed} contains all fixed parameters estimates, standard errors, CIs, wald test statistic and corresponding pvalue when possible; \code{sdres} the estimated residual error; \code{VarEA} a matrix containing the estimated random effects covariance matrix of the eight random effects: four for each marker with a general correlation structure between them; \code{optpar} the optimal parameters maximizing the log-likelihood; \code{covariate} the covariate declared in the function call; \code{REadjust} the string indicating how random effects structure is handled as declared in the function call, \code{invhessian} the covariance matrix containing all the standard errors and correlations of the parameter estimates;
#' @export
#'
#' @examples
bircpme <- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, adapt = FALSE, param = NULL){
  
  # errors handling
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (covariate != "NULL" | REadjust != "no") stop("It has not been implemented yet. Sorry for the inconvenience.")
  if (adapt == TRUE) stop("Adaptive gaussian quadrature not yet implemented")
  
  # data handling
  scorevar1 = longdata[,all.vars(formu)[1]]
  scorevar2 = longdata[,all.vars(formu)[2]]
  timevar = longdata[,all.vars(formu)[3]]
  groupvar = longdata[,all.vars(formu)[4]]
  if (covariate != "NULL") adjustvar = longdata[,covariate]
  ngroupvar = rep(seq(length(unique(groupvar))), by(longdata,groupvar,function(x){return(dim(x)[[1]])}))
  longdata <- cbind(longdata, ngroupvar)
  
  # estimation of univariate models and initialization
  print("I need to estimate both univariate models first...")
  rcpmeObj1 <- rcpme(longdata, as.formula(paste(all.vars(formu)[1], "~", all.vars(formu)[3], "|", all.vars(formu)[4])), covariate = covariate, REadjust = REadjust, gamma = gamma, nbnodes = 10)
  rcpmeObj2 <- rcpme(longdata, as.formula(paste(all.vars(formu)[2], "~", all.vars(formu)[3], "|", all.vars(formu)[4])), covariate = covariate, REadjust = REadjust, gamma = gamma, nbnodes = 10)
  
  if (is.null(param)){
    param <- as.numeric(c(rcpmeObj1$optpar, rcpmeObj2$optpar, rep(0.1,10)))
  }
  
  if (length(param) != 34){
    stop("Initial parameters vector must be a vector of size 34.")
  }
  
  # nodes and weights
  ghcoeff <- gauherm(nbnodes, 2)
  nodes <- ghcoeff$x
  weights <- ghcoeff$w
  # non adaptive
  if (adapt == FALSE){
    nodes <- sqrt(2) * nodes
    weights <- weights / pi
    newnodes = NULL; newweights = NULL;
  }
  # adaptive
  if (adapt == TRUE){
    RE1 <- REestimate(rcpmeObj1, longdata, var = TRUE); RE2 <- REestimate(rcpmeObj2, longdata, var = TRUE);
    REs <- lapply(mapply(FUN = function(a,b){return(list(c(a$par[4], b$par[4], a$var[4,4], b$var[4,4])))}, RE1, RE2), function(x){return(list("par"=x[c(1,2)], "var"=matrix(c(sqrt(x[3]),0,0,sqrt(x[4])), byrow=TRUE, nrow=2)))})
    newnodes <- mapply(function(a){return(list(t(apply(nodes,1,function(x){return(sqrt(2)*a[[2]]%*%x+a[[1]])}))))},REs)
    newweights <- lapply(REs,function(a){return(weights*det(sqrt(2)*a[[2]]))*2*apply(nodes,1,function(x){return(exp(t(x)%*%x))})})
  }
  
  # optimization
  # bilvsblNC(param, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, grp=ngroupvar, weights=weights,  nodes=nodes, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust)
  print("I can begin the optimization. Please be aware that it can take some time to run.")
  
  # opt <- bilvsblNC(param, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust)
  # opt <- optim(param, bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, method="BFGS", control = list(trace=TRUE, REPORT = 10))
  opt <- optim(param, bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, method="BFGS", hessian = TRUE)
  
  # OUT : fixed parameters ===========================================================================================
  if (covariate == "NULL"){
    invhessian <- solve(opt$hessian)
    tab <- cbind(opt$par[c(1,2,3,4)],sqrt(diag(invhessian))[c(1,2,3,4)])
    tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2])
    tab <- cbind(tab, opt$par[c(1,2,3,4)+12], sqrt(diag(invhessian))[c(1,2,3,4)+12])
    tab <- cbind(tab, tab[,5]-1.96*tab[,6],tab[,5]+1.96*tab[,6])
    covs <- c(invhessian[1,13], invhessian[2,14], invhessian[3,15], invhessian[4,16])
    tab <- cbind(tab, sqrt((tab[,1]-tab[,5])**2/(tab[,2]**2+tab[,6]**2-2*covs)), 1-pchisq(sqrt((tab[,1]-tab[,5])**2/(tab[,2]**2+tab[,6]**2-2*covs)), df=1))
    rownames(tab) <- c("beta0", "beta1", "beta2", "mutau")
    colnames(tab) <-  c(paste(all.vars(formu)[1],":", c("par", "se(par)", "ICinf", "ICsup"),sep = ""), paste(all.vars(formu)[2],":", c("par", "se(par)", "ICinf", "ICsup"),sep = ""), "Wald stat.", "pvalue")
  }

  # OUT : variance parameters ========================================================================================
  if (REadjust == "no" | covariate == "NULL"){
    sdres <- abs(opt$par[c(5,17)]); names(sdres) = c(all.vars(formu)[1],all.vars(formu)[2]);
    U0_1 <- opt$par[6]; U01_1 <- opt$par[7]; U02_1 <- opt$par[8]; U1_1 <- opt$par[9]; U12_1 <- opt$par[10]; U2_1 <- opt$par[11]; Utau_1 <- opt$par[12];
    U0_2 <- opt$par[18]; U01_2 <- opt$par[19]; U02_2 <- opt$par[20]; U1_2 <- opt$par[21]; U12_2 <- opt$par[22]; U2_2 <- opt$par[23]; Utau_2 <- opt$par[24];
    Utau_12 <- opt$par[25]; U1_12 <- opt$par[26]; U2_12 <- opt$par[27]; U3_12 <- opt$par[28]; U4_12 <- opt$par[29]; U5_12 <- opt$par[30];
    U6_12 <- opt$par[31]; U7_12 <- opt$par[32]; U8_12 <- opt$par[33]; U9_12 <- opt$par[34];
    U <- matrix(c(U0_1,U01_1,U02_1,0,U1_12,U2_12,U3_12,0,0,U1_1,U12_1,0,U4_12,U5_12,U6_12,0,0,0,U2_1,0,U7_12,U8_12,U9_12,0,0,0,0,Utau_1,0,0,0,Utau_12,0,0,0,0,U0_2,U01_2,U02_2,0,0,0,0,0,0,U1_2,U12_2,0,0,0,0,0,0,0,U2_2,0,0,0,0,0,0,0,0,Utau_2),byrow=TRUE, nrow=8)
    B <- t(U) %*% U
  }
  
  # return(opt)
  # return(list("Loglik" = opt$value, optpar= opt$par, conv = opt$convergence))
  return(list("Loglik" = opt$value, "fixed" = round(tab,3), "sdres"=sdres, "VarEA" = B, optpar= opt$par, "covariate" = covariate, "REadjust" = REadjust, "invhessian"= invhessian, "conv" = opt$convergence, "init" = param))
}

#' Individual prediction based on a \code{rcpme} model
#'
#' @param rcpmeObj An object which contains the output of a \code{rcpme} function call on which individual predictions are to be computed. 
#'
#' @return A list with each individual prediction
#' @export
#'
#' @examples
IndPred <- function(rcpmeObj){
  
  formu <- rcpmeObj$call$formu
  timevar = all.vars(formu)[2]
  groupvar = all.vars(formu)[3]
  
  longdata <- get(paste(rcpmeObj$call$longdata))
  REesti <- REestimate(rcpmeObj, longdata)
  betas <- rcpmeObj$fixed[,1]
  
  if (rcpmeObj$covariate == "NULL"){
    betas <- as.numeric(rcpmeObj$fixed[,1])
    Ypred <- mapply(function(a,b){return(CPmodel(betas[1]+a[1], betas[2]+a[2], betas[3]+a[3], betas[4]+a[4], t = b[,timevar], gamma = 0.1))}, REesti, by(longdata,longdata[,groupvar],function(x){return(x)}))
  }
  
  if (rcpmeObj$covariate != "NULL"){
    betas <- as.numeric(rcpmeObj$fixed[,1])
    covInd = rcpmeObj$covariate
    Ypred <- mapply(function(a,b){return(CPmodel(betas[1]+betas[2]*b[,covInd][1]+a[1], betas[3]+betas[4]*b[,covInd][1]+a[2], betas[5]+betas[6]*b[,covInd][1]+a[3], betas[7]+betas[8]*b[,covInd][1]+a[4], t = b[,timevar], gamma = 0.1))}, REesti, by(longdata,longdata[,groupvar],function(x){return(x)}))
  }
  return(Ypred) 
}
