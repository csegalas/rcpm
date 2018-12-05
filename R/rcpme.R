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
geneData <- function(n = 100, hyp = "null", params, pNA = 0, DO = 0, vis = c(0, 3, 6, 9, 12, 15, 18, 21), beta2alea=FALSE){
  
  nbvis = length(vis)
  donnees <- matrix(NA, nrow = n * nbvis, ncol = ifelse(hyp=="altaleabiv", 4, 3))
  colnames(donnees) <- if (hyp=="altaleabiv"){c("ID", "date", "score1", "score2")} else {c("ID", "date", "score")}
  
  if (hyp=="null"){
    if (length(params) != 6){
      stop("Under the null, params must be of length 6.")
    }
    else {
      Beta0 <- params[1]
      Beta1 <- params[2]
      sigma <- params[3]
      sigmab0 <- params[4]
      sigmab1 <- params[5]
      sigmab01 <- params[6] * sigmab0 * sigmab1
      B = matrix(c(sigmab0**2, sigmab01, sigmab01, sigmab1**2), nrow = 2)
      bis <- matrix(NA, nrow = n, ncol = 2)
      for (i in seq(n)){
        # génération des effets aléatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = Beta0 + bis[i,1] + (Beta1 + bis[i,2]) * vis[j] + eps
        }
      }
    }
  }  else if (hyp=="alt") {
    if (length(params) != 9){
      stop("Under the alternative, params must be of length 9.")
    }
    else {
      Beta0 <- params[1]
      Beta1 <- params[2]
      Beta2 <- params[3]
      mutau <- params[4]
      sigma <- params[5]
      sigmab0 <- params[6]
      sigmab1 <- params[7]
      sigmabtau <- params[8]
      sigmab01 <- params[9] * sigmab0 * sigmab1
      B = matrix(c(sigmab0**2, sigmab01, 0, sigmab01, sigmab1**2, 0, 0, 0, sigmabtau**2), nrow = 3)
      bis <- matrix(NA, nrow = n, ncol = 3)
      for (i in seq(n)){
        # génération des effets aléatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0 + bis[i,1], Beta1 + bis[i,2], Beta2, mutau + bis[i,3], vis[j], 0.1) + eps
        }
      }
    }
  } else if (hyp=="altalea"){
    if (length(params) != 12){
      stop("Under the alternative with beta_2 random, params must be of length 12.")
    }
    else {
      Beta0 <- params[1]
      Beta1 <- params[2]
      Beta2 <- params[3]
      mutau <- params[4]
      sigma <- params[5]
      sigmab0 <- params[6]
      sigmab1 <- params[7]
      sigmab2 <- params[8]
      sigmabtau <- params[9]
      sigmab01 <- params[10] * sigmab0 * sigmab1
      sigmab02 <- params[11] * sigmab0 * sigmab2
      sigmab12 <- params[12] * sigmab2 * sigmab1
      B = matrix(c(sigmab0**2, sigmab01, sigmab02, 0, sigmab01, sigmab1**2, sigmab12, 0, sigmab02, sigmab12, sigmab2**2, 0, 0, 0, 0, sigmabtau**2), nrow = 4)
      bis <- matrix(NA, nrow = n, ncol = 4)
      for (i in seq(n)){
        # génération des effets aléatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0,0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0 + bis[i,1], Beta1 + bis[i,2], Beta2 + bis[i,3], mutau + bis[i,4], vis[j], 0.1) + eps
        }
      }
    }
  } else if (hyp=="altaleabiv"){
    if (length(params) != 34){
      stop("Under the bivariate alternative, params must be of length 34.")
    }
    
    else {
      Beta0_1 <- params[1]; Beta1_1 <- params[2]; Beta2_1 <- params[3]; mutau_1 <- params[4]; Vareps_1 <- params[5]; B0_1 <- params[6]; B01_1 <- params[7]; B02_1 <- params[8];
      B1_1 <- params[9]; B12_1 <- params[10]; B2_1 <- params[11]; Btau_1 <- params[12];
      Beta0_2 <- params[13]; Beta1_2 <- params[14]; Beta2_2 <- params[15]; mutau_2 <- params[16]; Vareps_2 <- params[17]; B0_2 <- params[18]; B01_2 <- params[19]; B02_2 <- params[20];
      B1_2 <- params[21]; B12_2 <- params[22]; B2_2 <- params[23]; Btau_2 <- params[24];
      Btau_12 <- params[25]; B1_12 <- params[26]; B2_12 <- params[27]; B3_12 <- params[28]; B4_12 <- params[29]; B5_12 <- params[30];
      B6_12 <- params[31]; B7_12 <- params[32]; B8_12 <- params[33]; B9_12 <- params[34];
      
      B <- matrix(c(B0_1, B01_1, B02_1, 0, B1_12, B2_12, B3_12, 0,
                    B01_1, B1_1, B12_1, 0, B4_12, B5_12, B6_12, 0,
                    B02_1, B12_1, B2_1, 0, B7_12, B8_12, B9_12, 0,
                    0, 0, 0, Btau_1, 0, 0, 0, Btau_12,
                    B1_12, B4_12, B7_12, 0, B0_2, B01_2, B02_2, 0,
                    B2_12, B5_12, B8_12, 0, B01_2, B1_1, B12_2, 0,
                    B3_12, B6_12, B9_12, 0, B02_2, B12_2, B2_2, 0,
                    0, 0, 0, Btau_12, 0, 0, 0, Btau_2), byrow = TRUE, nrow = 8)
      
      bis <- matrix(NA, nrow = n, ncol = 8);
      for (i in seq(n)){
        bis[i,] <- rmvnorm(1, mean = rep(0,8), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps1 = rnorm(1,0,sqrt(Vareps_1)); eps2 <- rnorm(1,0,sqrt(Vareps_2));
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0_1 + bis[i,1], Beta1_1 + bis[i,2], Beta2_1 + bis[i,3], mutau_1 + bis[i,4], vis[j], 0.1) + eps1
          donnees[(nbvis*(i-1))+j,4] = CPmodel(Beta0_2 + bis[i,5], Beta1_2 + bis[i,6], Beta2_2 + bis[i,7], mutau_2 + bis[i,8], vis[j], 0.1) + eps2
        }
      }
      
    }
  }
  
  donnees <- as.data.frame(donnees)
  if (pNA!=0 & DO !=0){
    stop("pNA and DO mustn\'t be both positive")
  }
  if (pNA<0 | pNA>1 | DO<0 | DO>1){
    stop("pNA and DO takes values in [0;1]")
  }
  
  # censures
  if (pNA !=0){
    censures <- ifelse(c(t(cbind(rep(0,n),matrix(rbinom((nbvis-1)*n,1,pNA),nrow = n, ncol = nbvis-1))))==1,NA,0)
    donnees$score <- ifelse(is.na(censures),NA,donnees$score)
  }
  
  # drop-out
  if (DO !=0){
    geo <- rgeom(n, DO)+1
    dropout_times <- ifelse(geo>7,NA,geo)
    for (i in seq(n)){ 
      if (!is.na(dropout_times[i])){
        donnees$score[donnees$ID==i][-(1:(dropout_times[i]))] = NA
      }
    }
  }
  
  return(donnees)
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
                       method = "ML", control = nlme::lmeControl(opt = "optim"))  
      
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
                       method = "ML", control = nlme::lmeControl(opt = "optim"))  
      
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

#' Random Effects estimate
#'
#' @param rcpmeObj An object which contains the output of a \code{rcpme} function call.
#' @param var A boolean indicating if the covariance matrices of random effects estimation should be returned. Default to \code{FALSE}.
#'
#' @return A list with each item containing the subject specific random effects estimation and its covariance matrix if \code{var = TRUE}.
#' @export
#'
#' @examples
REestimate <- function(rcpmeObj, longdata, var = FALSE, onlytau = FALSE){
  
  formu <- rcpmeObj$formula
  scorevar = all.vars(formu)[1]
  timevar = all.vars(formu)[2]
  groupvar = all.vars(formu)[3]
  
  if (onlytau == FALSE){
    re <- rep(0,4)
    if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar);return("par"=opt$b);} else {return("par"=re);}}))}
    if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar);return(list("par"=opt$b,"var"=matrix(c(opt$v[c(1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10)]), nrow = 4, byrow=TRUE)));} else {return(list("par"=re, "var"=diag(4)))}}))}
  }
  
  else {
    re <- 0
    if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar);return("par"=opt$b);} else {return("par"=re);}}))}
    if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar);return(list("par"=opt$b,"var"=opt$v));} else {return(list("par"=0, "var"=1))}}))}
    # if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){opt<-optim(par = re, fn = IndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, method = "BFGS");return("par"=opt$par)}))}
    # if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){opt<-optim(par = re, fn = IndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, hessian = TRUE, method = "BFGS");return(list("par"=opt$par,"var"=solve(opt$hessian)))}))}
    
  }
  
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
bircpme <- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, adapt = FALSE, param = NULL, nproc = 1){
  
  # errors handling
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (covariate != "NULL" | REadjust != "no") stop("It has not been implemented yet. Sorry for the inconvenience.")
  if (!is.null(param)){
    if (length(param) != 34){ stop("Initial parameters vector must be a vector of size 34.")}
    B <- matrix(c(param[c(6,7,8)], 0, param[c(26,27,28)], 0, 
                  param[c(7,9,10)], 0, param[c(29,30,31)], 0, 
                  param[c(8,10,11)], 0, param[c(32,33,34)], 0, 
                  rep(0,3), param[12], rep(0,3), param[25],
                  param[c(26,29,32)], 0, param[c(18,19,20)], 0,
                  param[c(27,30,33)], 0, param[c(19,21,22)], 0,
                  param[c(28,31,34)], 0, param[c(20,22,23)], 0,
                  rep(0,3), param[25], rep(0,3), param[24]), byrow = TRUE, nrow = 8)
    if (class(try(chol(B), silent = TRUE)) == "try-error") stop("Input parameters doesn't definite a positive definite covariance matrix. You may remove the input param option so that the algorithm will automatically chose initial parameters.")
    param <- VarToChol(param);remove(B);
  }

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
  rcpmeObj1 <- rcpme(longdata, as.formula(paste(all.vars(formu)[1], "~", all.vars(formu)[3], "|", all.vars(formu)[4])), covariate = covariate, REadjust = REadjust, gamma = gamma, nbnodes = 20)
  rcpmeObj2 <- rcpme(longdata, as.formula(paste(all.vars(formu)[2], "~", all.vars(formu)[3], "|", all.vars(formu)[4])), covariate = covariate, REadjust = REadjust, gamma = gamma, nbnodes = 20)
  
  if (is.null(param)){
    param <- as.numeric(c(rcpmeObj1$optpar, rcpmeObj2$optpar, rep(0,10))) # en sortie de rcpme j'ai les params de Chol pour le moment donc j'y touche pas
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
    RE1 <- REestimate(rcpmeObj1, longdata, var = TRUE, onlytau = TRUE); RE2 <- REestimate(rcpmeObj2, longdata, var = TRUE, onlytau = TRUE);
    REs <- lapply(mapply(FUN = function(a,b){return(list(c(a$par, b$par, a$var, b$var)))}, RE1, RE2), function(x){return(list("par"=x[c(1,2)], "var"=matrix(c(sqrt(x[3]),0,0,sqrt(x[4])), byrow=TRUE, nrow=2)))})
    newnodes <- mapply(function(a){return(list(t(apply(nodes,1,function(x){return(sqrt(2)*a[[2]]%*%x+a[[1]])}))))},REs)
    newweights <- lapply(REs,function(a){return(weights*2*det(a[[2]])*apply(nodes,1,function(x){return(exp(t(x)%*%x))}))})
  }
  
  # optimization
  print("I can begin the optimization. Please be aware that it can take some time to run.")
  
  if (nproc == 1){
    opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust)
  }
  else {
    opt <- marqLevAlgParallel::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, nproc = nproc)
  }
  

  # OUT : fixed parameters ===========================================================================================
  if (covariate == "NULL"){
    invhessian <- CholToVarCovMatrix(opt) # JE PEUX LAISSER OU REMETTRE LE INVHESSIAN OBTENU A PARTIR DE OPT$V CAR ON SE FICHE DES IC SUR LES PARAMS DE VARIANCE (SAUF DANS LES SIMUS) COMME CA PAS DE PROBLEME DE SOUPLESSE SI PLUS D'EA DANS LA DEF DU MODELE ::

    # invhessian <- diag(34)
    # invhessian[upper.tri(invhessian, diag=TRUE)] <- opt$v
    # invhessian <- invhessian + t(invhessian) - diag(diag(invhessian)) # si je ne remets que cette partie par contre il faut enlever le invhessian en sortie car il est faux !!

    hats <- CholToVar(opt$b)

    tab <- cbind(hats[c(1,2,3,4)],sqrt(diag(invhessian)[c(1,2,3,4)]))
    tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2])
    tab <- cbind(tab, hats[c(1,2,3,4)+12],sqrt(diag(invhessian)[c(1,2,3,4)+12]))
    tab <- cbind(tab, tab[,5]-1.96*tab[,6],tab[,5]+1.96*tab[,6])
    covs <- c(invhessian[1,13], invhessian[2,14], invhessian[3,15], invhessian[4,16])
    tab <- cbind(tab, sqrt((tab[,1]-tab[,5])**2/(tab[,2]**2+tab[,6]**2-2*covs)), 1-pchisq(sqrt((tab[,1]-tab[,5])**2/(tab[,2]**2+tab[,6]**2-2*covs)), df=1))
    rownames(tab) <- c("beta0", "beta1", "beta2", "mutau")
    colnames(tab) <-  c(paste(all.vars(formu)[1],":", c("par", "se(par)", "ICinf", "ICsup"),sep = ""), paste(all.vars(formu)[2],":", c("par", "se(par)", "ICinf", "ICsup"),sep = ""), "Wald stat.", "pvalue")
  }

  # OUT : variance parameters ========================================================================================
  if (REadjust == "no" | covariate == "NULL"){
    sdres <- sqrt(hats[c(5,17)]); names(sdres) = c(all.vars(formu)[1],all.vars(formu)[2]);
    B <- matrix(c(hats[c(6,7,8)], 0, hats[c(26,27,28)], 0,
                  hats[c(7,9,10)], 0, hats[c(29,30,31)], 0,
                  hats[c(8,10,11)], 0, hats[c(32,33,34)], 0,
                  rep(0,3), hats[12], rep(0,3), hats[25],
                  hats[c(26,29,32)], 0, hats[c(18,19,20)], 0,
                  hats[c(27,30,33)], 0, hats[c(19,21,22)], 0,
                  hats[c(28,31,34)], 0, hats[c(20,22,23)], 0,
                  rep(0,3), hats[25], rep(0,3), hats[24]), byrow = TRUE, nrow = 8)
  }

  # return(opt)
  # return(list("Loglik" = opt$value, optpar= opt$par, conv = opt$convergence))
  return(list("Loglik" = opt$fn.value, "fixed" = round(tab,3), "sdres" = sdres, "VarEA" = B, optpar= hats, "covariate" = covariate, "REadjust" = REadjust, "invhessian"= invhessian, "conv" = opt$istop, "init" = CholToVar(param), "niter" = opt$iter))
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


VarToChol <- function(param){
  # this function computes a set of cholesky type parameters from a set of variance type parameters so that user
  # can input variance parameters directly
  
  B <- matrix(c(param[c(6,7,8)], 0, param[c(26,27,28)], 0, 
                param[c(7,9,10)], 0, param[c(29,30,31)], 0, 
                param[c(8,10,11)], 0, param[c(32,33,34)], 0, 
                rep(0,3), param[12], rep(0,3), param[25],
                param[c(26,29,32)], 0, param[c(18,19,20)], 0,
                param[c(27,30,33)], 0, param[c(19,21,22)], 0,
                param[c(28,31,34)], 0, param[c(20,22,23)], 0,
                rep(0,3), param[25], rep(0,3), param[24]), byrow = TRUE, nrow = 8)
  
  U <- chol(B)
  
  return(c(param[1:4], sqrt(param[5]), U[1,1], U[1,2], U[1,3], U[2,2], U[2,3], U[3,3], U[4,4],
           param[13:16], sqrt(param[17]), U[5,5], U[5,6], U[5,7], U[6,6], U[6,7], U[7,7], U[8,8],
           U[4,8], U[1,5], U[1,6], U[1,7], U[2,5], U[2,6], U[2,7], U[3,5], U[3,6], U[3,7]))
}
CholToVar <- function(x){
  
  Uepsa <- x[5]; U0a <- x[6]; U01a <- x[7]; U02a <- x[8]; U1a <- x[9]; U12a <- x[10]; U2a <- x[11]; Utaua <- x[12];
  Uepsb <- x[17]; U0b <- x[18]; U01b <- x[19]; U02b <- x[20]; U1b <- x[21]; U12b <- x[22]; U2b <- x[23]; Utaub <- x[24];
  Utau <- x[25]; U1 <- x[26]; U2 <- x[27]; U3 <- x[28]; U4 <- x[29]; U5 <- x[30]; U6 <- x[31]; U7 <- x[32]; U8 <- x[33]; U9 <- x[34];
  
  return(c(x[1:4], Uepsa**2, U0a**2, U01a * U0a, U02a * U0a, U01a**2+U1a**2, U02a * U01a + U12a * U1a, U02a**2 + U12a**2 + U2a**2, Utaua**2,
           x[13:16], Uepsb**2, U1**2 + U4**2 + U7**2 + U0b**2, U1*U2 + U4*U5 + U7*U8 + U0b*U01b, U1*U3 + U4*U6 + U7*U9 + U0b*U02b,
           U2**2+U5**2+U8**2+U01b**2+U1b**2, U3*U2 + U6*U5 + U9*U8 + U01b*U02b + U1b*U12b, U3**2 + U6**2 + U9**2 + U02b**2 + U12b**2 + U2b**2,
           Utau**2 + Utaub**2, Utaua*Utau, U0a*U1, U0a*U2, U0a*U3, U01a*U1 + U1a*U4, U01a*U2 + U1a*U5, U01a*U3 + U1a*U6, 
           U02a*U1 + U12a*U4 + U2a*U7, U02a*U2 + U12a*U5 + U2a * U8, U02a*U3 + U12a*U6 + U2a*U9))
}
CholToVarCovMatrix <- function(sim){
  
  x <- sim$b
  
  Uepsa <- x[5]; U0a <- x[6]; U01a <- x[7]; U02a <- x[8]; U1a <- x[9]; U12a <- x[10]; U2a <- x[11]; Utaua <- x[12];
  Uepsb <- x[17]; U0b <- x[18]; U01b <- x[19]; U02b <- x[20]; U1b <- x[21]; U12b <- x[22]; U2b <- x[23]; Utaub <- x[24];
  Utau <- x[25]; U1 <- x[26]; U2 <- x[27]; U3 <- x[28]; U4 <- x[29]; U5 <- x[30]; U6 <- x[31]; U7 <- x[32]; U8 <- x[33]; U9 <- x[34];
  
  vi <- diag(34)
  vi[upper.tri(vi, diag=TRUE)] <- sim$v
  vi <- vi + t(vi) - diag(diag(vi))
  
  Di <- matrix(c(1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 2*Uepsa, rep(0,34), 2*U0a, rep(0,33), U01a, U0a, rep(0,32), U02a, 0, 
                 U0a, rep(0,32), 2*U01a, 0, 2*U1a, rep(0,31), U02a, U01a, U12a, U1a, rep(0,31), 2*U02a, 0, 2*U12a, 2*U2a, rep(0,34),
                 2*Utaua, rep(0,34), 1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 2*Uepsb, rep(0,34), 2*U0b, rep(0,7),
                 2*U1, 0, 0, 2*U4, 0, 0, 2*U7, rep(0,19), U01b, U0b, rep(0,6), U2, U1, 0, U5, U4, 0, U8, U7,
                 rep(0,18), U02b, 0, U0b, rep(0,5), U3, 0, U1, U6, 0, U4, U9, 0, U7, rep(0,18), 2*U01b, 0, 2*U1b, rep(0,5), 2*U2,
                 0, 0, 2*U5, 0, 0, 2*U8, rep(0,19), U02b, U01b, U12b, U1b, rep(0,4), U3, U2, 0, U6, U5, 0, U9, U8, rep(0,19),
                 2*U02b, 0, 2*U12b, 2*U2b, rep(0,4), 2*U3, 0, 0, 2*U6, 0, 0, 2*U9, rep(0,23), 2*Utaub, 2*Utau, rep(0,20), Utau, rep(0,12),
                 Utaua, rep(0,14), U1, rep(0,19), U0a, rep(0,13), U2, rep(0,20), U0a, rep(0,12), U3, rep(0,21), U0a, rep(0,12),
                 U1,0,U4, rep(0,16), U01a,0,0,U1a, rep(0,11), U2, 0, U5, rep(0,17), U01a, 0, 0, U1a, rep(0,10), U3, 0, U6, rep(0,18),
                 U01a, 0, 0, U1a, rep(0,10), U1, 0, U4, U7, rep(0,14), U02a, 0, 0, U12a, 0, 0, U2a, rep(0,9), U2, 0, U5, U8, rep(0,15), 
                 U02a, 0, 0, U12a, 0, 0, U2a, rep(0,8), U3, 0, U6, U9, rep(0,16), U02a, 0, 0, U12a, 0, 0, U2a), byrow = TRUE, nrow = 34)
  
  nvi <- Di %*% vi %*% t(Di)
  
  return(nvi)
}