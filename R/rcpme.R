#' Random Change Point Mixed Model
#'
#' @param longdata A longitudinal dataset containing all variables used in the formula \code{formu}
#' @param formu A formula object describing which variables are to be used. The formula has to be of the following form \code{markervar ~ scorevar | groupvar} for the function to work.
#' @param covariate An optional string indicating a binary covariate to add on the fixed effects, i.e. intercept, mean slope, difference of slopes and changepoint date. The parameter \code{REadjust} indicates how this covariate influences the random effects variance structure. Default to NULL, i.e. no covariates.
#' @param REadjust An optional string value indicating how the random effects variance structure depends on \code{covariate}. "no" means that the structure doesn't depend upon \code{covariate}. "prop" indicates that the random effects variance structure is proportional according to \code{covariate} value. "yes" indicates that there is two different random effects variance structures, i.e. one for each level of \code{covariate}. Default to "no".
#' @param gamma A numeric parameter indicating how smooth the trajectory is on the changepoint date. It should be small according to the time variable scale. Default to 0.1.
#' @param nbnodes A numeric parameter indicating how many nodes are to be used for the gaussian quadrature for numerical integration. Default to 10.
#' @param param An optional vector parameter that contains initial parameter for the optimization of the log-likelihood. Default to NULL.
#' @param model An optional string indicating which formulation of the random changepoint exists. The first model, `test`, is \eqn{Y_ij = \beta_{0i} + \beta_{1i}t_{ij} + \beta_{2i}\sqrt{(t_{ij}-\tau_i)^2+\gamma} + \epsilon_{ij}} used by the `testRCPMM` function. The second is `bw` for the Bacon-Watts formulation of the model \eqn{Y_ij = \beta_{0i} + \beta_{1i}(t_{ij}-\tau_i) + \beta_{2i}\sqrt{(t_{ij}-\tau_i)^2+\gamma} + \epsilon_{ij}}. The third option is `isplines` for the I-spline model. When used for estimation purpose, you should either `bw` or `isplines` which has clear interpretability properties. Default to `bw`
#' @param link An optional string indicating which link function is to be used. This link function is used to deal with non-gaussian data. With `link=splines` the model estimates an appropriate I-spline link function `g` so that `g(scorevar)` is a gaussian variable. If data is already gaussian, you can chose `link=linear` so that no link function will be estimated. Default to `linear`.
#' @param statut An optional string indicating a binary variable from which two class are considered: a linear class for subjects with \code{statut=0} and a random changepoint class for subjects with \code{statut=1}. Default to NULL.
#' @param membership an optional string indicating whether the proportion of the controls belonging to the class with an accelerating decline. This proportion can take values from zero to one, with respectively 0 meaning non of the controls have an alternative trajectory while 1 means that all controls have an accelerated trajectory. "NULL" means that the proportion will be estimated using the semi-latent class.
#' @param two-means a variable indicating whether the random changepoint should be estimated using two means. 
#' @return The output contains several objects : \code{call} is the function call; \code{loglik} is the value of the log-likelihood at the optimum; \code{formula} is the formula describing which variables are used in the model; \code{fixed} contains all fixed parameters estimates, standard errors, CIs, wald test statistic and corresponding pvalue when possible; \code{sdres} the estimated residual error; \code{VarEA} a 4x4 matrix or a list of 4x4 matrices - if there is some covariate for example - containing the estimated random effects covariance matrix; \code{optpar} the optimal parameters maximizing the log-likelihood; \code{covariate} the covariate declared in the function call; \code{REadjust} the string indicating how random effects structure is handled as declared in the function call, \code{invhessian} the covariance matrix containing all the standard errors and correlations of the parameter estimates; \code{conv} an index of successful convergance, equals to 1 if success; \code{init} the initial values vector; \code{model} the model used during estimation; \code{gamma} the value of gamma used during estimation; \code{link} the link function used during estimation.
#' @export
#'
#' @examples

library(marqLevAlg)
print('succesfully loaded')
rcpme<- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, param = NULL, model = "test", link = "linear", statut = NULL, latent = FALSE, classprob = NULL, two_means = FALSE, membership = NULL, lambda = 0, only_cases = FALSE, maxiter = 500, epsa = 1e-04, nproc = 1, verbose = FALSE) {
  
  if(!is.null(membership)){
    if((membership <= 0) | (membership >= 1)) stop("membership is not well defined")
  }
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (REadjust == "prop") stop("It has not been implemented yet. Sorry for the inconvenience...")
  lparam = 12 +2*(model == "isplines") + 4*(link=="splines") + (covariate != "NULL")*(4*(REadjust=="no")+ 5*(REadjust=="prop")+11*(REadjust=="yes"))
  if (!is.null(param) & (length(param) != lparam)) stop(paste("Initial parameters vector must be a vector of size ", lparam, ".", sep = ""))
  rk0 = 3 - (model == "isplines")
  rk1 = 4 + (link == "linear") - (model == "isplines")
  rk2 = rk1 + 7 + (covariate !="NULL")*(4*(REadjust=="no") + 5*(REadjust=="prop") + 11*(REadjust=="yes"))
  lgtclassprob = length(all.vars(classprob))
  
  # =============================================
  
  if(!is.null(statut)){
    # Extract the variable names from the formula
    score_var_name = all.vars(formu)[1]
    time_var_name = all.vars(formu)[2]
    group_var_name = all.vars(formu)[3]
    
    # Create subsets for control and cases
    longdata1 <- longdata[longdata$statut == 0,]
    longdata2 <- longdata[longdata$statut == 1,]
    
    # Extract score, time, and group variables for controls
    scorevar <- longdata1[[score_var_name]]
    timevar <- longdata1[[time_var_name]]
    groupvar <- longdata1[[group_var_name]]
    
    # Calculate the number of repetitions for each group for controls
    ngroupvar <- rep(seq_along(unique(groupvar)), times = table(groupvar))
    
    # Bind the new variable to the data frame for controls
    longdata <- cbind(longdata1, ngroupvar)
    
    # Apply the transformation (assumed function datatrans available)
    objtrans <- datatrans(scorevar, ngroupvar, link)
    
    # Repeat the process for cases
    scorevar2 <- longdata2[[score_var_name]]
    timevar2 <- longdata2[[time_var_name]]
    groupvar2 <- longdata2[[group_var_name]]
    ngroupvar2 <- rep(seq_along(unique(groupvar2)), times = table(groupvar2))
    longdata2 <- cbind(longdata2, ngroupvar2)
    objtrans2 <- datatrans(scorevar2, ngroupvar2, link)
  }
  
  else {
    scorevar = longdata[,all.vars(formu)[1]]
    timevar = longdata[,all.vars(formu)[2]]
    groupvar = longdata[,all.vars(formu)[3]]
    ngroupvar = rep(seq(length(unique(groupvar))), lapply(split(longdata, groupvar),function(x){return(dim(x)[[1]])})) 
    longdata <- cbind(longdata, ngroupvar)
    objtrans <- datatrans(scorevar, ngroupvar, link)
  }
  
  if (covariate != "NULL") {
    adjustvar = longdata[,covariate]
  } 
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
      # lmm <- lcmm(fixed = score ~ date,
      #             random = ~ date, 
      #             subject = all.vars(formu)[3], 
      #             data = longdata, ng = 1, link="3-quant-splines")
      
      
      param <- c(lmm$coefficients$fixed[1:2], -0.5, median(timevar), nlme::VarCorr(lmm)[c(3,1),2],0,0,1,0,1,1)
      param <- as.numeric(param)
    }
    
    if (link == "splines"){
      param <- c(param[-5], rep(1,5))
      # param <- c(param, rep(1,5))
    }
    
    if (model == "isplines"){
      param <- c(param[-3], rep(2,3))
    }
    if (two_means) {
      param <- c(param, 15, 10)
    }
    
    if (!is.null(classprob) & latent == TRUE){
      param <- c(param, rep(0, lgtclassprob+1))
    }
  }
  #
  
  if(is.null(statut)){
    opt <- marqLevAlg(b=param,fn=lvsblNCgenR, minimize = FALSE, data=split(longdata, longdata[,"ngroupvar"]),nq=nbnodes,grp=ngroupvar,weights=weights, nodes=nodes, scorevar = all.vars(formu)[1], timevar = all.vars(formu)[2], covariate = covariate, REadjust = REadjust, model = model, link = link, objtrans = objtrans, gamma = gamma, loglik = TRUE, maxiter = maxiter, epsa = epsa, epsb = epsa, epsd = epsa, nproc = nproc, print.info = verbose)
  }
  else {
    # optimisation du melange
    opt <- marqLevAlg(b=param,fn=lvsblclass_penalized , minimize = FALSE, data1=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}),data2=by(longdata2,longdata2[,"ngroupvar2"],function(x){return(x)}),nq=nbnodes,grp=ngroupvar,grp2=ngroupvar2,weights=weights, nodes=nodes, scorevar = all.vars(formu)[1], timevar = all.vars(formu)[2], covariate = covariate, REadjust = REadjust, model = model, link = link, objtrans = objtrans, objtrans2 = objtrans2, gamma = gamma, latent = latent, classprob = classprob, two_means = two_means, membership = membership, lambda = lambda, only_cases = only_cases, maxiter = maxiter, epsa = epsa, epsb = epsa, epsd = epsa, print.info = verbose)
  }
  
  # OUT : fixed parameters ===========================================================================================
  
  if (covariate != "NULL"){
    invhessian <- diag(lparam)
    
    invhessian[upper.tri(invhessian, diag=TRUE)] <- opt$v
    invhessian <- invhessian + t(invhessian) - diag(diag(invhessian)) 
    tab <- cbind(opt$b[c(1,rk2+1,2,rk2+2,3,rk2+3,4,rk2+4)],sqrt(diag(invhessian))[c(1,rk2+1,2,rk2+2,3,rk2+3,4,rk2+4)])
    tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2], tab[,1]/tab[,2],1-pchisq(tab[,1]**2/tab[,2]**2,df=1)); tab[c(5,7),c(5,6)] = NA;
    rownames(tab) <- c("beta0", paste("beta0*",covariate),"beta1", paste("beta1*",covariate),"beta2", paste("beta2*",covariate),"mutau", paste("mutau*",covariate))
    colnames(tab) <- c("par", "se(par)", "ICinf", "ICsup", "Wald stat.", "pvalue")
  }
  
  else {
    invhessian <- diag(lparam)
    invhessian[upper.tri(invhessian, diag=TRUE)] <- opt$v
    invhessian <- invhessian + t(invhessian) - diag(diag(invhessian))
    if (model == "isplines"){
      tab <- cbind(c(opt$b[c(1,2,rk0+1)]),c(sqrt(diag(invhessian))[c(1,2,rk0+1)]))
      tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2], tab[,1]/tab[,2], 1-pchisq(tab[,1]**2/tab[,2]**2,df=1))
      tab[3,c(5,6)] = NA
      rownames(tab) <- c("beta0", "beta1", "mutau")
    }
    if (model != "isplines"){
      tab <- cbind(opt$b[c(1,2,3,4)],sqrt(diag(invhessian))[c(1,2,3,4)])
      tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2], tab[,1]/tab[,2], 1-pchisq(tab[,1]**2/tab[,2]**2,df=1))
      tab[c(3,4),c(5,6)] = NA
      rownames(tab) <- c("beta0", "beta1", "beta2", "mutau")
    }
    colnames(tab) <- c("par", "se(par)", "ICinf", "ICsup", "Wald stat.", "pvalue")
  }
  
  # OUT : variance parameters ===========================================================================================
  
  if (REadjust == "no" | covariate == "NULL"){
    if (link == "linear") seps <- abs(opt$b[rk1]) else seps = 1;
    U0 <- opt$b[rk1+1]; U01 <- opt$b[rk1+2]; U02 <- opt$b[rk1+3]; U1 <- opt$b[rk1+4]; U12 <- opt$b[rk1+5]; U2 <- opt$b[rk1+6]; stau <- abs(opt$b[rk1+7]);
    U <- matrix(c(U0,U01,U02,0,U1,U12,0,0,U2), byrow = TRUE, nrow=3); 
    B <- t(U) %*% U;
    vars <- data.frame(c(seps, sqrt(diag(B)), stau, B[1,2:3], B[2,3]))
    rownames(vars) <- c("sdres", "sd0", "sd1", "sd2", "sdtau", "cov01", "cov02", "cov12")
    colnames(vars) <- c("par")
    VarEA <- matrix(0, nrow=4,ncol=4); VarEA[1:3,1:3] <- B; VarEA[4,4] <- stau**2;
  }
  
  if (REadjust == "yes"){###ATTENTION LINEARE SPLINES EPSILON !!
    if (link == "linear") seps <- abs(opt$b[rk1]) else seps = 1; 
    U0 <- opt$b[rk1+1]; U01 <- opt$b[rk1+2]; U02 <- opt$b[rk1+3]; U1 <- opt$b[rk1+4]; U12 <- opt$b[rk1+5]; U2 <- opt$b[rk1+6]; utau <- abs(opt$b[rk1+7]);
    U <- matrix(c(U0,U01,U02,0,U1,U12,0,0,U2),byrow = TRUE, nrow=3); 
    B <- t(U) %*% U;
    V0 <- opt$b[rk1+1]+opt$b[rk1+12]; V01 <- opt$b[rk1+2]+opt$b[rk1+13]; V02 <- opt$b[rk1+3]+opt$b[rk1+14]; V1 <- opt$b[rk1+4]+opt$b[rk1+15]; V12 <- opt$b[rk1+5]+opt$b[rk1+16]; V2 <- opt$b[rK1+6]+opt$b[rk1+17]; vtau <- abs(opt$b[rk1+7]+opt$b[rk1+18]);
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
    if (link == "linear") seps <- abs(opt$b[rk1]); 
    U0 <- opt$b[rk1+1]; U01 <- opt$b[rk1+2]; U02 <- opt$b[rk1+3]; U1 <- opt$b[rk1+4]; U12 <- opt$b[rk1+5]; U2 <- opt$b[rk1+6]; utau <- abs(opt$b[rk1+7]);
    U <- matrix(c(U0,U01,U02,0,U1,U12,0,0,U2),byrow = TRUE,nrow=3); 
    B <- t(U) %*% U;
    coeff <- abs(1+param[rk1+12])
    vars <- data.frame(c(seps, sqrt(diag(B)), utau, B[1,2:3], B[2,3], coeff))
    rownames(vars) <- c("sdres", "sd0", "sd1", "sd2", "sdtau", "cov01", "cov02", "cov12", "prop.")
    colnames(vars) <- c("par")
  }
  
  return(list("call" = as.list(match.call()), "Loglik" = opt$fn.value, "formula" = formu, "fixed" = round(tab,3), "sdres"=seps, "VarEA" = VarEA, optpar= opt$b, "covariate" = covariate, "REadjust" = REadjust, "invhessian" = invhessian, "conv" = opt$istop, "init" = param, "model" = model, "statut" = statut, "gamma" = gamma, "link" = link))
}