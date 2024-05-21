

semi_latent_regularized <- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, param = NULL, model = "test", link = "linear", statut = NULL, classprob = NULL, lambda = 0) {
  
  # prepare the parameters
  
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (REadjust == "prop") stop("It has not been implemented yet. Sorry for the inconvenience...")
  lparam = 12 +2*(model == "isplines") + 4*(link=="splines") + (covariate != "NULL")*(4*(REadjust=="no")+ 5*(REadjust=="prop")+11*(REadjust=="yes"))
  if (!is.null(param) & (length(param) != lparam)) stop(paste("Initial parameters vector must be a vector of size ", lparam, ".", sep = ""))
  rk0 = 3 - (model == "isplines")
  rk1 = 4 + (link == "linear") - (model == "isplines")
  rk2 = rk1 + 7 + (covariate !="NULL")*(4*(REadjust=="no") + 5*(REadjust=="prop") + 11*(REadjust=="yes"))
  lgtclassprob = length(all.vars(classprob))
  
  #########
  
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
  objtrans <- datatrans(scorevar, ngroupvar, link)
  
  # Repeat the process for cases
  scorevar2 <- longdata2[[score_var_name]]
  timevar2 <- longdata2[[time_var_name]]
  groupvar2 <- longdata2[[group_var_name]]
  ngroupvar2 <- rep(seq_along(unique(groupvar2)), times = table(groupvar2))
  longdata2 <- cbind(longdata2, ngroupvar2)
  objtrans2 <- datatrans(scorevar2, ngroupvar2, link)
  
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
      
      param <- c(lmm$coefficients$fixed[1:2], -0.5, median(timevar), nlme::VarCorr(lmm)[c(3,1),2],0,0,1,0,1,1)
      param <- as.numeric(param)
    }
    
    if (link == "splines"){
      param <- c(param[-5], rep(1,5))
    }
    
    if (model == "isplines"){
      param <- c(param[-3], rep(2,3))
    }
    
    if (!is.null(classprob)){
      param <- c(param, rep(0, lgtclassprob+1))
    }
  }
  
  opt <- marqLevAlg(b=param,fn=lvsblclass, minimize = FALSE, data1=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}),data2=by(longdata2,longdata2[,"ngroupvar2"],function(x){return(x)}),nq=nbnodes,grp=ngroupvar,grp2=ngroupvar2,weights=weights, nodes=nodes, scorevar = all.vars(formu)[1], timevar = all.vars(formu)[2], covariate = covariate, REadjust = REadjust, model = model, link = link, objtrans = objtrans, objtrans2 = objtrans2, gamma = gamma, latent = latent, classprob = classprob)

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
  
  
  

  