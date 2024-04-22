#' Bivariate Random Change Point Mixed Model
#'
#' @param longdata A dataframe containing the variables used in the formula  \code{formu}
#' @param formu A formula object describing which variables are to be used. The formula has to be of the following form \code{markervar1 + markervar2 ~ scorevar | grouvpar} for the function to work.
#' @param covariate An optional string indicating a binary covariate to add on the fixed effects, i.e. intercept, mean slope, difference of slopes and changepoint date. The parameter \code{REadjust} indicates how this covariate influences the random effects variance structure. Default to NULL, i.e. no covariates.
#' @param REadjust An optional string value indicating how the random effects variance structure depends on \code{covariate}. "no" means that the structure doesn't depend upon \code{covariate}. "prop" indicates that the random effects variance structure is proportionnal according to \code{covariate} value. "yes" indicates that there is two different random effects variance structures, i.e. one for each level of \code{covariate}. Default to "no".
#' @param gamma A numeric parameter indicating how smooth the trajectory is on the changepoint date. Default to 0.1.
#' @param nbnodes A numeric parameter indicating how many nodes are to be used for the gaussian quadrature for numerical integration. Default to 10.
#' @param adapt A boolean indicating whether adaptive gaussian quadrature should be used for numerical integration. Default to FALSE.
#' @param param An optional vector parameter that contains initial parameter for the optimization of the log-likelihood. Default to NULL.
#' @param nproc An optional integer specifying the number of processors for parallelisation of the optimization algorithm. Default to 1.
#' @param model An optional string indicating which formulation of the random changepoint exists. The first is `test` which is used by the `testRCPMM` function, the second is `bw` for the Bacon-Watts formulation of the model, the third is `isplines` for the I-spline model. When used for estimation purpose, you should either `bw` or `isplines` which has better interpretability properties. Default to `bw`
#' @param link1 An optional string indicating which link function is to be used for the first marker. This link function is used to deal with non-gaussian data. With `link=splines` the model estimates an appropriate I-spline link function `g` so that `g(scorevar)` is a gaussian variable. If data is already gaussian, you can chose `link=linear` so that no link function will be estimated. Default to `linear`.
#' @param link2 Same as \code{link1} but for the second marker. Default to `linear`.
#' @param twostep An optional boolean to specify if a two-step pseudo adaptive Gaussian quadrature should be used. Currently not working. Default to `FALSE`.
#' 
#' @return The output contains several objects : \code{loglik} is the value of the log-likelihood at the optimum; \code{fixed} contains all fixed parameters estimates, standard errors, CIs, wald test statistic and corresponding pvalue when possible; \code{sdres} the estimated residual error; \code{VarEA} a matrix containing the estimated random effects covariance matrix of the eight random effects: four for each marker with a general correlation structure between them; \code{optpar} the optimal parameters maximizing the log-likelihood; \code{covariate} the covariate declared in the function call; \code{REadjust} the string indicating how random effects structure is handled as declared in the function call, \code{invhessian} the covariance matrix containing all the standard errors and correlations of the parameter estimates;
#' @export
#'
#' @examples
bircpme <- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, adapt = FALSE, param = NULL, nproc = 1, model = "test", link1 = "linear", link2 = "linear", twostep = FALSE){
  
  # errors handling
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (covariate != "NULL" | REadjust != "no") stop("It has not been implemented yet. Sorry for the inconvenience.")
  if (!(link1 %in% c("linear", "splines")) | !(link2 %in% c("linear", "splines"))) stop("Link must be either linear or splines.")
  rk1 = 5 - (link1 == "splines"); rk2 = rk1 + 12 - (link2 == "splines"); rk3 = rk2 + 17; rk4 = rk3 + (link1 == "splines")*5;
  if (!is.null(param)){
    lparam = 34 + (link1 == "splines")*4 + (link2 == "splines")*4
    if (length(param) != lparam){ stop(paste("Initial parameters vector must be a vector of size ", lparam, " not ", length(param), ".", sep = ""))}
    B <- matrix(c(param[c(rk1+1,rk1+2,rk1+3)], 0, param[c(rk2+9,rk2+10,rk2+11)], 0, 
                  param[c(rk1+2,rk1+4,rk1+5)], 0, param[c(rk2+12,rk2+13,rk2+14)], 0, 
                  param[c(rk1+3,rk1+5,rk1+6)], 0, param[c(rk2+15,rk2+16,rk2+17)], 0, 
                  rep(0,3), param[rk1+7], rep(0,3), param[rk2+8],
                  param[c(rk2+9,rk2+12,rk2+15)], 0, param[c(rk2+1,rk2+2,rk2+3)], 0,
                  param[c(rk2+10,rk2+13,rk2+16)], 0, param[c(rk2+2,rk2+4,rk2+5)], 0,
                  param[c(rk2+11,rk2+14,rk2+17)], 0, param[c(rk2+3,rk2+5,rk2+6)], 0,
                  rep(0,3), param[rk2+8], rep(0,3), param[rk2+7]), byrow = TRUE, nrow = 8)
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
  rcpmeObj1 <- rcpme(longdata, as.formula(paste(all.vars(formu)[1], "~", all.vars(formu)[3], "|", all.vars(formu)[4])), covariate = covariate, REadjust = REadjust, gamma = gamma, nbnodes = 20, model = model, link = link1, statut = NULL)
  rcpmeObj2 <- rcpme(longdata, as.formula(paste(all.vars(formu)[2], "~", all.vars(formu)[3], "|", all.vars(formu)[4])), covariate = covariate, REadjust = REadjust, gamma = gamma, nbnodes = 20, model = model, link = link2, statut = NULL)
  if ((rcpmeObj1$conv > 1) | (rcpmeObj2$conv > 1)) warning("One of the two univariate model did not converge...")
  
  # initialization
  if (is.null(param)){
    if ((link1 == "linear") & (link2 == "linear")) param <- as.numeric(c(rcpmeObj1$optpar, rcpmeObj2$optpar, rep(0,10)))
    if ((link1 == "splines") & (link2 == "splines")) param <- as.numeric(c(rcpmeObj1$optpar[1:11], rcpmeObj2$optpar[1:11], rep(0,10), rcpmeObj1$optpar[12:16], rcpmeObj2$optpar[12:16]))
    if ((link1 == "splines") & (link2 == "linear")) param <- as.numeric(c(rcpmeObj1$optpar[1:11], rcpmeObj2$optpar, rep(0,10), rcpmeObj1$optpar[12:16]))
    if ((link1 == "linear") & (link2 == "splines")) param <- as.numeric(c(rcpmeObj1$optpar, rcpmeObj2$optpar[1:11], rep(0,10), rcpmeObj2$optpar[12:16]))
  }
  
  # data transformation
  objtrans1 <- datatrans(scorevar1, ngroupvar, link1)
  objtrans2 <- datatrans(scorevar2, ngroupvar, link2)
  
  # nodes and weights
  ghcoeff <- gauherm(nbnodes, 2)
  nodes <- ghcoeff$x
  weights <- ghcoeff$w
  
  # non adaptive
  if (adapt == FALSE){
    nodes <- sqrt(2) * nodes
    weights <- weights / pi
    newnodes = NULL; newweights = NULL;
    
    # optimization
    print("I can begin the optimization. Please be aware that it can take some time to run.")
    
    if (nproc == 1){
      opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, model = model, link1 = link1, link2 = link2, objtrans1 = objtrans1, objtrans2 = objtrans2, gamma = gamma)
    }
    else {
      opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, model = model, link1 = link1, link2 = link2, objtrans1 = objtrans1, objtrans2 = objtrans2, gamma = gamma, nproc = nproc)
    }
    
  }
  # adaptive
  if (adapt == TRUE){
    RE1 <- REestimate(rcpmeObj1, longdata, var = TRUE, onlytau = TRUE); RE2 <- REestimate(rcpmeObj2, longdata, var = TRUE, onlytau = TRUE);
    REs <- lapply(mapply(FUN = function(a,b){return(list(c(a$par, b$par, a$var, b$var)))}, RE1, RE2), function(x){return(list("par"=x[c(1,2)], "var"=matrix(c(sqrt(x[3]),0,0,sqrt(x[4])), byrow=TRUE, nrow=2)))})
    newnodes <- mapply(function(a){return(list(t(apply(nodes,1,function(x){return(sqrt(2)*a[[2]]%*%x+a[[1]])}))))},REs)
    newweights <- lapply(REs,function(a){return(weights*2*det(a[[2]])*apply(nodes,1,function(x){return(exp(t(x)%*%x))}))})
    
    # 1st optimization
    print("I can begin the optimization. Please be aware that it can take some time to run.")
    
    if (nproc == 1){
      opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, model = model, link1 = link1, link2 = link2, objtrans1 = objtrans1, objtrans2 = objtrans2, gamma = gamma)
    }
    else {
      opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, model = model, link1 = link1, link2 = link2, objtrans1 = objtrans1, objtrans2 = objtrans2, gamma = gamma, nproc = nproc)
    }
    
    if (twostep == TRUE){
      # updating gh nodes and weights
      print("I am updating nodes and weights of gaussian quadrature.")
      hats <- CholToVar(opt$b, link1 = link1, link2 = link2)
      B <- matrix(c(hats[c(rk1+1,rk1+2,rk1+3)], 0, hats[c(rk2+9,rk2+10,rk2+11)], 0,
                    hats[c(rk1+2,rk1+4,rk1+5)], 0, hats[c(rk2+12,rk2+13,rk2+14)], 0,
                    hats[c(rk1+3,rk1+5,rk1+6)], 0, hats[c(rk2+15,rk2+16,rk2+17)], 0,
                    rep(0,3), hats[rk1+7], rep(0,3), hats[rk2+8],
                    hats[c(rk2+9,rk2+12,rk2+15)], 0, hats[c(rk2+1,rk2+2,rk2+3)], 0,
                    hats[c(rk2+10,rk2+13,rk2+16)], 0, hats[c(rk2+2,rk2+4,rk2+5)], 0,
                    hats[c(rk2+11,rk2+14,rk2+17)], 0, hats[c(rk2+3,rk2+5,rk2+6)], 0,
                    rep(0,3), hats[rk2+8], rep(0,3), hats[rk2+7]), byrow = TRUE, nrow = 8)
      paramSpl <- NULL
      if (link1 == "splines"){paramSpl <- c(paramSpl, opt$b[(rk3+1):(rk3+5)])}
      if (link2 == "splines"){paramSpl <- c(paramSpl, opt$b[(rk4+1):(rk4+5)])}
      sdres <- c(1,1)
      if (link1 == "linear") sdres[1] = sqrt(hats[rk1])
      if (link2 == "linear") sdres[2] = sqrt(hats[rk2])
      esti <- list("call" = NULL, "Loglik" = NULL, "formula" = formu, "fixed" = NULL, "sdres"=sdres, "VarEA" = B, "optpar"= hats, "covariate" = NULL, "REadjust" = NULL, "invhessian" = NULL, "conv" = NULL, "init" = NULL, "niter" = NULL, "model" = model, "gamma" = gamma, "link1" = link1, "link2" = link2, "paramSpl" = paramSpl)
      remove(paramSpl, sdres)
      RE <- REestimate(esti, longdata, var = TRUE, onlytau = TRUE);
      RE <- lapply(RE, function(x){return(list("par"=x$par, "var"=chol(matrix(c(x$var[1], x$var[2], x$var[2], x$var[3]), nrow=2))))})
      ghcoeff <- gauherm(10, 2)
      nodes <- ghcoeff$x
      weights <- ghcoeff$w
      newnodes <- mapply(function(a){return(list(t(apply(nodes,1,function(x){return(sqrt(2)*a[[2]]%*%x+a[[1]])}))))},RE)
      newweights <- lapply(RE,function(a){return(weights*2*det(a[[2]])*apply(nodes,1,function(x){return(exp(t(x)%*%x))}))})
      param <- opt$b
      
      # 2nd and last optimization
      print("I begin the last optimization.")
      if (nproc == 1){
        opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, model = model, link1 = link1, link2 = link2, objtrans1 = objtrans1, objtrans2 = objtrans2, gamma = gamma)
        
      }
      else {
        opt <- marqLevAlg::marqLevAlg(b=param, fn=bilvsblNC, data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}), nq=nbnodes, adapt = adapt, grp=ngroupvar, weights=weights,  nodes=nodes, newnodes = newnodes, newweights = newweights, scorevar1 = all.vars(formu)[1], scorevar2 = all.vars(formu)[2], timevar = all.vars(formu)[3], covariate = covariate, REadjust = REadjust, model = model, link1 = link1, link2 = link2, objtrans1 = objtrans1, objtrans2 = objtrans2, gamma = gamma, nproc = nproc)
      }
    }
  }
  
  
  
  
  # OUT : fixed parameters ===========================================================================================
  if (covariate == "NULL"){
    invhessian <- CholToVarCovMatrix(opt, link1 = link1, link2 = link2) # JE PEUX LAISSER OU REMETTRE LE INVHESSIAN OBTENU A PARTIR DE OPT$V CAR ON SE FICHE DES IC SUR LES PARAMS DE VARIANCE (SAUF DANS LES SIMUS) COMME CA PAS DE PROBLEME DE SOUPLESSE SI PLUS D'EA DANS LA DEF DU MODELE ::
    
    # invhessian <- diag(34)
    # invhessian[upper.tri(invhessian, diag=TRUE)] <- opt$v
    # invhessian <- invhessian + t(invhessian) - diag(diag(invhessian)) # si je ne remets que cette partie par contre il faut enlever le invhessian en sortie car il est faux !!
    
    hats <- CholToVar(opt$b, link1 = link1, link2 = link2)
    
    tab <- cbind(hats[c(1,2,3,4)],sqrt(diag(invhessian)[c(1,2,3,4)]))
    tab <- cbind(tab, tab[,1]-1.96*tab[,2],tab[,1]+1.96*tab[,2])
    tab <- cbind(tab, hats[c(1,2,3,4)+12-(link1=="splines")],sqrt(diag(invhessian)[c(1,2,3,4)+12-(link1=="splines")]))
    tab <- cbind(tab, tab[,5]-1.96*tab[,6],tab[,5]+1.96*tab[,6])
    covs <- c(invhessian[1,1+12-(link1=="splines")], invhessian[2,2+12-(link1=="splines")], invhessian[3,3+12-(link1=="splines")], invhessian[4,4+12-(link1=="splines")])
    tab <- cbind(tab, sqrt((tab[,1]-tab[,5])**2/(tab[,2]**2+tab[,6]**2-2*covs)), 1-pchisq(sqrt((tab[,1]-tab[,5])**2/(tab[,2]**2+tab[,6]**2-2*covs)), df=1))
    rownames(tab) <- c("beta0", "beta1", "beta2", "mutau")
    colnames(tab) <-  c(paste(all.vars(formu)[1],":", c("par", "se(par)", "ICinf", "ICsup"),sep = ""), paste(all.vars(formu)[2],":", c("par", "se(par)", "ICinf", "ICsup"),sep = ""), "Wald stat.", "pvalue")
  }
  
  # OUT : variance parameters ========================================================================================
  if (REadjust == "no" | covariate == "NULL"){
    
    sdres <- c(1,1)
    if (link1 == "linear") sdres[1] = sqrt(hats[rk1])
    if (link2 == "linear") sdres[2] = sqrt(hats[rk2])
    names(sdres) = c(all.vars(formu)[1],all.vars(formu)[2]);
    B <- matrix(c(hats[c(rk1+1,rk1+2,rk1+3)], 0, hats[c(rk2+9,rk2+10,rk2+11)], 0,
                  hats[c(rk1+2,rk1+4,rk1+5)], 0, hats[c(rk2+12,rk2+13,rk2+14)], 0,
                  hats[c(rk1+3,rk1+5,rk1+6)], 0, hats[c(rk2+15,rk2+16,rk2+17)], 0,
                  rep(0,3), hats[rk1+7], rep(0,3), hats[rk2+8],
                  hats[c(rk2+9,rk2+12,rk2+15)], 0, hats[c(rk2+1,rk2+2,rk2+3)], 0,
                  hats[c(rk2+10,rk2+13,rk2+16)], 0, hats[c(rk2+2,rk2+4,rk2+5)], 0,
                  hats[c(rk2+11,rk2+14,rk2+17)], 0, hats[c(rk2+3,rk2+5,rk2+6)], 0,
                  rep(0,3), hats[rk2+8], rep(0,3), hats[rk2+7]), byrow = TRUE, nrow = 8)
  }
  
  paramSpl <- NULL
  if (link1 == "splines"){paramSpl <- c(paramSpl, opt$b[(rk3+1):(rk3+5)])}
  if (link2 == "splines"){paramSpl <- c(paramSpl, opt$b[(rk4+1):(rk4+5)])}
  
  return(list("call" = as.list(match.call()), "Loglik" = opt$fn.value, "formula" = formu, "fixed" = round(tab,3), "sdres"=sdres, "VarEA" = B, optpar= hats, "covariate" = covariate, "REadjust" = REadjust, "invhessian" = invhessian, "conv" = opt$istop, "init" = CholToVar(param), "niter" = opt$iter, "model" = model, "gamma" = gamma, "link1" = link1, "link2" = link2, "paramSpl" = paramSpl))
  
}
