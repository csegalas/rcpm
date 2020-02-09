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

CPmodel <- function(a,b,c,d,t,gamma,model="test"){
  if (model=="test") return(a + b*t + c*sqrt((t-d)**2+gamma))
  if (model=="bw") return(a + b*(t-d) + c*sqrt((t-d)**2+gamma))
}

geneData <- function(n = 100, hyp = "null", params, pNA = 0, DO = 0, vis = c(0, 3, 6, 9, 12, 15, 18, 21), beta2alea=FALSE, model = "test", gamma = 0.1){
  
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
        # generation des effets aleatoires pour l'individu i
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
        # generation des effets aleatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0 + bis[i,1], Beta1 + bis[i,2], Beta2, mutau + bis[i,3], vis[j], gamma, model = model) + eps
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
        # generation des effets aleatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0,0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0 + bis[i,1], Beta1 + bis[i,2], Beta2 + bis[i,3], mutau + bis[i,4], vis[j], gamma, model = model) + eps
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
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0_1 + bis[i,1], Beta1_1 + bis[i,2], Beta2_1 + bis[i,3], mutau_1 + bis[i,4], vis[j], gamma, model = model) + eps1
          donnees[(nbvis*(i-1))+j,4] = CPmodel(Beta0_2 + bis[i,5], Beta1_2 + bis[i,6], Beta2_2 + bis[i,7], mutau_2 + bis[i,8], vis[j], gamma, model = model) + eps2
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

geneDataSpl  <- function(n = 100, params, pNA = 0, DO = 0, vis = c(0, 3, 6, 9, 12, 15, 18, 21), gamma = 0.1, pcas){
  
  nbvis = length(vis)
  donnees <- matrix(NA, nrow = n * nbvis, ncol = 4)
  colnames(donnees) <- {c("ID", "date", "score", "statut")}
  
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
      # cas ou temoin
      statut = rbinom(1, 1, pcas)
      donnees[((nbvis*(i-1))+1):(nbvis*i),4] <- statut 
      
      # generation des effets aleatoires pour l'individu i
      bis[i,] <- rmvnorm(1, mean = c(0,0,0,0), sigma = B)
      donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
      donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
      for (j in seq(nbvis)){
        eps = rnorm(1,0,sigma)
        donnees[(nbvis*(i-1))+j,3] = Beta0 + bis[i,1] + (Beta1 + bis[i,2])*vis[j] + (statut==1)*(Beta2 + bis[i,3])*(sum(4*ispline(vis[j] - (mutau + bis[i,4]),0,20,5))) + eps
      }
    }
  }
  
  donnees <- as.data.frame(donnees)
  if ((pNA!=0) & (DO !=0)){
    stop("pNA and DO mustn\'t be both positive")
  }
  if ((pNA<0) | (pNA>1) | (DO<0) | (DO>1)){
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


plottrans <- function(rcpmeObj, longdata){ # plotte la transfo crude / gaussian
  varsformu = all.vars(rcpmeObj$formula)
  
  if (length(varsformu) == 3){  # univ
    par(mfrow=c(1,1))
    scorevar = varsformu[1]
    isOut = iSpline(sort(longdata[,scorevar][!is.na(longdata[,scorevar])]), knots=quantile(longdata[,scorevar], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    plot(sort(longdata[,scorevar][!is.na(longdata[,scorevar])]), isOut %*% tail(rcpmeObj$optpar,5)**2, type = "l", xlab = paste("crude", scorevar), ylab = paste("transformed", scorevar))
  }
  
  if (length(varsformu) == 4){ # biv
    par(mfrow=c(1,2))
    scorevar1 = varsformu[1]
    scorevar2 = varsformu[2]
    isOut1 = iSpline(sort(longdata[,scorevar1][!is.na(longdata[,scorevar1])]), knots=quantile(longdata[,scorevar1], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    isOut2 = iSpline(sort(longdata[,scorevar2][!is.na(longdata[,scorevar2])]), knots=quantile(longdata[,scorevar2], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    plot(sort(longdata[,scorevar1][!is.na(longdata[,scorevar1])]), isOut1 %*% rcpmeObj$paramSpl[1:5]**2, type = "l", xlab = paste("crude", scorevar1), ylab = paste("transformed", scorevar1))
    plot(sort(longdata[,scorevar2][!is.na(longdata[,scorevar2])]), isOut2 %*% rcpmeObj$paramSpl[6:10]**2, type = "l", xlab = paste("crude", scorevar2), ylab = paste("transformed", scorevar2))
    par(mfrow=c(1,1))
  }
}

datatrans <- function(Y, ngroupvar, trans){
  if (trans == "linear"){
    return(by(Y, ngroupvar, function(x) return(list("isOut" = NULL, "isOutDeriv" = NULL))))
  }
  if (trans == "splines"){
    isOut = iSpline(Y, knots=quantile(Y, probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    # isOut = iSpline(Y[!is.na(Y)], knots=c(5,10), Boundary.knots = c(-10,30), degree=2, derivs=0, intercept = T)
    isOutDeriv = deriv(isOut)
    return(by(cbind(isOut, isOutDeriv), ngroupvar, function(x){return(list("isOut"=na.omit(as.matrix(x[,1:5])), "isOutDeriv"=na.omit(as.matrix(x[,6:10]))))}))
  }
}

transY <- function(rcpmeObj, longdata){ # calcule les scores transformes a partir de l'esti Spl d'un modele
  varsformu = all.vars(rcpmeObj$formula)
  if (length(varsformu) == 3){  # univ
    scorevar = varsformu[1]
    isOut = iSpline(longdata[,scorevar], knots=quantile(longdata[,scorevar], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    tY = isOut %*% tail(rcpmeObj$optpar,5)**2
  }
  if (length(varsformu) == 4){ # biv
    scorevar1 = varsformu[1]
    scorevar2 = varsformu[2]
    isOut1 = iSpline(longdata[,scorevar1], knots=quantile(longdata[,scorevar1], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    isOut2 = iSpline(longdata[,scorevar2], knots=quantile(longdata[,scorevar2], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
    tY1 = isOut1 %*% rcpmeObj$paramSpl[1:5]**2
    tY2 = isOut2 %*% rcpmeObj$paramSpl[6:10]**2
    tY = cbind(tY1, tY2)
  }
  return(tY)
}

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
#'
#' @return The output contains several objects : \code{call} is the function call; \code{loglik} is the value of the log-likelihood at the optimum; \code{formula} is the formula describing which variables are used in the model; \code{fixed} contains all fixed parameters estimates, standard errors, CIs, wald test statistic and corresponding pvalue when possible; \code{sdres} the estimated residual error; \code{VarEA} a 4x4 matrix or a list of 4x4 matrices - if there is some covariate for example - containing the estimated random effects covariance matrix; \code{optpar} the optimal parameters maximizing the log-likelihood; \code{covariate} the covariate declared in the function call; \code{REadjust} the string indicating how random effects structure is handled as declared in the function call, \code{invhessian} the covariance matrix containing all the standard errors and correlations of the parameter estimates; \code{conv} an index of successful convergance, equals to 1 if success; \code{init} the initial values vector; \code{model} the model used during estimation; \code{gamma} the value of gamma used during estimation; \code{link} the link function used during estimation.
#' @export
#'
#' @examples
rcpme <- function(longdata, formu, covariate = "NULL", REadjust = "no", gamma = 0.1, nbnodes = 10, param = NULL, model = "test", link = "linear", statut = NULL){
  
  if (covariate == "NULL" & REadjust != "no") stop("Need a covariate to adjust random effects variance structure.")
  if (REadjust == "prop") stop("It has not been implemented yet. Sorry for the inconvenience...")
  lparam = 12 +2*(model == "isplines") + 4*(link=="splines") + (covariate != "NULL")*(4*(REadjust=="no")+ 5*(REadjust=="prop")+11*(REadjust=="yes"))
  if (!is.null(param) & (length(param) != lparam)) stop(paste("Initial parameters vector must be a vector of size ", lparam, ".", sep = ""))
  rk0 = 3 - (model == "isplines")
  rk1 = 4 + (link == "linear") - (model == "isplines")
  rk2 = rk1 + 7 + (covariate !="NULL")*(4*(REadjust=="no") + 5*(REadjust=="prop") + 11*(REadjust=="yes"))
  
  # =============================================
  
  if(!is.null(statut)){
    # on construit la base control et la base cas
    longdata1 = longdata[longdata[statut] == 0,]
    longdata2 = longdata[longdata[statut] == 1,]
    
    # on cree les objets pr les controls
    scorevar = longdata1[,all.vars(formu)[1]]
    timevar = longdata1[,all.vars(formu)[2]]
    groupvar = longdata1[,all.vars(formu)[3]]
    ngroupvar = rep(seq(length(unique(groupvar))), by(longdata1,groupvar,function(x){return(dim(x)[[1]])}))
    longdata <- cbind(longdata1, ngroupvar)
    objtrans <- datatrans(scorevar, ngroupvar, link)
    
    # et pour les cas
    scorevar2 = longdata2[,all.vars(formu)[1]]
    timevar2 = longdata2[,all.vars(formu)[2]]
    groupvar2 = longdata2[,all.vars(formu)[3]]
    ngroupvar2 = rep(seq(length(unique(groupvar2))), by(longdata2,groupvar2,function(x){return(dim(x)[[1]])}))
    longdata2 <- cbind(longdata2, ngroupvar2)
    objtrans2 <- datatrans(scorevar2, ngroupvar2, link)
  }
  
  if(is.null(statut)){
    scorevar = longdata[,all.vars(formu)[1]]
    timevar = longdata[,all.vars(formu)[2]]
    groupvar = longdata[,all.vars(formu)[3]]
    ngroupvar = rep(seq(length(unique(groupvar))), by(longdata,groupvar,function(x){return(dim(x)[[1]])}))
    longdata <- cbind(longdata, ngroupvar)
    objtrans <- datatrans(scorevar, ngroupvar, link)
  }

  if (covariate != "NULL") adjustvar = longdata[,covariate]
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
  }


  if(is.null(statut)){
    opt <- marqLevAlg(b=param,fn=lvsblNCgen,data=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}),nq=nbnodes,grp=ngroupvar,weights=weights, nodes=nodes, scorevar = all.vars(formu)[1], timevar = all.vars(formu)[2], covariate = covariate, REadjust = REadjust, model = model, link = link, objtrans = objtrans, gamma = gamma)
  }
  
  if(!is.null(statut)){
    # optimisation du melange
    opt <- marqLevAlg(b=param,fn=lvsblclass,data1=by(longdata,longdata[,"ngroupvar"],function(x){return(x)}),data2=by(longdata2,longdata2[,"ngroupvar2"],function(x){return(x)}),nq=nbnodes,grp=ngroupvar,grp2=ngroupvar2,weights=weights, nodes=nodes, scorevar = all.vars(formu)[1], timevar = all.vars(formu)[2], covariate = covariate, REadjust = REadjust, model = model, link = link, objtrans = objtrans, objtrans2 = objtrans2, gamma = gamma)
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
  
  return(list("call" = as.list(match.call()), "Loglik" = opt$fn.value, "formula" = formu, "fixed" = round(tab,3), "sdres"=seps, "VarEA" = VarEA, optpar= opt$b, "covariate" = covariate, "REadjust" = REadjust, "invhessian" = invhessian, "conv" = opt$istop, "init" = param, "model" = model, "gamma" = gamma, "link" = link))
}


lvsblclass <- function(param, data1, data2, nq, grp, grp2, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, objtrans2, gamma, latent = FALSE){

  rk1 = 4 + (link == "linear") - (model == "isplines")

  # on vire les parametres inutiles ds l'esti du modele lineaire
  if (link == "linear"){
    param2 <- param[c(1,2, rk1, rk1+1, rk1+2, rk1+4)]
  }
  if (link == "isplines"){ # a finir celui ci...
    param2 <- param[c(1,2,rk1,rk1+1,rk1+2,rk1+3)]
  }
  
  loglik1 <- lvsblNCgen(param, data2, nq, grp2, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans2, gamma)
  
  #probi = 0
  #if (latent = TRUE){
  #  probi = ()
  #}
  
  # loglik2 <- lvsbllin(param2, data1, nq, grp, weights, nodes, scorevar, timevar, link, objtrans)
  loglik2 <- lvsbllin(param2, data1, grp, scorevar, timevar, link, objtrans)
  out <- loglik1 + loglik2
  # 
  # if (latent = TRUE){
  #   probi <- 
  #   loglik3 <- lvsblNCgen(param, data1, nq, grp, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma)
  #   out = loglik1 + ()
  # }
  
  return(out)
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
  model <- rcpmeObj$model
  if (is.null(model)){model <- "test"} # pour les vieilles sorties sans "model"
  gamma <- rcpmeObj$gamma
  if (is.null(gamma)){gamma <- 0.1} # pour les vieilles sorties sans "gamma"
  varsformu = all.vars(formu)
  
  if (length(varsformu) == 3){ # if univariate
    
    scorevar = varsformu[1]
    timevar = varsformu[2]
    groupvar = varsformu[3]
    
    link <- rcpmeObj$link
    if (is.null(gamma)){link <- "linear"} # pour les vieilles sorties sans "link"
    
    if (link == "splines"){
      isOut = iSpline(longdata[,scorevar][!is.na(longdata[,scorevar])], knots=quantile(longdata[,scorevar], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
      longdata[,scorevar][!is.na(longdata[,scorevar])] <- isOut %*% rcpmeObj$optpar[12:16]**2
    }   
    
    if (onlytau == FALSE){
      re <- rep(0,4)
      if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop == 1) {return("par"=opt$b);} else {return("par"=re);}} else {return("par"=re);}}))}
      if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop == 1) {return(list("par"=opt$b,"var"=matrix(c(opt$v[c(1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10)]), nrow = 4, byrow=TRUE)));} else {return(list("par"=re, "var"=diag(4)))}} else {return(list("par"=re, "var"=diag(4)))}}))}
    }
    
    else {
      re <- 0
      if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop ==1) {return("par"=opt$b);} else {return("par"=re);}} else {return("par"=re);}}))}
      if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop == 1) {return(list("par"=opt$b,"var"=opt$v));} else {return(list("par"=0, "var"=1))};} else {return(list("par"=0, "var"=1))}}))}
    }
  }
  
  if (length(varsformu) == 4){ # if bivariate
    scorevar1 = varsformu[1]
    scorevar2 = varsformu[2]
    timevar = varsformu[3]
    groupvar = varsformu[4]
    link1 = rcpmeObj$link1
    link2 = rcpmeObj$link2

    if (link1 == "splines"){
      isOut1 = iSpline(longdata[,scorevar1][!is.na(longdata[,scorevar1])], knots=quantile(longdata[,scorevar1], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
      longdata[,scorevar1][!is.na(longdata[,scorevar1])] <- isOut1 %*% rcpmeObj$paramSpl[1:5]**2
    }   

    if (link2 == "splines"){
      isOut2 = iSpline(longdata[,scorevar2][!is.na(longdata[,scorevar2])], knots=quantile(longdata[,scorevar2], probs=c(1/3,2/3), na.rm=T), degree=2, derivs=0, intercept = T)
      longdata[,scorevar2][!is.na(longdata[,scorevar2])] <- isOut2 %*% rcpmeObj$paramSpl[6:10]**2
    }   
  
    if (onlytau == FALSE){
      if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- rep(0,8)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2, blinding= TRUE);
          if(opt$istop == 1) {return("par"=opt$b);} 
          else {return("par"=re);}
        }
        if(((sum(!is.na(x[,scorevar1]))==0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- rep(0,4)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2, blinding= TRUE);
          if(opt$istop == 1) {return("par"=c(rep(0,4),opt$b));} 
          else {return("par"=rep(0,8));}
        }
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))==0)) {
          re <- rep(0,4)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2, blinding= TRUE);
          if(opt$istop == 1) {return("par"=c(opt$b,rep(0,4)));} 
          else {return("par"=rep(0,8));}
        }
        else {return("par"=rep(0,8));}
      }))}
      
      if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- rep(0,8)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return(list("par"=opt$b,"var"=matrix(c(opt$v[c(1,2,4,7,11,16,22,29,2,3,5,8,12,17,23,30,4,5,6,9,13,18,24,31,7,8,9,10,14,19,25,32,11,12,13,14,15,20,26,33,16,17,18,19,20,21,27,34,22,23,24,25,26,27,28,35,29,30,31,32,33,34,35,36)]), nrow = 8, byrow=TRUE)));} 
          else {return(list("par"=re, "var"=diag(8)))}
      }
        if(((sum(!is.na(x[,scorevar1]))==0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- rep(0,4)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return(list("par"=c(rep(0,4),opt$b), "var" = rbind(cbind(diag(4), matrix(rep(0,16),nrow=4)), cbind(matrix(rep(0,16),nrow=4), matrix(opt$v[c(1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10)], nrow = 4, byrow = TRUE)))));} 
          else {return(list("par"=rep(0,8), "var"=diag(8)))} 
      }
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))==0)) {
          re <- rep(0,4)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return(list("par"=c(opt$b,rep(0,4)), "var" = rbind(cbind(matrix(opt$v[c(1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10)], nrow = 4, byrow = TRUE), matrix(rep(0,16),nrow=4)), cbind(matrix(rep(0,16),nrow=4), diag(4)))));} 
          else {return(list("par"=rep(0,8), "var"=diag(8)))} 
    }
        else {return(list("par"=rep(0,8), "var"=diag(8)))} 
      }))}
    
    }
    
    else {
      re <- c(0,0)  
      if (var == FALSE) {return(by(longdata,longdata[,groupvar],function(x){
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- rep(0,2)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return("par"=opt$b);} 
          else {return("par"=c(0,0))};}
        if(((sum(!is.na(x[,scorevar1]))==0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- 0
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return("par"=c(0,opt$b));} 
          else {return("par"=c(0,0))};}
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))==0)) {
          re <- 0
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return("par"=c(opt$b,0));} 
          else {return("par"=c(0,0))};}
        else {return("par"=c(0,0))}}))}
      
      if (var == TRUE) {return(by(longdata,longdata[,groupvar],function(x){
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- rep(0,2)
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return(list("par"=opt$b,"var"=opt$v));} 
          else {return(list("par"=c(0,0), "var"=diag(2)))};}
        if(((sum(!is.na(x[,scorevar1]))==0)) & (sum(!is.na(x[,scorevar2]))!=0)) {
          re <- 0
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return(list("par"=c(0,opt$b),"var"=matrix(c(1,0,0,opt$v), nrow=2, byrow=TRUE)));} 
          else {return(list("par"=c(0,0), "var"=diag(2)))};}
        if(((sum(!is.na(x[,scorevar1]))!=0)) & (sum(!is.na(x[,scorevar2]))==0)) {
          re <- 0
          opt<-marqLevAlg(b = re, fn = BivIndRePostDis2, rcpmeObj = rcpmeObj, data = x, scorevar1 = scorevar1, scorevar2 = scorevar2, timevar = timevar, model = model, gamma = gamma, link1 = link1, link2 = link2);
          if(opt$istop == 1) {return(list("par"=c(opt$b,0),"var"=matrix(c(opt$v,0,0,1), nrow =2, byrow=TRUE)));} 
          else {return(list("par"=c(0,0), "var"=diag(2)))};}
        else {return(list("par"=c(0,0), "var"=diag(2)))}}))}
    }
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
  model <- rcpmeObj$model
  if (is.null(model)){model <- "test"}
  gamma <- rcpmeObj$gamma
  if (is.null(gamma)){gamma <- 0.1}
  
  longdata <- get(paste(rcpmeObj$call$longdata))
  REesti <- REestimate(rcpmeObj, longdata)
  betas <- rcpmeObj$fixed[,1]
  
  if (rcpmeObj$covariate == "NULL"){
    betas <- as.numeric(rcpmeObj$fixed[,1])
    Ypred <- mapply(function(a,b){return(CPmodel(betas[1]+a[1], betas[2]+a[2], betas[3]+a[3], betas[4]+a[4], t = b[,timevar], gamma = gamma, model = model))}, REesti, by(longdata,longdata[,groupvar],function(x){return(x)}))
  }
  
  if (rcpmeObj$covariate != "NULL"){
    betas <- as.numeric(rcpmeObj$fixed[,1])
    covInd = rcpmeObj$covariate
    Ypred <- mapply(function(a,b){return(CPmodel(betas[1]+betas[2]*b[,covInd][1]+a[1], betas[3]+betas[4]*b[,covInd][1]+a[2], betas[5]+betas[6]*b[,covInd][1]+a[3], betas[7]+betas[8]*b[,covInd][1]+a[4], t = b[,timevar], gamma = gamma, model = model))}, REesti, by(longdata,longdata[,groupvar],function(x){return(x)}))
  }
  return(Ypred) 
}

BivIndPred <- function(rcpmeObj){
  
  formu <- rcpmeObj$call$formu
  scorevar1 = all.vars(formu)[1]
  scorevar2 = all.vars(formu)[2]
  timevar = all.vars(formu)[3]
  groupvar = all.vars(formu)[4]
  model <- rcpmeObj$model
  gamma <- rcpmeObj$gamma
  
  longdata <- get(paste(rcpmeObj$call$longdata))
  REesti <- REestimate(rcpmeObj, longdata)
  betas1 <- as.numeric(rcpmeObj$fixed[,1])
  betas2 <- as.numeric(rcpmeObj$fixed[,5])
  
  Ypred <- mapply(function(a,b){return(list(scorevar1=CPmodel(betas1[1]+a[1], betas1[2]+a[2], betas1[3]+a[3], betas1[4]+a[4], t = b[,timevar], gamma = gamma, model = model), scorevar2 = CPmodel(betas2[1]+a[5], betas2[2]+a[6], betas2[3]+a[7], betas2[4]+a[8], t = b[,timevar], gamma = gamma, model = model)))}, REesti, by(longdata,longdata[,groupvar],function(x){return(x)}))
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
CholToVar <- function(x, link1 = "linear", link2 = "linear"){
  
  rk1 = 5 - (link1 == "splines"); rk2 = rk1 + 12 - (link2 == "splines"); rk3 = rk2 + 17; rk4 = rk3 + (link1 == "splines")*5;
  Uepsa = NULL; if (link1 == "linear") {Uepsa <- x[rk1];}
  Uepsb = NULL; if (link2 == "linear") {Uepsb <- x[rk2];}
  
  U0a <- x[rk1+1]; U01a <- x[rk1+2]; U02a <- x[rk1+3]; U1a <- x[rk1+4]; U12a <- x[rk1+5]; U2a <- x[rk1+6]; Utaua <- x[rk1+7];
  U0b <- x[rk2+1]; U01b <- x[rk2+2]; U02b <- x[rk2+3]; U1b <- x[rk2+4]; U12b <- x[rk2+5]; U2b <- x[rk2+6]; Utaub <- x[rk2+7];
  Utau <- x[rk2+8]; U1 <- x[rk2+9]; U2 <- x[rk2+10]; U3 <- x[rk2+11]; U4 <- x[rk2+12]; U5 <- x[rk2+13]; U6 <- x[rk2+14]; U7 <- x[rk2+15]; U8 <- x[rk2+16]; U9 <- x[rk2+17];
 
  paramSpl1 <- NULL; paramSpl2 <- NULL;
  # if (link1 == "splines") {paramSpl1 <- x[(rk3+1):(rk3+5)]}  
  # if (link2 == "splines") {paramSpl2 <- x[(rk4+1):(rk4+5)]}
  
  return(c(x[1:4], Uepsa**2, U0a**2, U01a * U0a, U02a * U0a, U01a**2+U1a**2, U02a * U01a + U12a * U1a, U02a**2 + U12a**2 + U2a**2, Utaua**2,
           x[(1:4)+12-(link1=="splines")], Uepsb**2, U1**2 + U4**2 + U7**2 + U0b**2, U1*U2 + U4*U5 + U7*U8 + U0b*U01b, U1*U3 + U4*U6 + U7*U9 + U0b*U02b,
           U2**2+U5**2+U8**2+U01b**2+U1b**2, U3*U2 + U6*U5 + U9*U8 + U01b*U02b + U1b*U12b, U3**2 + U6**2 + U9**2 + U02b**2 + U12b**2 + U2b**2,
           Utau**2 + Utaub**2, Utaua*Utau, U0a*U1, U0a*U2, U0a*U3, U01a*U1 + U1a*U4, U01a*U2 + U1a*U5, U01a*U3 + U1a*U6, 
           U02a*U1 + U12a*U4 + U2a*U7, U02a*U2 + U12a*U5 + U2a * U8, U02a*U3 + U12a*U6 + U2a*U9, paramSpl1, paramSpl2))
}
CholToVarCovMatrix <- function(sim, link1 = "linear", link2 = "linear"){
  
  x <- sim$b
  
  rk1 = 5 - (link1 == "splines"); rk2 = rk1 + 12 - (link2 == "splines"); rk3 = rk2 + 17; rk4 = rk3 + (link1 == "splines")*5;
  lparam = 34 + (link1 == "splines")*4 + (link2 == "splines")*4
  
  Uepsa = 1; if (link1 == "linear") {Uepsa <- x[rk1];}
  Uepsb = 1; if (link2 == "linear") {Uepsb <- x[rk2];}
  
  U0a <- x[rk1+1]; U01a <- x[rk1+2]; U02a <- x[rk1+3]; U1a <- x[rk1+4]; U12a <- x[rk1+5]; U2a <- x[rk1+6]; Utaua <- x[rk1+7];
  U0b <- x[rk2+1]; U01b <- x[rk2+2]; U02b <- x[rk2+3]; U1b <- x[rk2+4]; U12b <- x[rk2+5]; U2b <- x[rk2+6]; Utaub <- x[rk2+7];
  Utau <- x[rk2+8]; U1 <- x[rk2+9]; U2 <- x[rk2+10]; U3 <- x[rk2+11]; U4 <- x[rk2+12]; U5 <- x[rk2+13]; U6 <- x[rk2+14]; U7 <- x[rk2+15]; U8 <- x[rk2+16]; U9 <- x[rk2+17];
  
  vi <- diag(lparam)
  vi[upper.tri(vi, diag=TRUE)] <- sim$v
  vi <- vi + t(vi) - diag(diag(vi))
  vi <- vi[1:rk3,1:rk3]
  
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
  
  if ((link1 == "splines") & (link2 == "linear")){Di <- Di[-5,-5]}
  if ((link1 == "linear") & (link2 == "splines")){Di <- Di[-17,-17]}
  if ((link1 == "splines") & (link2 == "splines")){Di <- Di[-c(5,17),-c(5,17)]}
  
  nvi <- Di %*% vi %*% t(Di)
  
  return(nvi)
}
