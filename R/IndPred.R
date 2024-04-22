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