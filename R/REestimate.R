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
      if (var == FALSE) {return(lapply(split(longdata,longdata[,groupvar]),function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis, minimize = FALSE, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link); if(opt$istop == 1) {return("par"=opt$b);} else {return("par"=re);}} else {return("par"=re);}}))}
      if (var == TRUE) {return(lapply(split(longdata,longdata[,groupvar]),function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis, minimize = FALSE, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop == 1) {return(list("par"=opt$b,"var"=matrix(c(opt$v[c(1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10)]), nrow = 4, byrow=TRUE)));} else {return(list("par"=re, "var"=diag(4)))}} else {return(list("par"=re, "var"=diag(4)))}}))}
    }
    
    else {
      re <- 0
      if (var == FALSE) {return(lapply(split(longdata,longdata[,groupvar]),function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis2, minimize = FALSE, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop ==1) {return("par"=opt$b);} else {return("par"=re);}} else {return("par"=re);}}))}
      if (var == TRUE) {return(lapply(split(longdata,longdata[,groupvar]),function(x){if(sum(!is.na(x[,scorevar])!=0)) {opt<-marqLevAlg(b = re, fn = IndRePostDis2, minimize = FALSE, rcpmeObj = rcpmeObj, data = x, scorevar = scorevar, timevar = timevar, model = model, gamma = gamma, link = link);if(opt$istop == 1) {return(list("par"=opt$b,"var"=opt$v));} else {return(list("par"=0, "var"=1))};} else {return(list("par"=0, "var"=1))}}))}
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