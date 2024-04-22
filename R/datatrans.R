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