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