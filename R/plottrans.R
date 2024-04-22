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