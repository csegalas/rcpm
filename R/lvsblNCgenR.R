lvsblNCgenR <- function(param,data,nq,grp,weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma, loglik){
  
  temp <- lvsblNCgen(param,data,nq,grp,weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma, loglik)
  
  if (loglik == TRUE){
    out = sum(temp)
  }
  if (loglik == FALSE){
    out = sum(log(temp))
  }
  return(out)  
}
