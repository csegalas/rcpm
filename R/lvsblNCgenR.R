lvsblNCgenR <- function(param,data,nq,grp,weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma, loglik) {
  
  temp <- lvsblNCgen(param,data,nq,grp,weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma, loglik, two_means = FALSE)
  
  if (loglik == TRUE){
    out = sum(temp) }
  else {
    out = sum(log(temp))
  }
  return(out)  
}
