# function for individual random effects estimation
estisranef <- function(data, init, estis, estiVarEA, formu){
  
  opt <- optim(par = init, loglikmode, estis = estis, estiVarEA = estiVarEA, data = data, scorevar = data[,all.vars(formu)[1]], timevar = data[,all.vars(formu)[2]], hessian = TRUE, method=c("L-BFGS-B"))
  
  ranef <- opt$par
  ranefvar <- solve(opt$hessian)
  
  return(list(ranef,matrix(ranefvar, byrow=TRUE, nrow=2)))
}
