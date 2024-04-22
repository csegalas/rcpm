# function for individual random effects estimation
loglikmode <- function(param, estis, estiVarEA, data, scorevar, timevar){
  
  score = scorevar[!is.na(scorevar)]
  date = timevar[!is.na(scorevar)]
  
  b0 <- param[1]
  b1 <- param[2]
  
  Beta0 <- estis[1]
  Beta1 <- estis[2]
  sigma <- estis[3]
  
  MATmu <- Beta0 + b0 + (Beta1 + b1) * date
  
  return(-sum(dnorm(score, mean = MATmu, sigma, log = TRUE)) - mvtnorm::dmvnorm(c(b0,b1), c(0,0), estiVarEA, log = TRUE))
}