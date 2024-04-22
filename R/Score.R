Score <- function(tau,estis,estiVarEA,longdata,newnodes,newweights,nbnodes,nd,groupvar, MATmus, formu, pert){
  mutau <- tau[1]
  sigmatau <- tau[2]
  
  Beta0 <- estis[1]
  Beta1 <- estis[2]
  sigma <- estis[3]
  
  # U0 <- estis[4]
  # U01 <- estis[5]
  # U1 <- estis[6]
  # U <- matrix(c(U0,U01,0,0,U1,0,0,0,sigmatau), byrow = TRUE,nrow = 3)
  # 
  # B <- t(U) %*% U
  
  B <- rbind(cbind(estiVarEA,c(0,0)),c(0,0,sigmatau**2))
  
  return(-tst(by(longdata, groupvar, ScoreInd, Beta0, Beta1, sigma, mutau, sigmatau, B, MATmus, nbnodes, nd, newnodes, newweights, all.vars(formu)[1],all.vars(formu)[2],"ngroupvar"), pert)) # rcpp
  # return(-tst(by(data, grp, ScoreInd, sigma, mutau, sigmatau, B, MATmus, nq, nd, newnodes, newweights)))
}
