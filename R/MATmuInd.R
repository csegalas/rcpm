MATmuInd <- function(data, Beta0, Beta1, sigma, newnodes, nbnodes, nd=3, formu){
  
  scorevar = data[!is.na(data[,all.vars(formu)[1]]),all.vars(formu)[1]]
  timevar = data[!is.na(data[,all.vars(formu)[1]]),all.vars(formu)[2]]
  lgt = length(scorevar)
  
  if (lgt==0){
    return(0)
  } else {
    idind <- data[,"ngroupvar"][1]
    MATmu <- (Beta0 + newnodes[[idind]][,1]) %*% t(rep(1,lgt)) + (Beta1+newnodes[[idind]][,2]) %*% t(timevar)
    diffY_MATmu <- (rep(1,nbnodes**nd) %*% t(scorevar)) - MATmu
    expNorm <- exp(-0.5*(diffY_MATmu/sigma)**2)
    if (lgt == 1){
      # prodIndDiff <- <- rep(1,nq**nd) modifie pour avoir type matrix
      prodIndDiff <- matrix(rep(1,nbnodes**nd),nrow = nbnodes**nd, ncol = 1)
    } else if (lgt == 2){
      prodIndDiff <- sapply(1:lgt,function(x){return(expNorm[,-x])})
    } else { 
      prodIndDiff <- sapply(1:lgt,function(x){return(apply(expNorm[,-x],1,prod))})
    }
  } 
  return(list(MATmu, diffY_MATmu, expNorm, prodIndDiff))
}