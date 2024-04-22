MATmu <- function(estis, longdata, newnodes, nbnodes, nd=3, groupvar, formu){
  
  Beta0 <- estis[1]
  Beta1 <- estis[2]
  sigma <- estis[3]
  
  # return(by(longdata, groupvar, MATmuInd, Beta0, Beta1, sigma, newnodes, nbnodes, nd, formu))
  return(lapply(split(longdata, groupvar), function(x) return(MATmuInd(x, Beta0, Beta1, sigma, newnodes, nbnodes, nd, formu))))
}