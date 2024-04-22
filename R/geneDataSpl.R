geneDataSpl  <- function(n = 100, params, pNA = 0, DO = 0, vis = c(0, 3, 6, 9, 12, 15, 18, 21), gamma = 0.1, pcas, pi = 0){
  
  nbvis = length(vis)
  donnees <- matrix(NA, nrow = n * nbvis, ncol = 4)
  colnames(donnees) <- {c("ID", "date", "score", "statut")}
  
  if (length(params) != 12){
    stop("Under the alternative with beta_2 random, params must be of length 12.")
  }
  else {
    Beta0 <- params[1]
    Beta1 <- params[2]
    Beta2 <- params[3]
    mutau <- params[4]
    sigma <- params[5]
    sigmab0 <- params[6]
    sigmab1 <- params[7]
    sigmab2 <- params[8]
    sigmabtau <- params[9]
    sigmab01 <- params[10] * sigmab0 * sigmab1
    sigmab02 <- params[11] * sigmab0 * sigmab2
    sigmab12 <- params[12] * sigmab2 * sigmab1
    B = matrix(c(sigmab0**2, sigmab01, sigmab02, 0, sigmab01, sigmab1**2, sigmab12, 0, sigmab02, sigmab12, sigmab2**2, 0, 0, 0, 0, sigmabtau**2), nrow = 4)
    bis <- matrix(NA, nrow = n, ncol = 4)
    for (i in seq(n)){
      # cas ou temoin
      statut = rbinom(1, 1, pcas)
      donnees[((nbvis*(i-1))+1):(nbvis*i),4] <- statut 
      
      if (statut == 1){
        delta = 1
      } else if (statut == 0 & pi == 0){
        delta = 0
      } else if (statut == 0 & pi > 0){
        delta = rbinom(1, 1, pi)
      }
      
      # generation des effets aleatoires pour l'individu i
      bis[i,] <- rmvnorm(1, mean = c(0,0,0,0), sigma = B)
      donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
      donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
      for (j in seq(nbvis)){
        eps = rnorm(1,0,sigma)
        donnees[(nbvis*(i-1))+j,3] = Beta0 + bis[i,1] + (Beta1 + bis[i,2])*vis[j] + (delta)*(Beta2 + bis[i,3])*(sum(c(4,4,4)*ispline2(vis[j] - (mutau + bis[i,4]),0,20,5,3,FALSE))) + eps
      }
    }
  }
  
  donnees <- as.data.frame(donnees)
  if ((pNA!=0) & (DO !=0)){
    stop("pNA and DO mustn\'t be both positive")
  }
  if ((pNA<0) | (pNA>1) | (DO<0) | (DO>1)){
    stop("pNA and DO takes values in [0;1]")
  }
  
  # censures
  if (pNA !=0){
    censures <- ifelse(c(t(cbind(rep(0,n),matrix(rbinom((nbvis-1)*n,1,pNA),nrow = n, ncol = nbvis-1))))==1,NA,0)
    donnees$score <- ifelse(is.na(censures),NA,donnees$score)
  }
  
  # drop-out
  if (DO !=0){
    geo <- rgeom(n, DO)+1
    dropout_times <- ifelse(geo>7,NA,geo)
    for (i in seq(n)){ 
      if (!is.na(dropout_times[i])){
        donnees$score[donnees$ID==i][-(1:(dropout_times[i]))] = NA
      }
    }
  }
  
  return(donnees)
}
