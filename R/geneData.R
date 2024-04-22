geneData <- function(n = 100, hyp = "null", params, pNA = 0, DO = 0, vis = c(0, 3, 6, 9, 12, 15, 18, 21), beta2alea=FALSE, model = "test", gamma = 0.1){
  
  nbvis = length(vis)
  donnees <- matrix(NA, nrow = n * nbvis, ncol = ifelse(hyp=="altaleabiv", 4, 3))
  colnames(donnees) <- if (hyp=="altaleabiv"){c("ID", "date", "score1", "score2")} else {c("ID", "date", "score")}
  
  if (hyp=="null"){
    if (length(params) != 6){
      stop("Under the null, params must be of length 6.")
    }
    else {
      Beta0 <- params[1]
      Beta1 <- params[2]
      sigma <- params[3]
      sigmab0 <- params[4]
      sigmab1 <- params[5]
      sigmab01 <- params[6] * sigmab0 * sigmab1
      B = matrix(c(sigmab0**2, sigmab01, sigmab01, sigmab1**2), nrow = 2)
      bis <- matrix(NA, nrow = n, ncol = 2)
      for (i in seq(n)){
        # generation des effets aleatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = Beta0 + bis[i,1] + (Beta1 + bis[i,2]) * vis[j] + eps
        }
      }
    }
  }  else if (hyp=="alt") {
    if (length(params) != 9){
      stop("Under the alternative, params must be of length 9.")
    }
    else {
      Beta0 <- params[1]
      Beta1 <- params[2]
      Beta2 <- params[3]
      mutau <- params[4]
      sigma <- params[5]
      sigmab0 <- params[6]
      sigmab1 <- params[7]
      sigmabtau <- params[8]
      sigmab01 <- params[9] * sigmab0 * sigmab1
      B = matrix(c(sigmab0**2, sigmab01, 0, sigmab01, sigmab1**2, 0, 0, 0, sigmabtau**2), nrow = 3)
      bis <- matrix(NA, nrow = n, ncol = 3)
      for (i in seq(n)){
        # generation des effets aleatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0 + bis[i,1], Beta1 + bis[i,2], Beta2, mutau + bis[i,3], vis[j], gamma, model = model) + eps
        }
      }
    }
  } else if (hyp=="altalea"){
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
        # generation des effets aleatoires pour l'individu i
        bis[i,] <- rmvnorm(1, mean = c(0,0,0,0), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps = rnorm(1,0,sigma)
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0 + bis[i,1], Beta1 + bis[i,2], Beta2 + bis[i,3], mutau + bis[i,4], vis[j], gamma, model = model) + eps
        }
      }
    }
  } else if (hyp=="altaleabiv"){
    if (length(params) != 34){
      stop("Under the bivariate alternative, params must be of length 34.")
    }
    
    else {
      Beta0_1 <- params[1]; Beta1_1 <- params[2]; Beta2_1 <- params[3]; mutau_1 <- params[4]; Vareps_1 <- params[5]; B0_1 <- params[6]; B01_1 <- params[7]; B02_1 <- params[8];
      B1_1 <- params[9]; B12_1 <- params[10]; B2_1 <- params[11]; Btau_1 <- params[12];
      Beta0_2 <- params[13]; Beta1_2 <- params[14]; Beta2_2 <- params[15]; mutau_2 <- params[16]; Vareps_2 <- params[17]; B0_2 <- params[18]; B01_2 <- params[19]; B02_2 <- params[20];
      B1_2 <- params[21]; B12_2 <- params[22]; B2_2 <- params[23]; Btau_2 <- params[24];
      Btau_12 <- params[25]; B1_12 <- params[26]; B2_12 <- params[27]; B3_12 <- params[28]; B4_12 <- params[29]; B5_12 <- params[30];
      B6_12 <- params[31]; B7_12 <- params[32]; B8_12 <- params[33]; B9_12 <- params[34];
      
      B <- matrix(c(B0_1, B01_1, B02_1, 0, B1_12, B2_12, B3_12, 0,
                    B01_1, B1_1, B12_1, 0, B4_12, B5_12, B6_12, 0,
                    B02_1, B12_1, B2_1, 0, B7_12, B8_12, B9_12, 0,
                    0, 0, 0, Btau_1, 0, 0, 0, Btau_12,
                    B1_12, B4_12, B7_12, 0, B0_2, B01_2, B02_2, 0,
                    B2_12, B5_12, B8_12, 0, B01_2, B1_1, B12_2, 0,
                    B3_12, B6_12, B9_12, 0, B02_2, B12_2, B2_2, 0,
                    0, 0, 0, Btau_12, 0, 0, 0, Btau_2), byrow = TRUE, nrow = 8)
      
      bis <- matrix(NA, nrow = n, ncol = 8);
      for (i in seq(n)){
        bis[i,] <- rmvnorm(1, mean = rep(0,8), sigma = B)
        donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i
        donnees[((nbvis*(i-1))+1):(nbvis*i),2] <- vis
        for (j in seq(nbvis)){
          eps1 = rnorm(1,0,sqrt(Vareps_1)); eps2 <- rnorm(1,0,sqrt(Vareps_2));
          donnees[(nbvis*(i-1))+j,3] = CPmodel(Beta0_1 + bis[i,1], Beta1_1 + bis[i,2], Beta2_1 + bis[i,3], mutau_1 + bis[i,4], vis[j], gamma, model = model) + eps1
          donnees[(nbvis*(i-1))+j,4] = CPmodel(Beta0_2 + bis[i,5], Beta1_2 + bis[i,6], Beta2_2 + bis[i,7], mutau_2 + bis[i,8], vis[j], gamma, model = model) + eps2
        }
      }
      
    }
  }
  
  donnees <- as.data.frame(donnees)
  if (pNA!=0 & DO !=0){
    stop("pNA and DO mustn\'t be both positive")
  }
  if (pNA<0 | pNA>1 | DO<0 | DO>1){
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
