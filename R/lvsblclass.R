lvsblclass <- function(param, data1, data2, nq, grp, grp2, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, objtrans2, gamma, latent, classprob){
  
  rk0 = 3 - (model == "isplines")
  rk1 = 4 + (link == "linear") - (model == "isplines")
  rk4 = rk1 + 7 + (covariate !="NULL")*(4*(REadjust=="no") + 5*(REadjust=="prop") + 11*(REadjust=="yes")) + 5*(link=="isplines") + 3*(model=="isplines")
  
  # on vire les parametres inutiles ds l'esti du modele lineaire
  if (link == "linear"){
    paramlin <- param[c(1,2, rk1, rk1+1, rk1+2, rk1+4)]
  }
  
  if (link == "isplines"){ # a finir celui ci...
    paramlin <- param[c(1,2,rk1,rk1+1,rk1+2,rk1+3)]
  }
  
  
  loglik1 <- sum(lvsblNCgen(param, data2, nq, grp2, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans2, gamma, loglik = TRUE))
  
  if (latent == FALSE){
    loglik2 <- sum(lvsbllin(paramlin, data1, grp, scorevar, timevar, link, objtrans, loglik = TRUE))
    out = loglik1 + loglik2
  }
  
  if (latent == TRUE){
    lgtclassprob = length(all.vars(classprob))
    if (lgtclassprob > 0){
      # extract class membership parameters
      classmbpar <- as.matrix(param[seq(rk4+2,rk4+lgtclassprob+1)], ncol = lgtclassprob)
      # extract class membership dep var
      classmbvar <- matrix(unlist(lapply(data1, function(x) return(x[1,all.vars(classprob)]))), ncol = lgtclassprob, byrow = TRUE)
      prob <- exp(param[rk4+1] + classmbvar %*% classmbpar)/(1+exp(param[rk4+1] + classmbvar %*% classmbpar))
    } else {
      prob <- exp(param[rk4+1])/(1+exp(param[rk4+1]))
    }
    lik3 <- lvsblNCgen(param, data1, nq, grp, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma, loglik = FALSE)
    lik4 <- lvsbllin(paramlin, data1, grp, scorevar, timevar, link, objtrans, loglik = FALSE)
    out = loglik1 + sum(log((1-prob)*lik4 + prob*lik3))
  }
  
  return(out)
}

