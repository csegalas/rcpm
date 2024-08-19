"L'amour pour principe, l'ordre pour base, et le progrès pour but; 
tel est, 
d'après ce long discours préliminaire, 
le caractère fondamental du régime définitif que le positivisme vient inaugurer
-Auguste Comte"
lvsblclass_penalized <- function(param, data1, data2, nq, grp, grp2, weights, nodes, scorevar, timevar, covariate, age_of_diagnosis, REadjust, model, link, objtrans, objtrans2, gamma, latent, classprob, two_means, intercept, membership, lambda, only_cases = FALSE){
  
  rk0 = 3 - (model == "isplines")
  rk1 = 4 + (link == "linear") - (model == "isplines")
  rk2 = rk1 + 7
  
  # on vire les parametres inutiles ds l'esti du modele lineaire
  if (link == "linear"){
    if(intercept){
      if(covariate != 'NULL'){
        paramlin <- param[c(1,2, rk1, rk1+1, rk1+2, rk1+4,rk2, rk2+1, rk2+2, rk2+3, rk4)]
      }
      else{
        paramlin <- param[c(1,2, rk1, rk1+1, rk1+2, rk1+4, rk4)]
      }
    }
    else {
      if(covariate != 'NULL'){
        paramlin <- param[c(1,2, rk1, rk1+1, rk1+2, rk1+4,rk2, rk2+1, rk2+2, rk2+3)]
        print(paramlin)
        print('params: ')
        print(param)
      }
      else{
        paramlin <- param[c(1,2, rk1, rk1+1, rk1+2, rk1+4, rk4)]
      }
    }
  }
  
  if (link == "isplines"){ 
    if(intercept){
      paramlin <- param[c(1,2,rk1,rk1+1,rk1+2,rk1+3,rk2, rk2+1, rk2+2, rk2+3, rk4)]

    }
    else {
      paramlin <- param[c(1,2,rk1,rk1+1,rk1+2,rk1+3,rk2, rk2+1, rk2+2, rk2+3)]
    }
  }
  
  
  loglik1 <- sum(lvsblNCgen(param, data2, nq, grp2, weights, nodes, scorevar, timevar, covariate, age_of_diagnosis, REadjust, model, link, objtrans2, gamma, loglik = TRUE, two_means = FALSE, intercept = FALSE))
  
  if(only_cases) {
    out <- loglik1 
  }
  else {
    if (latent == FALSE){
      loglik2 <- sum(lvsbllin(paramlin, data1, grp, scorevar, timevar, link, objtrans, loglik = TRUE, intercept = intercept, covariate))
      out = loglik1 + loglik2
    }
    if (latent == TRUE){
      if(is.null(membership)){
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
      }
      else {
        prob <- membership
      }
      lik3 <- lvsblNCgen(param, data1, nq, grp, weights, nodes, scorevar, timevar, covariate, age_of_diagnosis, REadjust, model, link, objtrans, gamma, loglik = FALSE, two_means = two_means, intercept = intercept)
      lik4 <- lvsbllin(paramlin, data1, grp, scorevar, timevar, link, objtrans, loglik = FALSE, intercept = intercept, covariate)
      if(is.null(membership)){
        out = loglik1 + sum(log((1-prob)*lik4 + prob*lik3)) - lambda * param[rk4+1] 
      }
      else{
        out = loglik1 + sum(log((1-prob)*lik4 + prob*lik3))
      }
    }
    
  }
  
  return(out)
}
