# Function to select parameters based on link and model type
getParamSlice <- function(param, rk1, link) {
  if (link == "linear") {
    return(param[c(1, 2, rk1, rk1+1, rk1+2, rk1+4)])
  } else {  # Assuming the only other option is "isplines"
    return(param[c(1, 2, rk1, rk1+1, rk1+2, rk1+3)])
  }
}

# Function to compute class membership probabilities
computeClassProb <- function(param, rk4, data, classprob) {
  vars <- all.vars(classprob)
  lgtclassprob <- length(vars)
  if (lgtclassprob > 0) {
    # Precompute indices once
    param_indices <- seq(rk4 + 2, rk4 + lgtclassprob + 1)
    classmbpar <- matrix(param[param_indices], ncol = lgtclassprob)
    
    # Optimizing data extraction:
    # Extract all needed variables at once into a matrix format if feasible
    classmbvar <- sapply(data[vars], `[`, 1)  # Assuming data is a list or data frame
    classmbvar <- matrix(classmbvar, ncol = lgtclassprob, byrow = TRUE)
    
    # Compute the logistic function in a vectorized manner
    linear_predictor <- param[rk4 + 1] + classmbvar %*% classmbpar
    return(exp(linear_predictor) / (1 + exp(linear_predictor)))
  } else {
    return(exp(param[rk4 + 1]) / (1 + exp(param[rk4 + 1])))
  }
}

lvsblclass <- function(param, data1, data2, nq, grp, grp2, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, objtrans2, gamma, latent, classprob){
  
  rk0 = 3 - (model == "isplines")
  rk1 = 4 + (link == "linear") - (model == "isplines")
  rk4 = rk1 + 7 + (covariate != "NULL")*(4*(REadjust == "no") + 5*(REadjust == "prop") + 11*(REadjust == "yes")) + 5*(link == "isplines") + 3*(model == "isplines")
  
  paramlin <- getParamSlice(param, rk1, link)
  
  loglik1 <- sum(lvsblNCgen(param, data2, nq, grp2, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans2, gamma, loglik = TRUE))
  
  if (!latent) {
    loglik2 <- sum(lvsbllin(paramlin, data1, grp, scorevar, timevar, link, objtrans, loglik = TRUE))
    out = loglik1 + loglik2
  } else {
    prob <- computeClassProb(param, rk4, data1, classprob)
    lik3 <- lvsblNCgen(param, data1, nq, grp, weights, nodes, scorevar, timevar, covariate, REadjust, model, link, objtrans, gamma, loglik = FALSE)
    lik4 <- lvsbllin(paramlin, data1, grp, scorevar, timevar, link, objtrans, loglik = FALSE)
    out = loglik1 + sum(log((1 - prob) * lik4 + prob * lik3))
  }
  
  return(out)
}

