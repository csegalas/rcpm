BivIndPred <- function(rcpmeObj){
  
  formu <- rcpmeObj$call$formu
  scorevar1 = all.vars(formu)[1]
  scorevar2 = all.vars(formu)[2]
  timevar = all.vars(formu)[3]
  groupvar = all.vars(formu)[4]
  model <- rcpmeObj$model
  gamma <- rcpmeObj$gamma
  
  longdata <- get(paste(rcpmeObj$call$longdata))
  REesti <- REestimate(rcpmeObj, longdata)
  betas1 <- as.numeric(rcpmeObj$fixed[,1])
  betas2 <- as.numeric(rcpmeObj$fixed[,5])
  
  Ypred <- mapply(function(a,b){return(list(scorevar1=CPmodel(betas1[1]+a[1], betas1[2]+a[2], betas1[3]+a[3], betas1[4]+a[4], t = b[,timevar], gamma = gamma, model = model), scorevar2 = CPmodel(betas2[1]+a[5], betas2[2]+a[6], betas2[3]+a[7], betas2[4]+a[8], t = b[,timevar], gamma = gamma, model = model)))}, REesti, by(longdata,longdata[,groupvar],function(x){return(x)}))
  return(Ypred) 
}