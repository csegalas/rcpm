CPmodel <- function(a, b, c, d, t, gamma, model="test", statut=NULL, etas = NULL){
  
  if (model == "test") return(a + b*t + c*sqrt((t-d)**2+gamma))
  if (model == "bw") return(a + b*(t-d) + c*sqrt((t-d)**2+gamma))
  if (model == "isplines") {
    ispline_t <- matrix(NA, nrow = 3, ncol = length(t))
    for (i in seq(t)) ispline_t[,i] <- ispline(t[i] - d, 0, 20, 5)
    return(a + b*t + statut*c*(etas^2 %*% ispline_t))
  }
}