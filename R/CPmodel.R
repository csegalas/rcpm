CPmodel <- function(a,b,c,d,t,gamma,model="test"){
  if (model=="test") return(a + b*t + c*sqrt((t-d)**2+gamma))
  if (model=="bw") return(a + b*(t-d) + c*sqrt((t-d)**2+gamma))
}