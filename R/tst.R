# functions for the test statistic
tst <- function(x, pert){
  return(sum(x * pert)**2/(sum(x**2)))
}