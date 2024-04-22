gauherm <- function(n,m){
  grid <- gauher(n)
  idx <- as.matrix(expand.grid(rep(list(1:n),m)))
  points <- matrix(grid$x[idx],nrow(idx),m)
  weights <- apply(matrix(grid$w[idx],nrow(idx),m),1,prod)
  
  return(list(x = points, w = weights))
}