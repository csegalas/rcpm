

truncated_normal_quadrature <- function(N, a = NULL, b = NULL, mu, sigma){
  
  if(!is.null(a) & !is.null(b)){
    if( a > b) stop('A can not be bigger than B')
  }
  if(is.null(a) & is.null(b)){
    moment = moments_normal(2*N+1, mu, sigma)
  }
  else if(!is.null(a) & is.null(b)) {
    moment = moments_truncated_normal_a(2*N+1, mu, sigma, a)
  }
  else if(is.null(a) & !is.null(b)) {
    moment = moments_truncated_normal_b(2*N+1, mu, sigma, b)
  }
  else {
    moment = moments_truncated_normal_ab(2*N+1, mu, sigma, a, b)
  } 
  result <- moment_method(N, moment)
  return(result)
  
}