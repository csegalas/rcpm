VarToChol <- function(param){
  # this function computes a set of cholesky type parameters from a set of variance type parameters so that user
  # can input variance parameters directly
  
  B <- matrix(c(param[c(6,7,8)], 0, param[c(26,27,28)], 0, 
                param[c(7,9,10)], 0, param[c(29,30,31)], 0, 
                param[c(8,10,11)], 0, param[c(32,33,34)], 0, 
                rep(0,3), param[12], rep(0,3), param[25],
                param[c(26,29,32)], 0, param[c(18,19,20)], 0,
                param[c(27,30,33)], 0, param[c(19,21,22)], 0,
                param[c(28,31,34)], 0, param[c(20,22,23)], 0,
                rep(0,3), param[25], rep(0,3), param[24]), byrow = TRUE, nrow = 8)
  
  U <- chol(B)
  
  return(c(param[1:4], sqrt(param[5]), U[1,1], U[1,2], U[1,3], U[2,2], U[2,3], U[3,3], U[4,4],
           param[13:16], sqrt(param[17]), U[5,5], U[5,6], U[5,7], U[6,6], U[6,7], U[7,7], U[8,8],
           U[4,8], U[1,5], U[1,6], U[1,7], U[2,5], U[2,6], U[2,7], U[3,5], U[3,6], U[3,7]))
}
