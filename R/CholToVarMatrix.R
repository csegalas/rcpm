CholToVarCovMatrix <- function(sim, link1 = "linear", link2 = "linear"){
  
  x <- sim$b
  
  rk1 = 5 - (link1 == "splines"); rk2 = rk1 + 12 - (link2 == "splines"); rk3 = rk2 + 17; rk4 = rk3 + (link1 == "splines")*5;
  lparam = 34 + (link1 == "splines")*4 + (link2 == "splines")*4
  
  Uepsa = 1; if (link1 == "linear") {Uepsa <- x[rk1];}
  Uepsb = 1; if (link2 == "linear") {Uepsb <- x[rk2];}
  
  U0a <- x[rk1+1]; U01a <- x[rk1+2]; U02a <- x[rk1+3]; U1a <- x[rk1+4]; U12a <- x[rk1+5]; U2a <- x[rk1+6]; Utaua <- x[rk1+7];
  U0b <- x[rk2+1]; U01b <- x[rk2+2]; U02b <- x[rk2+3]; U1b <- x[rk2+4]; U12b <- x[rk2+5]; U2b <- x[rk2+6]; Utaub <- x[rk2+7];
  Utau <- x[rk2+8]; U1 <- x[rk2+9]; U2 <- x[rk2+10]; U3 <- x[rk2+11]; U4 <- x[rk2+12]; U5 <- x[rk2+13]; U6 <- x[rk2+14]; U7 <- x[rk2+15]; U8 <- x[rk2+16]; U9 <- x[rk2+17];
  
  vi <- diag(lparam)
  vi[upper.tri(vi, diag=TRUE)] <- sim$v
  vi <- vi + t(vi) - diag(diag(vi))
  vi <- vi[1:rk3,1:rk3]
  
  Di <- matrix(c(1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 2*Uepsa, rep(0,34), 2*U0a, rep(0,33), U01a, U0a, rep(0,32), U02a, 0, 
                 U0a, rep(0,32), 2*U01a, 0, 2*U1a, rep(0,31), U02a, U01a, U12a, U1a, rep(0,31), 2*U02a, 0, 2*U12a, 2*U2a, rep(0,34),
                 2*Utaua, rep(0,34), 1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 1, rep(0,34), 2*Uepsb, rep(0,34), 2*U0b, rep(0,7),
                 2*U1, 0, 0, 2*U4, 0, 0, 2*U7, rep(0,19), U01b, U0b, rep(0,6), U2, U1, 0, U5, U4, 0, U8, U7,
                 rep(0,18), U02b, 0, U0b, rep(0,5), U3, 0, U1, U6, 0, U4, U9, 0, U7, rep(0,18), 2*U01b, 0, 2*U1b, rep(0,5), 2*U2,
                 0, 0, 2*U5, 0, 0, 2*U8, rep(0,19), U02b, U01b, U12b, U1b, rep(0,4), U3, U2, 0, U6, U5, 0, U9, U8, rep(0,19),
                 2*U02b, 0, 2*U12b, 2*U2b, rep(0,4), 2*U3, 0, 0, 2*U6, 0, 0, 2*U9, rep(0,23), 2*Utaub, 2*Utau, rep(0,20), Utau, rep(0,12),
                 Utaua, rep(0,14), U1, rep(0,19), U0a, rep(0,13), U2, rep(0,20), U0a, rep(0,12), U3, rep(0,21), U0a, rep(0,12),
                 U1,0,U4, rep(0,16), U01a,0,0,U1a, rep(0,11), U2, 0, U5, rep(0,17), U01a, 0, 0, U1a, rep(0,10), U3, 0, U6, rep(0,18),
                 U01a, 0, 0, U1a, rep(0,10), U1, 0, U4, U7, rep(0,14), U02a, 0, 0, U12a, 0, 0, U2a, rep(0,9), U2, 0, U5, U8, rep(0,15), 
                 U02a, 0, 0, U12a, 0, 0, U2a, rep(0,8), U3, 0, U6, U9, rep(0,16), U02a, 0, 0, U12a, 0, 0, U2a), byrow = TRUE, nrow = 34)
  
  if ((link1 == "splines") & (link2 == "linear")){Di <- Di[-5,-5]}
  if ((link1 == "linear") & (link2 == "splines")){Di <- Di[-17,-17]}
  if ((link1 == "splines") & (link2 == "splines")){Di <- Di[-c(5,17),-c(5,17)]}
  
  nvi <- Di %*% vi %*% t(Di)
  
  return(nvi)
}
