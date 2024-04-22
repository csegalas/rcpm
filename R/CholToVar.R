CholToVar <- function(x, link1 = "linear", link2 = "linear"){
  
  rk1 = 5 - (link1 == "splines"); rk2 = rk1 + 12 - (link2 == "splines"); rk3 = rk2 + 17; rk4 = rk3 + (link1 == "splines")*5;
  Uepsa = NULL; if (link1 == "linear") {Uepsa <- x[rk1];}
  Uepsb = NULL; if (link2 == "linear") {Uepsb <- x[rk2];}
  
  U0a <- x[rk1+1]; U01a <- x[rk1+2]; U02a <- x[rk1+3]; U1a <- x[rk1+4]; U12a <- x[rk1+5]; U2a <- x[rk1+6]; Utaua <- x[rk1+7];
  U0b <- x[rk2+1]; U01b <- x[rk2+2]; U02b <- x[rk2+3]; U1b <- x[rk2+4]; U12b <- x[rk2+5]; U2b <- x[rk2+6]; Utaub <- x[rk2+7];
  Utau <- x[rk2+8]; U1 <- x[rk2+9]; U2 <- x[rk2+10]; U3 <- x[rk2+11]; U4 <- x[rk2+12]; U5 <- x[rk2+13]; U6 <- x[rk2+14]; U7 <- x[rk2+15]; U8 <- x[rk2+16]; U9 <- x[rk2+17];
  
  paramSpl1 <- NULL; paramSpl2 <- NULL;
  # if (link1 == "splines") {paramSpl1 <- x[(rk3+1):(rk3+5)]}  
  # if (link2 == "splines") {paramSpl2 <- x[(rk4+1):(rk4+5)]}
  
  return(c(x[1:4], Uepsa**2, U0a**2, U01a * U0a, U02a * U0a, U01a**2+U1a**2, U02a * U01a + U12a * U1a, U02a**2 + U12a**2 + U2a**2, Utaua**2,
           x[(1:4)+12-(link1=="splines")], Uepsb**2, U1**2 + U4**2 + U7**2 + U0b**2, U1*U2 + U4*U5 + U7*U8 + U0b*U01b, U1*U3 + U4*U6 + U7*U9 + U0b*U02b,
           U2**2+U5**2+U8**2+U01b**2+U1b**2, U3*U2 + U6*U5 + U9*U8 + U01b*U02b + U1b*U12b, U3**2 + U6**2 + U9**2 + U02b**2 + U12b**2 + U2b**2,
           Utau**2 + Utaub**2, Utaua*Utau, U0a*U1, U0a*U2, U0a*U3, U01a*U1 + U1a*U4, U01a*U2 + U1a*U5, U01a*U3 + U1a*U6, 
           U02a*U1 + U12a*U4 + U2a*U7, U02a*U2 + U12a*U5 + U2a * U8, U02a*U3 + U12a*U6 + U2a*U9, paramSpl1, paramSpl2))
}
