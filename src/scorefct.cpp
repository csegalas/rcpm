#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dn(double x, double mu, double sd) {
  return pow(2*M_PI,-0.5) * pow(sd,-1) * exp(-0.5*pow((x-mu)/sd,2));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrmarma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd = false) {
  const double log2pi = std::log(2.0 * M_PI);
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma,"upper"))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}



// mspline <- function(x, tmin, tmax, tint){
//   t <- c(rep(tmin,3),tint,rep(tmax,3))
//   mspl <- rep(0,4)
//   mspl[1] <- (3*(t[4]-x)^2/((t[4]-t[1])*(t[4]-t[2])*(t[4]-t[3])))*(x<t[4])
//   mspl[2] <- (3*(x-t[2])*(t[4]-x)/((t[5]-t[2])*(t[4]-t[2])*(t[4]-t[3])) + 3*(t[5]-x)*(x-t[3])/((t[5]-t[2])*(t[5]-t[3])*(t[4]-t[3])))*(x<t[4]) + (3*(t[5]-x)^2/((t[5]-t[2])*(t[5]-t[3])*(t[5]-t[4])))*(x>=t[4])
//   mspl[3] <- (3*(x-t[3])^2/((t[6]-t[3])*(t[5]-t[3])*(t[4]-t[3])))*(x<t[4]) + (3*(x-t[3])*(t[5]-x)/((t[6]-t[3])*(t[5]-t[3])*(t[5]-t[4])) + 3*(t[6]-x)*(x-t[4])/((t[6]-t[3])*(t[6]-t[4])*(t[5]-t[4])))*(x>=t[4])
//   mspl[4] <- (3*(x-t[4])^2/((t[7]-t[4])*(t[6]-t[4])*(t[5]-t[4])))*(x>=t[4])
//   
//   return(mspl)
// }

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::colvec mspline(double x, double tmin, double tmax, double tint){
  arma::vec t = arma::zeros<arma::vec>(7);
  t(0) = tmin; t(1) = tmin; t(2) = tmin; t(3) = tint;
  t(4) = tmax; t(5) = tmax; t(6) = tmax;
  
  arma::colvec mspl = arma::zeros<arma::colvec>(4);
  
  mspl(0) = (3*pow((t(3)-x),2)/((t(3)-t(0))*(t(3)-t(1))*(t(3)-t(2))))*(x<t(3));
  mspl(1) = (3*(x-t(1))*(t(3)-x)/((t(4)-t(1))*(t(3)-t(1))*(t(3)-t(2))) + 3*(t(4)-x)*(x-t(2))/((t(4)-t(1))*(t(4)-t(2))*(t(3)-t(2))))*(x<t(3)) + (3*pow((t(4)-x),2)/((t(4)-t(1))*(t(4)-t(2))*(t(4)-t(3))))*(x>=t(3));
  mspl(2) = (3*pow((x-t(2)),2)/((t(5)-t(2))*(t(4)-t(2))*(t(3)-t(2))))*(x<t(3)) + (3*(x-t(2))*(t(4)-x)/((t(5)-t(2))*(t(4)-t(2))*(t(4)-t(3))) + 3*(t(5)-x)*(x-t(3))/((t(5)-t(2))*(t(5)-t(3))*(t(4)-t(3))))*(x>=t(3));
  mspl(3) = (3*pow((x-t(3)),2)/((t(6)-t(3))*(t(5)-t(3))*(t(4)-t(3))))*(x>=t(3));
  
  return mspl;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::colvec ispline(double x, double tmin, double tmax, double tint){
  arma::vec t = arma::zeros<arma::vec>(7);
  t(0) = tmin; t(1) = tmin; t(2) = tmin; t(3) = tint;
  t(4) = tmax; t(5) = tmax; t(6) = tmax;
  
  arma::colvec ispl = arma::zeros<arma::colvec>(4);
  arma::colvec mspl = mspline(x,tmin,tmax,tint);
  
  if (x >= tmin & x <= tmax){
    ispl(0) = (((x-t(0))*mspl(0) + (t(4)-x)*mspl(1) + (x-t(1))*mspl(1) + (t(5)-x)*mspl(2) + (x-t(2))*mspl(2) + (t(6)-x)*mspl(3))/3)*(x<t(3)) + 1*(x>=t(3));
    ispl(1) = (((x-t(1))*mspl(1) + (t(5)-x)*mspl(2) + (x-t(2))*mspl(2) + (t(6)-x)*mspl(3))/3) + ((x-t(3))*mspl(3)/3)*(x>=t(3));
    ispl(2) = ((x-t(2))*mspl(2) + (t(6)-x)*mspl(3))/3 + ((x-t(3))*mspl(3)/3)*(x>=t(3));
    ispl(3) = ((x-t(3))*mspl(3)/3)*(x>=t(3));
  }
  
  IntegerVector idx = IntegerVector::create(1,2,3);
  
  arma::colvec ispl2 = arma::zeros<arma::colvec>(3);
  ispl2(0) = ispl(1); ispl2(1) = ispl(2); ispl2(2) = ispl(3);
  
  return ispl2;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dmvnrmarma1d(arma::rowvec x, arma::rowvec mean, arma::mat sigma, bool logd = false) {
  const double log2pi = std::log(2.0 * M_PI);
  // int n = x.n_rows;
  int xdim = x.n_cols;
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma,"upper"))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  arma::vec z = rooti * arma::trans(x - mean) ;
  out = constants - 0.5 * arma::sum(z%z) + rootisum;
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double ScoreInd(DataFrame data, double Beta0, double Beta1, double sigma, double mutau, double sigmatau, arma::mat B, List MATmus, int nbnodes, int nd, List newnodes, List newweights, String scorevar, String timevar, String groupvar){
  
  NumericVector score = data[scorevar];
  LogicalVector indNoNA = !is_na(score);
  NumericVector scoreNoNA = score[indNoNA];
  NumericVector time = data[timevar];
  NumericVector timeNoNA = time[indNoNA];
  
  int lgt = scoreNoNA.length();
  if (lgt == 0){
    return 0;
  }
  else {
    NumericVector idsind = data[groupvar];
    int idind = idsind(0)-1;
    List MATmusInd = as<List>(MATmus[idind]);
    NumericMatrix MATmuInd = as<NumericMatrix>(MATmusInd[0]);
    NumericMatrix diffY_MATmuInd = as<NumericMatrix>(MATmusInd[1]);
    NumericMatrix expNormInd = as<NumericMatrix>(MATmusInd[2]);
    NumericMatrix prodIndDiffInd = as<NumericMatrix>(MATmusInd[3]);
    arma::mat noeuds = as<arma::mat>(newnodes[idind]);
    NumericVector poids = as<NumericVector>(newweights[idind]);
    arma::vec normRE = dmvnrmarma(noeuds, arma::zeros(1,3), B);
    double res1 = 0;
    double res2 = 0;
    for (int k = 0; k < pow(nbnodes,nd); k++){
      // On parcourt les noeuds
      double int1 = 1;
      double int2 = 0;
      for (int j = 0; j < lgt; j++){
        // On parcourt les mesures
        double mu = MATmuInd(k,j);
        double diff = diffY_MATmuInd(k,j);
        double expo = expNormInd(k,j);
        double prod = prodIndDiffInd(k,j);
        int1 *= dn(scoreNoNA(j),mu,sigma);
        int2 += expo * pow(sigma,-2) * diff * prod  * pow(pow(timeNoNA(j)-mutau-noeuds(k,2),2)+0.1,0.5);
      }
      res1 = res1 + int1 * poids(k) * normRE(k);
      res2 = res2 + int2 * pow(pow(2*M_PI,0.5)*sigma,-lgt) * poids(k) * normRE(k);
    }
    double res = 0;
    res = res2/res1;
    return res;
  }
}

// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// double lvsbl(NumericVector param, List data, int nq, NumericVector grp, NumericVector weights, NumericVector nodes){
//   double Beta0=param(0); double Beta1=param(1); double Beta2=param(2);
//   double tau=param(3); double Utau=std::exp(param(11)); double sigma=std::exp(param(4));
//   
//   arma::mat U = arma::zeros<arma::mat>(3,3);
//   U(0,0) = std::exp(param(5)); U(0,1) = param(6); U(0,2)=param(7); U(1,1)=std::exp(param(8)); U(1,2)=param(9); U(2,2)=std::exp(param(10));
//   arma::mat B = trans(U) * U;
//   
//   int N = max(grp);
//   double out = 0;
//   for (int i = 0; i<N; i++){
//     DataFrame datai = data[i];
//     
//     NumericVector score = datai["score"];
//     LogicalVector indNoNA = !is_na(score);
//     arma::colvec scoreNoNA = as<arma::colvec>(score[indNoNA]);
//     NumericVector time = datai["date"];
//     arma::colvec timeNoNA = as<arma::colvec>(time[indNoNA]);
//     
//     int lgt = scoreNoNA.n_elem;
//     double res = 0;
//     
//     if (lgt != 0){
//       arma::colvec Betas = arma::zeros<arma::colvec>(3);
//       Betas(0) = Beta0; Betas(1) = Beta1; Betas(2)=Beta2;
//       arma::mat Zk = arma::ones<arma::mat>(lgt,3);
//       Zk.col(1) = timeNoNA;
//       arma::colvec muk = arma::zeros<arma::colvec>(lgt);
//       arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
//       for (int k = 0; k < nq; ++k){
//         Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+Utau*nodes(k)),2)+0.1,0.5);
//         muk = Zk * Betas;
//         Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt); 
//         double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
//         res = res + (f * weights(k));
//       }
//       out = out + log(res);
//     }
//   }
//   
//   return out;
// }

// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// double lvsblNC(NumericVector param, List data, int nq, NumericVector grp, NumericVector weights, NumericVector nodes){
//   double Beta0=param(0); double Beta1=param(1); double Beta2=param(2);
//   double tau=param(3); double Utau=param(11); double sigma=param(4);
//   
//   arma::mat U = arma::zeros<arma::mat>(3,3);
//   U(0,0) = param(5); U(0,1) = param(6); U(0,2)=param(7); U(1,1)=param(8); U(1,2)=param(9); U(2,2)=param(10);
//   arma::mat B = trans(U) * U;
//   
//   int N = max(grp);
//   double out = 0;
//   for (int i = 0; i<N; i++){
//     DataFrame datai = data[i];
//     
//     NumericVector score = datai["score"];
//     LogicalVector indNoNA = !is_na(score);
//     arma::colvec scoreNoNA = as<arma::colvec>(score[indNoNA]);
//     NumericVector time = datai["date"];
//     arma::colvec timeNoNA = as<arma::colvec>(time[indNoNA]);
//     
//     int lgt = scoreNoNA.n_elem;
//     double res = 0;
//     
//     if (lgt != 0){
//       arma::colvec Betas = arma::zeros<arma::colvec>(3);
//       Betas(0) = Beta0; Betas(1) = Beta1; Betas(2)=Beta2;
//       arma::mat Zk = arma::ones<arma::mat>(lgt,3);
//       Zk.col(1) = timeNoNA;
//       arma::colvec muk = arma::zeros<arma::colvec>(lgt);
//       arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
//       for (int k = 0; k < nq; ++k){
//         Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+Utau*nodes(k)),2)+0.1,0.5);
//         muk = Zk * Betas;
//         Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt); 
//         double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
//         res = res + (f * weights(k));
//       }
//       out = out + log(res);
//     }
//   }
//   
//   return -out;
// }

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat transfY(arma::colvec Y, String link, arma::colvec param, List objtransi){
  int lgt = Y.n_rows;
  arma::mat out = arma::zeros<arma::mat>(lgt,2);
  if (link == "linear"){
    out.col(0) = Y; out.col(1) = arma::ones<arma::colvec>(lgt);
  }
  if (link == "splines"){
    arma::mat isOuti = as<arma::mat>(objtransi[0]);
    arma::mat isOutDerivi = as<arma::mat>(objtransi[1]);
    // double normparam = pow(param(0),2) + pow(param(1),2) + pow(param(2),2) + pow(param(3),2) + pow(param(4),2);
    out.col(0) = isOuti * param; out.col(1) = isOutDerivi * param;
    //out.col(0) = isOuti * pow(param,2)/normparam; out.col(1) = isOutDerivi * pow(param,2)/normparam;
  }
  return(out);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double lvsblNCgen(NumericVector param, List data, int nq, NumericVector grp, NumericVector weights, NumericVector nodes, String scorevar, String timevar, String covariate, String REadjust, String model, String link, List objtrans, double gamma){
  
  //// extraction des parametres selon les cas (link, cov/REadjust, model)
  // CP trajectory parameters
  double Beta0=param(0); double Beta1=param(1);
  int rk0 = 2; double Beta2 = param(rk0); 
  if (model == "isplines") {rk0 = 1; Beta2 = -1;}
  double tau=param(rk0+1);
  arma::colvec Betas = arma::zeros<arma::colvec>(3);
  Betas(0) = Beta0; Betas(1) = Beta1; Betas(2) = Beta2;
  
  // variance parameters
  double sigma = 1; int rk1 = rk0+1; 
  if (link == "linear") {rk1 = rk0+2; sigma = param(rk1);} 
  // double sigma = param(4); int rk1 = 4;
  arma::mat U = arma::zeros<arma::mat>(3,3);
  arma::mat B = arma::zeros<arma::mat>(3,3);
  U(0,0) = param(rk1+1); U(0,1) = param(rk1+2); U(0,2) = param(rk1+3); 
  U(1,1) = param(rk1+4); U(1,2) = param(rk1+5); U(2,2) = param(rk1+6);
  B = trans(U) * U; double Utau=param(rk1+7);
  
  // link parameters
  int rk2 = rk1 + 7;
  if ((covariate != "NULL") & (REadjust == "no")) {rk2 = rk1 + 11;}
  if ((covariate != "NULL") & (REadjust == "yes")) {rk2 = rk1 + 18;}
  if ((covariate != "NULL") & (REadjust == "prop")) {rk2 = rk1 + 12;}
  arma::colvec paramtrans = arma::zeros<arma::colvec>(5);
  
  int rk3 = rk2; // isplines index
  if (link == "splines"){
    rk3 = rk2 + 5;
    IntegerVector idx = IntegerVector::create(rk2+1, rk2+2, rk2+3, rk2+4, rk2+5);
    paramtrans = pow(as<arma::rowvec>(param[idx]),2);
  }
  
  // ispline model parameters
  arma::rowvec paramspl = arma::ones<arma::rowvec>(3);
  if (model == "isplines"){
    IntegerVector idx = IntegerVector::create(rk3+1, rk3+2, rk3+3);
    paramspl = pow(as<arma::rowvec>(param[idx]),2);
  }
 
  int N = max(grp);
  double out = 0;
  for (int i = 0; i<N; i++){
    
    DataFrame datai = as<DataFrame>(data[i]);
    NumericVector score = datai[scorevar];
    LogicalVector indNoNA = !is_na(score);
    arma::colvec scoreNoNA = as<arma::colvec>(score[indNoNA]);
    
    int lgt = scoreNoNA.n_elem;
    double res = 0;

    if (lgt > 0){

      NumericVector time = datai[timevar];
      arma::colvec timeNoNA = as<arma::colvec>(time[indNoNA]);
      List objtransi = as<List>(objtrans[i]);
      arma::mat tY = transfY(scoreNoNA, link, paramtrans, objtransi);

      if (covariate == "NULL"){
        // if no covariate
        arma::mat Zk = arma::ones<arma::mat>(lgt,3);
        if ((model == "test") | (model == "isplines")){Zk.col(1) = timeNoNA;}
        arma::colvec muk = arma::zeros<arma::colvec>(lgt);
        arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
        for (int k = 0; k < nq; ++k){
          if (model == "bw"){Zk.col(1) = timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+Utau*nodes(k));}
          if ((model == "bw") | (model == "test")) {Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+Utau*nodes(k)),2)+gamma,0.5);}
          if (model == "isplines") {
            arma::mat basespl = arma::zeros<arma::mat>(3,lgt);
            for (int l = 0; l<lgt; l++){
              basespl.col(l) = ispline(timeNoNA(l) - tau - Utau*nodes(k), 0, 20, 5);
            }
            Zk.col(2) = trans(paramspl * basespl);
            }
          muk = Zk * Betas;
          Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt);
          double f = dmvnrmarma1d(trans(tY.col(0)), trans(muk), Vk);
          res = res + (f * weights(k));
        }
        out = out + log(res) + sum(log(tY.col(1)));
      }

      if (covariate != "NULL"){
        // if covariate : need to take it into account
        NumericVector adjust = datai[covariate]; double adjusti = adjust(0);
        Betas(0) = Beta0 + param(rk1+8)*adjusti; Betas(1) = Beta1 + param(rk1+9)*adjusti; Betas(2)=Beta2+param(rk1+10)*adjusti;


        if (REadjust == "prop"){
          U(0,0) = (1 + param(rk1+12)*adjusti) * param(rk1+1); U(0,1) = (1 + param(rk1+12)*adjusti) * param(rk1+2); U(0,2)=(1 + param(rk1+12)*adjusti) * param(rk1+3);
          U(1,1)=(1 + param(rk1+12)*adjusti) * param(rk1+4); U(1,2)=(1 + param(rk1+12)*adjusti) * param(rk1+5); U(2,2)=(1 + param(rk1+12)*adjusti) * param(rk1+6);
          B = trans(U) * U;
          Utau = (1 + param(rk1+12)*adjusti) * param(rk1+7);
        }

        if (REadjust == "yes"){
          U(0,0) = param(rk1+1)+param(rk1+12)*adjusti; U(0,1) = param(rk1+2)+param(rk1+13)*adjusti; U(0,2)=param(rk1+3)+param(rk1+14)*adjusti;
          U(1,1)=param(rk1+4)+param(rk1+15)*adjusti; U(1,2)=param(rk1+5)+param(rk1+16)*adjusti; U(2,2)=param(rk1+6)+param(rk1+17)*adjusti;
          B = trans(U) * U;
          Utau = param(rk1+7) + param(rk1+18)*adjusti;
        }

        arma::mat Zk = arma::ones<arma::mat>(lgt,3);
        if (model == "test") {Zk.col(1) = timeNoNA;}
        arma::colvec muk = arma::zeros<arma::colvec>(lgt);
        arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
        for (int k = 0; k < nq; ++k){
          if (model == "bw") {Zk.col(1) = timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+param(rk1+11)*adjusti+Utau*nodes(k));}
          Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+param(rk1+11)*adjusti+Utau*nodes(k)),2)+gamma,0.5);
          muk = Zk * Betas;  Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt);
          double f = dmvnrmarma1d(trans(tY.col(0)), trans(muk), Vk);
          res = res + (f * weights(k));
        }
        out = out + log(res) + sum(log(tY.col(1)));
      }
    }
  }
  return out;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double lvsbllin(NumericVector param, List data, int nq, NumericVector grp, NumericVector weights, NumericVector nodes, String scorevar, String timevar, String link, List objtrans){
  
  //// extraction des parametres selon les cas (link) et (covariate & REadjust)
  // CP trajectory parameters
  double Beta0=param(0); double Beta1=param(1);
  arma::colvec Betas = arma::zeros<arma::colvec>(2);
  Betas(0) = Beta0; Betas(1) = Beta1;
  
  // variance parameters
  double sigma = 1; int rk1 = 1; 
  if (link == "linear") {rk1 = 2; sigma = param(rk1);} 
  arma::mat U = arma::zeros<arma::mat>(2,2);
  arma::mat B = arma::zeros<arma::mat>(2,2);
  U(0,0) = param(rk1+1); U(0,1) = param(rk1+2); U(1,1) = param(rk1+3);
  B = trans(U) * U;
  
  // link parameters
  int rk2 = rk1 + 3;
  arma::colvec paramtrans = arma::zeros<arma::colvec>(5);
  if (link == "splines"){
    IntegerVector idx = IntegerVector::create(rk2+1, rk2+2, rk2+3, rk2+4, rk2+5);
    paramtrans = pow(as<arma::rowvec>(param[idx]),2);
  }
  
  int N = max(grp);
  double out = 0;
  for (int i = 0; i<N; i++){
    
    DataFrame datai = as<DataFrame>(data[i]);
    NumericVector score = datai[scorevar];
    LogicalVector indNoNA = !is_na(score);
    arma::colvec scoreNoNA = as<arma::colvec>(score[indNoNA]);
    
    int lgt = scoreNoNA.n_elem;
    double res = 0;
    
    if (lgt > 0){
      
      NumericVector time = datai[timevar];
      arma::colvec timeNoNA = as<arma::colvec>(time[indNoNA]);
      List objtransi = as<List>(objtrans[i]);
      arma::mat tY = transfY(scoreNoNA, link, paramtrans, objtransi);
      
      // if no covariate
      arma::mat Zk = arma::ones<arma::mat>(lgt,2); Zk.col(1) = timeNoNA;
      arma::colvec muk = arma::zeros<arma::colvec>(lgt);
      arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
      for (int k = 0; k < nq; ++k){
        muk = Zk * Betas;
        Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt);
        double f = dmvnrmarma1d(trans(tY.col(0)), trans(muk), Vk);
        res = res + (f * weights(k));
      }
      out = out + log(res) + sum(log(tY.col(1)));
    }
  }
  return out;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double bilvsblNC(NumericVector param, List data, int nq, bool adapt, NumericVector grp, NumericVector weights, NumericMatrix nodes, List newnodes, List newweights, String scorevar1, String scorevar2, String timevar, String covariate, String REadjust, String model, String link1, String link2, List objtrans1, List objtrans2, double gamma){
  
  // index
  int rk1 = 4 - (link1 == "splines"); 
  int rk2 = rk1 + 12 - (link2 == "splines"); 
  int rk3 = rk2 + 17; 
  int rk4 = rk3 + (link1 == "splines")*5;
  
  // CP trajectory parameters
  double Beta0_1=param(0); double Beta1_1=param(1); double Beta2_1=param(2);
  double tau_1=param(3); double Utau_1=param(rk1+7);
  double Beta0_2=param(rk1+8); double Beta1_2=param(rk1+9); double Beta2_2=param(rk1+10);
  double tau_2=param(rk1+11); double Utau_2=param(rk2+7);
  double sigma_1 = 1; double sigma_2 = 1;
  if (link1 == "linear") {sigma_1 = param(4);} 
  if (link2 == "linear") {sigma_2 = param(rk1+12);} 
  
  // transformation parameters
  arma::colvec paramtrans1 = arma::zeros<arma::colvec>(5);
  arma::colvec paramtrans2 = arma::zeros<arma::colvec>(5);
  if (link1 == "splines"){
    IntegerVector idx1 = IntegerVector::create(rk3+1, rk3+2, rk3+3, rk3+4, rk3+5);
    paramtrans1 = pow(as<arma::colvec>(param[idx1]),2);
  }
  if (link2 == "splines"){
    IntegerVector idx2 = IntegerVector::create(rk4+1, rk4+2, rk4+3, rk4+4, rk4+5);
    paramtrans2 = pow(as<arma::colvec>(param[idx2]),2);
  }
  
  // non adaptive
  arma::mat nnodes(pow(nq,2),2);
  NumericVector nweights(pow(nq,2));
  if (adapt == FALSE){
    for (int i = 0; i<pow(nq,2); i++){
      nnodes(i,0) = nodes(i,0) * Utau_1 + nodes(i,1) * param(rk2+8);
      nnodes(i,1) = nodes(i,1) * Utau_2;
      nweights = weights;
    }
  }
  
  // adaptative
  arma::rowvec mure(2); arma::mat Bre(2,2);
  if (adapt == TRUE){
    arma::mat Ure(2,2);
    Ure(0,0) = Utau_1; Ure(1,1) = Utau_2; Ure(0,1) = param(rk2+8); Bre = trans(Ure) * Ure;
    mure = arma::zeros<arma::rowvec>(2);
  }
  
  int N = max(grp);
  double out = 0;
  for (int i = 0; i<N; i++){
    DataFrame datai = as<DataFrame>(data[i]);
    
    NumericVector score1 = datai[scorevar1];
    NumericVector score2 = datai[scorevar2];
    LogicalVector indNoNA1 = !is_na(score1);
    LogicalVector indNoNA2 = !is_na(score2);
    arma::colvec scoreNoNA1 = as<arma::colvec>(score1[indNoNA1]);
    arma::colvec scoreNoNA2 = as<arma::colvec>(score2[indNoNA2]);
    NumericVector time = datai[timevar];
    arma::colvec timeNoNA1 = as<arma::colvec>(time[indNoNA1]);
    arma::colvec timeNoNA2 = as<arma::colvec>(time[indNoNA2]);
    
    // transformation objects
    List objtrans1i = as<List>(objtrans1[i]);
    arma::mat tY1 = transfY(scoreNoNA1, link1, paramtrans1, objtrans1i);
    List objtrans2i = as<List>(objtrans2[i]);
    arma::mat tY2 = transfY(scoreNoNA2, link2, paramtrans2, objtrans2i);
    
    int lgt1 = scoreNoNA1.n_elem; int lgt2 = scoreNoNA2.n_elem; int lgt = lgt1 + lgt2;
    arma::colvec scoreNoNA = join_cols(tY1.col(0), tY2.col(0));
    double res = 0;
    
    if (lgt != 0){
    
      arma::mat U = arma::zeros<arma::mat>(6,6);
      arma::mat B = arma::zeros<arma::mat>(6,6);
      
      if (covariate == "NULL"){
        
        U(0,0) = param(rk1+1); U(0,1) = param(rk1+2); U(0,2) = param(rk1+3); U(1,1) = param(rk1+4); U(1,2) = param(rk1+5); U(2,2) = param(rk1+6);
        U(3,3) = param(rk2+1); U(3,4) = param(rk2+2); U(3,5) = param(rk2+3); U(4,4) = param(rk2+4); U(4,5) = param(rk2+5); U(5,5) = param(rk2+6);
        U(0,3) = param(rk2+9); U(0,4) = param(rk2+10); U(0,5) = param(rk2+11); U(1,3)=param(rk2+12); U(1,4)=param(rk2+13); U(1,5) = param(rk2+14);
        U(2,3) = param(rk2+15); U(2,4) = param(rk2+16); U(2,5) = param(rk2+17);
        
        B = trans(U) * U;
        
        arma::colvec Betas = arma::zeros<arma::colvec>(6);
        Betas(0) = Beta0_1; Betas(1) = Beta1_1; Betas(2)=Beta2_1;
        Betas(3) = Beta0_2; Betas(4) = Beta1_2; Betas(5)=Beta2_2;
        
        arma::mat Zk = arma::zeros<arma::mat>(lgt,6);
        arma::mat Zk01 = arma::zeros<arma::mat>(lgt1,3);
        arma::mat Zk02 = arma::zeros<arma::mat>(lgt2,3);
        arma::mat Zk1 = arma::ones<arma::mat>(lgt1,3);
        arma::mat Zk2 = arma::ones<arma::mat>(lgt2,3);
        if (model == "test"){Zk1.col(1) = timeNoNA1; Zk2.col(1) = timeNoNA2;}
        arma::colvec muk = arma::zeros<arma::colvec>(lgt);
        arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
        arma::mat Sigma1 = pow(sigma_1,2) * arma::eye<arma::mat>(lgt1,lgt1);
        arma::mat Sigma2 = pow(sigma_2,2) * arma::eye<arma::mat>(lgt2,lgt2);
        arma::mat Sigma0 = arma::zeros<arma::mat>(lgt1,lgt2);
        arma::mat Sigma = join_cols(join_rows(Sigma1,Sigma0), join_rows(trans(Sigma0),Sigma2));
        arma::vec qre = arma::ones<arma::vec>(pow(nq,2));
        
        if (adapt == TRUE){
          nnodes = as<arma::mat>(newnodes[i]);
          nweights = as<NumericVector>(newweights[i]);
          qre = dmvnrmarma(nnodes, mure, Bre);
        }
        
        for (int k = 0; k < pow(nq,2); ++k){
          if (model == "bw"){
            Zk1.col(1) = timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(tau_1+nnodes(k,0)); 
            Zk2.col(1) = timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(tau_2+nnodes(k,1));}
            Zk1.col(2) = pow(pow(timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(tau_1+nnodes(k,0)),2)+0.1,0.5);
            Zk2.col(2) = pow(pow(timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(tau_2+nnodes(k,1)),2)+0.1,0.5);
            Zk = join_cols(join_rows(Zk1,Zk01),join_rows(Zk02,Zk2));
            muk = Zk * Betas; Vk = (Zk * B) * trans(Zk) + Sigma;
            double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
            res = res + (f * qre(k) * nweights(k));
        }
        out = out + log(res) + sum(log(tY1.col(1))) + sum(log(tY2.col(1)));
      }
      
      // if (covariate != "NULL"){
      //   NumericVector adjust = datai[covariate];
      //   double adjusti = adjust(0);
      //   
      //   if (REadjust == "no"){
      //     U(0,0) = param(5); U(0,1) = param(6); U(0,2)=param(7); U(1,1)=param(8); U(1,2)=param(9); U(2,2)=param(10); B = trans(U) * U;
      //   }
      //   
      //   if (REadjust == "prop"){
      //     U(0,0) = (1 + param(16)*adjusti) * param(5); U(0,1) = (1 + param(16)*adjusti) * param(6); U(0,2)=(1 + param(16)*adjusti) * param(7);
      //     U(1,1)=(1 + param(16)*adjusti) * param(8); U(1,2)=(1 + param(16)*adjusti) * param(9); U(2,2)=(1 + param(16)*adjusti) * param(10);
      //     B = trans(U) * U;
      //     Utau = (1 + param(16)*adjusti) * Utau;
      //   }
      //   
      //   if (REadjust == "yes"){
      //     U(0,0) = param(5)+param(16)*adjusti; U(0,1) = param(6)+param(17)*adjusti; U(0,2)=param(7)+param(18)*adjusti;
      //     U(1,1)=param(8)+param(19)*adjusti; U(1,2)=param(9)+param(20)*adjusti; U(2,2)=param(10)+param(21)*adjusti;
      //     B = trans(U) * U;
      //     Utau = param(11) + param(22)*adjusti;
      //   }
      //   
      //   arma::colvec Betas = arma::zeros<arma::colvec>(3);
      //   Betas(0) = Beta0 + param(12)*adjusti; Betas(1) = Beta1 + param(13)*adjusti; Betas(2)=Beta2+param(14)*adjusti;
      //   
      //   arma::mat Zk = arma::ones<arma::mat>(lgt,3);
      //   Zk.col(1) = timeNoNA;
      //   arma::colvec muk = arma::zeros<arma::colvec>(lgt);
      //   arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
      //   for (int k = 0; k < nq; ++k){
      //     Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+param(15)*adjusti+Utau*nodes(k)),2)+0.1,0.5);
      //     muk = Zk * Betas;
      //     Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt); 
      //     double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
      //     res = res + (f * weights(k));
      //   }
      //   out = out + log(res);
      // }
    }
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double IndRePostDis(arma::rowvec re, DataFrame data, List rcpmeObj, String scorevar, String timevar, String model, double gamma, String link){
  
  NumericVector scoreAll = data[scorevar];
  NumericVector timeAll = data[timevar];
  LogicalVector indNoNA = !is_na(scoreAll);
  arma::colvec score = as<arma::colvec>(scoreAll[indNoNA]);
  arma::colvec time = as<arma::colvec>(timeAll[indNoNA]);
  int lgt = score.n_elem;
  
  if (lgt == 0){
    return 0;
  }
  
  String covariate =  as<String>(rcpmeObj[7]);
  NumericVector mus(lgt);
  arma::mat estiVarEA(4,4);
  
  if (covariate == "NULL"){
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    IntegerVector idx = IntegerVector::create(0, 1, 2, 3);
    NumericVector betas = par[idx];
    for (int i = 0; i<lgt; ++i){
      if (model == "test"){
        mus(i) = par(0) + re(0) + (par(1) + re(1)) * time(i) + (par(2)+re(2)) * pow(pow(time(i) - par(3) - re(3),2)+gamma,0.5);
      }
      if (model == "bw"){
        mus(i) = par(0) + re(0) + (par(1) + re(1)) * (time(i) - par(3) - re(3)) + (par(2)+re(2)) * pow(pow(time(i) - par(3) - re(3),2)+gamma,0.5);
      }
    }
  }

  NumericVector adjust;
  double adjusti;
  if (covariate != "NULL"){
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    int rk1 = 3 + (link == "linear"); 
    IntegerVector idx = IntegerVector::create(0, rk1+8, 1, rk1+9, 2, rk1+10, 3, rk1+11);
    NumericVector betas = par[idx];
    adjust = data[covariate];
    adjusti = adjust(0);
    for (int i =0; i<lgt; i++){
      if (model == "test"){
        mus(i) = betas(0) + betas(1)*adjusti + re(0) + (betas(2) + betas(3)*adjusti + re(1)) * time(i) + (betas(4) + betas(5)*adjusti + re(2)) * pow(pow(time(i) - betas(6) - betas(7)*adjusti - re(3),2)+gamma,0.5);
      }
      if (model == "bw"){
        mus(i) = betas(0) + betas(1)*adjusti + re(0) + (betas(2) + betas(3)*adjusti + re(1)) * (time(i) - betas(6) - betas(7)*adjusti - re(3)) + (betas(4) + betas(5)*adjusti + re(2)) * pow(pow(time(i) - betas(6) - betas(7)*adjusti - re(3),2)+gamma,0.5);
      }
    }
  }

  double sdres = as<double>(rcpmeObj[4]); // return 1 if splines
  String REadjust = as<String>(rcpmeObj[8]);

  if ((REadjust == "no") | (covariate == "NULL")){
    estiVarEA = as<arma::mat>(rcpmeObj[5]);
  }

  if (REadjust == "yes"){ // CAUTION only fit for 0/1 covariate
    List estiVarEAs = as<List>(rcpmeObj[5]);
    estiVarEA = as<arma::mat>(estiVarEAs[adjusti]);
  }

  arma::rowvec mure = arma::zeros<arma::rowvec>(4);
  bool logtrue = true;
  double out = dmvnrmarma1d(re,mure,estiVarEA,logtrue);
  for (int i = 0; i<lgt; i++){
    out = out - 0.5*log(2*M_PI) - log(sdres) - (0.5/pow(sdres,2))*pow(score(i)-mus(i),2);
  }
  
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double IndRePostDis2(double re, DataFrame data, List rcpmeObj, String scorevar, String timevar, String model, double gamma, String link){
  
  // extracting data
  NumericVector scoreAll = data[scorevar];
  NumericVector timeAll = data[timevar];
  LogicalVector indNoNA = !is_na(scoreAll);
  arma::colvec score = as<arma::colvec>(scoreAll[indNoNA]);
  arma::colvec time = as<arma::colvec>(timeAll[indNoNA]);
  int lgt = score.n_elem;
  
  // extracting estimated parameters
  arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
  arma::mat varB = estiVarEA.submat(0,0,2,2);
  NumericVector par = as<NumericVector>(rcpmeObj[6]);
  IntegerVector idx = IntegerVector::create(0, 1, 2);
  arma::colvec betas = as<arma::colvec>(par[idx]);
  double sdtau = sqrt(estiVarEA(3,3)); double sdres = as<double>(rcpmeObj[4]);
  
  // distribution of Y_i | tau_i and tau_i
  arma::colvec mus(lgt);
  arma::mat Zk = arma::ones<arma::mat>(lgt,3);
  if (model == "test") {Zk.col(1) = time;}
  if (model == "bw") {Zk.col(1) = time - arma::ones<arma::colvec>(lgt)*(par[3]+re);}
  Zk.col(2) = pow(pow(time - arma::ones<arma::colvec>(lgt)*(par[3]+re),2)+gamma,0.5);
  mus = Zk * betas; arma::mat B = Zk * varB * trans(Zk) + pow(sdres,2)*arma::eye<arma::mat>(lgt,lgt);
  bool logtrue = true; double out = dmvnrmarma1d(trans(score),trans(mus),B,logtrue);
  // out = out - 0.5*(log(2*M_PI) + pow(re,2));
  out = out - 0.5*(log(2*M_PI)+pow(re/sdtau,2)) - log(sdtau);
  
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double BivIndRePostDis(arma::rowvec re, DataFrame data, List rcpmeObj, String scorevar1, String scorevar2, String timevar, String model, double gamma, String link1, String link2){
  
  double out = 0;
  
  // index
  int rk1 = 4 - (link1 == "splines");
  int rk2 = rk1 + 12 - (link2 == "splines");
  int rk3 = rk2 + 17;
  int rk4 = rk3 + (link1 == "splines")*5;
  
  // extracting data
  NumericVector score1 = data[scorevar1]; NumericVector score2 = data[scorevar2];
  LogicalVector indNoNA1 = !is_na(score1); LogicalVector indNoNA2 = !is_na(score2);
  arma::colvec scoreNoNA1 = as<arma::colvec>(score1[indNoNA1]); arma::colvec scoreNoNA2 = as<arma::colvec>(score2[indNoNA2]);
  NumericVector time = data[timevar];
  arma::colvec timeNoNA1 = as<arma::colvec>(time[indNoNA1]); arma::colvec timeNoNA2 = as<arma::colvec>(time[indNoNA2]);
  
  int lgt1 = scoreNoNA1.n_elem; int lgt2 = scoreNoNA2.n_elem; int lgt = lgt1 + lgt2;
  arma::colvec scoreNoNA = join_cols(scoreNoNA1, scoreNoNA2);
  if (lgt2 == 0 & lgt1 != 0){
    // extracting estimated parameters
    arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
    arma::mat B1 = estiVarEA.submat(0,0,3,3);
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    arma::colvec betas = arma::zeros<arma::colvec>(3);
    
    betas(0) = par(0) + re(0); betas(1) = par(1) + re(1);
    betas(2) = par(2) + re(2);
    arma::colvec sigma2 = as<arma::colvec>(rcpmeObj[4]);
    arma::mat Sigma1 = pow(sigma2(0),2) * arma::eye<arma::mat>(lgt1,lgt1);
    
    
    // distribution of Y_i | tau_i and tau_i
    arma::mat Z1 = arma::ones<arma::mat>(lgt1,3);
    if (model == "test"){Z1.col(1) = timeNoNA1;}
    if (model == "bw"){Z1.col(1) = timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(3));}
    Z1.col(2) = pow(pow(timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(3)),2)+gamma,0.5);
    arma::colvec mu = arma::zeros<arma::colvec>(lgt1); mu = Z1 * betas;

    arma::rowvec mure = arma::zeros<arma::rowvec>(4);
    bool logtrue = true; double out = dmvnrmarma1d(trans(scoreNoNA1), trans(mu), Sigma1, logtrue);
    out = out + dmvnrmarma1d(re, mure, B1, logtrue);
    
    return(out);
  }
  if (lgt1 == 0 & lgt2 != 0){
    // extracting estimated parameters
    arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
    arma::mat B2 = estiVarEA.submat(4,4,7,7);
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    arma::colvec betas = arma::zeros<arma::colvec>(3);
    
    betas(0) = par(rk1+8) + re(0); betas(1) = par(rk1+9) + re(1); 
    betas(2) = par(rk1+10) + re(2);
    arma::colvec sigma2 = as<arma::colvec>(rcpmeObj[4]);
    arma::mat Sigma2 = pow(sigma2(1),2) * arma::eye<arma::mat>(lgt2,lgt2);
    
    // distribution of Y_i | tau_i and tau_i
    arma::mat Z2 = arma::ones<arma::mat>(lgt2,3);
    if (model == "test"){ Z2.col(1) = timeNoNA2;}
    if (model == "bw"){Z2.col(1) = timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(3));}
    Z2.col(2) = pow(pow(timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(3)),2)+gamma,0.5);
    arma::colvec mu = arma::zeros<arma::colvec>(lgt2); mu = Z2 * betas;

    arma::rowvec mure = arma::zeros<arma::rowvec>(4);
    bool logtrue = true; double out = dmvnrmarma1d(trans(scoreNoNA2), trans(mu), Sigma2, logtrue);
    out = out + dmvnrmarma1d(re, mure, B2, logtrue);
    
    return(out);
  }
  if ((lgt1 != 0) & (lgt2 != 0)){
    // extracting estimated parameters
    arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    arma::colvec betas = arma::zeros<arma::colvec>(6);
    
    betas(0) = par(0) + re(0); betas(1) = par(1) + re(1);
    betas(2) = par(2) + re(2); betas(3) = par(rk1+8) + re(4);
    betas(4) = par(rk1+9) + re(5); betas(5) = par(rk1+10) + re(6);
    arma::colvec sigma2 = as<arma::colvec>(rcpmeObj[4]);
    arma::mat Sigma1 = pow(sigma2(0),2) * arma::eye<arma::mat>(lgt1,lgt1);
    arma::mat Sigma2 = pow(sigma2(1),2) * arma::eye<arma::mat>(lgt2,lgt2);
    arma::mat Sigma0 = arma::zeros<arma::mat>(lgt1,lgt2);
    arma::mat Sigma = join_cols(join_rows(Sigma1,Sigma0), join_rows(trans(Sigma0),Sigma2));
    
    
    // distribution of Y_i | tau_i and tau_i
    arma::mat Z = arma::zeros<arma::mat>(lgt,6);
    arma::mat Z01 = arma::zeros<arma::mat>(lgt1,3);
    arma::mat Z02 = arma::zeros<arma::mat>(lgt2,3);
    arma::mat Z1 = arma::ones<arma::mat>(lgt1,3);
    arma::mat Z2 = arma::ones<arma::mat>(lgt2,3);
    if (model == "test"){Z1.col(1) = timeNoNA1; Z2.col(1) = timeNoNA2;}
    if (model == "bw"){Z1.col(1) = timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(3)); Z2.col(1) = timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(7));}
    Z1.col(2) = pow(pow(timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(3)),2)+gamma,0.5);
    Z2.col(2) = pow(pow(timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(7)),2)+gamma,0.5);
    Z = join_cols(join_rows(Z1,Z01),join_rows(Z02,Z2));
    arma::colvec mu = arma::zeros<arma::colvec>(lgt); mu = Z * betas;

    arma::rowvec mure = arma::zeros<arma::rowvec>(8);
    bool logtrue = true; double out = dmvnrmarma1d(trans(scoreNoNA), trans(mu), Sigma, logtrue);
    out = out + dmvnrmarma1d(re, mure, estiVarEA, logtrue);
    
    return(out);
  }
  
  return(out);
}

// [[Rcpp::depends("RcppArmadillo")]] A ADAPTER
// [[Rcpp::export]]
double BivIndRePostDis2(arma::colvec re, DataFrame data, List rcpmeObj, String scorevar1, String scorevar2, String timevar, String model, double gamma, String link1, String link2){
  
  double out = 0;
  
  // index
  int rk1 = 4 - (link1 == "splines");
  int rk2 = rk1 + 12 - (link2 == "splines");
  int rk3 = rk2 + 17;
  int rk4 = rk3 + (link1 == "splines")*5;
  
  // extracting data
  NumericVector score1 = data[scorevar1]; NumericVector score2 = data[scorevar2];
  LogicalVector indNoNA1 = !is_na(score1); LogicalVector indNoNA2 = !is_na(score2);
  arma::colvec scoreNoNA1 = as<arma::colvec>(score1[indNoNA1]); arma::colvec scoreNoNA2 = as<arma::colvec>(score2[indNoNA2]);
  NumericVector time = data[timevar];
  arma::colvec timeNoNA1 = as<arma::colvec>(time[indNoNA1]); arma::colvec timeNoNA2 = as<arma::colvec>(time[indNoNA2]);
  
  int lgt1 = scoreNoNA1.n_elem; int lgt2 = scoreNoNA2.n_elem; int lgt = lgt1 + lgt2;
  arma::colvec scoreNoNA = join_cols(scoreNoNA1, scoreNoNA2);

  if ((lgt1 == 0) & (lgt2 != 0)){
    // extracting estimated parameters
    arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
    arma::mat B2 = estiVarEA.submat(4,4,6,6);
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    IntegerVector idx = IntegerVector::create(rk1+8, rk1+9, rk1+10);
    arma::colvec betas = as<arma::colvec>(par[idx]);
    arma::colvec sigma2 = as<arma::colvec>(rcpmeObj[4]);
    arma::mat Sigma2 = pow(sigma2(1),2) * arma::eye<arma::mat>(lgt2,lgt2);
    
    double mure = 0;
    double Bre = sqrt(estiVarEA(7,7));
    
    // distribution of Y_i | tau_i and tau_i
    arma::mat Z2 = arma::ones<arma::mat>(lgt2,3);
    if (model == "test"){Z2.col(1) = timeNoNA2;}
    if (model == "bw"){Z2.col(1) = timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(0));}
    Z2.col(2) = pow(pow(timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(0)),2)+gamma,0.5);
    arma::colvec mu = arma::zeros<arma::colvec>(lgt2); mu = Z2 * betas;
    arma::mat V = arma::zeros<arma::mat>(lgt2,lgt2); V = (Z2 * B2) * trans(Z2) + Sigma2;
    bool logtrue = true; double out = dmvnrmarma1d(trans(scoreNoNA2), trans(mu), V, logtrue);
    out = out + log(dn(re(0), mure, Bre));
    
    return out;
  }
  if ((lgt1 != 0) & (lgt2 == 0)){
    // extracting estimated parameters
    arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
    arma::mat B1 = estiVarEA.submat(0,0,2,2);
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    IntegerVector idx = IntegerVector::create(0, 1, 2);
    arma::colvec betas = as<arma::colvec>(par[idx]);
    arma::colvec sigma2 = as<arma::colvec>(rcpmeObj[4]);
    arma::mat Sigma1 = pow(sigma2(0),2) * arma::eye<arma::mat>(lgt1,lgt1);
    
    double mure = 0;
    double Bre = sqrt(estiVarEA(3,3)); 
    
    // distribution of Y_i | tau_i and tau_i
    arma::mat Z1 = arma::ones<arma::mat>(lgt1,3);
    if (model == "test"){Z1.col(1) = timeNoNA1;}
    if (model == "bw"){Z1.col(1) = timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(0)); }
    Z1.col(2) = pow(pow(timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(0)),2)+gamma,0.5);
    
    arma::colvec mu = arma::zeros<arma::colvec>(lgt1); mu = Z1 * betas;
    arma::mat V = arma::zeros<arma::mat>(lgt1,lgt1); V = (Z1 * B1) * trans(Z1) + Sigma1;
    
    bool logtrue = true; double out = dmvnrmarma1d(trans(scoreNoNA1), trans(mu), V, logtrue);
    out = out + log(dn(re(0), mure, Bre));
    
    return out;
  }
  if ((lgt1 != 0) & (lgt2 != 0)){
    // extracting estimated parameters
    arma::mat estiVarEA = as<arma::mat>(rcpmeObj[5]);
    arma::mat B12 = estiVarEA.submat(0,4,2,6); arma::mat B21 = estiVarEA.submat(4,0,6,2);
    arma::mat B1 = estiVarEA.submat(0,0,2,2); arma::mat B2 = estiVarEA.submat(4,4,6,6);
    arma::mat B = arma::zeros<arma::mat>(6,6); B = join_cols(join_rows(B1,B21),join_rows(B12,B2));
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    IntegerVector idx = IntegerVector::create(0, 1, 2, rk1+8, rk1+9, rk1+10);
    arma::colvec betas = as<arma::colvec>(par[idx]);
    arma::colvec sigma2 = as<arma::colvec>(rcpmeObj[4]);
    double sigma2_1 = pow(sigma2(0),2); double sigma2_2 = pow(sigma2(1),2);
    arma::mat Sigma1 = pow(sigma2_1,2) * arma::eye<arma::mat>(lgt1,lgt1);
    arma::mat Sigma2 = pow(sigma2_2,2) * arma::eye<arma::mat>(lgt2,lgt2);
    arma::mat Sigma0 = arma::zeros<arma::mat>(lgt1,lgt2);
    arma::mat Sigma = join_cols(join_rows(Sigma1,Sigma0), join_rows(trans(Sigma0),Sigma2));
    
    arma::mat mure = arma::zeros<arma::colvec>(2);
    arma::mat Bre = arma::zeros<arma::mat>(2,2); Bre(0,0) = estiVarEA(3,3); 
    Bre(0,1) = estiVarEA(3,7); Bre(1,0) = estiVarEA(7,3); Bre(1,1) = estiVarEA(7,7);
    
    // distribution of Y_i | tau_i and tau_i
    arma::mat Z = arma::zeros<arma::mat>(lgt,6);
    arma::mat Z01 = arma::zeros<arma::mat>(lgt1,3);
    arma::mat Z02 = arma::zeros<arma::mat>(lgt2,3);
    arma::mat Z1 = arma::ones<arma::mat>(lgt1,3);
    arma::mat Z2 = arma::ones<arma::mat>(lgt2,3);
    if (model == "test"){Z1.col(1) = timeNoNA1; Z2.col(1) = timeNoNA2;}
    if (model == "bw"){Z1.col(1) = timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(0)); Z2.col(1) = timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(1));}
    Z1.col(2) = pow(pow(timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(par[3]+re(0)),2)+gamma,0.5);
    Z2.col(2) = pow(pow(timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(par[rk1+11]+re(1)),2)+gamma,0.5);
    Z = join_cols(join_rows(Z1,Z01),join_rows(Z02,Z2));
    
    arma::colvec mu = arma::zeros<arma::colvec>(lgt); mu = Z * betas;
    arma::mat V = arma::zeros<arma::mat>(lgt,lgt); V = (Z * B) * trans(Z) + Sigma;
    
    bool logtrue = true; double out = dmvnrmarma1d(trans(scoreNoNA), trans(mu), V, logtrue);
    out = out + dmvnrmarma1d(trans(re), trans(mure), Bre, logtrue);
    
    return out;
  }

  return out;
}