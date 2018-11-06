#include <RcppArmadillo.h>
using namespace Rcpp;


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

// VERSION PAQUID
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
double lvsblNCgen(NumericVector param, List data, int nq, NumericVector grp, NumericVector weights, NumericVector nodes, String scorevar, String timevar, String covariate, String REadjust){
  double Beta0=param(0); double Beta1=param(1); double Beta2=param(2);
  double tau=param(3); double Utau=param(11); double sigma=param(4);
  
  int N = max(grp);
  double out = 0;
  for (int i = 0; i<N; i++){
    DataFrame datai = as<DataFrame>(data[i]);
    
    NumericVector score = datai[scorevar];
    LogicalVector indNoNA = !is_na(score);
    arma::colvec scoreNoNA = as<arma::colvec>(score[indNoNA]);
    NumericVector time = datai[timevar];
    arma::colvec timeNoNA = as<arma::colvec>(time[indNoNA]);

    int lgt = scoreNoNA.n_elem;
    double res = 0;
    
    if (lgt != 0){

      arma::mat U = arma::zeros<arma::mat>(3,3);
      arma::mat B = arma::zeros<arma::mat>(3,3);
      
      if (covariate == "NULL"){
        U(0,0) = param(5); U(0,1) = param(6); U(0,2)=param(7); U(1,1)=param(8); U(1,2)=param(9); U(2,2)=param(10);
        B = trans(U) * U;
        
        arma::colvec Betas = arma::zeros<arma::colvec>(3);
        Betas(0) = Beta0; Betas(1) = Beta1; Betas(2)=Beta2;
        
        arma::mat Zk = arma::ones<arma::mat>(lgt,3);
        Zk.col(1) = timeNoNA;
        arma::colvec muk = arma::zeros<arma::colvec>(lgt);
        arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
        for (int k = 0; k < nq; ++k){
          Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+Utau*nodes(k)),2)+0.1,0.5);
          muk = Zk * Betas;
          Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt); 
          double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
          res = res + (f * weights(k));
        }
        out = out + log(res);
      }
      
      if (covariate != "NULL"){
        NumericVector adjust = datai[covariate];
        double adjusti = adjust(0);
        
        if (REadjust == "no"){
          U(0,0) = param(5); U(0,1) = param(6); U(0,2)=param(7); U(1,1)=param(8); U(1,2)=param(9); U(2,2)=param(10); B = trans(U) * U;
        }
        
        if (REadjust == "prop"){
          U(0,0) = (1 + param(16)*adjusti) * param(5); U(0,1) = (1 + param(16)*adjusti) * param(6); U(0,2)=(1 + param(16)*adjusti) * param(7);
          U(1,1)=(1 + param(16)*adjusti) * param(8); U(1,2)=(1 + param(16)*adjusti) * param(9); U(2,2)=(1 + param(16)*adjusti) * param(10);
          B = trans(U) * U;
          Utau = (1 + param(16)*adjusti) * Utau;
        }
        
        if (REadjust == "yes"){
          U(0,0) = param(5)+param(16)*adjusti; U(0,1) = param(6)+param(17)*adjusti; U(0,2)=param(7)+param(18)*adjusti;
          U(1,1)=param(8)+param(19)*adjusti; U(1,2)=param(9)+param(20)*adjusti; U(2,2)=param(10)+param(21)*adjusti;
          B = trans(U) * U;
          Utau = param(11) + param(22)*adjusti;
        }
        
        arma::colvec Betas = arma::zeros<arma::colvec>(3);
        Betas(0) = Beta0 + param(12)*adjusti; Betas(1) = Beta1 + param(13)*adjusti; Betas(2)=Beta2+param(14)*adjusti;
        
        arma::mat Zk = arma::ones<arma::mat>(lgt,3);
        Zk.col(1) = timeNoNA;
        arma::colvec muk = arma::zeros<arma::colvec>(lgt);
        arma::mat Vk = arma::zeros<arma::mat>(lgt,lgt);
        for (int k = 0; k < nq; ++k){
          Zk.col(2) = pow(pow(timeNoNA - arma::ones<arma::colvec>(lgt)*(tau+param(15)*adjusti+Utau*nodes(k)),2)+0.1,0.5);
          muk = Zk * Betas;  Vk = (Zk * B) * trans(Zk) + pow(sigma,2) * arma::eye<arma::mat>(lgt,lgt); 
          double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
          res = res + (f * weights(k));
        }
        out = out + log(res);
      }
    }
  }
  return -out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double bilvsblNC(NumericVector param, List data, int nq, bool adapt, NumericVector grp, NumericVector weights, NumericMatrix nodes, List newnodes, List newweights, String scorevar1, String scorevar2, String timevar, String covariate, String REadjust){
  double Beta0_1=param(0); double Beta1_1=param(1); double Beta2_1=param(2);
  double tau_1=param(3); double Utau_1=param(11); double sigma_1=param(4);
  
  double Beta0_2=param(12); double Beta1_2=param(13); double Beta2_2=param(14);
  double tau_2=param(15); double Utau_2=param(23); double sigma_2=param(16);

 // non adaptive
 arma::mat nnodes(pow(nq,2),2);
 NumericVector nweights(pow(nq,2));
  if (adapt == FALSE){
    for (int i = 0; i<pow(nq,2); i++){
      nnodes(i,0) = nodes(i,0) * Utau_1 + nodes(i,1) * param(24);
      nnodes(i,1) = nodes(i,1) * Utau_2;
      nweights = weights;
    }
  }
  
  // adaptative
  arma::rowvec mure(2); arma::mat Bre(2,2);
  if (adapt == TRUE){
    arma::mat Ure(2,2);
    Ure(0,0) = Utau_1; Ure(1,1) = Utau_2; Ure(0,1) = param(24); Bre = trans(Ure) * Ure;
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
    
    int lgt1 = scoreNoNA1.n_elem; int lgt2 = scoreNoNA2.n_elem; int lgt = lgt1 + lgt2;
    arma::colvec scoreNoNA = join_cols(scoreNoNA1, scoreNoNA2);
    double res = 0;
    
    if (lgt != 0){
      
      arma::mat U = arma::zeros<arma::mat>(6,6);
      arma::mat B = arma::zeros<arma::mat>(6,6);
      
      if (covariate == "NULL"){

        U(0,0) = param(5); U(0,1) = param(6); U(0,2)=param(7); U(1,1)=param(8); U(1,2)=param(9); U(2,2)=param(10);
        U(3,3) = param(17); U(3,4) = param(18); U(3,5) = param(19); U(4,4) = param(20); U(4,5) = param(21); U(5,5) = param(22);
        U(0,3) = param(25); U(0,4) = param(26); U(0,5) = param(27); U(1,3)=param(28); U(1,4)=param(29); U(1,5) = param(30);
        U(2,3) = param(31); U(2,4) = param(32); U(2,5) = param(33);
        
        B = trans(U) * U;
        
        arma::colvec Betas = arma::zeros<arma::colvec>(6);
        Betas(0) = Beta0_1; Betas(1) = Beta1_1; Betas(2)=Beta2_1;
        Betas(3) = Beta0_2; Betas(4) = Beta1_2; Betas(5)=Beta2_2;
        
        arma::mat Zk = arma::zeros<arma::mat>(lgt,6);
        arma::mat Zk01 = arma::zeros<arma::mat>(lgt1,3);
        arma::mat Zk02 = arma::zeros<arma::mat>(lgt2,3);
        arma::mat Zk1 = arma::ones<arma::mat>(lgt1,3);
        arma::mat Zk2 = arma::ones<arma::mat>(lgt2,3);
        Zk1.col(1) = timeNoNA1; Zk2.col(1) = timeNoNA2;
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
          Zk1.col(2) = pow(pow(timeNoNA1 - arma::ones<arma::colvec>(lgt1)*(tau_1+nnodes(k,0)),2)+0.1,0.5);
          Zk2.col(2) = pow(pow(timeNoNA2 - arma::ones<arma::colvec>(lgt2)*(tau_2+nnodes(k,1)),2)+0.1,0.5);
          Zk = join_cols(join_rows(Zk1,Zk01),join_rows(Zk02,Zk2));
          muk = Zk * Betas; Vk = (Zk * B) * trans(Zk) + Sigma;
          double f = dmvnrmarma1d(trans(scoreNoNA),trans(muk),Vk);
          res = res + (f * qre(k) * nweights(k));
        }
        out = out + log(res);
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
double IndRePostDis(arma::rowvec re, DataFrame data, List rcpmeObj, String scorevar, String timevar){
  
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
      mus(i) = par(0) + re(0) + (par(1) + re(1)) * time(i) + (par(2)+re(2)) * pow(pow(time(i) - par(3) - re(3),2)+0.1,0.5);
    }
  }

  double adjusti;
  if (covariate != "NULL"){
    NumericVector par = as<NumericVector>(rcpmeObj[6]);
    IntegerVector idx = IntegerVector::create(0, 1, 2, 3, 4, 5, 6, 7);
    NumericVector betas = par[idx];
    NumericVector adjust = data[covariate];
    adjusti = adjust(0);
    for (int i =0; i<lgt; i++){
      mus(i) = par(0) + par(1)*adjusti + re(0) + (par(2) + par(3)*adjusti + re(1)) * time(i) + (par(4) + par(5)*adjusti + re(2)) * pow(pow(time(i) - par(6) - par(7)*adjusti - re(3),2)+0.1,0.5);
    }
  }

  double sdres = as<double>(rcpmeObj[4]);
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
double IndRePostDis2(double re, DataFrame data, List rcpmeObj, String scorevar, String timevar){
  
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
  double sdtau = sqrt(estiVarEA(3,3)); double sdres = par[4];
  
  // distribution of Y_i | tau_i and tau_i
  arma::colvec mus(lgt);
  arma::mat Zk = arma::ones<arma::mat>(lgt,3); Zk.col(1) = time;
  Zk.col(2) = pow(pow(time - arma::ones<arma::colvec>(lgt)*(par[3]+re),2)+0.1,0.5);
  mus = Zk * betas; arma::mat B = Zk * varB * trans(Zk) + pow(sdres,2)*arma::eye<arma::mat>(lgt,lgt);
  bool logtrue = true; double out = dmvnrmarma1d(trans(score),trans(mus),B,logtrue);
  // out = out - 0.5*(log(2*M_PI) + pow(re,2));
  out = out - 0.5*(log(2*M_PI)+pow(re/sdtau,2)) - log(sdtau);
  
  return out;
}