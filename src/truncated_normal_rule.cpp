# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

List jacobi_eigenvalue(int n, NumericMatrix A, int it_max);
List moment_method(int n, NumericVector moment);
double r8_factorial2 ( int n );
double r8_choose ( int n, int k );
NumericVector moments_normal ( int m, double mu, double sigma );
NumericVector moments_truncated_normal_ab ( int m, double mu, double sigma,
  double a, double b );
NumericVector moments_truncated_normal_a ( int m, double mu, double sigma,
  double a );
NumericVector moments_truncated_normal_b ( int m, double mu, double sigma,
  double b );
double normal_01_cdf ( double x );
double normal_01_pdf ( double x );
double r8_factorial ( int n );
double r8_huge ();
double r8_mop (int i);
NumericMatrix r8mat_cholesky_factor_upper (int n, NumericMatrix a);
NumericMatrix r8mat_copy_new ( int m, int n, double a1[] );
NumericVector r8mat_diag_get_vector ( int n, double a[], double v[] );
NumericMatrix r8mat_identity  ( int n, double a[] );
double truncated_normal_ab_moment ( int order, double mu, double s, double a,
  double b );
double truncated_normal_a_moment ( int order, double mu, double s, double a );
double truncated_normal_b_moment ( int order, double mu, double s, double b );

//****************************************************************************80

// [[Rcpp::export]]
List jacobi_eigenvalue(int n, NumericMatrix A, int it_max) {
  NumericMatrix V(n, n);
  NumericVector d(n);
  NumericVector bw(n);
  NumericVector zw(n);
  
  // Initialize V to the identity matrix and d to the diagonal of A
  for (int i = 0; i < n; i++) {
    V(i, i) = 1.0;
    d[i] = A(i, i);
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  
  int it_num = 0;
  int rot_num = 0;
  
  for (it_num = 0; it_num < it_max; it_num++) {
    double thresh = 0.0;
    
    // Calculate the off-diagonal norm (threshold)
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        thresh += A(i, j) * A(i, j);
      }
    }
    thresh = sqrt(thresh);
    
    if (thresh == 0.0) break;
    
    for (int p = 0; p < n; p++) {
      for (int q = p + 1; q < n; q++) {
        double gapq = 10.0 * fabs(A(p, q));
        double termp = gapq + fabs(d[p]);
        double termq = gapq + fabs(d[q]);
        
        if (4 * it_num > it_max && termp == fabs(d[p]) && termq == fabs(d[q])) {
          A(p, q) = 0.0;
        } else if (fabs(A(p, q)) > thresh) {
          double t = (d[q] - d[p]) / (2.0 * A(p, q));
          double c = 1.0 / sqrt(1 + t * t);
          double s = t * c;
          double tau = s / (1.0 + c);
          double h = t * A(p, q);
          
          zw[p] -= h;
          zw[q] += h;
          d[p] -= h;
          d[q] += h;
          
          A(p, q) = 0.0;
          
          for (int j = 0; j < p; j++) {
            double g = A(j, p);
            double h = A(j, q);
            A(j, p) = g - s * (h + g * tau);
            A(j, q) = h + s * (g - h * tau);
          }
          for (int j = p + 1; j < q; j++) {
            double g = A(p, j);
            double h = A(j, q);
            A(p, j) = g - s * (h + g * tau);
            A(j, q) = h + s * (g - h * tau);
          }
          for (int j = q + 1; j < n; j++) {
            double g = A(p, j);
            double h = A(q, j);
            A(p, j) = g - s * (h + g * tau);
            A(q, j) = h + s * (g - h * tau);
          }
          
          // Update the eigenvector matrix
          for (int j = 0; j < n; j++) {
            double g = V(j, p);
            double h = V(j, q);
            V(j, p) = g - s * (h + g * tau);
            V(j, q) = h + s * (g - h * tau);
          }
          rot_num++;
        }
      }
    }
    
    for (int i = 0; i < n; i++) {
      bw[i] += zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
  
  return List::create(Named("values") = d, Named("vectors") = V,
                      Named("iterations") = it_num, Named("rotations") = rot_num);
}
//****************************************************************************80
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List moment_method(int n, Rcpp::NumericVector moment) {
  // Define the (n+1)x(n+1) Hankel matrix from the moments
  Function chol("chol");   
  Function eigen("eigen");
  Rcpp::NumericMatrix H(n + 1, n + 1);
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      H(i, j) = moment[i + j];
    }
  }
  
  // Compute the Cholesky factorization of the Hankel matrix
  Rcpp::NumericMatrix R = chol(H);
  
  // Compute alpha and beta coefficients from the Cholesky factor
  Rcpp::NumericVector alpha(n);
  Rcpp::NumericVector beta(n - 1);
  alpha[0] = R(1, 0) / R(0, 0);
  for (int i = 1; i < n; i++) {
    alpha[i] = R(i, i+1) / R(i, i) - R(i-1, i) / R(i-1, i-1);
  }
  for (int i = 0; i < n - 1; i++) {
    beta[i] = R(i+1, i+1) / R(i, i);
  }
  
  // Construct the Jacobi matrix and perform custom jacobi eigenvalue decomposition
  NumericMatrix jacobi(n, n);
  for (int i = 0; i < n; i++) {
    jacobi(i, i) = alpha[i];
    if (i < n - 1) {
      jacobi(i, i + 1) = beta[i];
      jacobi(i + 1, i) = beta[i];
    }
  }
  
  List eigenResults = eigen(jacobi);  // Assume it_max = 100 for example
  
  NumericVector eigenvalues = eigenResults["values"];
  NumericMatrix eigenvectors = eigenResults["vectors"];
  
  // Calculate weights using the first moment and the first column of eigenvectors
  NumericVector w(n);
  for (int i = 0; i < n; i++) {
    w[i] = moment[0] * pow(eigenvectors(0, i), 2);
  }
  
  // Return the points and weights
  return Rcpp::List::create(Rcpp::Named("points") = eigenvalues,
                            Rcpp::Named("weights") = w, Rcpp::Named("R") = R, Rcpp::Named("jacobi") = jacobi, Rcpp::Named("H")= H);
}
//****************************************************************************80
// Function to compute the double factorial
double r8_factorial2(int n) {
  double value = 1;
  while (n > 1) {
    value *= n;
    n -= 2;
  }
  return value;
}

// Function to compute binomial coefficients
double r8_choose(int n, int k) {
  if (k > n) {
    return 0;
  }
  if (k == 0 || k == n) {
    return 1;
  }
  double num = 1;
  for (int i = 1; i <= k; ++i) {
    num *= (n + 1 - i);
  }
  double denom = 1;
  for (int i = 1; i <= k; ++i) {
    denom *= i;
  }
  return num / denom;
}

// [[Rcpp::export]]
NumericVector moments_normal(int m, double mu, double sigma) {
  NumericVector w(m);
  int j_hi;
  double t;
  
  for (int k = 0; k < m; k++) {
    t = 0.0;
    j_hi = k / 2;
    for (int j = 0; j <= j_hi; j++) {
      t += r8_choose(k, 2 * j) * r8_factorial2(2 * j - 1) 
      * pow(sigma, 2 * j) * pow(mu, k - 2 * j);
    }
    w[k] = t;
  }
  
  return w;
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector moments_truncated_normal_ab(int m, double mu, double sigma, double a, double b) {
  if (a >= b) {
    stop("Lower truncation limit A must be less than upper truncation limit B.");
  }
  
  NumericVector w(m);
  
  for (int order = 0; order < m; order++) {
    w[order] = truncated_normal_ab_moment(order, mu, sigma, a, b);
  }
  
  return w;
}
//****************************************************************************80
// [[Rcpp::export]]
NumericVector moments_truncated_normal_a(int m, double mu, double sigma, double a) {
  NumericVector w(m);
  
  for (int order = 0; order < m; order++) {
    w[order] = truncated_normal_a_moment(order, mu, sigma, a);
  }
  
  return w;
}
//****************************************************************************80

// [[Rcpp::export]]
NumericVector moments_truncated_normal_b(int m, double mu, double sigma, double b) {
  NumericVector w(m);
  
  for (int order = 0; order < m; order++) {
    w[order] = truncated_normal_b_moment(order, mu, sigma, b);
  }
  
  return w;
}
//****************************************************************************80

double normal_01_pdf(double x) {
  return R::dnorm(x, 0.0, 1.0, false);
}
//****************************************************************************80

double normal_01_cdf(double x) {
  return R::pnorm(x, 0.0, 1.0, true, false);
}
//****************************************************************************80

//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}

//****************************************************************************80

// [[Rcpp::export]]
double r8_huge() {
  return 1.0E+30;  // Define a huge number that is not considered infinite.
}
// [[Rcpp::export]]
double r8_mop(int i) {
  // Return -1 if i is odd, 1 if i is even
  return (i % 2 == 0) ? 1.0 : -1.0;
}
// [[Rcpp::export]]
NumericMatrix r8mat_cholesky_factor_upper(int n, NumericMatrix a) {
  // Create a copy of the matrix 'a' to store the Cholesky factors
  NumericMatrix c(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      c(i, j) = a(i, j);
    }
  }
  
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++) {
      c(j, i) = 0.0;
    }
    for (int i = j; i < n; i++) {
      double sum = c(i, j);
      for (int k = 0; k < j; k++) {
        sum -= c(j, k) * c(i, k);
      }
      if (i == j) {
        if (sum <= 0.0) {
          stop("Matrix is not positive definite.");
        }
        c(j, i) = sqrt(sum);
      } else {
        if (c(j, j) != 0.0) {
          c(j, i) = sum / c(j, j);
        } else {
          c(j, i) = 0.0;
        }
      }
    }
  }
  
  return c;
}
//****************************************************************************80

// [[Rcpp::export]]
NumericMatrix r8mat_copy_new(int m, int n, NumericMatrix a1) {
  NumericMatrix a2(m, n);
  
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      a2(i, j) = a1(i, j);
    }
  }
  
  return a2;
}
//****************************************************************************80
// [[Rcpp::export]]
NumericVector r8mat_diag_get_vector(int n, NumericMatrix a) {
  NumericVector v(n);
  
  for (int i = 0; i < n; i++) {
    v[i] = a(i, i);  // Accessing the diagonal elements directly
  }
  
  return v;
}
//****************************************************************************80
// [[Rcpp::export]]
NumericMatrix r8mat_identity(int n) {
  NumericMatrix a(n, n);
  
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      if (i == j) {
        a(i, j) = 1.0;  // Set diagonal elements to 1
      } else {
        a(i, j) = 0.0;  // Set off-diagonal elements to 0
      }
    }
  }
  
  return a;
}


//****************************************************************************80
// [[Rcpp::export]]
double truncated_normal_ab_moment(int order, double mu, double s, double a, double b) {
  if (order < 0) {
    stop("ORDER must be non-negative.");
  }
  if (s <= 0.0) {
    stop("Standard deviation S must be positive.");
  }
  if (b <= a) {
    stop("Upper bound B must be greater than lower bound A.");
  }
  
  double a_h = (a - mu) / s;
  double a_pdf = normal_01_pdf(a_h);
  double a_cdf = normal_01_cdf(a_h);
  
  if (a_cdf == 0.0) {
    stop("CDF at A is zero, leading to division by zero.");
  }
  
  double b_h = (b - mu) / s;
  double b_pdf = normal_01_pdf(b_h);
  double b_cdf = normal_01_cdf(b_h);
  
  if (b_cdf == 0.0) {
    stop("CDF at B is zero, leading to division by zero.");
  }
  
  double moment = 0.0;
  double irm2 = 0.0;
  double irm1 = 0.0;
  double ir;
  
  for (int r = 0; r <= order; r++) {
    if (r == 0) {
      ir = 1.0;
    } else if (r == 1) {
      ir = - (b_pdf - a_pdf) / (b_cdf - a_cdf);
    } else {
      ir = (double)(r - 1) * irm2 - (pow(b_h, r - 1) * b_pdf - pow(a_h, r - 1) * a_pdf) / (b_cdf - a_cdf);
    }
    
    moment += R::choose(static_cast<double>(order), static_cast<double>(r)) * pow(mu, order - r) * pow(s, r) * ir;
    
    irm2 = irm1;
    irm1 = ir;
  }
  
  return moment;
}
// [[Rcpp::export]]
double truncated_normal_a_moment(int order, double mu, double sigma, double a) {
  if (order < 0) {
    stop("Order must be non-negative.");
  }
  
  double moment = r8_mop(order) * truncated_normal_b_moment(order, -mu, sigma, -a);
  return moment;
}
//****************************************************************************80

// [[Rcpp::export]]
double truncated_normal_b_moment(int order, double mu, double sigma, double b) {
  if (order < 0) {
    Rcpp::stop("ORDER must be non-negative.");
  }
  
  // Compute the standardized upper limit
  double h = (b - mu) / sigma;
  double h_pdf = R::dnorm(h, 0.0, 1.0, false); // Normal PDF
  double h_cdf = R::pnorm(h, 0.0, 1.0, true, false); // Normal CDF
  
  if (h_cdf == 0.0) {
    Rcpp::stop("CDF of standardized upper limit is 0, which may lead to division by zero.");
  }
  
  double f = h_pdf / h_cdf;
  double moment = 0.0;
  double irm2 = 0.0;
  double irm1 = 0.0;
  double ir = 0.0;
  
  for (int r = 0; r <= order; r++) {
    if (r == 0) {
      ir = 1.0;
    } else if (r == 1) {
      ir = -f;
    } else {
      ir = -pow(h, r - 1) * f + static_cast<double>(r - 1) * irm2;
    }
    
    moment += R::choose(static_cast<double>(order), static_cast<double>(r))
      * pow(mu, order - r) * pow(sigma, r) * ir;
    
    irm2 = irm1;
    irm1 = ir;
  }
  
  return moment;
}