//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppDist.h>
//#include <truncnorm.h>

using namespace Rcpp;


// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function that computes pi(U; beta)

// Input: beta = (beta0, beta1): Vector of parameters
//        U_i: Covariate (Scalar)

// [[Rcpp::export]] 
double pi_U_beta_cpp(NumericVector beta,
                     double U_i){
  
  // Extract beta0 and beta1 from beta
  double beta0_pi = beta(0);
  double beta1_pi = beta(1);
  // Compute the linear predictor
  double num0 = beta0_pi + beta1_pi * U_i;
  // Compute the numerator of pi(U; beta)
  double num = exp(num0);
  // Compute the denominator of pi(U; beta)
  double den = 1 + num;
  // return "num/den"
  double val = num / den;
  return(val);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function that computes the pdf of a normal random variable

// Input: x: scalar
//        mean: scalar
//        sd: scalar

// [[Rcpp::export]] 
double dnorm_cpp(double x,
                 double mean,
                 double sd){
  double pi = 3.14159265358979323846264338328; 
  double val = 1 / ( sqrt(2 * pi) * sd) * 
    exp(-0.5 * ( (x - mean) / sd) * ( (x - mean) / sd));
  return(val);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write an "rnorm" function in C++

// Input: n: number of realizations to generate
//        mu: expected value of the random variables
//        sd: standard deviation of the random variables

// [[Rcpp::export]]
NumericVector rnorm_cpp(double n,
                        double mu,
                        double sd){
  
  // Initialize a vector to store the results
  NumericVector Z(n); // initialize
  
  // Run a for loop over 1:n
  for (int i = 0 ; i < n ; ++i){
    
    // Use "rnorm" to generate observations from the standard normal distribution
    NumericVector val_Z = rnorm(1); // Have to generate it as a vector for it to work
    double val_Z2 = val_Z(0); // This part is needed for some reason
    
    
    // multiply "val_Z" by "sd" and add in "mu"
    double val = val_Z2 * sd + mu;
    
    // Update "Z[i]"
    Z(i) = val;
  }
  
  // return "Z"
  return(Z);
}

// ---------------------------------------------------------------------------------------------------------------------------------------------


// [[Rcpp::export]]
NumericVector Test_rtruncnorm(double x,
                              double mean,
                              double sd,
                              double a,
                              double b){
  
  NumericVector val = rtruncnorm(x, mean, sd, a, b);
  return(val);
}

// VERIFIED!!

// [[Rcpp::export]]
double Test_dtruncnorm(double x,
                       double a,
                       double b){
  
  double val1 = dnorm_cpp(x,
                          0,
                          1);
  double val2 = R::pnorm(a, 
                         0, 
                         1, 
                         TRUE, 
                         FALSE);
  double val3 = R::pnorm(b, 
                         0, 
                         1, 
                         TRUE, 
                         FALSE);
  double val = val1 / (val3-val2);
  return(val);
}

// VERIFIED!!

// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function that computes f(u | X <= a) I(X* = a)
// Input: u: Value evaluating the pdf (Scalar)
//        Xstar: Observed covariate value (Scalar)
//        LLOD: Lower limit of detection (Scalar)
//        mu_U: Mean of U and X (Scalar)
//        sigma2_U: Variance of U (Scalar)
//        sigma2_X: Variance of X (Scalar)
//        MeasurementErrorLogi: Do we want to account for measurement error (Logical)

// [[Rcpp::export]]
double ProbUGivenXstarLessLLOD_func(double u,
                                    double LLOD,
                                    double mu_U,
                                    double mu_X,
                                    double sigma2_U,
                                    double sigma2_X,
                                    bool MeasurementErrorLogi = true){
  
  double val = 1; // Initialize
  
  
  // Proceed into cases depending if we want to account for measurement error or not
  if (MeasurementErrorLogi == true){
    
    // Compute f(u | Xstar <= LLOD)
    double val1 = dnorm_cpp(u,
                            mu_U,
                            sqrt(sigma2_U));
    double val2 = R::pnorm(LLOD, 
                           u - mu_U + mu_X, 
                           sqrt(sigma2_X - sigma2_U), 
                           TRUE, 
                           FALSE);
    double val3 = R::pnorm(LLOD, 
                           mu_X, 
                           sqrt(sigma2_X), 
                           TRUE, 
                           FALSE);
    val = val1 * val2 / val3;
    
    
  } else{
    
    // Compute f(x | Xstar <= LLOD)
    double val1 = dnorm_cpp(u,
                            mu_X,
                            sqrt(sigma2_X));
    double val2 = 1;
    double val3 = R::pnorm(LLOD, 
                           mu_X, 
                           sqrt(sigma2_X), 
                           TRUE, 
                           FALSE);
    val = val1 * val2 / val3;
    
  }
  
  return(val);
  
}

// VERIFIED!!

// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function that generates variables from either 
//      (i) a normal distribution with mean mu and variance sigma2
//      (ii) a truncated normal distribution with mean mu and variance sigma2 (pre-truncation)
// Input: X_i: Scalar covariate
//        LLOD: Lower limit of detection (Scalar)
//        mu_U: Mean of U (Scalar)
//        mu_X: Mean of X (Scalar)
//        sigma2_U: Variance of U (Scalar)
//        sigma2_X: Variance of X (Scalar)
//        MeasurementErrorLogi: Do we want to account for measurement error (Logical)
//        B: Number of iterations for the Monte Carlo approximation (Scalar; default is 10000)

// [[Rcpp::export]]
NumericVector GenerateRVs_func(double Xstar,
                               double LLOD,
                               double mu_U,
                               double mu_X,
                               double sigma2_U,
                               double sigma2_X,
                               int B = 10000,
                               bool MeasurementErrorLogi = true){
  
  NumericVector val(B); // Initialize
  
  // Break into cases depending if we want to account for measurement error or not
  if (MeasurementErrorLogi == true){
    
    // Generate random variables from [U|X]: normal with mean and variance defined below
    double mu_U_X = mu_U + sigma2_U / sigma2_X * (Xstar - mu_X);
    double sigma2_U_X = sigma2_U - (sigma2_U * sigma2_U / sigma2_X);
    
    // Generate B observations from [U|X]
    val = rnorm(B, mu_U_X, sqrt(sigma2_U_X));
    
  } else{ // "MeasurementErrorLogi == T"
    
    // Generate B observations from a truncated normal distribution
    //    with mean mu and variance sigma2 (pre-truncation)
    val = rtruncnorm(B, 
                     mu_X, 
                     sqrt(sigma2_X), 
                     R_NegInf, 
                     LLOD);
  }
  return(val);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function for the likelihood function corresponding to the observed 
//    data likelihood function

// Input: beta: Vector of parameters
//        Y_vec: Response (Vector)
//        Xstar_vec: Noisy Covariate (Vector)
//        Omega_vec: IPWs (Vector)
//        LLOD: Lower limit of detection (Scalar)
//        mu_U: Mean of U (Scalar)
//        mu_X: Mean of X (Scalar)
//        sigma2_U: Variance of U (Scalar)
//        sigma2_X: Variance of X (Scalar)
//        B: Number of iterations for the Monte Carlo approximation (Scalar; default is 10000)
//        MeasurementErrorLogi: Do we want to account for measurement error (Logical)

// [[Rcpp::export]]
NumericVector LogLikFcn_cpp(NumericVector beta,
                            NumericVector Y_vec,
                            NumericVector Xstar_vec,
                            NumericVector Omega_vec,
                            double LLOD,
                            double mu_U,
                            double mu_X,
                            double sigma2_U,
                            double sigma2_X,
                            int B = 10000,
                            bool MeasurementErrorLogi = true){
  
  // Get the number of units
  int n = Y_vec.size();
  
  // Run a for loop over each unit, 
  //    Generate observations from [U | X; alpha]
  //    Take a Monte Carlo expectation
  NumericVector val(n); // Initialize
  
  for (int i = 0 ; i < n ; ++i){ 
    
    // Extract elements from "Y_vec", "Xstar_vec", and "Omega_vec"
    double Y_i = Y_vec(i);
    double Xstar_i = Xstar_vec(i);
    double IPW_i = Omega_vec(i);
    
    // Break into cases depending if we want to account for measurement error or not
    if (MeasurementErrorLogi == true){
      
      // Proceed to generate B random variables from [U | Xstar]
      NumericVector U_i_B = GenerateRVs_func(Xstar_i,
                                             LLOD,
                                             mu_U,
                                             mu_X,
                                             sigma2_U,
                                             sigma2_X,
                                             B,
                                             true); // generated from normal distribution [U|X]
      
      // Run a loop to compute [Y_i | U_{ib}]
      double val_Z = 0; // initialize
      
      for (int b = 0 ; b < B ; ++b){ 
        
        // Compute pi(beta; U_i[b])
        double U_ib = U_i_B(b);
        double prob = pi_U_beta_cpp(beta, U_ib);
        
        // Need to compute f(U_{ib} | Xstar_i) / g(U_{ib} | Xstar_i)
        double omega_b = 1; // Initialize; correct if "Xstar_i > LLOD"
        // THIS CONDITION WILL NEVER BE SATISFIED!
        if (Xstar_i <= LLOD){
          
          // Compute f(U_{ib} | Xstar_i)
          double fb_val = ProbUGivenXstarLessLLOD_func(U_ib,
                                                       LLOD,
                                                       mu_U,
                                                       mu_X,
                                                       sigma2_U,
                                                       sigma2_X,
                                                       true);
          
          // Compute g(U_{ib} | Xstar_i)
          double mu_U_X = mu_U + sigma2_U / sigma2_X * (Xstar_i - mu_X);
          double sigma2_U_X = sigma2_U - (sigma2_U * sigma2_U / sigma2_X);
          double gb_val = dnorm_cpp(U_ib,
                                    mu_U_X,
                                    sqrt(sigma2_U_X));
          
          // Compute f(U_{ib} | Xstar_i) / g(U_{ib} | Xstar_i)
          omega_b = fb_val / gb_val;
        }
        
        // Compute [Y_i | U_{ib}] and update "val_Z"
        val_Z += ( (Y_i * prob) + ( (1-Y_i) * (1-prob) ) ) * omega_b;
      }
      
      // Update "val_Z"
      val_Z = val_Z / B;
      
      // Store "IPW_i * log(val_Z)"
      val(i) = IPW_i * log(val_Z);
    }
    
    else{ // "MeasurementErrorLogi == T"
      
      // Break into cases if "Xstar_i > LLOD" or not
      if (Xstar_i > LLOD){
        
        // No need to do any integration, proceed to compute contribution to the log likelihood function
        double prob = pi_U_beta_cpp(beta, Xstar_i);
        double val_Z = (Y_i * prob) + ( (1-Y_i) * (1-prob) );
        
        // Store "IPW_i * log(val_Z)"
        val(i) = IPW_i * log(val_Z);
        
      } else{
        
        // Proceed to generate B random variables from [X | Xstar]
        NumericVector X_i_B = GenerateRVs_func(Xstar_i,
                                               LLOD,
                                               mu_U,
                                               mu_X,
                                               sigma2_U,
                                               sigma2_X,
                                               B,
                                               false); // generated from truncated normal distribution
        
        // Run a loop to compute [Y_i | X_{ib}]
        double val_Z = 0; // initialize
        
        for (int b = 0 ; b < B ; ++b){ 
          
          // Compute pi(beta; X_i[b])
          double X_ib = X_i_B(b);
          double prob_b = pi_U_beta_cpp(beta, X_ib);
          
          double omega_b = 1;
          
          
          // Compute [Y_i | U_{ib}] and update "val_Z"
          val_Z += ( (Y_i * prob_b) + ( (1-Y_i) * (1-prob_b) ) ) * omega_b;
        }
        
        // Update "val_Z"
        val_Z = val_Z / B;
        
        // Store "IPW_i * log(val_Z)"
        val(i) = IPW_i * log(val_Z);
      }
    }
  } // end of for loop (1:n)
  
  
  return(val);
}



// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function for the score function corresponding to the observed 
//    data likelihood function

// Input: beta: Vector of parameters
//        Y_vec: Response (Vector)
//        Xstar_vec: Noisy Covariate (Vector)
//        Omega_vec: IPWs (Vector)
//        LLOD: Lower limit of detection (Scalar)
//        mu_U: Mean of U (Scalar)
//        mu_X: Mean of X (Scalar)
//        sigma2_U: Variance of U (Scalar)
//        sigma2_X: Variance of X (Scalar)
//        B: Number of iterations for the Monte Carlo approximation (Scalar; default is 10000)
//        MeasurementErrorLogi: Do we want to account for measurement error (Logical)


// [[Rcpp::export]]
NumericMatrix ScoreFunction_cpp(NumericVector beta,
                                NumericVector Y_vec,
                                NumericVector Xstar_vec,
                                NumericVector Omega_vec,
                                double LLOD,
                                double mu_U,
                                double mu_X,
                                double sigma2_U,
                                double sigma2_X,
                                int B = 10000,
                                bool MeasurementErrorLogi = true){
  
  // Get the number of units
  int n = Y_vec.size();
  
  // Run a for loop over each unit, 
  //    Generate observations from [U | X; alpha]
  //    Take a Monte Carlo expectation
  NumericMatrix val(n,2); // Initialize
  
  for (int i = 0 ; i < n ; ++i){ 
    
    // Extract elements from "Y", "Xstar", and "Omega_vec"
    double Y_i = Y_vec(i);
    double Xstar_i = Xstar_vec(i);
    double IPW_i = Omega_vec(i);
    
    // Break into cases depending if we want to account for measurement error or not
    if (MeasurementErrorLogi == true){
      
      // Proceed to generate B random variables from [U | Xstar]
      NumericVector U_i_B = GenerateRVs_func(Xstar_i,
                                             mu_U,
                                             mu_X,
                                             sigma2_U,
                                             sigma2_X,
                                             B); // generated from normal distribution [U|X]
      
      // Run a loop to compute [Y_i | U_{ib}]
      double val_Z = 0; // initialize
      
      for (int b = 0 ; b < B ; ++b){ 
        
        // Compute pi(beta; U_i[b])
        double U_ib = U_i_B(b);
        double prob = 0; // Initialize
        prob = pi_U_beta_cpp(beta, U_ib);
        
        
        // Compute [Y_i | U_{ib}] and update "val_Z"
        val_Z += ( (Y_i * prob) + ( (1-Y_i) * (1-prob) ) );
      }
      
      // Update "val_Z"
      val_Z = val_Z / B;
      
      // Run another for loop to get components we want
      double val1 = 0; // initialize
      double val2 = 0; // initialize
      
      for (int b = 0 ; b < B ; ++b){ 
        
        // Compute pi(beta; U_i[b])
        double U_ib = U_i_B(b);
        double prob = 0; // Initialize
        prob = pi_U_beta_cpp(beta, U_ib);
        
        
        // Get the two components of the estimating function
        val1 += (2*Y_i-1) * prob * (1 - prob) / val_Z;
        val2 += (2*Y_i-1) * prob * (1 - prob) * U_ib / val_Z;
      }
      
      // Store "val1 * IPW_i" and "val2 * IPW_i"
      val(i,0) = IPW_i * val1 / B;
      val(i,1) = IPW_i * val2 / B;
      
    } else{ // "MeasurementErrorLogi == T"
      
      double val1 = 0; // initialize
      double val2 = 0; // initialize
      
      // No need to do any integration, proceed to compute contribution to the score function
      double prob = 0; // Initialize
      prob = pi_U_beta_cpp(beta, Xstar_i);
      
      double val_Z = (Y_i * prob) + ( (1-Y_i) * (1-prob) );
      
      // Get the two components of the estimating function
      val1 = (2*Y_i-1) * prob * (1 - prob) / val_Z;
      val2 = (2*Y_i-1) * prob * (1 - prob) * Xstar_i / val_Z;
      // val1 = (Y_i - prob);
      // val2 = (Y_i - prob) * Xstar_i;
      
      // Store val1 and val2
      val(i,0) = IPW_i * val1; 
      val(i,1) = IPW_i * val2;
      
    }
  } // end of for loop (1:n)
  return(val);
}