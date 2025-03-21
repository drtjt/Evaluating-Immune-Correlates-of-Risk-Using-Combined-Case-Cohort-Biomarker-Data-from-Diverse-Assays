# ------------------------------------------------------------------------------
# Simulation Study for 
#   "Meta-Analysis Approaches for Evaluating Immune Correlates of Risk 
#     Using Combined Case-Cohort Biomarker Data from Diverse Assays"
#
# Date Created: March 21, 2025
#
# Author: Trevor Thomson
#         tthomson@fredhutch.org
# ------------------------------------------------------------------------------


rm(list =ls())
set.seed(5347547)

library(MASS)
library(pbapply)
library(nleqslv)
library(Rcpp)
library(RcppArmadillo)
library(stringr)
library(numDeriv) # for "jacobian"
library(mvtnorm)
library(cubature) # for "cuhre"

# Specify the initial sizes of the datasets
N0 = 250
N1 = 15000
N2 = 15000

# Specify the size of the first phase sampling under the case-cohort study design
n0 = N0 # (no case-cohort sampling)
n1 = 0.08 * N1
n2 = 0.08 * N2

# Specify beta
Beta = c(-2.5, -1.5)

# Specify the mean and variance of U
mu_U_vec = c(2.8, 2.6, 1.9)
sigma2_U_vec = c(0.4, 0.5, 0.2)

# Specify the variance of measurement errors across cohorts
MeasurementErrorVariance_Setting = 1
# NEED TO CHECK!
if (MeasurementErrorVariance_Setting == 1){
  sigma2_eps_1_vec = sigma2_U_vec * (1.05-1) # Ratio of sigma2_X / sigma2_U = 1.05
  sigma2_eps_2_vec = sigma2_U_vec * (1.05-1) # Ratio of sigma2_X / sigma2_U = 1.05
}
if (MeasurementErrorVariance_Setting == 2){
  sigma2_eps_1_vec = sigma2_U_vec * (1.25-1) # Ratio of sigma2_X / sigma2_U = 1.25
  sigma2_eps_2_vec = sigma2_U_vec * (1.25-1) # Ratio of sigma2_X / sigma2_U = 1.25
}

# Specify the value of delta
delta_vec = c(0.6, 0.6, 0.6)

# Have a setting for the LLOD
Setting_LLOD = 3
if (Setting_LLOD == 1){
  LLOD_X1 = c(-Inf, -Inf, -Inf)
  LLOD_X2 = c(-Inf, -Inf, -Inf)
}
if (Setting_LLOD == 2){
  LLOD_X1 = c(qnorm(0.10,
                    mean = mu_U_vec[1],
                    sd = sqrt(sigma2_U_vec[1] + sigma2_eps_1_vec[1])),
              qnorm(0.10,
                    mean = mu_U_vec[2],
                    sd = sqrt(sigma2_U_vec[2] + sigma2_eps_1_vec[2])),
              qnorm(0.10,
                    mean = mu_U_vec[3],
                    sd = sqrt(sigma2_U_vec[3] + sigma2_eps_1_vec[3])))
  LLOD_X2 = c(qnorm(0.10,
                    mean = mu_U_vec[1] + delta_vec[1],
                    sd = sqrt(sigma2_U_vec[1] + sigma2_eps_2_vec[1])),
              qnorm(0.10,
                    mean = mu_U_vec[2] + delta_vec[2],
                    sd = sqrt(sigma2_U_vec[2] + sigma2_eps_2_vec[2])),
              qnorm(0.10,
                    mean = mu_U_vec[3] + + delta_vec[3],
                    sd = sqrt(sigma2_U_vec[3] + sigma2_eps_2_vec[3])))
}
if (Setting_LLOD == 3){
  LLOD_X1 = c(qnorm(0.20,
                    mean = mu_U_vec[1],
                    sd = sqrt(sigma2_U_vec[1] + sigma2_eps_1_vec[1])),
              qnorm(0.20,
                    mean = mu_U_vec[2],
                    sd = sqrt(sigma2_U_vec[2] + sigma2_eps_1_vec[2])),
              qnorm(0.20,
                    mean = mu_U_vec[3],
                    sd = sqrt(sigma2_U_vec[3] + sigma2_eps_1_vec[3])))
  LLOD_X2 = c(qnorm(0.20,
                    mean = mu_U_vec[1] + delta_vec[1],
                    sd = sqrt(sigma2_U_vec[1] + sigma2_eps_2_vec[1])),
              qnorm(0.20,
                    mean = mu_U_vec[2] + delta_vec[2],
                    sd = sqrt(sigma2_U_vec[2] + sigma2_eps_2_vec[2])),
              qnorm(0.20,
                    mean = mu_U_vec[3] + + delta_vec[3],
                    sd = sqrt(sigma2_U_vec[3] + sigma2_eps_2_vec[3])))
}

# Number of simulations
M = 1000


# ---------------------------------------------------------------------------
#           TABLE OF CONTENTS
#
# (0) Define Functions
#
#     (0.1) Estimate Mean and Variance from Normal Distribution with a LLOD
#
#     (0.2) Probability of X1 Given X2star
#
#     (0.3) Score Function Wrapper
#
#     (0.4) Beta Estimation & Variance Wrapper
#
#     (0.5) Data Analysis Wrapper
#
#     (0.6) ProbYGivenX
#
# (1) Initialize vectors and matrices
#
# ------------------------------------------------------------------------
#                     RUN SIMULATION
# ------------------------------------------------------------------------
#
# (2) Generate Covariates and the Response for each Cohort
#
# (3) Estimate Parameters with Population 1
#
# (4) Estimate Parameters with Population 2
#
# (5) Estimate Parameters with Populations 1 & 2
#
# ------------------------------------------------------------------------
#                     END OF SIMULATION
# ------------------------------------------------------------------------
#
# (5) Save / Load Results
#
# (6) Summarize Results
# ---------------------------------------------------------------------------

# ... Source C++ functions
sourceCpp('Cpp_Functions_GitHub.cpp')

f = function(x){
  val = 1 + exp(-x)
  val2 = 1 / val
  return(val2)
} 

# ---------------------------------------------------------------------------
# (0) Define Functions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#   (0.1) Estimate Mean and Variance from Normal Distribution with a LLOD
# ---------------------------------------------------------------------------

# Let X ~ N(mu, sigma^2). Suppose we observe W = max(X, a), where a is known.
#   Given a random sample of {W_1, ..., W_n}, the ith unit's contribution to the
#   likelihood function is 
#   phi(W_i; mu, sigma^2)^I(W_i > a) * Phi(W_i; mu, sigma^2)^I(W_i = a)
#   where phi(x;mu,sigma^2) and Phi(x;mu,sigma^2) are the PDF and CDF evaluated at x
#   of the normal distribution with mean mu and variance sigma^2, respectively.

# --- INPUT:
# W: Data of left-censored observations
# mu: mean parameter
# sigma2: variance parameter
# LLOD: Lower limit of detection (Default is "-Inf")

# --- OUTPUT:
# U(delta): value of estimating function evaluated at delta 
#       (Scalar if "Sum == T"; Vector of length n otherwise)

LogLik_LeftCensored_Normal = function(W,
                                      mu, 
                                      sigma2,
                                      LLOD = -Inf){
  
  # Compute log(dnorm(W[i], mu, sigma2)) for each i; we will update this value for the units with a value of the LLOD
  valA = dnorm(W, 
               mean = mu, 
               sd = sqrt(sigma2), 
               log = T) 
  
  # Compute log(pnorm(LLOD, mean = mu, sd = sqrt(sigma2)))
  valB = pnorm(LLOD, 
               mean = mu, 
               sd = sqrt(sigma2), 
               log.p = T)
  
  # Get the indices that report the LLOD
  ind = which(W <= LLOD) # "<=" for numerical stability
  
  # Update "valA[ind]" to "valB"
  if (length(ind) > 0) valA[ind] = valB
  
  # return "sum(valA)"
  val = sum(valA)
  return(val)
}


# Let X ~ N(mu, sigma^2). Suppose we observe W = max(X, a), where a is known.
#   The function below corresponds to the score function

# --- INPUT:
# W: Data of left-censored observations
# mu: mean parameter
# sigma2: variance parameter
# LLOD: Lower limit of detection (Default is "-Inf")
# Sum: Do we want to take the sum of each component of the estimating function?
#       (Logical; default is "T")

# --- OUTPUT:
# U(delta): value of estimating function evaluated at delta 
#       (Scalar if "Sum == T"; Vector of length n otherwise)


ScoreFcn_LeftCensored_Normal = function(W,
                                        mu, 
                                        sigma2,
                                        LLOD = -Inf,
                                        Sum = T){
  
  # For units with "W[i]  > LLOD", their contribution to the score function is...
  valA = (W - mu) / sigma2
  valB = 0.5/sigma2 * ( (W-mu)^2 / sigma2 - 1)
  
  # For units with "W[i] = LLOD", their contribution to the score function is... 
  integrand1 = function(x){
    dnorm(x, mean = mu, sd = sqrt(sigma2)) * (x - mu) / sigma2
  }
  valA.2 = cuhre(f=integrand1, 
                 nComp = 1,
                 lowerLimit = -Inf, 
                 upperLimit = LLOD)$integral
  valA.2 = valA.2 / pnorm(LLOD, mean = mu, sd = sqrt(sigma2))
  
  integrand2 = function(x){
    0.5 * dnorm(x, mean = mu, sd = sqrt(sigma2)) / sigma2 * ( (x - mu)^2/sigma2 - 1 )
  }
  valB.2 = cuhre(f=integrand2, 
                 nComp = 1,
                 lowerLimit = -Inf, 
                 upperLimit = LLOD)$integral
  valB.2 = valB.2 / pnorm(LLOD, mean = mu, sd = sqrt(sigma2))
  
  # Get the indices that report the LLOD
  ind = which(W <= LLOD) # "<=" for numerical stability
  # Update "valA" and "valB"
  if (length(ind) > 0){
    valA[ind] = valA.2
    valB[ind] = valB.2
  }
  val = cbind(valA, valB)
  if (Sum == T) val = as.numeric(colSums(val))
  return(val)
}


# ---------------------------------------------------------------------------
# (0.2) Probability of X1 Given X2star
# ---------------------------------------------------------------------------

# Following the derivation in the Appendix of my notes, we will compute
#   P(x1 | X2 <= LLOD)
# Note that "X1" could be "U" or 
#           "X1" could be "X2"

# --- INPUT:
# x1: distribution evaluated at this value (Scalar)
# LLOD_X2: Lower limit of detection (Scalar)
# mu_X1: Mean of X1 (Scalar)
# mu_X2: Mean of X2 (Scalar)
# sigma2_X1: Variance of X1 (Scalar)
# sigma2_X2: Variance of X2 (Scalar)
# sigma_X1X2: Covariance of X1 and X2 (Scalar)
# MeasurementErrorLogi = T: Indicator for accounting for X being a noisy measurement of U (Logical)

Distribution_X1GivenX2 = function(x1,
                                  LLOD_X2,
                                  mu_X1,
                                  mu_X2,
                                  sigma2_X1,
                                  sigma2_X2,
                                  sigma_X1X2,
                                  MeasurementErrorLogi = T){
  
  # Compute the pdf of X1, where
  # X1 ~ N(mu_X1, sigma2_X1)
  val1 = dnorm(x1,
               mean = mu_X1,
               sd = sqrt(sigma2_X1))
  
  # Get the mean and variance for the second component
  m = (x1 * sigma_X1X2 - mu_X1 * sigma_X1X2 + mu_X2 * sigma2_X1) / sigma2_X1
  s = (sigma2_X1 * sigma2_X2 - sigma_X1X2^2) / sigma2_X1
  
  # Compute the second component (See the Appendix!) 
  # and update it if we want to account for measurement error
  val2 = 1 # Initialize
  if (MeasurementErrorLogi == T) val2 = pnorm(LLOD_X2,
                                              mean = m,
                                              sd = sqrt(s))
  
  # Compute the third component (See the Appendix!)
  val3 = pnorm(LLOD_X2,
               mean = mu_X2,
               sd = sqrt(sigma2_X2))
  val = val1 * val2 / val3
  return(val)
}


# ---------------------------------------------------------------------------
#   (0.3) Score Function Wrapper
# ---------------------------------------------------------------------------

# Write a wrapper function for the score function with a LLOD

# --- INPUT:
# Beta: Regression coefficients (Vector of length 2)
# Y_Data: Response (Vector)
# Xstar_Data: Covariate (Vector)
# Omega_Data: Weights (Vector)
# Strata: Stratification level a unit belongs to (Vector)
# LLOD_vec: Lower limit of detection (Vector)
# mu_U_vec: Mean of U (Vector)
# mu_X_vec: Mean of X (Vector)
# sigma2_U_vec: Variance of U (Vector)
# sigma_UX_vec: Covariance of U and X (Vector)
# sigma2_X_vec: Variance of X (Vector)
# B = 1500: Number of Monte Carlo samples (Scalar)
# MeasurementErrorLogi = T: Indicator for accounting for X being a noisy measurement of U (Logical)
# Sum: Do we want to take the sum of the estimating function?
#       (Logical; default is "T")

# --- OUTPUT:
# U: value of estimating function evaluated at Beta
#       (Vector if "Sum == T"; Matrix otherwise)

# NOTE: For calibrating dataset, the naive approach will have to be specified with "MeasurementErrorLogi == T"
#       but we'll need to be smart with the inputs

ScoreFunction_Wrapper = function(Beta,
                                 Y_Data,
                                 Xstar_Data,
                                 Omega_Data,
                                 Strata,
                                 LLOD_vec,
                                 mu_U_vec,
                                 mu_X_vec,
                                 sigma2_U_vec,
                                 sigma_UX_vec,
                                 sigma2_X_vec,
                                 Setting,
                                 B = 1500,
                                 MeasurementErrorLogi = T,
                                 Sum = T){
  # Set the seed 
  set.seed(seed0)
  
  # Get the unique levels of the stratification variable
  Unique_Strata = unique(Strata)
  
  # Run an lapply statement for each of the levels of the stratification variable
  val_list = lapply(Unique_Strata,
                    function(g){
                      
                      # Initialize a matrix to store contributions to the score function
                      val = matrix(0,
                                   nrow = sum(Strata == g),
                                   ncol = length(Beta))
                      
                      # Get the data corresponding to this group
                      Y_Data.g = Y_Data[which(Strata == g)]
                      Xstar_Data.g = Xstar_Data[which(Strata == g)]
                      Omega_Data.g = Omega_Data[which(Strata == g)]
                      LLOD_vec.g = LLOD_vec[which(Strata == g)]
                      mu_U_vec.g = mu_U_vec[which(Strata == g)]
                      mu_X_vec.g = mu_X_vec[which(Strata == g)]
                      sigma2_U_vec.g = sigma2_U_vec[which(Strata == g)]
                      sigma_UX_vec.g = sigma_UX_vec[which(Strata == g)]
                      sigma2_X_vec.g = sigma2_X_vec[which(Strata == g)]
                      
                      # Get the corresponding parameters for this group
                      LLOD.g =  LLOD_vec.g[1]
                      mu_U.g = mu_U_vec.g[1]
                      mu_X.g = mu_X_vec.g[1]
                      sigma2_U.g = sigma2_U_vec.g[1]
                      sigma_UX.g = sigma_UX_vec.g[1]
                      sigma2_X.g = sigma2_X_vec.g[1]
                      UL.g = ifelse(MeasurementErrorLogi == T, Inf, LLOD.g)
                      
                      # Get the response and covariate for those with Xstar > LLOD
                      Y_vec.g = Y_Data.g[which(Xstar_Data.g > LLOD.g)]
                      Xstar_vec.g = Xstar_Data.g[which(Xstar_Data.g > LLOD.g)]
                      Omega_vec.g = Omega_Data.g[which(Xstar_Data.g > LLOD.g)]
                      
                      # Compute the unit's contribution to the likelihood function for those with Xstar > LLOD
                      valA = ScoreFunction_cpp(beta = Beta,
                                               Y_vec = Y_vec.g,
                                               Xstar_vec = Xstar_vec.g,
                                               Omega_vec = Omega_vec.g,
                                               LLOD = LLOD.g,
                                               mu_U = mu_U.g,
                                               mu_X = mu_X.g,
                                               sigma2_U = sigma2_U.g,
                                               sigma2_X = sigma2_X.g,
                                               B = B, 
                                               MeasurementErrorLogi = MeasurementErrorLogi)
                      
                      # Update "val"
                      val[which(Xstar_Data.g > LLOD.g),] = valA
                      
                      # Compute the unit's contribution to the likelihood function for those with Xstar == LLOD
                      
                      # -----------------------
                      # Y.i == 1
                      # -----------------------
                      
                      # Compute $\int pi(u; \beta) dF(u|X \leq LLOD)$
                      Integrand_Den1 = function(x){
                        val1 = Distribution_X1GivenX2(x1 = x,
                                                      LLOD_X2 = LLOD.g,
                                                      mu_X1 = mu_U.g,
                                                      mu_X2 = mu_X.g,
                                                      sigma2_X1 = sigma2_U.g,
                                                      sigma2_X2 = sigma2_X.g,
                                                      sigma_X1X2 = sigma_UX.g,
                                                      MeasurementErrorLogi = MeasurementErrorLogi) # this is dF(u|X <= LLOD)
                        prob = pi_U_beta_cpp(beta = Beta,
                                             U_i = x)
                        val = prob * val1
                        return(val)
                      }
                      L1 = cuhre(f=Integrand_Den1, 
                                 nComp = 1,
                                 lowerLimit = -Inf, 
                                 upperLimit = UL.g )$integral
                      
                      # -----------------------
                      # Y.i == 0
                      # -----------------------
                      
                      # Compute $\int ( 1-pi(u; \beta) ) dF(u|X \leq LLOD)$
                      Integrand_Den0 = function(x){
                        val1 =  Distribution_X1GivenX2(x1 = x,
                                                       LLOD_X2 = LLOD.g,
                                                       mu_X1 = mu_U.g,
                                                       mu_X2 = mu_X.g,
                                                       sigma2_X1 = sigma2_U.g,
                                                       sigma2_X2 = sigma2_X.g,
                                                       sigma_X1X2 = sigma_UX.g,
                                                       MeasurementErrorLogi = MeasurementErrorLogi) # this is dF(u|X <= LLOD)
                        prob = pi_U_beta_cpp(beta = Beta,
                                             U_i = x)
                        val = (1-prob) * val1
                        return(val)
                      }
                      L0 = cuhre(f=Integrand_Den0, 
                                 nComp = 1,
                                 lowerLimit = -Inf, 
                                 upperLimit = UL.g )$integral
                      
                      # --------------------------------
                      # Compute each unit's contribution
                      #   to the score function
                      # --------------------------------
                      
                      # Compute $\int \pi(u;\beta) (1-\pi(u;\beta)) dF(u|X \leq LLOD)$
                      Integrand_Num1_beta0 = function(x){
                        val1 =  Distribution_X1GivenX2(x1 = x,
                                                       LLOD_X2 = LLOD.g,
                                                       mu_X1 = mu_U.g,
                                                       mu_X2 = mu_X.g,
                                                       sigma2_X1 = sigma2_U.g,
                                                       sigma2_X2 = sigma2_X.g,
                                                       sigma_X1X2 = sigma_UX.g,
                                                       MeasurementErrorLogi = MeasurementErrorLogi) # this is dF(u|X <= LLOD)
                        prob = pi_U_beta_cpp(beta = Beta,
                                             U_i = x)
                        val2 = prob * (1 - prob) * val1
                        return(val2)
                      }
                      Num1_beta0 = cuhre(f=Integrand_Num1_beta0, 
                                         nComp = 1,
                                         lowerLimit = -Inf, 
                                         upperLimit = UL.g )$integral
                      
                      # Compute $\int \pi(u;\beta) (1-\pi(u;\beta)) * u dF(u|X \leq LLOD)$
                      Integrand_Num1_beta1 = function(x){
                        val1 =  Distribution_X1GivenX2(x1 = x,
                                                       LLOD_X2 = LLOD.g,
                                                       mu_X1 = mu_U.g,
                                                       mu_X2 = mu_X.g,
                                                       sigma2_X1 = sigma2_U.g,
                                                       sigma2_X2 = sigma2_X.g,
                                                       sigma_X1X2 = sigma_UX.g,
                                                       MeasurementErrorLogi = MeasurementErrorLogi) # this is dF(u|X <= LLOD)
                        prob = pi_U_beta_cpp(beta = Beta,
                                             U_i = x)
                        val2 = prob * (1 - prob) * val1 * x
                        return(val2)
                      }
                      Num1_beta1 = cuhre(f=Integrand_Num1_beta1, 
                                         nComp = 1,
                                         lowerLimit = -Inf, 
                                         upperLimit = UL.g )$integral
                      
                      # Get the weight
                      IPW.g = 1 # Initialize
                      ind_g0 = which(Y_Data.g == 0)
                      if (length(ind_g0) > 0) IPW.g = Omega_Data.g[ind_g0[1]]
                      
                      # Determine which units have 
                      # "Xstar_Data.g <= LLOD.g" & "Y_vec.g == 0" 
                      # and update their contribution to the score function
                      ind0 = which(Xstar_Data.g <= LLOD.g & Y_Data.g == 0)
                      val_Y0_beta0 = -IPW.g * Num1_beta0 / L0 
                      val_Y0_beta1 = -IPW.g * Num1_beta1 / L0
                      val[ind0,] = c(val_Y0_beta0, val_Y0_beta1)
                      
                      # Determine which units have 
                      # "Xstar_Data.g <= LLOD.g" & "Y_vec.g == 1" 
                      # and update their contribution to the score function
                      ind1 = which(Xstar_Data.g <= LLOD.g & Y_Data.g == 1)
                      val_Y1_beta0 = Num1_beta0 / L1
                      val_Y1_beta1 = Num1_beta1 / L1
                      val[ind1,] = c(val_Y1_beta0, val_Y1_beta1)
                      
                      return(val)
                    })
  
  # Merge the lists together
  val = do.call(rbind,
                val_list)
  
  # Take the column sum if we want
  if (Sum == T) val = as.numeric(colSums(val))
  
  # return "val"
  return(val)
}



# ---------------------------------------------------------------------------
#   (0.4) Beta Estimation & Variance Wrapper
# ---------------------------------------------------------------------------

# Write a wrapper function to estimate beta (Commenting out estimating the variance)

# --- INPUT:
# Data: Dataframe that has the relevant information
# Setting: An integer that specifies different settings:
#     Setting == 1: No calibration taking place, only one population
#     Setting == 2: Yes Calibration taking place, only one population
#     Setting == 3: Yes Calibration taking place, with two populations
# alpha: An estimate of the nuisance parameters
# mu_U_vec_EF: Mean of U (Vector)
# mu_X_vec_EF: Mean of X (Vector)
# sigma2_U_vec_EF: Variance of U (Vector)
# sigma_UX_vec_EF: Covariance of U and X (Vector)
# sigma2_X_vec_EF: Variance of X (Vector)
# B = 1500: Number of Monte Carlo samples (Scalar)
# MeasurementErrorLogi = T: Indicator for accounting for X being a noisy measurement of U (Logical)
# seed0 = 13543645: An integer to specify the seed for the MCMLE
# start: Starting value for the estimation of beta
# index_beta: Indices corresponding to "Beta" in "Theta"
# Variance_Estimation = F: Indicator if we want to estimate the variance of "Theta"
#
# --- OUTPUT: A list containing
# mod: Model output of "nleqslv" for estimating beta


BetaEstimation_Wrapper = function(Data,
                                  Setting,
                                  alpha,
                                  mu_U_vec_EF,
                                  mu_X_vec_EF,
                                  sigma2_U_vec_EF,
                                  sigma_UX_vec_EF,
                                  sigma2_X_vec_EF,
                                  B = 1500,
                                  MeasurementErrorLogi,
                                  seed0 = 13543645,
                                  start,
                                  index_beta,
                                  Variance_Estimation = F){
  
  # Define the estimating function
  EF_func = function(Beta){
    set.seed(seed0)
    val = ScoreFunction_Wrapper(Beta = Beta,
                                Y_Data = Data$Y,
                                Xstar_Data = Data$Xstar,
                                Omega_Data = 1/Data$Prob,
                                Strata = Data$Strata,
                                LLOD_vec = Data$LLOD,
                                mu_U_vec = mu_U_vec_EF, 
                                mu_X_vec = mu_X_vec_EF, 
                                sigma2_U_vec = sigma2_U_vec_EF,
                                sigma_UX_vec = sigma_UX_vec_EF,
                                sigma2_X_vec = sigma2_X_vec_EF,
                                Setting = Setting,
                                B = B, 
                                MeasurementErrorLogi = MeasurementErrorLogi,
                                Sum = T)
    return(val)
  }
  
  # Test the estimating function
  # EF_func(Beta = start)
  
  # Solve the estimating equation
  mod = nleqslv(x = start,
                fn = EF_func,
                jacobian = T,
                control = list(trace = 0,
                               maxit = 500))
  
  # mod$x; start
  
  
  # ----------------------
  # Estimate the Variance
  # ----------------------
  
  A_mat = NA # Initialize
  A_sanitycheck = NA
  B_mat = NA
  B_sanitycheck = NA
  AB_mat = NA
  
  if (Variance_Estimation == T){
    
    # Specify the function to obtain "A_mat"
    A_func = function(param){
      
      # Set the seed and extract "Beta"
      set.seed(seed0)
      Beta = param[index_beta]
      
      # Go into cases for "Setting"
      if (Setting == 1){
        
        # Compute the estimating function for the nuisance parameter
        val_alpha = ScoreFcn_LeftCensored_Normal(W = Data$Xstar,
                                                 mu = param[1],
                                                 sigma2 = param[2],
                                                 LLOD = Data$LLOD[1])
        
        # Define vectors for parameters
        mu_U_vec = rep(param[1], dim(Data)[1])
        mu_X_vec = rep(param[1], dim(Data)[1])
        if (MeasurementErrorLogi == F){
          sigma2_U_vec = rep(param[2], dim(Data)[1])
        }
        if (MeasurementErrorLogi == T){
          sigma2_U_vec = sigma2_U_vec_EF
        }
        sigma_UX_vec = sigma_UX_vec_EF
        sigma2_X_vec = rep(param[2], dim(Data)[1])
      }
      if (Setting == 2){
        
        # Compute the estimating function for the nuisance parameter
        val_alpha = ScoreFcn_LeftCensored_Normal(W = Data$Xstar,
                                                 mu = param[1],
                                                 sigma2 = param[2],
                                                 LLOD = Data$LLOD[1])
        
        # Define vectors for parameters
        mu_U_vec = mu_U_vec_EF 
        mu_X_vec = rep(param[1], dim(Data)[1])
        sigma2_U_vec = sigma2_U_vec_EF
        sigma_UX_vec = sigma_UX_vec_EF 
        sigma2_X_vec = rep(param[2], dim(Data)[1])
      }
      
      # Compute the estimating function for the regression parameter
      val_beta = ScoreFunction_Wrapper(Beta = Beta,
                                       Y_Data = Data$Y,
                                       Xstar_Data = Data$Xstar,
                                       Omega_Data = 1/Data$Prob,
                                       Strata = Data$Strata,
                                       LLOD_vec = Data$LLOD,
                                       mu_U_vec = mu_U_vec,
                                       mu_X_vec = mu_X_vec,
                                       sigma2_U_vec = sigma2_U_vec,
                                       sigma_UX_vec = sigma_UX_vec,
                                       sigma2_X_vec = sigma2_X_vec,
                                       B = B,
                                       MeasurementErrorLogi = MeasurementErrorLogi,
                                       Sum = T)
      
      val = as.numeric(c(val_alpha, val_beta))
      return(val)
    }
    
    # sanity check: Should essentially be 0
    A_sanitycheck = A_func(param = c(alpha,
                                     mod$x))
    
    A_mat = jacobian(A_func,
                     x = c(alpha,
                           mod$x))
    
    B_func = function(param){
      
      # Set the seed and extract "Beta"
      set.seed(seed0)
      Beta = param[index_beta]
      
      # Go into cases for "Setting"
      if (Setting == 1){
        
        # Compute the estimating function for the nuisance parameter
        val_alpha = ScoreFcn_LeftCensored_Normal(W = Data$Xstar,
                                                 mu = param[1],
                                                 sigma2 = param[2],
                                                 LLOD = Data$LLOD[1],
                                                 Sum = F)
        # Define vectors for parameters
        mu_U_vec = rep(param[1], dim(Data)[1])
        mu_X_vec = rep(param[1], dim(Data)[1])
        if (MeasurementErrorLogi == F){
          sigma2_U_vec = rep(param[2], dim(Data)[1])
        }
        if (MeasurementErrorLogi == T){
          sigma2_U_vec = sigma2_U_vec_EF
        }
        sigma_UX_vec = sigma_UX_vec_EF
        sigma2_X_vec = rep(param[2], dim(Data)[1])
      }
      if (Setting == 2){
        
        # Compute the estimating function for the nuisance parameter
        val_alpha = ScoreFcn_LeftCensored_Normal(W = Data$Xstar,
                                                 mu = param[1],
                                                 sigma2 = param[2],
                                                 LLOD = Data$LLOD[1],
                                                 Sum = F)
        
        # Define vectors for parameters
        mu_U_vec = mu_U_vec_EF 
        mu_X_vec = rep(param[1], dim(Data)[1])
        sigma2_U_vec = sigma2_U_vec_EF 
        sigma_UX_vec = sigma_UX_vec_EF 
        sigma2_X_vec = rep(param[2], dim(Data)[1])
      }
      
      # Compute the estimating function for the regression parameter
      val_beta = ScoreFunction_Wrapper(Beta = Beta,
                                       Y_Data = Data$Y,
                                       Xstar_Data = Data$Xstar,
                                       Omega_Data = 1/Data$Prob,
                                       Strata = Data$Strata,
                                       LLOD_vec = Data$LLOD,
                                       mu_U_vec = mu_U_vec,
                                       mu_X_vec = mu_X_vec,
                                       sigma2_U_vec = sigma2_U_vec,
                                       sigma_UX_vec = sigma_UX_vec,
                                       sigma2_X_vec = sigma2_X_vec,
                                       B = B, 
                                       MeasurementErrorLogi = MeasurementErrorLogi,
                                       Sum = F)
      
      val = cbind(val_alpha, val_beta)
      return(val)
    }
    
    # sanity check: Should essentially be 0
    B_sanitycheck = colSums(B_func(param = c(alpha,
                                             mod$x)))
    
    B_mat = crossprod(B_func(c(alpha,
                               mod$x)))
    
    # Compute the sandwich variance estimator
    AB_mat = MASS::ginv(A_mat) %*% B_mat %*% t(MASS::ginv(A_mat))
    #Var_mat = diag(AB_mat)[index_beta]
  }
  
  # Return everything we want in a list
  val = list(mod = mod,
             A = A_mat,
             A_sanitycheck = A_sanitycheck,
             B = B_mat,
             B_sanitycheck = B_sanitycheck,
             AInv_B_AInv = AB_mat)
  return(val)
}


# -----------------------------------------------------------------------------
#   (0.5) Data Analysis Wrapper
# -----------------------------------------------------------------------------


# Write a wrapper function to fit all of the models

# --- INPUT:
# Data: Dataframe that has the relevant information
# Setting: An integer that specifies different settings:
#     Setting == 1: No calibration taking place, only one population (useful for population 1 only)
#     Setting == 2: Yes Calibration taking place, only one population (useful for population 2 only)
#     Setting == 3: Yes Calibration taking place, with two populations (useful for population 1 & 2)
# Var_MeasurementError: Estimated variance of measurement error
# Delta_Zero = F: Do we want delta to be zero, or do we want to estimate it?
# B = 1500: Number of Monte Carlo samples (Scalar)
# seed0 = 13543645: An integer to specify the seed for the MCMLE
#
# --- OUTPUT: A list containing
# mod_alpha: Model output for estimating "alpha"
# mod_beta_NoME: Model output for estimating "beta" by ignoring ME
# mod_beta = Beta_Estimation: Model output for proposed approach


Data_Analysis_Wrapper = function(Data,
                                 Setting,
                                 Var_MeasurementError,
                                 PairedSample_Estimates = c(NA,NA,NA,NA,NA),
                                 Delta_Zero = F,
                                 B = 1500,
                                 seed0 = 13543645){
  
  if (Setting %in% c(1,2)){
    # ----------------------------
    # Estimate the mean 
    #   and variance of log(ID50)
    # ----------------------------
    
    # Specify the corresponding score function
    EF_alpha = function(par){
      ScoreFcn_LeftCensored_Normal(W = Data$Xstar,
                                   mu = par[1],
                                   sigma2 = par[2],
                                   LLOD = Data$LLOD[1])
    }
    
    # Test out the function
    EF_alpha(par = c(mean(Data$Xstar),
                     var(Data$Xstar)))
    
    # Solve the estimating equation
    MLEsEF_alpha = nleqslv(x = c(mean(Data$Xstar),
                                 var(Data$Xstar)),
                           fn = EF_alpha,
                           jacobian = T,
                           control = list(trace = 0,
                                          maxit = 500))
    
    # MLEsEF_alpha$x; c(mean(Data$Xstar),
    #                   var(Data$Xstar)*(length(Data$Xstar)-1)/length(Data$Xstar))
  }
  if (Setting == 3){
    
    # Omit estimating alpha, we'll use the estimates from before
    MLEsEF_alpha = NA
  }
  
  
  # ---------------------------------------------------------------------
  # Depending on "Setting", define 
  #   "mu_U_vec", "mu_X_vec", "sigma2_U_vec", "sigma_UX_vec", "sigma2_X_vec"
  #   and "MeasurementErrorLogi" <- NO TYPO!
  # ---------------------------------------------------------------------
  if (Setting == 1){
    
    # No measurement error
    mu_U_vec_NoME = rep(MLEsEF_alpha$x[1],
                        dim(Data)[1]) 
    mu_X_vec_NoME = mu_U_vec_NoME
    sigma2_U_vec_NoME = rep(MLEsEF_alpha$x[2],
                            dim(Data)[1])
    sigma_UX_vec_NoME = rep(1,
                            dim(Data)[1])
    sigma2_X_vec_NoME = sigma2_U_vec_NoME 
    
    # Yes measurement error
    mu_U_vec_YesME = rep(MLEsEF_alpha$x[1],
                         dim(Data)[1]) 
    mu_X_vec_YesME = mu_U_vec_YesME
    sigma2_U_vec_YesME = rep(MLEsEF_alpha$x[2] - Var_MeasurementError,
                             dim(Data)[1])
    sigma_UX_vec_YesME = sigma2_U_vec_YesME
    sigma2_X_vec_YesME = rep(MLEsEF_alpha$x[2],
                             dim(Data)[1])
    MeasurementErrorLogi = F
    index_beta = c(3,4)
    
    alpha = MLEsEF_alpha$x
    Bridged_Mean = NA
  }
  if (Setting == 2){
    
    # Compute the bridged mean
    delta = ifelse(Delta_Zero == T, 
                   0,
                   PairedSample_Estimates[2] - PairedSample_Estimates[1])
    Bridged_Mean = MLEsEF_alpha$x[1] - delta
    
    # No measurement error
    mu_U_vec_NoME = rep(Bridged_Mean,
                        dim(Data)[1]) 
    mu_X_vec_NoME = rep(MLEsEF_alpha$x[1],
                        dim(Data)[1]) 
    
    Ratio_PairedSample_Variance = PairedSample_Estimates[5]/PairedSample_Estimates[3]
    if (Ratio_PairedSample_Variance > 1) Ratio_PairedSample_Variance = 1/Ratio_PairedSample_Variance
    Corr_PairedSample = PairedSample_Estimates[4] / sqrt(PairedSample_Estimates[3] * PairedSample_Estimates[5])
    
    sigma2_U_Est = MLEsEF_alpha$x[2] * Ratio_PairedSample_Variance
    sigma_UX_Est = Corr_PairedSample * sqrt(sigma2_U_Est) * sqrt(MLEsEF_alpha$x[2])
    
    sigma2_U_vec_NoME = rep(sigma2_U_Est, 
                            dim(Data)[1])
    sigma_UX_vec_NoME = rep(sigma_UX_Est,
                            dim(Data)[1])
    sigma2_X_vec_NoME = rep(MLEsEF_alpha$x[2],
                            dim(Data)[1])
    MeasurementErrorLogi = T
    
    # Yes measurement error
    mu_U_vec_YesME = rep(Bridged_Mean,
                         dim(Data)[1]) 
    mu_X_vec_YesME = rep(MLEsEF_alpha$x[1],
                         dim(Data)[1]) 
    sigma2_U_vec_YesME = rep(MLEsEF_alpha$x[2] - Var_MeasurementError,
                             dim(Data)[1])
    sigma_UX_vec_YesME = sigma2_U_vec_YesME
    sigma2_X_vec_YesME = rep(MLEsEF_alpha$x[2],
                             dim(Data)[1])
    index_beta = c(3,4)
    #Setting_Adjust = Setting + 0.5
    
    alpha = MLEsEF_alpha$x
    
  }
  if (Setting == 3){
    
    n1_use = dim(D1)[1]
    n2_use = dim(D2)[1]
    
    # We will only estimate parameters with measurement error
    # Obtain "mu_U_vec", "mu_X_vec", "sigma2_U_vec", "sigma2_X_vec"
    if (Delta_Zero == F){
      mu_U_vec_YesME = c(rep(Data_Analysis_1$mod_alpha$x[1],n1_use),
                         rep(Data_Analysis_2_BridgeYes$Bridged_Mean, n2_use)) 
    }
    if (Delta_Zero == T){
      mu_U_vec_YesME = c(rep(Data_Analysis_1$mod_alpha$x[1],n1_use),
                         rep(Data_Analysis_2_BridgeNo$Bridged_Mean, n2_use)) 
    }
    mu_X_vec_YesME = c(rep(Data_Analysis_1$mod_alpha$x[1],n1_use),
                       rep(Data_Analysis_2_BridgeYes$mod_alpha$x[1], n2_use))
    sigma2_U_vec_YesME = c(rep(Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2], n1_use),
                           rep(Data_Analysis_2_BridgeYes$mod_alpha$x[2] - sigma2_eps_2_vec[3], n2_use))
    sigma_UX_vec_YesME = sigma2_U_vec_YesME
    sigma2_X_vec_YesME = c(rep(Data_Analysis_1$mod_alpha$x[2], n1_use),
                           rep(Data_Analysis_2_BridgeYes$mod_alpha$x[2], n2_use))
    alpha = c(0,0) # This sepcification is irrelevant!
    index_beta = c(3,4)
    Bridged_Mean = NA
    
  }
  
  # Proceed to estimate beta!
  
  if (Setting == 1){
    # ----------------------------
    # Estimate beta without
    #   measurement error
    # ----------------------------
    
    Beta_Estimation_NoME = BetaEstimation_Wrapper(Data = Data,
                                                  Setting = Setting,
                                                  alpha = alpha,
                                                  mu_U_vec_EF = mu_U_vec_NoME,
                                                  mu_X_vec_EF = mu_X_vec_NoME,
                                                  sigma2_U_vec_EF = sigma2_U_vec_NoME,
                                                  sigma_UX_vec_EF = sigma_UX_vec_NoME,
                                                  sigma2_X_vec_EF = sigma2_X_vec_NoME,
                                                  B = B, # THIS SPECIFICATION MAY BE IRRELEVANT!
                                                  MeasurementErrorLogi =  MeasurementErrorLogi,
                                                  seed0 = seed0,
                                                  start = Beta,
                                                  index_beta = index_beta,
                                                  Variance_Estimation = F)
    Beta_Estimation_NoME$mod$x
    #sqrt(diag(Beta_Estimation_NoME$AInv_B_AInv))[-c(1,2)]
    
    
    # ----------------------------
    # Estimate beta with
    #  EVERYTHING ACCOUNTED FOR
    # ----------------------------
    
    Beta_Estimation = BetaEstimation_Wrapper(Data = Data,
                                             Setting = Setting,
                                             alpha = alpha,
                                             mu_U_vec_EF = mu_U_vec_YesME,
                                             mu_X_vec_EF = mu_X_vec_YesME,
                                             sigma2_U_vec_EF = sigma2_U_vec_YesME,
                                             sigma_UX_vec_EF = sigma_UX_vec_YesME,
                                             sigma2_X_vec_EF = sigma2_X_vec_YesME,
                                             B = B,
                                             MeasurementErrorLogi =  T,
                                             seed0 = seed0,
                                             start = Beta,
                                             index_beta = index_beta,
                                             Variance_Estimation = F)
    Beta_Estimation$mod$x
    #sqrt(diag(Beta_Estimation$AInv_B_AInv))[-c(1,2)]
  }
  if (Setting == 2){
    
    # Proceed to estimate beta (omit the naive estimate)
    Beta_Estimation_NoME = NA
    
    Beta_Estimation = BetaEstimation_Wrapper(Data = Data,
                                             Setting = Setting,
                                             alpha = alpha,
                                             mu_U_vec_EF = mu_U_vec_YesME,
                                             mu_X_vec_EF = mu_X_vec_YesME,
                                             sigma2_U_vec_EF = sigma2_U_vec_YesME,
                                             sigma_UX_vec_EF = sigma_UX_vec_YesME,
                                             sigma2_X_vec_EF = sigma2_X_vec_YesME,
                                             B = B,
                                             MeasurementErrorLogi =  T,
                                             seed0 = seed0,
                                             start = Beta,
                                             index_beta = index_beta,
                                             Variance_Estimation = F)
    Beta_Estimation$mod$x
    #sqrt(diag(Beta_Estimation$AInv_B_AInv))[-c(1,2)]
  }
  if (Setting == 3){
    
    # Proceed to estimate beta (omit the naive estimate)
    Beta_Estimation_NoME = NA
    
    Beta_Estimation = BetaEstimation_Wrapper(Data = Data,
                                             Setting = Setting,
                                             alpha = alpha,
                                             mu_U_vec_EF = mu_U_vec_YesME,
                                             mu_X_vec_EF = mu_X_vec_YesME,
                                             sigma2_U_vec_EF = sigma2_U_vec_YesME,
                                             sigma_UX_vec_EF = sigma_UX_vec_YesME,
                                             sigma2_X_vec_EF = sigma2_X_vec_YesME,
                                             B = B,
                                             MeasurementErrorLogi =  T,
                                             seed0 = seed0,
                                             start = Beta,
                                             index_beta = index_beta,
                                             Variance_Estimation = F)
  }
  
  return(list(Bridged_Mean = Bridged_Mean,
              mod_alpha = MLEsEF_alpha,
              mod_beta_NoME = Beta_Estimation_NoME,
              mod_beta = Beta_Estimation))
}

# ---------------------------------------------------------------------------
#   (0.6) ProbYGivenX
# ---------------------------------------------------------------------------

# Write a function that applies a Monte Carlo approximation
#   to $\pi(X;\beta) - \int \pi(u;\beta) dF(u|X)$
#   and obtain $\pi(X;\beta)$ over a grid of X-values

# --- INPUT:
# betahat: Estimated regression coefficients (Vector of length 2)
# X: Noisy covariate (Scalar)
# mu_U: Mean of true covariate (Scalar)
# mu_X: Mean of noisy covariate (Scalar)
# sigma2_U: Variance of true covariate (Scalar)
# sigma2_X: Variance of noisy covariate (Scalar)
# B: Number of iterations (Scalar, default is 50000)
# seed0: Seed to generate reproducible results (Scalar, default is 13543645)

# --- OUTPUT:
# [Y |x; betahat] (Vector)


ProbYGivenX = function(betahat,
                       X,
                       mu_U,
                       mu_X,
                       sigma2_U,
                       sigma2_X,
                       B = 10000,
                       seed0 = 13543645){
  
  # set the seed to obtain reproducible results
  if (!is.na(seed0)) set.seed(seed0)
  
  # Specify the mean and variance of [U | X]
  mu_U_X = mu_U + (sigma2_U / sigma2_X) * (X - mu_X)
  sigma2_U_X = sigma2_U - (sigma2_U^2 / sigma2_X)
  
  # Generate B normal random variables with mean "mu_U_X" and variance "sigma2_U_X"
  U.b = rnorm(B, mean = mu_U_X, sd = sqrt(sigma2_U_X))
  
  # Compute P(Y = 1 | U.b ; betahat)
  Prob = 1/(1+exp(-betahat[1] - betahat[2]*U.b))
  
  # Compute the mean and return it
  val = mean(Prob)
  
  return(val)
}


# ---------------------------------------------------------------------------
# (0.8) Regression Calibration Prediction Function
# ---------------------------------------------------------------------------

# Write a function that computes E(U | X*)

# --- INPUT:
# X: Covariates (Vector)
# LLOD: Lower Limit of Detection (Scalar)
# mu_U: Mean of true covariate (Scalar)
# mu_X: Mean of noisy covariate (Scalar)
# sigma2_U: Variance of true covariate (Scalar)
# sigma2_X: Variance of noisy covariate (Scalar)

# --- OUTPUT:
# E(U | X) for each unit (Vector)


RC_Predict_Function = function(X,
                               LLOD,
                               mu_U,
                               mu_X,
                               sigma2_U,
                               sigma2_X){
  
  # Extract the number of units
  n = length(X)
  
  # Proceed if there is a LLOD
  ExpectedValue_True = NA # Initialize
  if (LLOD > -Inf){
    
    # Compute E(U | X <= LLOD)
    Integrand = function(x){
      val1 = Distribution_X1GivenX2(x1 = x,
                                    LLOD_X2 = LLOD,
                                    mu_X1 = mu_U,
                                    mu_X2 = mu_X,
                                    sigma2_X1 = sigma2_U,
                                    sigma2_X2 = sigma2_X,
                                    sigma_X1X2 = sigma2_U,
                                    MeasurementErrorLogi = T)
      val = x * val1
      return(val)
    }
    ExpectedValue_True = cuhre(f = Integrand, 
                               nComp = 1,
                               lowerLimit = -Inf, 
                               upperLimit = Inf)
  }
  
  # Run an sapply statement for each i to compute E(U | X_i)
  val = sapply(1:n, function(i){
    
    X.i = X[i]
    
    val.i = mu_U + sigma2_U / sigma2_X * (X.i - mu_X) # Initialize; correct if "X.i > LLOD"
    
    # Update it if "X.i == LLOD"
    if (X.i <= LLOD) val.i = ExpectedValue_True$integral
    
    return(val.i)
  })
  
  return(val)
}


# ---------------------------------------------------------------------------
# (1) Initialize vectors and matrices
# ---------------------------------------------------------------------------

# Mean of Y1_D1
Mean_Y_D1_vec = rep(NA, M)
Mean_Y_D2_vec = rep(NA, M)

# Proportion LLOD
Prop_LLOD_X1_D0_vec = rep(NA, M)
Prop_LLOD_X2_D0_vec = rep(NA, M)
Prop_LLOD_X1_D1_vec = rep(NA, M)
Prop_LLOD_X2_D1_vec = rep(NA, M)
Prop_LLOD_X1_D2_vec = rep(NA, M)
Prop_LLOD_X2_D2_vec = rep(NA, M)

# Estimates of alpha and bridged mean
Est_alpha_D1_mat = matrix(NA, nrow = M, ncol = length(Beta))
Est_alpha_D2_mat = matrix(NA, nrow = M, ncol = length(Beta))
Est_BridgedMean_vec = rep(NA, M)

# Estimates of betahat
Betahat_TrueCovariate_D1_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_Naive_D1_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_RC_D1_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_WMCMLE_D1_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_WMCMLE_D2_BridgeYes_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_WMCMLE_D2_BridgeNo_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_RC_D12_BridgeYes_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_WMCMLE_D12_BridgeYes_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_RC_D12_BridgeNo_mat = matrix(NA, nrow = M, ncol = length(Beta))
Betahat_WMCMLE_D12_BridgeNo_mat = matrix(NA, nrow = M, ncol = length(Beta))

# Prob Y Given X
X1_grid = seq(0,6, length.out = 61)

# Estimates of pi(x; alpha, beta) for x in {x_1,...,x_{10}}
ProbYGivenX_TrueCovariate_D1_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_Naive_D1_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_RC_D1_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_WMCMLE_D1_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_RC_D12_BrigeYes_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_WMCMLE_D12_BridgeYes_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_RC_D12_BrigeNo_mat = matrix(NA, nrow = M, ncol = length(X1_grid))
ProbYGivenX_WMCMLE_D12_BridgeNo_mat = matrix(NA, nrow = M, ncol = length(X1_grid))

# Compute the true probability values for X = X1_grid
Prob_True = sapply(1:length(X1_grid), function(i){
  x = X1_grid[i]
  val = ProbYGivenX(betahat = Beta,
                    X = x,
                    mu_U = mu_U_vec[2],
                    mu_X = mu_U_vec[2],
                    sigma2_U = sigma2_U_vec[2],
                    sigma2_X = sigma2_U_vec[2] + sigma2_eps_1_vec[2],
                    B = 10000)
  return(val)
})

# ---------------------------------------------------------------------------
#                     RUN SIMULATION
# ---------------------------------------------------------------------------

for (m in 1:M){
  
  seed_initial = 5347547 + m
  set.seed(seed_initial)
  
  print(paste("Working on Simulation", m, "out of", M))
  
  # ---------------------------------------------------------------------------
  # (2) Generate Covariates and the Response for each Cohort
  # ---------------------------------------------------------------------------
  
  # -----------------------------------------
  #       Data "D0"
  # -----------------------------------------
  
  # Generate true biomarker value
  U_D0 = rnorm(N0, 
               mean = mu_U_vec[1],
               sd = sqrt(sigma2_U_vec[1]))
  
  # Generate measurement error
  MeasurementError_D0 = MASS::mvrnorm(n = N0, 
                                      mu = c(0, 0),
                                      Sigma = matrix(c(sigma2_eps_1_vec[1], 0,
                                                       0, sigma2_eps_2_vec[1]),
                                                     nrow = 2, ncol = length(Beta)) )
  
  # Generate noisy covariate values
  X1_D0 = U_D0 + MeasurementError_D0[,1]
  X2_D0 = delta_vec[1] + U_D0 + MeasurementError_D0[,2]
  
  Prop_LLOD_X1_D0_vec[m] = sum(X1_D0 <= LLOD_X1[1]) / length(X1_D0)
  Prop_LLOD_X2_D0_vec[m] = sum(X2_D0 <= LLOD_X2[1]) / length(X2_D0)
  
  # Compute max(Xk_D0, LLOD)
  X1star_D0 = pmax(X1_D0, LLOD_X1[1])
  X2star_D0 = pmax(X2_D0, LLOD_X2[1])
  
  
  # Generate "D0"
  D0 = data.frame(X1star = X1star_D0,
                  X2star = X2star_D0,
                  LLOD_X1 = LLOD_X1[1],
                  LLOD_X2 = LLOD_X2[1])
  
  
  # -----------------------------------------
  #       Data "D1"
  # -----------------------------------------
  
  # Generate true biomarker value
  U_D1 = rnorm(N1, 
               mean = mu_U_vec[2],
               sd = sqrt(sigma2_U_vec[2]))
  
  # Generate measurement error
  MeasurementError_D1 = MASS::mvrnorm(n = N1, 
                                      mu = c(0, 0),
                                      Sigma = matrix(c(sigma2_eps_1_vec[2], 0,
                                                       0, sigma2_eps_2_vec[2]),
                                                     nrow = 2, ncol = length(Beta)) )
  
  # Generate noisy covariate values
  X1_D1 = U_D1 + MeasurementError_D1[,1]
  X2_D1 = delta_vec[2] + U_D1 + MeasurementError_D1[,2]
  
  Prop_LLOD_X1_D1_vec[m] = sum(X1_D1 <= LLOD_X1[2]) / length(X1_D1)
  Prop_LLOD_X2_D1_vec[m] = sum(X2_D1 <= LLOD_X2[2]) / length(X2_D1)
  
  # Compute max(Xk_D1, LLOD)
  X1star_D1 = pmax(X1_D1, LLOD_X1[2])
  X2star_D1 = pmax(X2_D1, LLOD_X2[2])
  
  # Generate Y
  Y_D1 = rbinom(N1, size = 1, prob = f(Beta[1] + Beta[2] * U_D1))
  # mean(Y_D1)*100; var(Y_D1)*100
  
  # Proceed to sample units into the subcohort
  index_SubCohort_D1 = sample(1:N1,
                              size = n1)
  
  # Include all cases as well
  index_SubCohort_D1 = c(index_SubCohort_D1,
                         which(Y_D1 == 1))
  index_SubCohort_D1 = unique(index_SubCohort_D1)
  
  # Generate "D1"
  D1 = data.frame(Y = Y_D1,
                  Xstar = X1star_D1,
                  U = U_D1,
                  Strata = 1,
                  Prob = n1 / N1,
                  LLOD = LLOD_X1[2])
  D1 = D1[index_SubCohort_D1,]
  
  # Update the probability if Y == 1
  D1[which(D1$Y == 1),]$Prob = 1
  
  Mean_Y_D1_vec[m] = mean(D1$Y)
  
  
  # -----------------------------------------
  #       Data "D2"
  # -----------------------------------------
  
  # Generate true biomarker value
  U_D2 = rnorm(N2, 
               mean = mu_U_vec[3],
               sd = sqrt(sigma2_U_vec[3]))
  
  # Generate measurement error
  MeasurementError_D2 = MASS::mvrnorm(n = N2, 
                                      mu = c(0, 0),
                                      Sigma = matrix(c(sigma2_eps_1_vec[3], 0,
                                                       0, sigma2_eps_2_vec[3]),
                                                     nrow = 2, ncol = length(Beta)) )
  
  # Generate noisy covariate values
  X1_D2 = U_D2 + MeasurementError_D2[,1]
  X2_D2 = delta_vec[3] + U_D2 + MeasurementError_D2[,2]
  
  Prop_LLOD_X1_D2_vec[m] = sum(X1_D2 <= LLOD_X1[3]) / length(X1_D2)
  Prop_LLOD_X2_D2_vec[m] = sum(X2_D2 <= LLOD_X2[3]) / length(X2_D2)
  
  # Compute max(Xk_D2, LLOD)
  X1star_D2 = pmax(X1_D2, LLOD_X1[3])
  X2star_D2 = pmax(X2_D2, LLOD_X2[3])
 
  # Generate Y
  Y_D2 = rbinom(N2, size = 1, prob = f(Beta[1] + Beta[2] * U_D2))
  
  # Proceed to sample units into the subcohort
  index_SubCohort_D2 = sample(1:N2,
                              size = n2)
  
  # Include all cases as well
  index_SubCohort_D2 = c(index_SubCohort_D2,
                         which(Y_D2 == 1))
  index_SubCohort_D2 = unique(index_SubCohort_D2)
  
  # Generate "D2"
  D2 = data.frame(Y = Y_D2,
                  Xstar = X2star_D2,
                  U = U_D2,
                  Strata = 100,
                  Prob = n2 / N2,
                  LLOD = LLOD_X2[3])
  D2 = D2[index_SubCohort_D2,]
  
  # Update the probability if Y == 1
  D2[which(D2$Y == 1),]$Prob = 1
  
  Mean_Y_D2_vec[m] = mean(D2$Y)
  
  
  # ---------------------------------------------------------------------------
  # (3) Estimate Parameters with Population 1
  # ---------------------------------------------------------------------------
  
  # -------------------------
  #   Regression Parameters
  # -------------------------
  
  # Estimate parameters with the true covariate and store the estimated coefficients
  Data_Analysis_1_TrueCovariate = glm(Y ~ U,
                                      data = D1,
                                      family=binomial(link='logit'),
                                      weights = 1/D1$Prob)
  Betahat_TrueCovariate_D1_mat[m,] = as.numeric(Data_Analysis_1_TrueCovariate$coef)
  
  
  # Run the analysis with population 1
  seed0 = 13543645 + m
  Data_Analysis_1 = Data_Analysis_Wrapper(Data = D1,
                                          Setting = 1,
                                          Var_MeasurementError = sigma2_eps_1_vec[2],
                                          B = 1500)
  
  # Store the estimates of alpha and beta
  Est_alpha_D1_mat[m,] = Data_Analysis_1$mod_alpha$x
  Betahat_Naive_D1_mat[m,] = Data_Analysis_1$mod_beta_NoME$mod$x
  Betahat_WMCMLE_D1_mat[m,] = Data_Analysis_1$mod_beta$mod$x
  
  
  # Obtain the BLUP of "U"
  Uhat_D1 = RC_Predict_Function(X = D1$Xstar,
                                LLOD = D1$LLOD[1],
                                mu_U = Data_Analysis_1$mod_alpha$x[1],
                                mu_X = Data_Analysis_1$mod_alpha$x[1],
                                sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                                sigma2_X = Data_Analysis_1$mod_alpha$x[2])
  D1$Uhat = Uhat_D1
  
  # Fit a glm model with weights and store the coefficients
  Data_Analysis_1_RC = glm(Y ~ Uhat,
                           data = D1,
                           family=binomial(link='logit'),
                           weights = 1/D1$Prob)
  Betahat_RC_D1_mat[m,] = as.numeric(Data_Analysis_1_RC$coef)
  
  # -------------------------
  #   Risk Prediction
  # -------------------------
  
  # True Covariate
  Prob_TrueCovariate_1 = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = as.numeric(Data_Analysis_1_TrueCovariate$coef),
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_TrueCovariate_D1_mat[m,] = Prob_TrueCovariate_1
  
  # Naive Estimator
  Prob_Naive_1 = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    lin = Data_Analysis_1$mod_beta_NoME$mod$x[1] + Data_Analysis_1$mod_beta_NoME$mod$x[2]*x
    val = 1/(1+exp(-lin))
    return(val)
  })
  ProbYGivenX_Naive_D1_mat[m,] = Prob_Naive_1
  
  
  # Weighted Monte Carlo Maximum Likelihood Estimator
  Prob_WMCMLE_1 = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = Data_Analysis_1$mod_beta$mod$x,
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_WMCMLE_D1_mat[m,] = Prob_WMCMLE_1
  
  
  # Regression Calibration
  Prob_RC_1 = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = as.numeric(Data_Analysis_1_RC$coef),
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_RC_D1_mat[m,] = Prob_RC_1
  
  
  # ---------------------------------------------------------------------------
  # (4) Estimate Parameters with Population 2
  # ---------------------------------------------------------------------------
  
  # We will have a setting where we estimate delta = E(X2) - E(X1),
  #   and a setting where we set delta = 0
  
  # ---------------------------------
  #   Obtain Paired Sample Estimates
  # ---------------------------------
  
  # Specify the Objective Function
  NegLogLikFcn_PairedSample = function(mu_1,
                                       mu_2,
                                       sigma_11,
                                       sigma_12,
                                       sigma_22){
    
    X1star_D0 = D0$X1star
    X2star_D0 = D0$X2star
    
    # Define mu
    mu = c(mu_1, mu_2)
    
    # Define Sigma
    Sigma = matrix(c(sigma_11, sigma_12, sigma_12, sigma_22),
                   nrow = 2, ncol = length(Beta))
    
    # -------------------------------------
    # Case 1: X_{i1} > a_1 & X_{i2} > a_2
    # -------------------------------------
    Case1_index = which(X1star_D0 > D0$LLOD_X1[1] & X2star_D0 > D0$LLOD_X2[1])
    
    X_Case1 = cbind(X1star_D0[Case1_index],
                    X2star_D0[Case1_index])
    val1 = mvtnorm::dmvnorm(X_Case1, 
                            mean = mu, 
                            sigma = Sigma,
                            log = T)
    
    # -------------------------------------
    # Case 2: X_{i1} = a_1 & X_{i2} > a_2
    # -------------------------------------
    Case2_index = which(X1star_D0 == D0$LLOD_X1[1] & X2star_D0 > D0$LLOD_X2[1])
    
    X_Case2 = cbind(X1star_D0[Case2_index],
                    X2star_D0[Case2_index])
    mu_1_Case2 = mu_1 + sigma_12 / sigma_22 * (X_Case2[,2] - mu_2)
    sigma_11_Case2 = sigma_11 - sigma_12^2 / sigma_22
    
    val2 = log( dnorm(X_Case2[,2],
                      mean = mu_2,
                      sd = sqrt(sigma_22)) ) + log(pnorm(X_Case2[,1],
                                                         mean = mu_1_Case2,
                                                         sd = sqrt(sigma_11_Case2)) )
    
    # -------------------------------------
    # Case 3: X_{i1} > a_1 & X_{i2} == a_2
    # -------------------------------------
    Case3_index = which(X1star_D0 > D0$LLOD_X1[1] & X2star_D0 == D0$LLOD_X2[1])
    
    X_Case3 = cbind(X1star_D0[Case3_index],
                    X2star_D0[Case3_index])
    mu_2_Case3 = mu_2 + sigma_12 / sigma_11 * (X_Case2[,1] - mu_1)
    sigma_22_Case3 = sigma_22 - sigma_12^2 / sigma_11
    val3 = log(dnorm(X_Case3[,1],
                     mean = mu_1,
                     sd = sqrt(sigma_11))) + log(pnorm(X_Case3[,2],
                                                       mean = mu_2_Case3[1],
                                                       sd = sqrt(sigma_22_Case3)))
    
    # -------------------------------------
    # Case 4: X_{i1} == a_1 & X_{i2} == a_2
    # -------------------------------------
    
    Case4_index = which(X1star_D0 == D0$LLOD_X1[1] & X2star_D0 == D0$LLOD_X2[1])
    
    X_Case4 = cbind(X1star_D0[Case4_index],
                    X2star_D0[Case4_index])
    
    integrand = function(x){
      mu_x = mu_1 + sigma_12 / sigma_22 * (x - mu_2)
      val = dnorm(x,
                  mean = mu_2,
                  sd = sqrt(sigma_22)) * pnorm(D0$LLOD_X2[1],
                                               mean = mu_x,
                                               sd = sqrt(sigma_11_Case2))
      return(val)
    }
    Integral_val = cubature::cuhre(f=integrand,
                                   nComp = 1,
                                   lowerLimit= -Inf,
                                   upperLimit= D0$LLOD_X1[1] )$integral
    
    val4 = rep(log(Integral_val), length(Case4_index))
    
    # Take the sum and return the negative
    val = -sum(c(val1, val2, val3, val4))
    return(val)
  }
  
  
  # Setup the objective function
  OF_PairedSample = function(par){
    NegLogLikFcn_PairedSample(par[1],
                              par[2],
                              par[3],
                              par[4],
                              par[5])
  }
  
  # Test out the function
  start_PairedSample = c(mean(D0$X1star),
                         mean(D0$X2star),
                         var(D0$X1star),
                         cov(D0$X1star,
                             D0$X2star),
                         var(D0$X2star))
  OF_PairedSample(start_PairedSample)
  
  # Minimize the objective function
  MLE_PairedSample = suppressWarnings(nlm(f = OF_PairedSample,
                                          p = start_PairedSample))
  
  # -------------------------
  #   Regression Parameters
  # -------------------------
  
  # ------------------------------------------------------------------
  # Consider the case where we estimated delta from the paired sample
  # ------------------------------------------------------------------
  
  Data_Analysis_2_BridgeYes = Data_Analysis_Wrapper(Data = D2,
                                                    Setting = 2,
                                                    Delta_Zero = F,
                                                    Var_MeasurementError = sigma2_eps_2_vec[3],
                                                    PairedSample_Estimates = MLE_PairedSample$estimate,
                                                    B = 1500)
  
  # Store the value of the bridged mean
  Est_BridgedMean_vec[m] = Data_Analysis_2_BridgeYes$Bridged_Mean
  
  # Store the estimates of alpha and beta
  Est_alpha_D2_mat[m,] = Data_Analysis_2_BridgeYes$mod_alpha$x
  Betahat_WMCMLE_D2_BridgeYes_mat[m,] = Data_Analysis_2_BridgeYes$mod_beta$mod$x
  
  # --------------------------------------------------
  #   Consider the case where we set delta to zero
  # --------------------------------------------------
  
  Data_Analysis_2_BridgeNo = Data_Analysis_Wrapper(Data = D2,
                                                   Setting = 2,
                                                   Delta_Zero = T,
                                                   Var_MeasurementError = sigma2_eps_2_vec[3],
                                                   PairedSample_Estimates = MLE_PairedSample$estimate,
                                                   B = 1500)
  
  # Store the estimates of beta
  Betahat_WMCMLE_D2_BridgeNo_mat[m,] = Data_Analysis_2_BridgeNo$mod_beta$mod$x
  
  
  # ---------------------------------------------------------------------------
  # (5) Estimate Parameters with Populations 1 & 2
  # ---------------------------------------------------------------------------
  
  # We will have a setting where we estimate delta = E(X2) - E(X1),
  #   and a setting where we set delta = 0
  
  # Merge the datasets together
  D12 = data.frame(Y = c(D1$Y,
                         D2$Y),
                   Xstar = c(D1$Xstar,
                             D2$Xstar),
                   U = c(D1$U,
                         D2$U),
                   Strata = c(D1$Strata,
                              D2$Strata),
                   Sample1_ind = c(D1$Strata == 1,
                                   D2$Strata == 1)*1,
                   Sample2_ind = c(D1$Strata == 100,
                                   D2$Strata == 100)*1,
                   Prob = c(D1$Prob,
                            D2$Prob),
                   LLOD = c(D1$LLOD,
                            D2$LLOD))
  
  # ------------------------------------------------------------------
  # Consider the case where we estimated delta from the paired sample
  # ------------------------------------------------------------------
  
  # -------------------------
  #   Regression Parameters
  # -------------------------
  
  
  # Obtain the BLUP of "U" for population 2
  Uhat_BridgeYes_D2 = RC_Predict_Function(X = D2$Xstar,
                                          LLOD = D2$LLOD[1],
                                          mu_U = Data_Analysis_2_BridgeYes$Bridged_Mean,
                                          mu_X = Data_Analysis_2_BridgeYes$mod_alpha$x[1],
                                          sigma2_U = Data_Analysis_2_BridgeYes$mod_alpha$x[2] - sigma2_eps_2_vec[3],
                                          sigma2_X = Data_Analysis_2_BridgeYes$mod_alpha$x[2])
  
  D12$Uhat_BridgeYes = c(Uhat_D1, Uhat_BridgeYes_D2)
  
  # Fit a glm model with weights and store the coefficients
  Data_Analysis_12_RC_BridgeYes = glm(Y ~ Uhat_BridgeYes,
                                      data = D12,
                                      family=binomial(link='logit'),
                                      weights = 1/D12$Prob)
  Betahat_RC_D12_BridgeYes_mat[m,] = as.numeric(Data_Analysis_12_RC_BridgeYes$coef)
  
  
  # Directly solve the estimating equation
  Data_Analysis_12_BridgeYes_EE = Data_Analysis_Wrapper(Data = D12,
                                                        Setting = 3,
                                                        Var_MeasurementError = sigma2_eps_2_vec[3], # IRRELEVANT!
                                                        Delta_Zero = F,
                                                        PairedSample_Estimates = MLE_PairedSample$estimate, # IRRELEVANT!
                                                        B = 1500)
  # Store the estimate of beta
  Betahat_WMCMLE_D12_BridgeYes_mat[m,] = Data_Analysis_12_BridgeYes_EE$mod_beta$mod$x
  
  
  # -------------------------
  #   Risk Prediction
  # -------------------------
  
  # Regression Calibration
  Prob_RC_12_BridgeYes = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = as.numeric(Data_Analysis_12_RC_BridgeYes$coef),
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_RC_D12_BrigeYes_mat[m,] = Prob_RC_12_BridgeYes
  
  # WMCMLE
  Prob_WMCMLE_12_BridgeYes = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = Data_Analysis_12_BridgeYes_EE$mod_beta$mod$x,
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_WMCMLE_D12_BridgeYes_mat[m,] = Prob_WMCMLE_12_BridgeYes
  
  # --------------------------------------------------
  #   Consider the case where we set delta to zero
  # --------------------------------------------------
  
  
  # -------------------------
  #   Regression Parameters
  # -------------------------
  
  
  # Obtain the BLUP of "U" for population 2
  Uhat_BridgeNo_D2 = RC_Predict_Function(X = D2$Xstar,
                                         LLOD = D2$LLOD[1],
                                         mu_U = Data_Analysis_2_BridgeNo$Bridged_Mean,
                                         mu_X = Data_Analysis_2_BridgeNo$mod_alpha$x[1],
                                         sigma2_U = Data_Analysis_2_BridgeNo$mod_alpha$x[2] - sigma2_eps_2_vec[3],
                                         sigma2_X = Data_Analysis_2_BridgeNo$mod_alpha$x[2])
  
  D12$Uhat_BridgeNo = c(Uhat_D1, Uhat_BridgeNo_D2)
  
  # Fit a glm model with weights and store the coefficients
  Data_Analysis_12_RC_BridgeNo = glm(Y ~ Uhat_BridgeNo,
                                     data = D12,
                                     family=binomial(link='logit'),
                                     weights = 1/D12$Prob)
  Betahat_RC_D12_BridgeNo_mat[m,] = as.numeric(Data_Analysis_12_RC_BridgeNo$coef)
  
  
  # Directly solve the estimating equation
  Data_Analysis_12_BridgeNo_EE = Data_Analysis_Wrapper(Data = D12,
                                                       Setting = 3,
                                                       Var_MeasurementError = sigma2_eps_2_vec[3], # IRRELEVANT!
                                                       Delta_Zero = T,
                                                       PairedSample_Estimates = MLE_PairedSample$estimate, # IRRELEVANT!
                                                       B = 1500)
  # Store the estimate of beta
  Betahat_WMCMLE_D12_BridgeNo_mat[m,] = Data_Analysis_12_BridgeNo_EE$mod_beta$mod$x
  
  
  # -------------------------
  #   Risk Prediction
  # -------------------------
  
  # Regression Calibration
  Prob_RC_12_BridgeNo = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = as.numeric(Data_Analysis_12_RC_BridgeNo$coef),
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_RC_D12_BrigeYes_mat[m,] = Prob_RC_12_BridgeNo
  
  # WMCMLE
  Prob_WMCMLE_12_BridgeNo = sapply(1:length(X1_grid), function(i){
    x = X1_grid[i]
    set.seed(seed0)
    val = ProbYGivenX(betahat = Data_Analysis_12_BridgeNo_EE$mod_beta$mod$x,
                      X = x,
                      mu_U = Data_Analysis_1$mod_alpha$x[1],
                      mu_X = Data_Analysis_1$mod_alpha$x[1],
                      sigma2_U = Data_Analysis_1$mod_alpha$x[2] - sigma2_eps_1_vec[2],
                      sigma2_X = Data_Analysis_1$mod_alpha$x[2],
                      B = 10000,
                      seed0 = seed0)
    return(val)
  })
  ProbYGivenX_WMCMLE_D12_BridgeNo_mat[m,] = Prob_WMCMLE_12_BridgeNo
}
