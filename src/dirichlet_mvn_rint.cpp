#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // load namespace with multivariate distributions
  using namespace density;

  // DATA ----------------------------------------------------------------------

  DATA_MATRIX(Y2_ik);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(X2_ij);    // covariate matrix (a n-by-j matrix)
  DATA_IVECTOR(rfac2);    // vector of random factor levels
  DATA_INTEGER(n_rfac2);  // number of random factor levels
  
  //predictions
  DATA_MATRIX(pred_X2_ij);    // model matrix for predictions
  DATA_IVECTOR(pred_rfac2); // vector of predicted random intercepts


  // PARAMETERS ----------------------------------------------------------------

  PARAMETER_MATRIX(B2_jk); // parameter matrix
  PARAMETER_MATRIX(A2_hk);  // matrix of random intercepts (n_rfac2 x n_cat)
  PARAMETER(ln_sigma_A2); // among random intercept SD


  // DERIVED QUANTITIES --------------------------------------------------------

  int n_obs = Y2_ik.rows();         // number of observations
  int n_cat = Y2_ik.cols();         // number of categories
  int n_levels = pred_X2_ij.rows();   // number of covariates to make predictions on

  // Matrix for intermediate objects
  matrix<Type> Mu2_ik(n_obs, n_cat); // matrix of combined fixed/random eff

  // Covariance matrix for MVN random intercepts
  matrix<Type> cov_mat(n_cat, n_cat);
  for (int j = 0; j < n_cat; j++) {
    for (int jj = 0; jj < n_cat; jj++) {
      if (j == jj) {
        cov_mat(j, jj) = exp(ln_sigma_A2) * exp(ln_sigma_A2);
      } else {
        cov_mat(j, jj) = 0;
      }
    }
  }

  Type jll = 0; // initialize joint log-likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------
  
  matrix<Type> Mu2_ik_fx = X2_ij * B2_jk; // fixed effects

  for (int i = 0; i < n_obs; ++i) {
    for(int k = 0; k < n_cat; k++) {
      Mu2_ik(i, k) = Mu2_ik_fx(i, k) + A2_hk(rfac2(i), k);
    }
  }

  matrix<Type> gamma = exp(Mu2_ik.array()); // add random effect
  vector<Type> n_plus = Y2_ik.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  
  
  // LIKELIHOOD ----------------------------------------------------------------

  for(int i = 0; i <= (n_obs - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((Y2_ik(i, k) + gamma(i, k)));
      jll -= lgamma(gamma(i, k));
      jll -= lgamma((Y2_ik(i, k) + 1));
    }
  }

  Type jnll;
  jnll = -jll;

  
  // Probability of multivariate random intercepts
  for (int h = 0; h < n_rfac2; h++) {
    vector<Type> A2_hk_vec = A2_hk.row(h);
    MVNORM_t<Type> neg_log_dmvnorm(cov_mat);
    jnll += neg_log_dmvnorm(A2_hk_vec);
  }

  Type sigma_rfac2 = exp(ln_sigma_A2);
  ADREPORT(sigma_rfac2);
  

  // PREDICTIONS ---------------------------------------------------------------
  
  matrix<Type> pred_mus_fx(n_levels, n_cat);    //pred fixed effects on log scale
  matrix<Type> pred_mu(n_levels, n_cat);    //pred FE + RE on log scale
  matrix<Type> pred_gamma(n_levels, n_cat);  //transformed pred effects 
  vector<Type> pred_gamma_plus(n_levels);        
  vector<Type> pred_theta(n_levels); 
  matrix<Type> pred_pi(n_levels, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_levels); 
  matrix<Type> pred_pi_prop(n_levels, n_cat); // predicted counts as ppn.
  matrix<Type> logit_pred_pi_prop(n_levels, n_cat); 

  pred_mu_fx = pred_X2_ij * B2_jk; 

  for (int i = 0; i < n_levels; ++i) {
    for(int k = 0; k < n_cat; k++) {
      pred_mu(i, k) = pred_mu_fx(i, k) + A2_hk(pred_rfac2(i), k);
    }
  }
  pred_gamma = exp(pred_mu.array());
  pred_gamma_plus = pred_gamma.rowwise().sum();
  pred_theta = 1 / (pred_gamma_plus + 1);
  for(int m = 0; m < n_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_pi.rowwise().sum();
  for(int m = 0; m < n_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
      logit_pred_pi_prop(m, k) = logit(pred_pi_prop(m, k));
    }
  }

  ADREPORT(pred_mu);
  ADREPORT(logit_pred_pi_prop);
  
  // Return negative loglikelihood
  return jnll;
  
}
