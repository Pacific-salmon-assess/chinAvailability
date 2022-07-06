#include <TMB.hpp>

// List of matrices (necessary for smooths)
template <class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x){  // x = list passed from R
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{

  // DATA ----------------------------------------------------------------------
  
  // Abundance data
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(rfac1);
  DATA_INTEGER(n_rfac1);
  DATA_IVECTOR(b_smooth_start);
  DATA_STRUCT(Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix

  // Composition data
  DATA_MATRIX(Y2_ik);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(X2_ij);    // covariate matrix (a n-by-j matrix)
  DATA_IVECTOR(rfac2);    // vector of random factor levels
  DATA_INTEGER(n_rfac2);  // number of random factor levels
  
  // Abundance predictions
  DATA_MATRIX(pred_X1_ij); // matrix for FE predictions
  DATA_STRUCT(pred_Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(pred_Xs); // smoother linear effect matrix 

  // Composition predictions
  DATA_MATRIX(pred_X2_ij);    // model matrix for predictions

  // Conditionals
  DATA_INTEGER(random_walk);

  // PARAMETERS ----------------------------------------------------------------
  
  // Abundance parameters
  PARAMETER_VECTOR(b1_j); // fixed effects parameters
  PARAMETER(ln_phi);     // variance
  PARAMETER_VECTOR(bs);   // smoother linear effects
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  // random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER(ln_sigma_a1);
  PARAMETER_VECTOR(a1); // vector of random intercepts (RW)
  
  // Composition parameters
  PARAMETER_MATRIX(B2_jk); // parameter matrix
  PARAMETER(ln_sigma_A2); // among random intercept SD
  PARAMETER_VECTOR(A2_h);  // vector of random intercepts


  // DERIVED QUANTITIES --------------------------------------------------------

  int n1 = y1_i.size();
  int n2 = Y2_ik.rows();         // number of observations
  int n_cat = Y2_ik.cols();         // number of categories
  int n_predX = pred_X2_ij.rows(); // number of predictions  
  
  // Intermediate objects
  vector<Type> eta_i(n1);
  vector<Type> mu_i(n1); 
  matrix<Type> Mu2_ik(n2, n_cat); // matrix of combined fixed/random eff

  Type jnll = 0.0; // initialize joint negative log likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------

  // Abundance
  vector<Type> eta_fx_i = X1_ij * b1_j;
  // Smooths
  vector<Type> eta_smooth_i(X1_ij.rows());
  eta_smooth_i.setZero();

  for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
    vector<Type> beta_s(Zs(s).cols());
    beta_s.setZero();
    for (int j = 0; j < beta_s.size(); j++) {
      beta_s(j) = b_smooth(b_smooth_start(s) + j);
      jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
    }
    eta_smooth_i += Zs(s) * beta_s;
  }
  eta_smooth_i += Xs * bs;

  // Combine smooths and linear
  eta_i.setZero();
  for (int i = 0; i < n1; i++) {
    eta_i(i) = eta_fx_i(i) + eta_smooth_i(i); 
  }
  
  // Add random intercepts
  mu_i.setZero();
  for (int i = 0; i < n1; i++) {
    mu_i(i) = eta_i(i) + a1(rfac1(i));
  }

  // Composition (no random smooths)
  matrix<Type> Mu2_fx_ik = X2_ij * B2_jk; // fixed effects

  for (int i = 0; i < n2; ++i) {
    for(int k = 0; k < n_cat; k++) {
      Mu2_ik(i, k) = Mu2_fx_ik(i, k) + A2_h(rfac2(i));
    }
  }

  matrix<Type> Gamma = exp(Mu2_ik.array()); // add random effect
  vector<Type> n_plus = Y2_ik.rowwise().sum(); // row sum of response
  vector<Type> Gamma_plus = Gamma.rowwise().sum(); // row sum of gamma
 

  // ABUNDANCE LIKELIHOOD ------------------------------------------------------

  // Type s1, s2;
  vector<Type> s1(n1);
  vector<Type> s2(n1);
  for(int i = 0; i < n1; i++){
    s1(i) = mu_i(i); //mu
    s2(i) = 2.0 * (s1(i) - ln_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1(i), s2(i), true);
  }

  // Probability of abundance random coefficients (random walk)
  if (random_walk) {
    for (int k = 0; k < n_rfac1; k++){
      if (k == 0) {
        jnll -= dnorm(a1(k), Type(0.0), exp(ln_sigma_a1), true);  
      }
      if (k > 0) {
        jnll -= dnorm(a1(k), a1(k - 1), exp(ln_sigma_a1), true);
      }
    }
  } else {
    for (int k = 0; k < n_rfac1; k++){
      jnll -= dnorm(a1(k), Type(0.0), exp(ln_sigma_a1), true);  
    }
  }

  // COMPOSITION LIKELIHOOD ----------------------------------------------------

  Type jll = 0; // initialize joint log-likelihood

  for(int i = 0; i <= (n2 - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(Gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + Gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((Y2_ik(i, k) + Gamma(i, k)));
      jll -= lgamma(Gamma(i, k));
      jll -= lgamma((Y2_ik(i, k) + 1));
    }
  }

  jnll -= jll;
  
  // Probability of random intercepts
  if (random_walk) {
    for (int h = 0; h < n_rfac2; h++) {
      if (h == 0) {
        jnll -= dnorm(A2_h(h), Type(0.0), exp(ln_sigma_A2), true);  
      }
      if (h > 0) {
        jnll -= dnorm(A2_h(h), A2_h(h - 1), exp(ln_sigma_A2), true);
      }
    } 
  } else {
    for (int h = 0; h < n_rfac2; h++) {
      jnll -= dnorm(A2_h(h), Type(0.0), exp(ln_sigma_A2), true);
    }
  }

  Type sigma_rfac2 = exp(ln_sigma_A2);
  REPORT(sigma_rfac2);
  

  // PREDICTIONS ---------------------------------------------------------------

  // Predicted composition 
  matrix<Type> pred_Mu2(n_predX, n_cat);    //pred FE + RE on log scale
  matrix<Type> pred_Gamma(n_predX, n_cat);  //transformed pred effects 
  vector<Type> pred_Gamma_plus(n_predX);        
  vector<Type> pred_theta(n_predX); 
  matrix<Type> pred_Pi(n_predX, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_predX); 
  matrix<Type> pred_Pi_prop(n_predX, n_cat); // predicted counts as ppn.
  matrix<Type> logit_pred_Pi_prop(n_predX, n_cat); 

  pred_Mu2 = pred_X2_ij * B2_jk; 

  pred_Gamma = exp(pred_Mu2.array());
  pred_Gamma_plus = pred_Gamma.rowwise().sum();
  pred_theta = 1 / (pred_Gamma_plus + 1);
  for(int m = 0; m < n_predX; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_Pi(m, k) = pred_Gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_Pi.rowwise().sum();
  for(int m = 0; m < n_predX; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_Pi_prop(m, k) = pred_Pi(m, k) / pred_n_plus(m);
      logit_pred_Pi_prop(m, k) = logit(pred_Pi_prop(m, k));
    }
  }
 
  // Predicted abundance
  vector<Type> pred_mu1 = pred_X1_ij * b1_j;

  ADREPORT(pred_mu1);
  ADREPORT(logit_pred_Pi_prop);


  // smoothers
  vector<Type> pred_smooth_i(n_predX);
  pred_smooth_i.setZero();
  for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
    vector<Type> beta_s(pred_Zs(s).cols());
    beta_s.setZero();
    for (int j = 0; j < beta_s.size(); j++) {
      beta_s(j) = b_smooth(b_smooth_start(s) + j);
    }
    pred_smooth_i += pred_Zs(s) * beta_s;
  }
  pred_smooth_i += pred_Xs * bs;
  
  // combine fixed and smoothed predictions
  for (int m = 0; m < n_predX; m++) {
    pred_mu1(m) += pred_smooth_i(m);
  }
  // multiply abundance and composition predictions
  vector<Type> exp_pred_mu1 = exp(pred_mu1); // calculate real values for summing
  matrix<Type> pred_mu1_Pi(n_predX, n_cat);
  for (int m = 0; m < n_predX; m++) {
    for (int k = 0; k < n_cat; k++) {
      pred_mu1_Pi(m, k) = exp_pred_mu1(m) * pred_Pi_prop(m, k);
    }
  }
  matrix<Type> log_pred_mu1_Pi = log(pred_mu1_Pi.array());
  
  ADREPORT(log_pred_mu1_Pi);

  return jnll;
}
