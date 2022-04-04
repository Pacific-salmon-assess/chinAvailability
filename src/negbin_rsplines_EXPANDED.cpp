#include <TMB.hpp>

// List of matrices
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
  
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(rfac1);
  DATA_INTEGER(n_rfac1);
  DATA_IVECTOR(b_smooth_start);
  DATA_STRUCT(Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix

  //predictions
  DATA_MATRIX(pred_X1_ij); // matrix for FE predictions
  DATA_STRUCT(pred_Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(pred_Xs); // smoother linear effect matrix 
  DATA_IVECTOR(pred_rfac1);
  // vector of higher level aggregates used to generate predictions; length
  // is equal to the number of predictions made
  DATA_IVECTOR(pred_rfac_agg);
  DATA_IVECTOR(pred_rfac_agg_levels);
  

  // PARAMETERS ----------------------------------------------------------------
  
  PARAMETER_VECTOR(b1_j); // fixed effects parameters
  PARAMETER(ln_phi);     // variance
  PARAMETER_VECTOR(bs);   // smoother linear effects
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  // random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(a1);
  PARAMETER(ln_sigma_a1);
  

  // DERIVED -------------------------------------------------------------------

  int n1 = y1_i.size();
  int n_predX1 = pred_X1_ij.rows(); // number of finest scale predictions (abundance only)  
  int n_predX2 = pred_rfac_agg_levels.size();

  Type jnll = 0.0; // initialize joint negative log likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------

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
  vector<Type> eta_i(n1);
  eta_i.setZero();
  for (int i = 0; i < n1; i++) {
    eta_i(i) = eta_fx_i(i) + eta_smooth_i(i); 
  }
  
  // Add random intercepts
  vector<Type> mu_i(n1); 
  mu_i.setZero();
  for (int i = 0; i < n1; i++) {
    mu_i(i) = eta_i(i) + a1(rfac1(i));
  }


  // LIKELIHOOD ----------------------------------------------------------------

  // Type s1, s2;
  vector<Type> s1(n1);
  vector<Type> s2(n1);
  for(int i = 0; i < n1; i++){
    s1(i) = mu_i(i); //mu
    s2(i) = 2.0 * (s1(i) - ln_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1(i), s2(i), true);
  }

  // Report for residuals
  // ADREPORT(s1);
  // ADREPORT(s2);

  // Probability of random coefficients
  for(int k = 0; k < n_rfac1; k++){
    if (k == 0) {
      jnll -= dnorm(a1(k), Type(0.0), exp(ln_sigma_a1), true);  
    }
    if (k > 0) {
      jnll -= dnorm(a1(k), a1(k - 1), exp(ln_sigma_a1), true);
    }
  }


  // PREDICTIONS ---------------------------------------------------------------

  vector<Type> pred_mu1 = pred_X1_ij * b1_j;

  // smoothers
  vector<Type> pred_smooth_i(n_predX1);
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
  for (int i = 0; i < n_predX1; i++) {
    pred_mu1(i) += pred_smooth_i(i);
  }

  // add random intercepts 
  for (int i = 0; i < n_predX1; i++) {
    pred_mu1(i) += a1(pred_rfac1(i));
  }

  REPORT(pred_mu1);
  ADREPORT(pred_mu1);


  // Calculate predicted abundance based on higher level groupings
  
  vector<Type> pred_mu1_cumsum(n_predX2);
  vector<Type> ln_pred_mu1_cumsum(n_predX2);
  vector<Type> exp_pred_mu1 = exp(pred_mu1); // calculate real values for summing


  for (int i = 0; i < n_predX1; i++) {
    for (int m = 0; m < n_predX2; m++) {
      if (pred_rfac_agg(i) == pred_rfac_agg_levels(m)) {
        pred_mu1_cumsum(m) += exp_pred_mu1(i);
        ln_pred_mu1_cumsum(m) = log(pred_mu1_cumsum(m));
      }
    }
  }

  ADREPORT(pred_mu1_cumsum);
  ADREPORT(ln_pred_mu1_cumsum);


  return jnll;
}
