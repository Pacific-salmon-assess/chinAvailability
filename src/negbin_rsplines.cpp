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

  using namespace density;
  using namespace Eigen;

  // DATA ----------------------------------------------------------------------
  
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  // DATA_IVECTOR(factor1k_i);
  // DATA_INTEGER(nk1);
  DATA_IVECTOR(b_smooth_start);
  DATA_STRUCT(Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix

  //predictions
  DATA_MATRIX(pred_X1_ij); // matrix for FE predictions
  DATA_STRUCT(pred_Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(pred_Xs); // smoother linear effect matrix 
  // vector of higher level aggregates used to generate predictions; length
  // is equal to the number of predictions made
  // DATA_IVECTOR(pred_factor2k_h);
  // DATA_IVECTOR(pred_factor2k_levels);
  

  // PARAMETERS ----------------------------------------------------------------
  
  PARAMETER_VECTOR(b1_j); // fixed effects parameters
  PARAMETER(log_phi);     // variance
  PARAMETER_VECTOR(bs);   // smoother linear effects
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  // random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  // PARAMETER_VECTOR(z1_k);
  // PARAMETER(log_sigma_zk1);
  

  Type jnll = 0.0; // initialize joint negative log likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------

  vector<Type> eta_fixed_i = X1_ij * b1_j;

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
  int n1 = y1_i.size();
  vector<Type> eta_i(n1);
  eta_i.setZero();
  for (int i = 0; i < n1; i++) {
    eta_i(i) = eta_fixed_i(i) + eta_smooth_i(i); // + offset_i(i);
  }
  
  vector<Type> mu_i = exp(eta_i);


  // LIKELIHOOD ----------------------------------------------------------------

  // Type s1, s2;
  vector<Type> s1(n1);
  vector<Type> s2(n1);
  for(int i = 0; i < n1; i++){
    s1(i) = log(mu_i(i)); // + z1_k(factor1k_i(i)); //mu
    s2(i) = 2.0 * (s1(i) - log_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1(i), s2(i), true);
  }

  // Report for residuals
  // ADREPORT(s1);
  // ADREPORT(s2);

  // Probability of random coefficients
  // for(int k = 0; k < nk1; k++){
  //   if (k == 0) {
  //     jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);  
  //   }
  //   if (k > 0) {
  //     jnll -= dnorm(z1_k(k), z1_k(k - 1), exp(log_sigma_zk1), true);
  //   }
  // }


  // PREDICTIONS ---------------------------------------------------------------

  vector<Type> pred_fe = pred_X1_ij * b1_j;

  // smoothers
  vector<Type> pred_smooth_i(pred_X1_ij.rows());
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
  
  //combine fixed and smoothed predictions
  for (int i = 0; i < pred_X1_ij.rows(); i++) {
    pred_fe(i) += pred_smooth_i(i);
  }

  REPORT(pred_fe);
  ADREPORT(pred_fe);
  // ADREPORT(pred_abund);
  // ADREPORT(log_pred_abund);

  // Calculate predicted abundance based on higher level groupings
  // int n_preds = pred_factor2k_h.size();
  // int n_pred_levels = pred_factor2k_levels.size();
  // vector<Type> agg_pred_abund(n_pred_levels);
  // vector<Type> log_agg_pred_abund(n_pred_levels);

  // for (int h = 0; h < n_preds; h++) {
  //   for (int g = 0; g < n_pred_levels; g++) {
  //     if (pred_factor2k_h(h) == pred_factor2k_levels(g)) {
  //       agg_pred_abund(g) += pred_abund(h);
  //       log_agg_pred_abund(g) = log(agg_pred_abund(g));
  //     }
  //   }
  // }

  // ADREPORT(agg_pred_abund);
  // ADREPORT(log_agg_pred_abund);

  return jnll;
}
