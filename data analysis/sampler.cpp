// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

#include <RcppEigen.h>
#include <RcppGSL.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#include <random>
#include <vector>
#include <iostream>
#include <cmath>

// Function to calculate the value
double calculateM(const arma::mat& Z, const arma::mat& Minc){
  arma::vec nonZeroElements;
  int idx = 0;
  for (int i = 0; i < Z.n_rows; i++) {
    for (int j = 0; j < Z.n_cols; j++) {
      if (Minc(i, j) != 0 && Z(i, j) != 0) {
        nonZeroElements.resize(idx + 1); // Resize the vector to accommodate new element
        nonZeroElements(idx) = Z(i,j);
        idx++;
      }
    }
  }
  
  if (nonZeroElements.is_empty()) {
    // Handle the case where nonZeroElements is empty
    return 0.0; // Or any other appropriate value
  }
  
  arma::vec probs = {0.25, 0.5, 0.75};
  arma::vec quantiles = arma::quantile(nonZeroElements, probs);
  double median = quantiles(1);
  
  if (median == 0.0) {
    // Handle the case where the median is zero to avoid division by zero
    return 0.0; // Or any other appropriate value
  }
  
  return -std::log(0.5)/median;
}

arma::mat compute_b(int n, int q, const arma::mat& inv_A, const arma::vec& r) {
  
  arma::mat b(n, q, arma::fill::zeros);
  
  for (int d = 0; d < q; d++) {
    for (int h = 0; h<= d; h++) {
      b.col(d) += inv_A(d,h) * r.subvec(h * n, (h + 1) * n - 1);  
    }
  }
  
  return b;
  
}

std::vector<arma::mat> Tau_new(const arma::vec& rho, 
                               int n, 
                               int q, 
                               const std::vector<arma::mat>& W) {
  
  std::vector<arma::mat> Taur(q, arma::mat(n, n, arma::fill::zeros));
  
  for (int d = 0; d < q; d++) {
    
    const arma::mat& Winc = arma::trimatl(W[d]);
    const arma::vec ns = arma::sum(Winc, 1);
    
    Taur[d].diag() = (1 + (ns - 1) * (rho(d) * rho(d))) / (1 - (rho(d) * rho(d)));
    
    }
  
  return Taur;
  
}

std::vector<arma::mat> Q_new(const arma::vec& rho, 
                             int n, int q, 
                             const std::vector<arma::mat>& W, 
                             const std::vector<arma::mat>& Tau) {
  
  std::vector<arma::mat> Q(q, arma::mat(n, n, arma::fill::zeros));
  arma::mat I = arma::eye(n,n);
  
  for (int d = 0; d < q; d++) {
    
    const arma::mat& Winc = arma::trimatl(W[d]);
    const arma::vec ns = arma::sum(Winc, 1);
    
    arma::mat B_d = arma::zeros(n, n);
    
    // B_d is strictly lower triangular
    for (int i = 1; i < n; i++) {
      for (int j = 0; j < i; j++) {
        if (Winc(i, j) == 1) {
          B_d(i, j) = rho(d) / (1 + (ns(i) - 1) * (rho(d) * rho(d)));
        }
      }
    }
    
    arma::mat ImB_d =  I - B_d;
    
    Q[d] = ImB_d.t() * Tau[d] * ImB_d;
    
    Q[d] = 0.5 * (Q[d] + Q[d].t());
    
  }
  
  return Q;
}

arma::vec compute_inv_Vr_diag(int n,
                              int q,
                              const arma::mat& inv_A,
                              const std::vector<arma::mat>& Q) {

  arma::vec inv_Vr_diag(n * q, arma::fill::zeros);

  for (int d = 0; d < q; d++) {
    for (int h = d; h < q; h++) {
      inv_Vr_diag.subvec(d * n, (d + 1) * n - 1) +=  inv_A(h,d) * inv_A(h,d) * Q[h].diag();
    }
  }
  
  return inv_Vr_diag;
  
}

double ldmvnorm(int n, int q, 
                const arma::mat& b, const arma::mat& A,
                const std::vector<arma::mat>& Q, 
                const std::vector<arma::mat>& Tau) {
  
  double ldens = - n * arma::sum(arma::log(A.diag()));
  
  for(int d = 0; d < q; d++){
    ldens += 0.5 * arma::sum(arma::log(Tau[d].diag())) -
      0.5 * arma::dot(b.col(d), Q[d] * b.col(d));
  }
  
  return ldens;
  
}

// Function to compute stick-breaking weights

arma::vec makeprobs(arma::vec v) {
  
  int m = v.n_elem;
  arma::vec probs(m,arma::fill::zeros);
  
  probs(0) = v(0);
  probs(1) = std::exp(std::log(v(1)) + std::log(1 - v(0)));
  
  for (int i = 2; i < m - 1; i++) {
    probs(i) = std::exp(std::log(v(i)) + arma::sum(arma::log(1 - v.subvec(0, i - 1))));
  }
  
  probs(m - 1) = 1 - sum(probs.subvec(0, m - 2));
  
  return probs;
  
}

arma::ivec makeu(arma::vec F_r, arma::vec probs) {
  
  int m1 = F_r.n_elem;
  int m2 = probs.n_elem;
  
  arma::ivec u(m1, arma::fill::ones);
  
  for (int k = 0; k < m1; k++) {
    for (int l = 0; l < m2; l++) {
      if (l == 0 && F_r(k) <= probs(0)){
        u(k) = 1;
      } else if (l > 0 && arma::sum(probs.subvec(0, l - 1)) < F_r(k) && 
        F_r(k) <= arma::sum(probs.subvec(0, l))) {
        u(k) = l + 1;
      }
    }
  }
  
  return u;
  
}

int makeuk(double F_r, arma::vec probs) {
  
  int m2 = probs.n_elem;
  
  int uk = 1;
  
  for (int l = 0; l < m2; l++) {
    if (l == 0 && F_r <= probs(0)){
      uk = 1;
    } else if (l > 0 && arma::sum(probs.subvec(0, l - 1)) < F_r && 
      F_r <= arma::sum(probs.subvec(0, l))) {
      uk = l + 1;
    }
  }
  
  return uk;
  
}

// Function for numerical issues: inv_trans_par
double inv_trans_par(double x, double lb, double ub) {
  if (x < -745) {
    return (lb + ub * std::exp(-745)) / (1 + std::exp(-745));
  } else if (x > 16.81) {
    return (lb + ub * std::exp(16.81)) / (1 + std::exp(16.81));
  } else if (x >= -745 && x <= 16.81) {
    return (lb + ub * std::exp(x)) / (1 + std::exp(x));
  }
  return(0.0);
}

// Function to compute Jacobian matrix for updating A
double Jacob_A(const arma::mat& A) {
  double prod = 1.0;
  int q = A.n_cols;
  
  for (int i = 0; i < q; i++) {
    prod *= pow(A(i, i), q - i);
  }
  
  return pow(2, q) * prod;
}

// Function to replicate elements of a vector based on indices from another vector
arma::vec replicate_elements(const arma::vec& theta, const arma::ivec& u) {
  int nq = u.n_elem;
  arma::vec result(nq, arma::fill::zeros);
  
  for (int i = 0; i < nq; i++) {
    int index = u(i) - 1; // Adjust index to 0-based
    if (index >= 0 && index < theta.n_elem) {
      result(i) = theta(index);
    } else {
      // Handle out-of-bounds indices if needed
      result(i) = NA_REAL; // Placeholder value
    }
  }
  
  return result;
}

// [[Rcpp::export]]

// Function to perform MADAGAR algorithm

List MADAGAR(arma::vec y, arma::mat X, arma::mat Z1, arma::mat Z2, arma::mat Z3, arma::vec E, std::string cvrts,
             int q, arma::mat Winc, arma::mat Minc, double alpha, int n_atoms, 
             int runs, int burn, int thin) {
  
  List output(15);
  
  int nq = y.n_elem;
  int n = nq / q;
  int p = X.n_cols;
  
  int iterations = runs * thin + burn;
  
  arma::mat I = arma::eye(n,n);
  
  arma::mat keepbeta(runs, p);
  arma::mat keepphi(runs, nq);
  arma::mat keeptheta(runs, n_atoms);
  arma::imat keepu(runs, nq); 
  arma::mat keeprho(runs, q);
  arma::vec keeptaus(runs, arma::fill::zeros); 
  arma::mat keepv(runs, n_atoms);
  arma::mat keepr(runs, nq);
  arma::mat keepFr(runs, nq);
  arma::mat keepeta(runs, 3*q);
  arma::mat keepA(runs, q*(q + 1)/2);
  
  List keepW1(runs);
  List keepW2(runs);
  List keepW3(runs);
  List keepW4(runs);
  
  for (int i = 0; i < runs; i++) {
    keepW1[i] = arma::zeros<arma::mat>(n, n);
    keepW2[i] = arma::zeros<arma::mat>(n, n);
    keepW3[i] = arma::zeros<arma::mat>(n, n);
    keepW4[i] = arma::zeros<arma::mat>(n, n);
  }
  
  // initial values
  
  // Initialize theta with n_atoms elements, each set to 0.0
  arma::vec theta(n_atoms, arma::fill::zeros);
  
  // Initialize beta with p elements, each set to 0.0
  arma::vec beta(p, arma::fill::zeros);
  arma::vec X_beta = X * beta;
  
  double taub = 0.001;
  double taus = 1.0;
  
  double c = 2.0;
  double d2 = 1.0; // 1.0;
  
  // Initialize v with n_atoms elements, each set to 0.1
  arma::vec v = arma::vec(n_atoms).fill(0.1);
  arma::vec vv = arma::log(v) - arma::log(1 - v);
  
  // Initialize rho with q elements, each set to 0.7
  arma::vec rho = arma::vec(q).fill(0.1);
  arma::vec rhorho = arma::log(rho) - arma::log(1 - rho);
  
  double M1 = calculateM(Z1, Minc);
  double M2 = calculateM(Z2, Minc);
  double M3 = calculateM(Z3, Minc);
  
  arma::vec M1_vec = arma::vec(q).fill(M1);
  arma::vec M2_vec = arma::vec(q).fill(M2);
  arma::vec M3_vec = arma::vec(q).fill(M3);
  
  arma::vec M(3 * q);
  
  // Copy elements from M1_vec
  M.subvec(0, q - 1) = M1_vec;
  M.subvec(q, 2 * q - 1) = M2_vec;
  M.subvec(2 * q, 3 * q - 1) = M3_vec;
  
  // Initialize eta with 3*q elements, each set to 0.1
  arma::vec eta = arma::vec(3 * q).fill(0.1);
  arma::vec etaeta = arma::log(eta) - arma::log(M - eta);
  
  std::vector<arma::mat> W(q,Minc);
  
  std::vector<arma::mat> Tau = Tau_new(rho, n, q, W);

  std::vector<arma::mat> Q = Q_new(rho, n, q, W, Tau);

  arma::vec A_vec = arma::randn<arma::vec>(q * (q + 1) / 2, arma::distr_param(0.0, 0.015));

  arma::mat A(q, q, arma::fill::zeros);

  int idxA = 0;

  for (int d = 0; d < q; d++) {
    for (int h = 0; h <= d; h++) {
      if (h == d) {
        // If diagonal element, set to exp(rnorm(1))
        A(d, h) = std::exp(A_vec(idxA));
      } else {
        // If off-diagonal element, set to rnorm(1)
        A(d, h) = A_vec(idxA);
      }
      idxA++;
    }
  }

  arma::mat inv_A = arma::inv(arma::trimatl(A));
  
  arma::vec inv_Vr_diag = compute_inv_Vr_diag(n, q, inv_A, Q);
  
  arma::mat Sigma = A*A.t();
  arma::mat inv_Sigma = inv_A.t() * inv_A;
  
  arma::vec mu_r = arma::vec(nq, arma::fill::zeros);
  
  arma::mat Q1 = Q[0];
  arma::mat Q2 = Q[1];
  arma::mat Q3 = Q[2];
  arma::mat Q4 = Q[3];
  
  arma::mat invQ1 = inv_sympd(Q1,arma::inv_opts::allow_approx);
  arma::mat invQ2 = inv_sympd(Q2,arma::inv_opts::allow_approx);
  arma::mat invQ3 = inv_sympd(Q3,arma::inv_opts::allow_approx);
  arma::mat invQ4 = inv_sympd(Q4,arma::inv_opts::allow_approx);
  
  arma::mat kprod = arma::kron(A,I);
  
  // Initialize block diagonal matrix
  arma::mat invQ = arma::zeros(n * q, n * q);
  
  // Block 1
  invQ.submat(0, 0, n - 1, n - 1) = invQ1;
  // Block 2
  invQ.submat(n, n, 2 * n - 1, 2 * n - 1) = invQ2;
  // Block 3
  invQ.submat(2 * n, 2 * n, 3 * n - 1, 3 * n - 1) = invQ3;
  // Block 4
  invQ.submat(3 * n, 3 * n, 4 * n - 1, 4 * n - 1) = invQ4;
  
  arma::mat Vr = kprod * invQ * kprod.t();
  
  arma::vec r = arma::mvnrnd(mu_r, Vr, 1);
  
  arma::mat b = compute_b(n, q, inv_A, r);

  double r_dens = ldmvnorm(n, q, b, A, Q, Tau);

  arma::vec F_r = arma::normcdf(r, mu_r, arma::sqrt(Vr.diag()));
  
  Q1.reset();
  Q2.reset();
  Q3.reset();
  Q4.reset();

  invQ1.reset();
  invQ2.reset();
  invQ3.reset();
  invQ4.reset();
  
  kprod.reset();
  
  invQ.reset();
  
  Vr.reset();
  
  arma::vec probs = makeprobs(v);

  arma::ivec u = makeu(F_r,probs);
  
  double nu = 2.0;
  arma::mat R = 0.1 * arma::eye(q,q);
  arma::mat R_inv_Sigma = R * inv_Sigma;

  double lpA = -(nu + q + 1) / 2 * arma::log_det(Sigma).real() + std::log(Jacob_A(A)) -
    arma::sum(1/2  * nu * R_inv_Sigma.diag() - arma::log(A.diag()));

  // Define scaling factors
  double s_theta = 0.0;
  double s_beta = 0.0;
  arma::vec s_r = 0.0 * arma::ones<arma::vec>(nq);
  double s_vv = 0.0;
  double s_rhorho = 0.0;
  double s_etaeta = 0.0;
  double s_A = 0.0;

  // Copy inputs to matrices
  arma::vec m_theta = theta;
  arma::vec m_beta = beta;
  arma::vec m_r = r;
  arma::vec m_vv = vv;
  arma::vec m_rhorho = rhorho;
  arma::vec m_etaeta = etaeta;
  arma::vec m_A = A_vec;

  // // Create diagonal matrices
  arma::mat R_theta = arma::eye(n_atoms,n_atoms);
  arma::mat R_beta = arma::eye(p,p);
  arma::vec R_r = 0.01*arma::ones<arma::vec>(nq);
  arma::mat R_vv = arma::eye(n_atoms,n_atoms);
  arma::mat R_rhorho = arma::eye(q,q);
  arma::mat R_etaeta = arma::eye(q * 3, q * 3);
  arma::mat R_A = arma::eye(q * (q + 1) / 2, q * (q + 1) / 2);

  // Apply scaling
  arma::mat xi_theta = std::exp(s_theta) * R_theta;
  arma::mat xi_beta = std::exp(s_beta) * R_beta;
  arma::vec xi_r = arma::exp(s_r) % R_r;
  arma::mat xi_vv = std::exp(s_vv) * R_vv;
  arma::mat xi_rhorho = std::exp(s_rhorho) * R_rhorho;
  arma::mat xi_etaeta = std::exp(s_etaeta) * R_etaeta;
  arma::mat xi_A = std::exp(s_A) * R_A;

  // Initialize variables
  int accepttheta = 0;
  int acceptbeta = 0;
  arma::ivec acceptr(nq,arma::fill::zeros);
  int acceptv = 0;
  int acceptrho = 0;
  int accepteta = 0;
  int acceptA = 0;
  
  int g = 0;
  
  for (int iter = 0; iter < iterations; iter++) {
    // Print progress

    double progress = static_cast<double>(iter * 100) / iterations;
    if (progress >= 10.0 && progress <= 100.0 && fmod(progress, 10.0) == 0.0) {
      std::cout << "### Progress: " << progress << " %" << std::endl;
      std::cout << "### beta:     " << static_cast<double>(acceptbeta) / (iter) << std::endl;
      std::cout << "### theta:    " << static_cast<double>(accepttheta) / (iter) << std::endl;
      std::cout << "### r:        " << static_cast<double>(arma::sum(acceptr)/nq) / (iter) << std::endl;
      std::cout << "### v:        " << static_cast<double>(acceptv) / (iter) << std::endl;
      std::cout << "### rho:      " << static_cast<double>(acceptrho) / (iter) << std::endl;
      std::cout << "### eta:      " << static_cast<double>(accepteta) / (iter) << std::endl;
      std::cout << "### A:        " << static_cast<double>(acceptA) / (iter) << std::endl;
      std::cout << " " << std::endl;
    }

    double iter07 = pow(iter + 2, -0.7);

    /////////////////
    // Update beta //
    /////////////////

    arma::vec pro_beta = arma::mvnrnd(beta, xi_beta, 1);

    // Calculate X_pro_beta and X_beta using matrix multiplication
    arma::vec X_pro_beta = X * pro_beta;

    // Calculate MH_beta using vectorized operations
    double MH_beta = 0.0;

    for (int i = 0; i < p; i++){
      MH_beta += 0.5 * taub * (pow(beta(i),2.0) - pow(pro_beta(i),2.0));
    }

    for (int id = 0; id < nq; id++){
      MH_beta +=  y(id) * (X_pro_beta(id) - X_beta(id)) -
        E(id) * exp(theta(u(id) - 1)) * (exp(X_pro_beta(id)) - exp(X_beta(id)));
    }

    MH_beta = std::min(0.0, MH_beta);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_beta)) {
      beta = pro_beta;
      X_beta = X_pro_beta;
      acceptbeta++;
    }

    // Update s_beta
    s_beta += iter07 * (std::exp(MH_beta) - 0.234);
    
    // if (s_beta < - 4){
    //   s_beta = -4;
    // }

    // Calculate dbeta_2 and update m_beta in a single loop
    arma::vec dbeta_2 = beta - m_beta;
    m_beta += iter07 * dbeta_2;

    // Update R_beta and xi_beta using vectorized operations
    R_beta = (1 - iter07) * R_beta + iter07 * (dbeta_2 * dbeta_2.t());
    xi_beta = exp(s_beta) * R_beta;

    // Adjust diagonal of xi_beta using vectorized operations
    xi_beta.diag() += 1e-3 * xi_beta.diag();

    //////////////////
    // Update theta //
    //////////////////

    arma::vec pro_theta = arma::mvnrnd(theta, xi_theta, 1);

    // Calculate MH_theta
    double MH_theta = 0.0;

    for (int k = 0; k < n_atoms; k++){
      MH_theta += 0.5 * taus * (pow(theta(k),2.0) - pow(pro_theta(k),2.0));
    }

    for (int id = 0; id < nq; id++) {
      MH_theta += y(id) * (pro_theta(u(id) - 1) - theta(u(id) - 1)) -
        E(id) * std::exp(X_beta(id)) * 
        (std::exp(pro_theta(u(id) - 1)) - std::exp(theta(u(id) - 1)));
    }

    MH_theta = std::min(0.0, MH_theta);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_theta)) {
      theta = pro_theta;
      accepttheta++;
    }

    // Update s_theta
    s_theta += iter07 * (std::exp(MH_theta) - 0.234);
    
    // if (s_theta < - 4){
    //   s_theta = -4;
    // }

    // Calculate dtheta_2 and update m_theta in a single loop
    arma::vec dtheta_2 = theta - m_theta;
    m_theta += iter07 * dtheta_2;

    // Update R_theta and xi_theta using vectorized operations
    R_theta = (1 - iter07) * R_theta + iter07 * (dtheta_2 * dtheta_2.t());
    xi_theta = exp(s_theta) * R_theta;

    // Adjust diagonal of xi_theta using vectorized operations
    xi_theta.diag() += 1e-6 * xi_theta.diag();

    //////////////
    // Update r //
    //////////////

    arma::vec MH_r = arma::vec(nq, arma::fill::zeros);
    
    for (int id = 0; id < nq; id++) {
      
      arma::vec pro_r = r;
      arma::vec pro_Fr = F_r;
      arma::ivec pro_u = u;
      arma::mat pro_b = b;
      
      pro_r(id) = arma::randn<double>(arma::distr_param(r(id), std::sqrt(xi_r(id))));
      pro_Fr(id) = arma::normcdf(pro_r(id), 0.0, 1/std::sqrt(inv_Vr_diag(id)));
      pro_u(id) = makeuk(pro_Fr(id), probs);
      
      pro_b = compute_b(n, q, inv_A, pro_r);
      
      double pro_r_dens = ldmvnorm(n, q, pro_b, A, Q, Tau);
      
      MH_r(id) = pro_r_dens - r_dens +
        y(id) * (theta(pro_u(id) - 1) - theta(u(id) - 1)) -
        E(id) * std::exp(X_beta(id)) * 
        (std::exp(theta(pro_u(id) - 1)) - std::exp(theta(u(id) - 1)));
      
      MH_r(id) = std::min(0.0, MH_r(id));
      
      if (arma::randu() < std::exp(MH_r(id))) {
        r(id) = pro_r(id);
        F_r(id) = pro_Fr(id);
        u(id) = pro_u(id);
        b = pro_b;
        r_dens = pro_r_dens;
        acceptr(id)++;
      }
      
      // // Update s_r
      // s_r(id) += iter07 * (std::exp(MH_r(id)) - 0.234);
      // 
      // if (s_r(id) < -2) {
      //   s_r(id) = -2;
      // }
      // 
      // // Calculate dr_2 and update m_r in a single loop
      // double dr_2 = r(id) - m_r(id);
      // m_r(id) += iter07 * dr_2;
      // 
      // // Update R_r and xi_r using vectorized operations
      // R_r(id) = (1 - iter07) * R_r(id) + iter07 * dr_2 * dr_2;
      // xi_r(id) = std::exp(s_r(id)) * R_r(id);
      // 
      // // Adjust diagonal of xi_r using vectorized operations
      // xi_r(id) += 1e-3 * xi_r(id);
      
    }

    //////////////
    // Update v //
    //////////////

    arma::vec pro_vv = arma::mvnrnd(vv, xi_vv, 1);

    arma::vec pro_v = arma::vec(n_atoms, arma::fill::zeros);
    for (int k = 0; k < n_atoms; k++){
      pro_v(k) = inv_trans_par(pro_vv(k), 0.0, 1.0);
    }

    // Calculate probabilities
    arma::vec pro_probs = makeprobs(pro_v);

    // Calculate u
    arma::ivec pro_u = makeu(F_r, pro_probs);

    double MH_vv = arma::sum(alpha * (arma::log(1 - pro_v) - arma::log(1 - v)) + 
                             arma::log(pro_v) - arma::log(v));

    for (int id = 0; id < nq; id++) {
      MH_vv += y(id) * (theta(pro_u(id) - 1) - theta(u(id) - 1)) -
        E(id) * std::exp(X_beta(id)) * 
        (std::exp(theta(pro_u(id) - 1)) - std::exp(theta(u(id) - 1)));
    }

    MH_vv = std::min(0.0, MH_vv);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_vv)) {
      vv = pro_vv;
      v = pro_v;
      probs = pro_probs;
      u = pro_u;
      acceptv++;
    }

    // Update s_vv
    s_vv += iter07 * (std::exp(MH_vv) - 0.234);
    
    // if (s_vv < - 4){
    //   s_vv = -4;
    // }

    // Calculate dvv_2 and update m_vv in a single loop
    arma::vec dvv_2 = vv - m_vv;
    m_vv += iter07 * dvv_2;

    // Update R_vv and xi_vv using vectorized operations
    R_vv = (1 - iter07) * R_vv + iter07 * (dvv_2 * dvv_2.t());
    xi_vv = exp(s_vv) * R_vv;

    // Adjust diagonal of xi_vv using vectorized operations
    xi_vv.diag() += 1e-3 * xi_vv.diag();

    /////////////////
    // Update taus //
    /////////////////

    taus = arma::randg<double>(arma::distr_param(0.5 * n_atoms + c, 
                                                 1/(dot(theta,theta)/2 + d2)));

    ////////////////
    // Update rho //
    ////////////////

    arma::vec pro_rhorho = arma::mvnrnd(rhorho, xi_rhorho, 1);

    // Apply inverse transformation
    arma::vec pro_rho = arma::vec(q, arma::fill::zeros);
    for (int d = 0; d < q; d++){
      pro_rho(d) = inv_trans_par(pro_rhorho(d), 0.0, 1.0);
    }

    std::vector<arma::mat> pro_Tau = Tau_new(pro_rho, n, q, W);
    std::vector<arma::mat> pro_Q = Q_new(pro_rho, n, q, W, pro_Tau);

    arma::vec pro_inv_Vr_diag = compute_inv_Vr_diag(n, q, inv_A, pro_Q);

    double pro_r_dens = ldmvnorm(n, q, b, A, pro_Q, pro_Tau);

    double MH_rhorho = pro_r_dens - r_dens +
      arma::sum(arma::log(pro_rho) + arma::log(1 - pro_rho) - 
      arma::log(rho) - arma::log(1 - rho));

    MH_rhorho = std::min(0.0, MH_rhorho);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_rhorho)) {
      rhorho = pro_rhorho;
      rho = pro_rho;
      Tau = pro_Tau;
      Q = pro_Q;
      inv_Vr_diag = pro_inv_Vr_diag;
      r_dens = pro_r_dens;
      acceptrho++;
    }

    // Update s_rhorho
    s_rhorho += iter07 * (std::exp(MH_rhorho) - 0.234);
    
    // if (s_rhorho < -4){
    //   s_rhorho = -4;
    // }

    // Calculate drhorho_2 and update m_rhorho in a single loop
    arma::vec drhorho_2 = rhorho - m_rhorho;
    m_rhorho += iter07 * drhorho_2;

    // Update R_rhorho and xi_rhorho using vectorized operations
    R_rhorho = (1 - iter07) * R_rhorho + iter07 * (drhorho_2 * drhorho_2.t());
    xi_rhorho = exp(s_rhorho) * R_rhorho;

    // Adjust diagonal of xi_rhorho using vectorized operations
    xi_rhorho.diag() += 1e-3 * xi_rhorho.diag();

    ////////////////
    // Update eta //
    ////////////////
if (cvrts == "adj" || cvrts == "meanadj") {
    arma::vec pro_etaeta = arma::mvnrnd(etaeta, xi_etaeta, 1);

    arma::vec pro_eta1 = arma::vec(q, arma::fill::zeros);
    arma::vec pro_eta2 = arma::vec(q, arma::fill::zeros);
    arma::vec pro_eta3 = arma::vec(q, arma::fill::zeros);
    for (int d = 0; d < q; d++) {
      pro_eta1(d) = inv_trans_par(pro_etaeta(d), 0, M1);
      pro_eta2(d) = inv_trans_par(pro_etaeta(d + q), 0, M2);
      pro_eta3(d) = inv_trans_par(pro_etaeta(d + 2 * q), 0, M3);
    }

    // Concatenate pro_eta1, pro_eta2, and pro_eta3
    arma::vec pro_eta = arma::vec(3 * q, arma::fill::zeros);
    // Copy elements from pro_eta1, pro_eta2, and pro_eta3 to pro_eta
    pro_eta.subvec(0, q - 1) = pro_eta1;
    pro_eta.subvec(q, 2 * q - 1) = pro_eta2;
    pro_eta.subvec(2 * q, 3 * q - 1) = pro_eta3;

    // Initialize pro_W
    std::vector<arma::mat> pro_W(q);
    for (int d = 0; d < q; d++) {
      // Initialize as an n x n matrix filled with zeros
      pro_W[d] = arma::mat(n, n, arma::fill::zeros);
      for (int i = 0; i < Minc.n_rows; i++) {
        for (int j = 0; j < Minc.n_cols; j++) {
          if (Minc(i, j) == 1) {
            double value = std::exp(-Z1(i, j) * pro_eta1(d) -
                                    Z2(i, j) * pro_eta2(d) -
                                    Z3(i, j) * pro_eta3(d));
            pro_W[d](i, j) = (value >= 0.5) ? 1.0 : 0.0;
          }
        }
      }
    }

    pro_Tau = Tau_new(rho, n, q, pro_W);

    pro_Q = Q_new(rho, n, q, pro_W, pro_Tau);

    pro_inv_Vr_diag = compute_inv_Vr_diag(n, q, inv_A, pro_Q);

    pro_r_dens = ldmvnorm(n, q, b, A, pro_Q, pro_Tau);

    double MH_etaeta = pro_r_dens - r_dens +
      arma::sum(pro_etaeta - 2 * arma::log1p(arma::exp(pro_etaeta)) -
      etaeta + 2 * arma::log1p(arma::exp(etaeta)));

    MH_etaeta = std::min(0.0, MH_etaeta);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_etaeta)) {
      etaeta = pro_etaeta;
      eta = pro_eta;
      W = pro_W;
      Tau = pro_Tau;
      Q = pro_Q;
      inv_Vr_diag = pro_inv_Vr_diag;
      r_dens = pro_r_dens;
      accepteta++;
    }

    // Update s_etaeta
    s_etaeta += iter07 * (std::exp(MH_etaeta) - 0.234);

    // if (s_etaeta < - 4){
    //   s_etaeta = -4;
    // }

    // Calculate detaeta_2 and update m_etaeta in a single loop
    arma::vec detaeta_2 = etaeta - m_etaeta;
    m_etaeta += iter07 * detaeta_2;

    // Update R_etaeta and xi_etaeta using vectorized operations
    R_etaeta = (1 - iter07) * R_etaeta + iter07 * (detaeta_2 * detaeta_2.t());
    xi_etaeta = exp(s_etaeta) * R_etaeta;

    // Adjust diagonal of xi_etaeta using vectorized operations
    xi_etaeta.diag() += 1e-3 * xi_etaeta.diag();

  }

    //////////////
    // Update A //
    //////////////

    arma::vec pro_A_vec = arma::mvnrnd(A_vec, xi_A, 1);

    arma::mat pro_A(q, q, arma::fill::zeros);
    int idxproA = 0;
    for (int d = 0; d < q; d++) {
      for (int h = 0; h <= d; h++) {
        if (h == d) {
          pro_A(d, h) = std::exp(pro_A_vec(idxproA));
        } else {
          pro_A(d, h) = pro_A_vec(idxproA);
        }
        idxproA++;
      }
    }

    arma::mat inv_pro_A = arma::inv(arma::trimatl(pro_A));

    arma::mat pro_b = compute_b(n, q, inv_pro_A, r);
    
    pro_inv_Vr_diag = compute_inv_Vr_diag(n, q, inv_pro_A, Q);
    
    pro_r_dens = ldmvnorm(n, q, pro_b, pro_A, Q, Tau);

    arma::mat pro_Sigma = pro_A * pro_A.t();
    arma::mat inv_pro_Sigma = inv_pro_A.t() * inv_pro_A;

    arma::mat R_inv_pro_Sigma = R * inv_pro_Sigma;

    double pro_lpA = -(nu + q + 1) / 2 * arma::log_det(pro_Sigma).real() + 
      std::log(Jacob_A(pro_A)) -
      arma::sum(1/2  * nu * R_inv_pro_Sigma.diag() - arma::log(pro_A.diag()));

    double MH_A = pro_r_dens - r_dens + pro_lpA - lpA;

    MH_A = std::min(0.0, MH_A);

    if(arma::randu() < std::exp(MH_A)){
      A_vec = pro_A_vec;
      A = pro_A;
      inv_A = inv_pro_A;
      b = pro_b;
      inv_Vr_diag = pro_inv_Vr_diag;
      Sigma = pro_Sigma;
      inv_Sigma = inv_pro_Sigma;
      lpA = pro_lpA;
      r_dens = pro_r_dens;
      acceptA++;
    }

    // Update s_A
    s_A += iter07 * (std::exp(MH_A) - 0.234);

    if (s_A < -4){
      s_A = -4;
    }

    // Calculate dA_2 and update m_A in a single loop
    arma::vec dA_2 = A_vec - m_A;
    m_A += iter07 * dA_2;

    // Update R_A and xi_A using vectorized operations
    R_A = (1 - iter07) * R_A + iter07 * (dA_2 * dA_2.t());
    xi_A = exp(s_A) * R_A;

    // Adjust diagonal of xi_A using vectorized operations
    xi_A.diag() += 1e-3 * xi_A.diag();

    ////////////////////
    // record samples //
    ////////////////////

    if (iter >= burn & iter % thin == 0) {
      arma::vec phi = replicate_elements(theta,u);

      keepbeta.row(g) = beta.t();
      keepphi.row(g) = phi.t();
      keeptheta.row(g) = theta.t();
      keepu.row(g) = u.t();
      keeprho.row(g) = rho.t();
      keepv.row(g) = v.t();
      keepr.row(g) = r.t();
      keepFr.row(g) = F_r.t();
      keepeta.row(g) = eta.t();
      keeptaus[g] = taus;
      keepW1[g] = W[0];
      keepW2[g] = W[1];
      keepW3[g] = W[2];
      keepW4[g] = W[3];
      keepA.row(g) = A_vec.t();

      g++;
    }

  }

  output[0] = keepbeta;
  output[1] = keepphi;
  output[2] = keeptheta;
  output[3] = keepu;
  output[4] = keeprho;
  output[5] = keepv;
  output[6] = keepr;
  output[7] = keepFr;
  output[8] = keepeta;
  output[9] = keeptaus;
  output[10] = keepW1;
  output[11] = keepW2;
  output[12] = keepW3;
  output[13] = keepW4;
  output[14] = keepA;
  
  // Return something
  return output;
  
}
