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

std::vector<arma::mat> B_new(const arma::vec& rho, 
                             int n, 
                             int q, 
                             const std::vector<arma::mat>& W) {
  
  std::vector<arma::mat> B(q, arma::mat(n, n, arma::fill::zeros));
  arma::mat I = arma::eye(n,n);
  
  for (int d = 0; d < q; d++) {
    
    const arma::mat& Winc = arma::trimatl(W[d]);
    const arma::vec ns = arma::sum(Winc, 1);
    
    // B_d is strictly lower triangular
    for (int i = 1; i < n; i++) {
      for (int j = 0; j < i; j++) {
        if (Winc(i, j) == 1.0) {
          B[d](i, j) = rho(d) / (1 + (ns(i) - 1) * (rho(d) * rho(d)));
        }
      }
    }

  }
  
  return B;
}

std::vector<arma::mat> Q_new(int n, int q, 
                             const std::vector<arma::mat>& Tau, 
                             const std::vector<arma::mat>& B) {
  
  std::vector<arma::mat> Q(q, arma::mat(n, n, arma::fill::zeros));
  arma::mat I = arma::eye(n,n);
  
  for (int d = 0; d < q; d++) {
    
    arma::mat ImB_d =  I - B[d];
    
    Q[d] = ImB_d.t() * Tau[d] * ImB_d;
    
    Q[d] = 0.5 * (Q[d] + Q[d].t());
    
  }
  
  return Q;
}

std::vector<arma::mat> U_new(int n, int q, 
                             const std::vector<arma::mat>& Tau,
                             const std::vector<arma::mat>& B) {
  
  std::vector<arma::mat> U(q, arma::mat(n, n, arma::fill::zeros));
  arma::mat I = arma::eye(n,n);
  
  for (int d = 0; d < q; d++) {
    
    arma::mat Tausq_d = arma::diagmat(arma::sqrt(Tau[d].diag()));
    
    arma::mat ImB_d =  I - B[d];
    
    U[d] = ImB_d.t() * Tausq_d;
    
  }
  
  return U;
  
}

arma::vec compute_inv_Vr_diag(int n,
                              int q,
                              const std::vector<arma::mat>& Q) {
  
  arma::vec inv_Vr_diag(n * q, arma::fill::zeros);
  
  for (int d = 0; d < q; d++) {
    inv_Vr_diag.subvec(d * n, (d + 1) * n - 1) += Q[d].diag();
  }
  
  return inv_Vr_diag;
  
}

double ldmvnorm(int n, int q,
                const arma::vec& r,
                const std::vector<arma::mat>& U,
                const std::vector<arma::mat>& Tau,
                const arma::mat& L_dis) {
  
  double ldens = n/2* log(arma::det(L_dis));
  
  for(int d = 0; d < q; d++) {
    ldens += 0.5 * arma::sum(arma::log(Tau[d].diag())) - 0.5 * L_dis(d,d) * arma::dot(U[d].t() * r.subvec(d * n, (d + 1) * n - 1), U[d].t() * r.subvec(d * n, (d + 1) * n - 1)); 
    for(int h = 0; h < d; h++) {
      if(L_dis(d,h) == 1.0){
        ldens -= L_dis(d,h) * arma::dot(U[d].t() * r.subvec(d * n, (d + 1) * n - 1), U[h].t() * r.subvec(h * n, (h + 1) * n - 1));  
      }
    }
  }
  
  return ldens;
  
}

// Function to compute stick-breaking weights

arma::vec makeprobs(const arma::vec& v) {
  
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

arma::ivec makeu(const arma::vec& F_r, const arma::vec& probs) {
  
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

int makeuk(double F_r, const arma::vec& probs) {
  
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

List MADAGAR(const arma::vec& y, const arma::mat& X, const arma::mat& Z1,
             int q, const arma::mat& Minc, const arma::mat& W_dis,
             double alpha, int n_atoms,
             int runs, int burn, int thin) {
  
  List output(16);
  
  int nq = y.n_elem;
  int n = nq / q;
  int p = X.n_cols;
  
  int iterations = runs * thin + burn;
  
  arma::mat I = arma::eye(n,n);
  
  arma::mat keepbeta(runs, p);
  arma::mat keeptaud(runs, q);
  arma::mat keepphi(runs, nq);
  arma::mat keeptheta(runs, n_atoms);
  arma::imat keepu(runs, nq);
  arma::mat keeprho(runs, q);
  arma::vec keeptaus(runs, arma::fill::zeros);
  arma::mat keepv(runs, n_atoms);
  arma::mat keepr(runs, nq);
  arma::mat keepFr(runs, nq);
  arma::mat keepeta(runs, q);
  arma::vec keeprhodis(runs, arma::fill::zeros);
  
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
  
  arma::vec taud = arma::ones<arma::vec>(q);
  double taub = 0.0001;
  double taus = 1.0;
  
  double ce = 2.0;
  double de = 0.1;
  
  double c = 2.0;
  double d2 = 1.0; // 1.0;
  
  // Initialize v with n_atoms elements, each set to 0.1
  arma::vec v = arma::vec(n_atoms).fill(0.1);
  arma::vec vv = arma::log(v) - arma::log(1 - v);
  
  // Initialize rho with q elements, each set to 0.1
  arma::vec rho = arma::vec(q).fill(0.1);
  arma::vec rhorho = arma::log(rho) - arma::log(1 - rho);
  
  double M1 = calculateM(Z1, Minc);
  
  // Initialize eta with 3*q elements, each set to 0.1
  arma::vec eta = arma::vec(q).fill(0.1);
  arma::vec etaeta = arma::log(eta) - arma::log(M1 - eta);
  
  std::vector<arma::mat> W(q, Minc);
  
  std::vector<arma::mat> Tau = Tau_new(rho, n, q, W);
  std::vector<arma::mat> B = B_new(rho, n, q, W);
  std::vector<arma::mat> Q = Q_new(n, q, Tau, B);
  
  double rhodis = 0.1;
  double rhorhodis = std::log(1 + rhodis) - std::log(1 - rhodis);
  
  arma::mat D_dis = arma::diagmat(arma::sum(W_dis, 1));
  
  arma::mat L_dis = D_dis - rhodis * W_dis;
  
  std::vector<arma::mat> U = U_new(n, q, Tau, B);
  
  arma::vec inv_Vr_diag = compute_inv_Vr_diag(n, q, Q);
  
  arma::vec r = arma::mvnrnd(arma::vec(nq, arma::fill::zeros), arma::eye(nq, nq), 1);
  
  double r_dens = ldmvnorm(n, q, r, U, Tau, L_dis);
  
  arma::vec F_r = arma::normcdf(r, arma::vec(nq, arma::fill::zeros), arma::vec(nq, arma::fill::ones));
  
  arma::vec probs = makeprobs(v);
  
  arma::ivec u = makeu(F_r,probs);
  
  arma::vec phi = replicate_elements(theta,u);
  
  // Define scaling factors
  arma::vec s_r = 1.0 * arma::ones<arma::vec>(nq);
  double s_vv = 1.0;
  double s_rhorho = 1.0;
  double s_etaeta = 1.0;
  double s_rhorhodis = 1.0;
  
  // Copy inputs to matrices
  arma::vec m_r = r;
  arma::vec m_vv = vv;
  arma::vec m_rhorho = rhorho;
  arma::vec m_etaeta = etaeta;
  double m_rhorhodis = rhorhodis;
  
  // // Create diagonal matrices
  arma::vec R_r = 1.0 * arma::ones<arma::vec>(nq);
  arma::mat R_vv = 1.0 * arma::eye(n_atoms,n_atoms);
  arma::mat R_rhorho = 1.0 * arma::eye(q,q);
  arma::mat R_etaeta = 1.0 *  arma::eye(q, q);
  double R_rhorhodis = 1.0;
  
  // Apply scaling
  arma::vec xi_r = arma::exp(s_r) % R_r;
  arma::mat xi_vv = std::exp(s_vv) * R_vv;
  arma::mat xi_rhorho = std::exp(s_rhorho) * R_rhorho;
  arma::mat xi_etaeta = std::exp(s_etaeta) * R_etaeta;
  double xi_rhorhodis = std::exp(s_rhorhodis) * R_rhorhodis;
  
  // Initialize variables
  arma::ivec acceptr(nq,arma::fill::zeros);
  int acceptv = 0;
  int acceptrho = 0;
  int accepteta = 0;
  int acceptrhodis = 0;
  
  int g = 0;
  
  for (int iter = 0; iter < iterations; iter++) {
    // Print progress
    
    // std::cout << iter << std::endl;
    
    double progress = static_cast<double>(iter * 100) / iterations;
    if (progress >= 10.0 && progress <= 100.0 && fmod(progress, 10.0) == 0.0) {
      std::cout << "### Progress: " << progress << " %" << std::endl;
      std::cout << "### r:        " << static_cast<double>(arma::sum(acceptr)/nq) / (iter) << std::endl;
      std::cout << "### v:        " << static_cast<double>(acceptv) / (iter) << std::endl;
      std::cout << "### rho:      " << static_cast<double>(acceptrho) / (iter) << std::endl;
      std::cout << "### eta:      " << static_cast<double>(accepteta) / (iter) << std::endl;
      std::cout << "### rhodis:   " << static_cast<double>(acceptrhodis) / (iter) << std::endl;
      std::cout << " " << std::endl;
    }

    double iter07 = pow(iter + 2, -0.7);

    /////////////////
    // Update beta //
    /////////////////

    for (int d = 0; d < q; d++){

      arma::mat Sigmab = arma::inv_sympd(taud(d) *
        X.submat(d * n, d * 2, (d + 1) * n - 1, (d + 1) * 2 - 1).t() *
        X.submat(d * n, d * 2, (d + 1) * n - 1, (d + 1) * 2 - 1) +
        taub * arma::eye(2,2), arma::inv_opts::allow_approx);

      arma::vec mub =  Sigmab *
        taud(d) * X.submat(d * n, d * 2, (d + 1) * n - 1, (d + 1) * 2 - 1).t() *
        (y.subvec(d * n, (d + 1) * n - 1) - phi.subvec(d * n, (d + 1) * n - 1));

      beta.subvec(d * 2, (d + 1) * 2 - 1) = arma::mvnrnd(mub, Sigmab, 1);

      X_beta.subvec(d * n, (d + 1) * n - 1) =
        X.submat(d * n, d * 2, (d + 1) * n - 1, (d + 1) * 2 - 1) *
        beta.subvec(d * 2, (d + 1) * 2 - 1);

    }
    
    // std::cout << beta << std::endl;

    //////////////////
    // Update theta //
    //////////////////

    for (int k = 0; k < n_atoms; k++){

      double mm = 0.0;
      double MM = 0.0;

      for (int d = 0; d < q; d++){
        for (int i = 0; i < n; i++){
          if (u(d * n + i) == k + 1){
            mm += taud(d) * (y(d * n + i) - X_beta(d * n + i));
            MM += taud(d);
          }
        }
      }

      MM += taus;

      theta(k) = arma::randn<double>(arma::distr_param(mm/MM,1/std::sqrt(MM)));

    }

    // std::cout << theta << std::endl;
    
    phi = replicate_elements(theta, u);

    /////////////////
    // Update taud //
    /////////////////

    for (int d = 0; d < q; d++){
      taud(d) = arma::randg<double>(arma::distr_param(0.5 * n + ce,
                                    1/(dot(y.subvec(d * n, (d + 1) * n - 1) -
                                      X_beta.subvec(d * n, (d + 1) * n - 1) -
                                      phi.subvec(d * n, (d + 1) * n - 1),
                                      y.subvec(d * n, (d + 1) * n - 1) -
                                        X_beta.subvec(d * n, (d + 1) * n - 1) -
                                        phi.subvec(d * n, (d + 1) * n - 1))/2 + de)));
    }

    
    // std::cout << taud << std::endl;
    
    /////////////////
    // Update taus //
    /////////////////

    taus = arma::randg<double>(arma::distr_param(0.5 * n_atoms + c,
                                                 1/(dot(theta,theta)/2 + d2)));
    
    // std::cout << taus << std::endl;

    //////////////
    // Update r //
    //////////////

    arma::vec MH_r = arma::vec(nq, arma::fill::zeros);

    for (int id = 0; id < nq; id++) {

      arma::vec pro_r = r;
      arma::vec pro_Fr = F_r;
      arma::ivec pro_u = u;
      arma::vec pro_phi = phi;

      pro_r(id) = arma::randn<double>(arma::distr_param(r(id), std::sqrt(xi_r(id))));
      pro_Fr(id) = arma::normcdf(pro_r(id), 0.0, 1/std::sqrt(inv_Vr_diag(id)));
      pro_u(id) = makeuk(pro_Fr(id), probs);
      pro_phi(id) = theta(pro_u(id) - 1);

      double pro_r_dens = ldmvnorm(n, q, pro_r, U, Tau, L_dis);

      MH_r(id) = pro_r_dens - r_dens -
        0.5 * taud(id/n) *
        (std::pow(y(id) - X_beta(id) - theta(pro_u(id) - 1), 2) -
        std::pow(y(id) - X_beta(id) - theta(u(id) - 1), 2));

      MH_r(id) = std::min(0.0, MH_r(id));

      if (arma::randu() < std::exp(MH_r(id))) {
        r(id) = pro_r(id);
        F_r(id) = pro_Fr(id);
        u(id) = pro_u(id);
        phi(id) = pro_phi(id);
        r_dens = pro_r_dens;
        acceptr(id)++;
      }

      // Update s_r
      s_r(id) += iter07 * (std::exp(MH_r(id)) - 0.44);

      // Calculate dr_2 and update m_r in a single loop
      double dr_2 = r(id) - m_r(id);
      m_r(id) += iter07 * dr_2;

      // Update R_r and xi_r using vectorized operations
      R_r(id) = (1 - iter07) * R_r(id) + iter07 * dr_2 * dr_2;
      xi_r(id) = std::exp(s_r(id)) * R_r(id);

      // Adjust diagonal of xi_r using vectorized operations
      xi_r(id) += 1e-6 * xi_r(id);

    }
    
    // std::cout << "r" << std::endl;

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
    arma::vec pro_phi = replicate_elements(theta, pro_u);

    double MH_vv = arma::sum(alpha * (arma::log(1 - pro_v) - arma::log(1 - v)) +
                             arma::log(pro_v) - arma::log(v));

    for (int id = 0; id < nq; id++) {
      MH_vv -= 0.5 * taud(id/n) *
        (std::pow(y(id) - X_beta(id) - theta(pro_u(id) - 1), 2) -
        std::pow(y(id) - X_beta(id) - theta(u(id) - 1), 2));
    }

    MH_vv = std::min(0.0, MH_vv);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_vv)) {
      vv = pro_vv;
      v = pro_v;
      probs = pro_probs;
      u = pro_u;
      phi = pro_phi;
      acceptv++;
    }

    // Update s_vv
    s_vv += iter07 * (std::exp(MH_vv) - 0.234);

    // Calculate dvv_2 and update m_vv in a single loop
    arma::vec dvv_2 = vv - m_vv;
    m_vv += iter07 * dvv_2;

    // Update R_vv and xi_vv using vectorized operations
    R_vv = (1 - iter07) * R_vv + iter07 * (dvv_2 * dvv_2.t());
    xi_vv = exp(s_vv) * R_vv;

    // Adjust diagonal of xi_vv using vectorized operations
    xi_vv.diag() += 1e-6 * xi_vv.diag();
    
    // std::cout << v << std::endl;

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
    std::vector<arma::mat> pro_B = B_new(pro_rho, n, q, W);
    std::vector<arma::mat> pro_Q = Q_new(n, q, pro_Tau, pro_B);
    std::vector<arma::mat> pro_U = U_new(n, q, pro_Tau, pro_B);

    arma::vec pro_inv_Vr_diag = compute_inv_Vr_diag(n, q, pro_Q);

    double pro_r_dens = ldmvnorm(n, q, r, pro_U, pro_Tau, L_dis);

    double MH_rhorho = pro_r_dens - r_dens +
      arma::sum(arma::log(pro_rho) + arma::log(1 - pro_rho) -
      arma::log(rho) - arma::log(1 - rho));

    MH_rhorho = std::min(0.0, MH_rhorho);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_rhorho)) {
      rhorho = pro_rhorho;
      rho = pro_rho;
      Tau = pro_Tau;
      B = pro_B;
      Q = pro_Q;
      U = pro_U;
      inv_Vr_diag = pro_inv_Vr_diag;
      r_dens = pro_r_dens;
      acceptrho++;
    }

    // Update s_rhorho
    s_rhorho += iter07 * (std::exp(MH_rhorho) - 0.234);

    // Calculate drhorho_2 and update m_rhorho in a single loop
    arma::vec drhorho_2 = rhorho - m_rhorho;
    m_rhorho += iter07 * drhorho_2;

    // Update R_rhorho and xi_rhorho using vectorized operations
    R_rhorho = (1 - iter07) * R_rhorho + iter07 * (drhorho_2 * drhorho_2.t());
    xi_rhorho = exp(s_rhorho) * R_rhorho;

    // Adjust diagonal of xi_rhorho using vectorized operations
    xi_rhorho.diag() += 1e-6 * xi_rhorho.diag();
    
    //std::cout << MH_rhorho << std::endl;
    
    
    ////////////////
    // Update eta //
    ////////////////

    arma::vec pro_etaeta = arma::mvnrnd(etaeta, xi_etaeta, 1);

    arma::vec pro_eta = arma::vec(q, arma::fill::zeros);
    for (int d = 0; d < q; d++) {
      pro_eta(d) = inv_trans_par(pro_etaeta(d), 0, M1);
    }

    // Initialize pro_W
    std::vector<arma::mat> pro_W(q);
    for (int d = 0; d < q; d++) {
      pro_W[d] = arma::mat(n, n, arma::fill::zeros);
      for (int i = 0; i < Minc.n_rows; i++) {
        for (int j = 0; j < Minc.n_cols; j++) {
          if (Minc(i, j) == 1) {
            double value = std::exp(-Z1(i, j) * pro_eta(d));
            pro_W[d](i, j) = (value >= 0.5) ? 1.0 : 0.0;
          }
        }
      }
    }

    pro_Tau = Tau_new(rho, n, q, pro_W);
    pro_B = B_new(rho, n, q, pro_W);
    pro_Q = Q_new(n, q, pro_Tau, pro_B);
    pro_U = U_new(n, q, pro_Tau, pro_B);

    pro_inv_Vr_diag = compute_inv_Vr_diag(n, q, pro_Q);

    pro_r_dens = ldmvnorm(n, q, r, pro_U, pro_Tau, L_dis);

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
      B = pro_B;
      Q = pro_Q;
      U = pro_U;
      inv_Vr_diag = pro_inv_Vr_diag;
      r_dens = pro_r_dens;
      accepteta++;
    }

    // Update s_etaeta
    s_etaeta += iter07 * (std::exp(MH_etaeta) - 0.234);

    // Calculate detaeta_2 and update m_etaeta in a single loop
    arma::vec detaeta_2 = etaeta - m_etaeta;
    m_etaeta += iter07 * detaeta_2;

    // Update R_etaeta and xi_etaeta using vectorized operations
    R_etaeta = (1 - iter07) * R_etaeta + iter07 * (detaeta_2 * detaeta_2.t());
    xi_etaeta = exp(s_etaeta) * R_etaeta;

    // Adjust diagonal of xi_etaeta using vectorized operations
    xi_etaeta.diag() += 1e-6 * xi_etaeta.diag();

    // std::cout << eta << std::endl;
    
    ///////////////////
    // Update rhodis //
    ///////////////////

    double pro_rhorhodis = arma::randn(arma::distr_param(rhorhodis, std::sqrt(xi_rhorhodis)));

    // Apply inverse transformation

    double pro_rhodis = inv_trans_par(pro_rhorhodis, -1.0, 1.0);

    arma::mat pro_L_dis = D_dis - pro_rhodis * W_dis;

    pro_r_dens = ldmvnorm(n, q, r, U, Tau, pro_L_dis);

    double MH_rhorhodis = pro_r_dens - r_dens +
      pro_rhorhodis - 2 * std::log(1 + std::exp(pro_rhorhodis)) -
      rhorhodis + 2 * std::log(1 + std::exp(rhorhodis));

    MH_rhorhodis = std::min(0.0, MH_rhorhodis);

    if(arma::randu() < std::exp(MH_rhorhodis)){
      rhodis = pro_rhodis;
      rhorhodis = pro_rhorhodis;
      L_dis = pro_L_dis;
      r_dens = pro_r_dens;
      acceptrhodis++;
    }

    // Update s_rhorho
    s_rhorhodis += iter07 * (std::exp(MH_rhorhodis) - 0.234);

    // Calculate drhorho_2 and update m_rhorho in a single loop
    double drhorhodis_2 = rhorhodis - m_rhorhodis;
    m_rhorhodis += iter07 * drhorhodis_2;

    // Update R_rhorho and xi_rhorho using vectorized operations
    R_rhorhodis = (1 - iter07) * R_rhorhodis + iter07 * (drhorhodis_2 * drhorhodis_2);
    xi_rhorhodis = exp(s_rhorhodis) * R_rhorhodis;

    // Adjust diagonal of xi_rhorho using vectorized operations
    xi_rhorhodis += 1e-6 * xi_rhorhodis;
    
    // std::cout << rhodis << std::endl;

    ////////////////////
    // record samples //
    ////////////////////

    if (iter >= burn & iter % thin == 0) {

      keepbeta.row(g) = beta.t();
      keeptaud.row(g) = taud.t();
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
      keeprhodis[g] = rhodis;

      g++;

    }

  }

  output[0] = keepbeta;
  output[1] = keeptaud;
  output[2] = keepphi;
  output[3] = keeptheta;
  output[4] = keepu;
  output[5] = keeprho;
  output[6] = keepv;
  output[7] = keepr;
  output[8] = keepFr;
  output[9] = keepeta;
  output[10] = keeptaus;
  output[11] = keepW1;
  output[12] = keepW2;
  output[13] = keepW3;
  output[14] = keepW4;
  output[15] = keeprhodis;

  // Return something
  return output;

}