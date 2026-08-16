// [[Rcpp::depends(RcppArmadillo)]]

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <algorithm>
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

// The adjacency matrix is a deterministic function of eta for the two
// specifications that model adjacency.  Keeping this construction in one
// place prevents the initialized or proposed W from drifting away from eta.
std::vector<arma::mat> W_new_from_eta(
    const std::string& cvrts,
    int q,
    const arma::mat& Minc,
    const arma::mat& Z1,
    const arma::mat& Z2,
    const arma::mat& Z3,
    const arma::vec& eta) {
  std::vector<arma::mat> result(q, Minc);
  if (cvrts != "adj" && cvrts != "meanadj") {
    return result;
  }

  arma::vec eta1 = eta.subvec(0, q - 1);
  arma::vec eta2 = eta.subvec(q, 2 * q - 1);
  arma::vec eta3 = eta.subvec(2 * q, 3 * q - 1);
  for (int disease = 0; disease < q; disease++) {
    result[disease].zeros(Minc.n_rows, Minc.n_cols);
    for (arma::uword i = 0; i < Minc.n_rows; i++) {
      for (arma::uword j = 0; j < Minc.n_cols; j++) {
        if (Minc(i, j) == 1.0) {
          double similarity = std::exp(
            -Z1(i, j) * eta1(disease) -
            Z2(i, j) * eta2(disease) -
            Z3(i, j) * eta3(disease)
          );
          result[disease](i, j) = similarity >= 0.5 ? 1.0 : 0.0;
        }
      }
    }
  }
  return result;
}

int W_eta_mismatches(
    const std::vector<arma::mat>& W,
    const std::string& cvrts,
    int q,
    const arma::mat& Minc,
    const arma::mat& Z1,
    const arma::mat& Z2,
    const arma::mat& Z3,
    const arma::vec& eta) {
  std::vector<arma::mat> expected = W_new_from_eta(
    cvrts, q, Minc, Z1, Z2, Z3, eta
  );
  int mismatches = 0;
  for (int disease = 0; disease < q; disease++) {
    mismatches += static_cast<int>(arma::accu(W[disease] != expected[disease]));
  }
  return mismatches;
}

// Diagonals of the disease-specific DAGAR covariance matrices Q_d^{-1}.
// These are marginal variances, not reciprocals of diagonal precision terms.
std::vector<arma::vec> compute_Q_marginal_variances(
    const std::vector<arma::mat>& Q) {
  std::vector<arma::vec> result(Q.size());
  for (arma::uword disease = 0; disease < Q.size(); disease++) {
    arma::mat inverse_Q;
    bool ok = arma::inv_sympd(inverse_Q, Q[disease]);
    if (!ok || !inverse_Q.is_finite()) {
      stop(
        "Could not compute the marginal DAGAR covariance for disease %d.",
        static_cast<int>(disease + 1)
      );
    }
    result[disease] = inverse_Q.diag();
    if (arma::any(result[disease] <= 0.0)) {
      stop(
        "Non-positive marginal DAGAR variance for disease %d.",
        static_cast<int>(disease + 1)
      );
    }
  }
  return result;
}

// For r = (A kron I) epsilon with independent disease-specific
// epsilon_d ~ N(0, Q_d^{-1}), the marginal variance of r_d is
// sum_{h <= d} A_{dh}^2 diag(Q_h^{-1}).
arma::vec combine_marginal_Vr_diag(
    int n,
    int q,
    const arma::mat& A,
    const std::vector<arma::vec>& Q_marginal_variances) {
  arma::vec result(n * q, arma::fill::zeros);
  for (int disease = 0; disease < q; disease++) {
    arma::vec current(n, arma::fill::zeros);
    for (int component = 0; component <= disease; component++) {
      current +=
        A(disease, component) * A(disease, component) *
        Q_marginal_variances[component];
    }
    result.subvec(disease * n, (disease + 1) * n - 1) = current;
  }
  if (!result.is_finite() || arma::any(result <= 0.0)) {
    stop("The reconstructed marginal latent variances are invalid.");
  }
  return result;
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
  double logistic;
  if (x >= 0.0) {
    logistic = 1.0 / (1.0 + std::exp(-x));
  } else {
    double exp_x = std::exp(x);
    logistic = exp_x / (1.0 + exp_x);
  }
  if (logistic <= 0.0) {
    logistic = std::nextafter(0.0, 1.0);
  } else if (logistic >= 1.0) {
    logistic = std::nextafter(1.0, 0.0);
  }
  return lb + (ub - lb) * logistic;
}

// Numerically stable log(1 + exp(x)), used in transformed-parameter
// Jacobians for the optional joint rho--eta proposal.
double log1pexp_stable(double x) {
  if (x > 0.0) {
    return x + std::log1p(std::exp(-x));
  }
  return std::log1p(std::exp(x));
}

double finite_log_acceptance(double value) {
  return std::isfinite(value) ? std::min(0.0, value) : -arma::datum::inf;
}

// Draw a symmetric Gaussian random-walk proposal after explicitly
// symmetrising its adaptive covariance.  Tiny diagonal jitter is used only
// when roundoff has made that covariance fail Cholesky factorisation.  The
// same covariance applies in both proposal directions, so the MH ratio is
// unchanged.
arma::vec safe_random_walk(
    const arma::vec& current,
    const arma::mat& covariance,
    const std::string& block_name) {
  arma::mat symmetric = 0.5 * (covariance + covariance.t());
  if (!symmetric.is_finite()) {
    stop(
      "Non-finite adaptive covariance in the %s proposal.",
      block_name.c_str()
    );
  }

  arma::mat cholesky;
  double average_variance = arma::mean(arma::abs(symmetric.diag()));
  double jitter = std::max(1e-12, average_variance * 1e-12);
  bool ok = arma::chol(cholesky, symmetric, "lower");
  for (int attempt = 0; !ok && attempt < 8; attempt++) {
    arma::mat regularised = symmetric;
    regularised.diag() += jitter;
    ok = arma::chol(cholesky, regularised, "lower");
    jitter *= 10.0;
  }
  if (!ok) {
    stop(
      "Could not regularise the adaptive covariance in the %s proposal.",
      block_name.c_str()
    );
  }
  return current + cholesky * arma::randn<arma::vec>(current.n_elem);
}

arma::vec safe_zero_mean_gaussian(
    const arma::mat& covariance,
    const std::string& block_name) {
  arma::vec zero(covariance.n_rows, arma::fill::zeros);
  return safe_random_walk(zero, covariance, block_name);
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

// Stable log density contribution for the inverse-Wishart covariance prior in
// the unconstrained Cholesky parameterisation. For Sigma = A A',
// log|Sigma| = 2 sum log(A_ii), and the Jacobian can be evaluated directly on
// the log scale. This avoids determinant and product underflow when a diagonal
// Cholesky coordinate is very small.
double log_density_A(const arma::mat& A,
                     const arma::mat& inv_Sigma,
                     const arma::mat& R,
                     double nu) {
  int q = A.n_cols;
  arma::vec log_diagonal = arma::log(A.diag());
  if (!log_diagonal.is_finite() || !inv_Sigma.is_finite()) {
    return -arma::datum::inf;
  }

  double log_det_sigma = 2.0 * arma::sum(log_diagonal);
  double log_cholesky_jacobian = q * std::log(2.0);
  for (int row = 0; row < q; row++) {
    log_cholesky_jacobian += (q - row) * log_diagonal(row);
  }
  double trace_term = arma::trace(R * inv_Sigma);
  double result =
    -(nu + q + 1.0) / 2.0 * log_det_sigma +
    log_cholesky_jacobian -
    0.5 * nu * trace_term +
    arma::sum(log_diagonal);
  return std::isfinite(result) ? result : -arma::datum::inf;
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

// Difference in the Poisson log likelihood caused only by changing cluster
// allocations. Terms that do not depend on u cancel from the MH ratio.
double allocation_log_likelihood_difference(
    const arma::vec& y,
    const arma::vec& E,
    const arma::vec& X_beta,
    const arma::vec& theta,
    const arma::ivec& proposed_u,
    const arma::ivec& current_u) {
  double difference = 0.0;
  for (arma::uword id = 0; id < y.n_elem; id++) {
    double proposed_phi = theta(proposed_u(id) - 1);
    double current_phi = theta(current_u(id) - 1);
    difference +=
      y(id) * (proposed_phi - current_phi) -
      E(id) * std::exp(X_beta(id)) *
        (std::exp(proposed_phi) - std::exp(current_phi));
  }
  return difference;
}

// F_r is a deterministic marginal-Gaussian transform of r under the current
// graph/covariance state, and u is then deterministic given F_r and the stick
// weights. A retained draw outside this support is invalid. Fail immediately
// instead of silently writing an inconsistent state.
void assert_state_consistency(
    const arma::vec& r,
    const arma::vec& marginal_Vr_diag,
    const arma::vec& F_r,
    const arma::ivec& u,
    const arma::vec& probs,
    double tolerance,
    int chain_id,
    int iteration,
    double& maximum_F_gap_seen,
    int& maximum_u_mismatches_seen) {
  if (!marginal_Vr_diag.is_finite() ||
      arma::any(marginal_Vr_diag <= 0.0)) {
    stop(
      "Chain %d iteration %d has an invalid latent marginal variance.",
      chain_id,
      iteration
    );
  }
  arma::vec expected_F = arma::normcdf(
    r,
    arma::vec(r.n_elem, arma::fill::zeros),
    arma::sqrt(marginal_Vr_diag)
  );
  if (!expected_F.is_finite()) {
    stop(
      "Chain %d iteration %d produced non-finite reconstructed F_r.",
      chain_id,
      iteration
    );
  }
  double maximum_gap = arma::max(arma::abs(expected_F - F_r));
  arma::ivec expected_u = makeu(expected_F, probs);
  int mismatches = static_cast<int>(arma::accu(expected_u != u));
  maximum_F_gap_seen = std::max(maximum_F_gap_seen, maximum_gap);
  maximum_u_mismatches_seen = std::max(
    maximum_u_mismatches_seen,
    mismatches
  );
  if (maximum_gap > tolerance || mismatches > 0) {
    stop(
      "Chain %d iteration %d failed latent-state consistency: "
      "max |F_r - F_r(reconstructed)| = %.17g; u mismatches = %d.",
      chain_id,
      iteration,
      maximum_gap,
      mismatches
    );
  }
}

void assert_W_eta_consistency(
    const std::vector<arma::mat>& W,
    const std::string& cvrts,
    int q,
    const arma::mat& Minc,
    const arma::mat& Z1,
    const arma::mat& Z2,
    const arma::mat& Z3,
    const arma::vec& eta,
    int chain_id,
    int iteration,
    int& maximum_mismatches_seen) {
  int mismatches = W_eta_mismatches(
    W, cvrts, q, Minc, Z1, Z2, Z3, eta
  );
  maximum_mismatches_seen = std::max(maximum_mismatches_seen, mismatches);
  if (mismatches > 0) {
    stop(
      "Chain %d iteration %d failed graph-state consistency: "
      "W has %d entries inconsistent with eta.",
      chain_id,
      iteration,
      mismatches
    );
  }
}

// Restartable segment of the MADAGAR chain. Priors and latent
// Gaussian terms remain untempered; only the Poisson likelihood is raised to
// inverse_temperature. This makes complete-state swaps between independently
// updated replicas an exact parallel-tempering transition.
// [[Rcpp::export]]
List MADAGAR_tempered_segment(arma::vec y, arma::mat X, arma::mat Z1,
                           arma::mat Z2, arma::mat Z3, arma::vec E,
                           std::string cvrts, int q, arma::mat Winc,
                           arma::mat Minc, double alpha, int n_atoms,
                           int runs, int burn, int thin,
                           int chain_seed, int chain_id, int worker_pid,
                           std::string initialization_regime = "original",
                           std::string progress_log_path = "",
                           bool print_progress = true,
                           double state_consistency_tolerance = 1e-12,
                           double inverse_temperature = 1.0,
                           Nullable<List> supplied_state = R_NilValue,
                           Nullable<List> supplied_proposal = R_NilValue,
                           int adaptation_offset = 0,
                           bool adapt_proposals = true) {

  if (runs <= 0) {
    stop("runs must be positive.");
  }
  if (burn < 0) {
    stop("burn must be non-negative.");
  }
  if (thin <= 0) {
    stop("thin must be positive.");
  }
  if (chain_seed <= 0) {
    stop("chain_seed must be positive.");
  }
  if (chain_id <= 0) {
    stop("chain_id must be positive.");
  }
  if (q <= 0 || y.n_elem == 0 || y.n_elem % q != 0) {
    stop("length(y) must be positive and divisible by q.");
  }
  if (n_atoms < 2) {
    stop("n_atoms must be at least 2.");
  }
  if (!std::isfinite(alpha) || alpha <= 0.0) {
    stop("alpha must be positive and finite.");
  }
  if (!std::isfinite(state_consistency_tolerance) ||
      state_consistency_tolerance <= 0.0) {
    stop("state_consistency_tolerance must be positive and finite.");
  }
  if (!std::isfinite(inverse_temperature) ||
      inverse_temperature <= 0.0 || inverse_temperature > 1.0) {
    stop("inverse_temperature must be in (0, 1].");
  }
  if (adaptation_offset < 0) {
    stop("adaptation_offset must be non-negative.");
  }
  if (cvrts != "adj" && cvrts != "mean" && cvrts != "meanadj") {
    stop("cvrts must be one of 'adj', 'mean', or 'meanadj'.");
  }

  Function set_seed = Environment::base_env()["set.seed"];
  set_seed(chain_seed);

  int nq = y.n_elem;
  int n = nq / q;
  int p = X.n_cols;
  if (X.n_rows != static_cast<arma::uword>(nq) ||
      E.n_elem != static_cast<arma::uword>(nq)) {
    stop("X and E must have the same observation dimension as y.");
  }
  if (!y.is_finite() || !X.is_finite() || !E.is_finite() ||
      arma::any(y < 0.0) || arma::any(E <= 0.0)) {
    stop("y, X, and E must be finite; y must be non-negative and E positive.");
  }
  if (Minc.n_rows != static_cast<arma::uword>(n) ||
      Minc.n_cols != static_cast<arma::uword>(n) ||
      Winc.n_rows != static_cast<arma::uword>(n) ||
      Winc.n_cols != static_cast<arma::uword>(n) ||
      Z1.n_rows != static_cast<arma::uword>(n) ||
      Z1.n_cols != static_cast<arma::uword>(n) ||
      Z2.n_rows != static_cast<arma::uword>(n) ||
      Z2.n_cols != static_cast<arma::uword>(n) ||
      Z3.n_rows != static_cast<arma::uword>(n) ||
      Z3.n_cols != static_cast<arma::uword>(n)) {
    stop("Minc, Winc, Z1, Z2, and Z3 must all be n by n matrices.");
  }
  if (!Minc.is_finite() ||
      !arma::approx_equal(Minc, Minc.t(), "absdiff", 1e-12)) {
    stop("Minc must be finite and symmetric.");
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (Minc(i, j) == 1.0 &&
          (!std::isfinite(Z1(i, j)) ||
           !std::isfinite(Z2(i, j)) ||
           !std::isfinite(Z3(i, j)))) {
        stop("Adjacency covariates must be finite on every physical edge.");
      }
    }
  }

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

  std::vector<arma::mat> sumW(q);
  for (int disease = 0; disease < q; disease++) {
    sumW[disease] = arma::mat(n, n, arma::fill::zeros);
  }
  arma::imat keepWcardinality(runs, q, arma::fill::zeros);

  // initial values

  // Initialize theta with n_atoms elements, each set to 0.0
  arma::vec theta(n_atoms, arma::fill::zeros);

  // Initialize beta with p elements, each set to 0.0
  arma::vec beta(p, arma::fill::zeros);

  double taub = 0.001;
  double taus = 1.0;

  double c = 2.0;
  double d2 = 1.0; // 1.0;

  // Initialize v with n_atoms elements, each set to 0.1
  arma::vec v = arma::vec(n_atoms).fill(0.1);

  // Match sampler.cpp: initialize rho with q elements, each set to 0.1.
  arma::vec rho = arma::vec(q).fill(0.1);

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

  std::vector<arma::mat> W(q,Minc);

  arma::vec A_vec = arma::randn<arma::vec>(q * (q + 1) / 2, arma::distr_param(0.0, 0.015));

  if (initialization_regime == "original") {
    // Match the original sampler.cpp initialization.
  } else if (initialization_regime == "dense_adjacency") {
    rho.fill(0.85);
    if (cvrts == "adj" || cvrts == "meanadj") {
      eta = 0.05 * M;
    }
  } else if (initialization_regime == "sparse_adjacency") {
    rho.fill(0.20);
    if (cvrts == "adj" || cvrts == "meanadj") {
      eta = 0.85 * M;
    }
  } else if (initialization_regime == "random_prior_like") {
    theta = 0.25 * arma::randn<arma::vec>(n_atoms);
    beta = 0.10 * arma::randn<arma::vec>(p);
    taus = arma::randg<double>(arma::distr_param(c, 1.0 / d2));
    v = 0.05 + 0.90 * arma::randu<arma::vec>(n_atoms);
    rho = 0.05 + 0.90 * arma::randu<arma::vec>(q);
    if (cvrts == "adj" || cvrts == "meanadj") {
      eta = M % (0.05 + 0.90 * arma::randu<arma::vec>(3 * q));
    }
    A_vec = arma::randn<arma::vec>(
      q * (q + 1) / 2, arma::distr_param(0.0, 0.05)
    );
  } else if (initialization_regime == "strong_cross_dependence") {
    rho.fill(0.60);
    if (cvrts == "adj" || cvrts == "meanadj") {
      eta = 0.25 * M;
    }
    int index = 0;
    for (int row = 0; row < q; row++) {
      for (int column = 0; column <= row; column++) {
        A_vec(index) = (row == column) ? std::log(1.0) : 0.35;
        index++;
      }
    }
  } else if (initialization_regime == "weak_cross_dependence") {
    rho.fill(0.35);
    if (cvrts == "adj" || cvrts == "meanadj") {
      eta = 0.35 * M;
    }
    int index = 0;
    for (int row = 0; row < q; row++) {
      for (int column = 0; column <= row; column++) {
        A_vec(index) = (row == column) ? std::log(1.0) : 0.01;
        index++;
      }
    }
  } else {
    stop("Unknown initialization_regime.");
  }

  // Keep the regime-defining rho, eta, graph, and covariance choices while
  // dispersing the remaining posterior blocks as well. The original regime is
  // left exactly as in sampler.cpp, and random_prior_like already initializes
  // every block separately above.
  if (initialization_regime != "original" &&
      initialization_regime != "random_prior_like") {
    theta += 0.15 * arma::randn<arma::vec>(n_atoms);
    beta += 0.05 * arma::randn<arma::vec>(p);
    taus = std::exp(
      arma::randn<double>(arma::distr_param(0.0, 0.35))
    );
    v = arma::clamp(
      v + 0.08 * arma::randn<arma::vec>(n_atoms),
      0.02,
      0.40
    );
  }

  // A tempering driver swaps and restarts complete model states. All
  // deterministic quantities (W, Q, F_r, and u) are reconstructed below so
  // stale cached values can never enter a restarted transition.
  bool has_supplied_state = supplied_state.isNotNull();
  arma::vec supplied_r;
  if (has_supplied_state) {
    List restart_state = as<List>(supplied_state);
    const char* required_state_names[] = {
      "beta", "theta", "tau", "V", "rho", "eta", "A", "r"
    };
    for (const char* name : required_state_names) {
      if (!restart_state.containsElementNamed(name)) {
        stop("supplied_state is missing '%s'.", name);
      }
    }
    beta = as<arma::vec>(restart_state["beta"]);
    theta = as<arma::vec>(restart_state["theta"]);
    taus = as<double>(restart_state["tau"]);
    v = as<arma::vec>(restart_state["V"]);
    rho = as<arma::vec>(restart_state["rho"]);
    eta = as<arma::vec>(restart_state["eta"]);
    A_vec = as<arma::vec>(restart_state["A"]);
    supplied_r = as<arma::vec>(restart_state["r"]);
    if (beta.n_elem != static_cast<arma::uword>(p) ||
        theta.n_elem != static_cast<arma::uword>(n_atoms) ||
        v.n_elem != static_cast<arma::uword>(n_atoms) ||
        rho.n_elem != static_cast<arma::uword>(q) ||
        eta.n_elem != static_cast<arma::uword>(3 * q) ||
        A_vec.n_elem != static_cast<arma::uword>(q * (q + 1) / 2) ||
        supplied_r.n_elem != static_cast<arma::uword>(nq)) {
      stop("supplied_state has an incompatible parameter dimension.");
    }
  }

  if (!theta.is_finite() || !beta.is_finite() ||
      !v.is_finite() || !rho.is_finite() ||
      !eta.is_finite() || !A_vec.is_finite()) {
    stop("Initialization produced a non-finite parameter value.");
  }
  if (taus <= 0.0 || arma::any(v <= 0.0) || arma::any(v >= 1.0) ||
      arma::any(rho <= 0.0) || arma::any(rho >= 1.0)) {
    stop("Initialization violated a positivity or unit-interval constraint.");
  }

  if (cvrts == "adj" || cvrts == "meanadj") {
    if (M1 <= 0.0 || M2 <= 0.0 || M3 <= 0.0) {
      stop("Covariate dissimilarity bounds must be positive.");
    }
    for (int index = 0; index < 3 * q; index++) {
      eta(index) = std::max(
        1e-8,
        std::min(eta(index), M(index) - 1e-8)
      );
    }
    if (arma::any(eta <= 0.0) || arma::any(eta >= M)) {
      stop("Initialization violated an eta constraint.");
    }

    // This must also be done for the original-style regime: sampler.cpp used
    // eta = 0.1 but initially left W equal to Minc, which is not generally the
    // graph implied by eta.
    W = W_new_from_eta(cvrts, q, Minc, Z1, Z2, Z3, eta);
  }

  arma::vec vv = arma::log(v) - arma::log(1 - v);
  arma::vec rhorho = arma::log(rho) - arma::log(1 - rho);
  arma::vec etaeta = arma::log(eta) - arma::log(M - eta);
  arma::vec X_beta = X * beta;

  std::vector<arma::mat> Tau = Tau_new(rho, n, q, W);

  std::vector<arma::mat> Q = Q_new(rho, n, q, W, Tau);

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

  if (!A.is_finite() || arma::any(A.diag() <= 0.0)) {
    stop("Initialization produced an invalid Cholesky factor.");
  }

  arma::mat inv_A = arma::inv(arma::trimatl(A));

  arma::mat inv_Sigma = inv_A.t() * inv_A;

  arma::vec mu_r = arma::vec(nq, arma::fill::zeros);

  arma::mat kprod = arma::kron(A,I);

  // Initialize block diagonal matrix
  arma::mat invQ = arma::zeros(n * q, n * q);
  for (int disease = 0; disease < q; disease++) {
    arma::mat inverse_Q;
    bool inverse_ok = arma::inv_sympd(inverse_Q, Q[disease]);
    if (!inverse_ok || !inverse_Q.is_finite()) {
      stop(
        "Could not invert the initial DAGAR precision for disease %d.",
        disease + 1
      );
    }
    invQ.submat(
      disease * n,
      disease * n,
      (disease + 1) * n - 1,
      (disease + 1) * n - 1
    ) = inverse_Q;
  }

  arma::mat Vr = kprod * invQ * kprod.t();
  Vr = 0.5 * (Vr + Vr.t());

  std::vector<arma::vec> Q_marginal_variances =
    compute_Q_marginal_variances(Q);
  arma::vec marginal_Vr_diag = combine_marginal_Vr_diag(
    n, q, A, Q_marginal_variances
  );
  double initialization_variance_gap = arma::max(
    arma::abs(marginal_Vr_diag - Vr.diag())
  );
  if (!std::isfinite(initialization_variance_gap) ||
      initialization_variance_gap > 1e-8) {
    stop(
      "Marginal variance reconstruction failed at initialization: %.17g.",
      initialization_variance_gap
    );
  }

  arma::vec r = has_supplied_state ? supplied_r :
    safe_zero_mean_gaussian(Vr, "initial latent-r");

  arma::mat b = compute_b(n, q, inv_A, r);

  double r_dens = ldmvnorm(n, q, b, A, Q, Tau);

  arma::vec F_r = arma::normcdf(
    r,
    mu_r,
    arma::sqrt(marginal_Vr_diag)
  );

  kprod.reset();

  invQ.reset();

  Vr.reset();

  arma::vec probs = makeprobs(v);

  arma::ivec u = makeu(F_r,probs);

  if (!r.is_finite() || !F_r.is_finite() ||
      arma::any(u < 1) || arma::any(u > n_atoms)) {
    stop("Initialization produced an invalid latent state or allocation.");
  }
  for (int disease = 0; disease < q; disease++) {
    if (!W[disease].is_finite() ||
        arma::any(arma::vectorise(W[disease]) < 0.0) ||
        arma::any(arma::vectorise(W[disease]) > 1.0)) {
      stop("Initialization produced an invalid adjacency state.");
    }
  }

  List initial_W(q);
  for (int disease = 0; disease < q; disease++) {
    initial_W[disease] = W[disease];
  }

  List initial_values = List::create(
    _["chain_id"] = chain_id,
    _["seed"] = chain_seed,
    _["initialization_regime"] = initialization_regime,
    _["beta"] = beta,
    _["theta"] = theta,
    _["phi"] = replicate_elements(theta, u),
    _["tau"] = taus,
    _["V"] = v,
    _["rho"] = rho,
    _["eta"] = eta,
    _["W"] = initial_W,
    _["A"] = A_vec,
    _["r"] = r,
    _["F_r"] = F_r,
    _["u"] = u
  );

  double nu = 2.0;
  arma::mat R = 0.1 * arma::eye(q,q);
  double lpA = log_density_A(A, inv_Sigma, R, nu);

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

  // Proposal adaptation belongs to a temperature, not to the model state.
  // The R driver therefore keeps this object at its temperature when model
  // states are swapped and passes it back on the next segment.
  if (supplied_proposal.isNotNull()) {
    List proposal = as<List>(supplied_proposal);
    const char* required_proposal_names[] = {
      "s_theta", "s_beta", "s_r", "s_vv", "s_rhorho", "s_etaeta", "s_A",
      "m_theta", "m_beta", "m_r", "m_vv", "m_rhorho", "m_etaeta", "m_A",
      "R_theta", "R_beta", "R_r", "R_vv", "R_rhorho", "R_etaeta", "R_A"
    };
    for (const char* name : required_proposal_names) {
      if (!proposal.containsElementNamed(name)) {
        stop("supplied_proposal is missing '%s'.", name);
      }
    }
    s_theta = as<double>(proposal["s_theta"]);
    s_beta = as<double>(proposal["s_beta"]);
    s_r = as<arma::vec>(proposal["s_r"]);
    s_vv = as<double>(proposal["s_vv"]);
    s_rhorho = as<double>(proposal["s_rhorho"]);
    s_etaeta = as<double>(proposal["s_etaeta"]);
    s_A = as<double>(proposal["s_A"]);
    m_theta = as<arma::vec>(proposal["m_theta"]);
    m_beta = as<arma::vec>(proposal["m_beta"]);
    m_r = as<arma::vec>(proposal["m_r"]);
    m_vv = as<arma::vec>(proposal["m_vv"]);
    m_rhorho = as<arma::vec>(proposal["m_rhorho"]);
    m_etaeta = as<arma::vec>(proposal["m_etaeta"]);
    m_A = as<arma::vec>(proposal["m_A"]);
    R_theta = as<arma::mat>(proposal["R_theta"]);
    R_beta = as<arma::mat>(proposal["R_beta"]);
    R_r = as<arma::vec>(proposal["R_r"]);
    R_vv = as<arma::mat>(proposal["R_vv"]);
    R_rhorho = as<arma::mat>(proposal["R_rhorho"]);
    R_etaeta = as<arma::mat>(proposal["R_etaeta"]);
    R_A = as<arma::mat>(proposal["R_A"]);
    xi_theta = std::exp(s_theta) * R_theta;
    xi_beta = std::exp(s_beta) * R_beta;
    if (s_r.n_elem != static_cast<arma::uword>(nq) ||
        m_r.n_elem != static_cast<arma::uword>(nq) ||
        R_r.n_elem != static_cast<arma::uword>(nq) ||
        !s_r.is_finite() || !m_r.is_finite() || !R_r.is_finite() ||
        arma::any(R_r <= 0.0)) {
      stop("The supplied scalar-r proposal state is invalid.");
    }
    xi_r = arma::exp(s_r) % R_r;
    for (int id = 0; id < nq; id++) {
      xi_r(id) = std::max(1e-12, xi_r(id));
    }
    xi_vv = std::exp(s_vv) * R_vv;
    xi_rhorho = std::exp(s_rhorho) * R_rhorho;
    xi_etaeta = std::exp(s_etaeta) * R_etaeta;
    xi_A = std::exp(s_A) * R_A;
  }

  // Initialize variables
  int accepttheta = 0;
  int acceptbeta = 0;
  arma::ivec acceptr(nq,arma::fill::zeros);
  int acceptv = 0;
  int acceptrho = 0;
  int accepteta = 0;
  int acceptA = 0;
  double proposed_rho_u_changes = 0.0;
  double accepted_rho_u_changes = 0.0;
  double proposed_eta_u_changes = 0.0;
  double accepted_eta_u_changes = 0.0;
  double proposed_A_u_changes = 0.0;
  double accepted_A_u_changes = 0.0;
  int invariant_checks = 0;
  double maximum_F_gap_seen = 0.0;
  int maximum_u_mismatches_seen = 0;
  int maximum_W_eta_mismatches_seen = 0;

  int g = 0;

  std::ofstream progress_log;
  bool write_progress_log = progress_log_path.size() > 0;
  if (write_progress_log) {
    progress_log.open(progress_log_path.c_str(), std::ios::out | std::ios::app);
    if (!progress_log.is_open()) {
      stop("Could not open progress log file.");
    }
  }

  std::time_t start_time = std::time(nullptr);
  int next_progress_percent = 5;

  assert_state_consistency(
    r,
    marginal_Vr_diag,
    F_r,
    u,
    probs,
    state_consistency_tolerance,
    chain_id,
    -1,
    maximum_F_gap_seen,
    maximum_u_mismatches_seen
  );
  assert_W_eta_consistency(
    W, cvrts, q, Minc, Z1, Z2, Z3, eta,
    chain_id, -1, maximum_W_eta_mismatches_seen
  );
  invariant_checks++;

  for (int iter = 0; iter < iterations; iter++) {
    int completed_iterations = iter + 1;
    int progress_percent =
      static_cast<int>(
        std::floor(
          100.0 * completed_iterations /
            static_cast<double>(iterations)
        )
      );
    if (print_progress &&
        progress_percent >= next_progress_percent) {
      std::time_t now = std::time(nullptr);
      double elapsed_seconds = std::difftime(now, start_time);
      double remaining_seconds =
        elapsed_seconds *
        (iterations - completed_iterations) /
        static_cast<double>(completed_iterations);
      std::string phase = adapt_proposals ? "warmup" : "sampling";

      std::ostringstream message;
      message
        << "### Chain " << chain_id
        << " | seed " << chain_seed
        << " | worker PID " << worker_pid
        << " | regime " << initialization_regime
        << " | phase " << phase
        << " | iteration " << completed_iterations
        << " / " << iterations
        << " (" << progress_percent << "%)"
        << " | elapsed_sec " << elapsed_seconds
        << " | remaining_sec " << remaining_seconds
        << "\n"
        << "### acceptance beta "
        << static_cast<double>(acceptbeta) /
          std::max(1, completed_iterations)
        << " | theta "
        << static_cast<double>(accepttheta) /
          std::max(1, completed_iterations)
        << " | r "
        << static_cast<double>(arma::sum(acceptr)) /
          static_cast<double>(nq) /
          std::max(1, completed_iterations)
        << " | v "
        << static_cast<double>(acceptv) /
          std::max(1, completed_iterations)
        << " | rho "
        << static_cast<double>(acceptrho) /
          std::max(1, completed_iterations)
        << " | eta "
        << static_cast<double>(accepteta) /
          std::max(1, completed_iterations)
        << " | A "
        << static_cast<double>(acceptA) /
          std::max(1, completed_iterations)
        << "\n";

      Rcpp::Rcout << message.str();
      if (write_progress_log) {
        progress_log << message.str();
        progress_log.flush();
      }

      while (next_progress_percent <= progress_percent) {
        next_progress_percent += 5;
      }
    }

    double iter07 = pow(adaptation_offset + iter + 2, -0.7);

    /////////////////
    // Update beta //
    /////////////////

    arma::vec pro_beta = safe_random_walk(beta, xi_beta, "beta");

    // Calculate X_pro_beta and X_beta using matrix multiplication
    arma::vec X_pro_beta = X * pro_beta;

    // Calculate MH_beta using vectorized operations
    double MH_beta = 0.0;

    for (int i = 0; i < p; i++){
      MH_beta += 0.5 * taub * (pow(beta(i),2.0) - pow(pro_beta(i),2.0));
    }

    for (int id = 0; id < nq; id++){
      MH_beta += inverse_temperature * (
        y(id) * (X_pro_beta(id) - X_beta(id)) -
        E(id) * exp(theta(u(id) - 1)) *
          (exp(X_pro_beta(id)) - exp(X_beta(id)))
      );
    }

    MH_beta = finite_log_acceptance(MH_beta);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_beta)) {
      beta = pro_beta;
      X_beta = X_pro_beta;
      acceptbeta++;
    }

    // Adapt only during warm-up; retained draws use a fixed transition kernel.
    if (adapt_proposals) {
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
    }

    //////////////////
    // Update theta //
    //////////////////

    arma::vec pro_theta = safe_random_walk(theta, xi_theta, "theta");

    // Calculate MH_theta
    double MH_theta = 0.0;

    for (int k = 0; k < n_atoms; k++){
      MH_theta += 0.5 * taus * (pow(theta(k),2.0) - pow(pro_theta(k),2.0));
    }

    for (int id = 0; id < nq; id++) {
      MH_theta += inverse_temperature * (
        y(id) * (pro_theta(u(id) - 1) - theta(u(id) - 1)) -
        E(id) * std::exp(X_beta(id)) *
          (std::exp(pro_theta(u(id) - 1)) - std::exp(theta(u(id) - 1)))
      );
    }

    MH_theta = finite_log_acceptance(MH_theta);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_theta)) {
      theta = pro_theta;
      accepttheta++;
    }

    if (adapt_proposals) {
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
    }

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
      pro_Fr(id) = arma::normcdf(
        pro_r(id),
        0.0,
        std::sqrt(marginal_Vr_diag(id))
      );
      pro_u(id) = makeuk(pro_Fr(id), probs);

      pro_b = compute_b(n, q, inv_A, pro_r);

      double pro_r_dens = ldmvnorm(n, q, pro_b, A, Q, Tau);

      MH_r(id) = pro_r_dens - r_dens + inverse_temperature * (
        y(id) * (theta(pro_u(id) - 1) - theta(u(id) - 1)) -
        E(id) * std::exp(X_beta(id)) *
          (std::exp(theta(pro_u(id) - 1)) - std::exp(theta(u(id) - 1)))
      );

      MH_r(id) = finite_log_acceptance(MH_r(id));

      if (arma::randu() < std::exp(MH_r(id))) {
        r(id) = pro_r(id);
        F_r(id) = pro_Fr(id);
        u(id) = pro_u(id);
        b = pro_b;
        r_dens = pro_r_dens;
        acceptr(id)++;
      }

      // Each r coordinate is proposed separately, so its Robbins--Monro
      // target is the conventional one-dimensional random-walk value 0.44,
      // rather than the 0.234 target used for multivariate blocks.  This
      // adaptation is active only during warm-up and its state is persisted
      // at the temperature level across the short tempering segments.
      if (adapt_proposals) {
        s_r(id) += iter07 * (std::exp(MH_r(id)) - 0.44);
        s_r(id) = std::max(-12.0, std::min(12.0, s_r(id)));

        double dr_2 = r(id) - m_r(id);
        m_r(id) += iter07 * dr_2;
        R_r(id) =
          (1.0 - iter07) * R_r(id) + iter07 * dr_2 * dr_2;
        if (!std::isfinite(R_r(id)) || R_r(id) <= 0.0) {
          stop("The adaptive scalar-r variance became invalid.");
        }
        xi_r(id) = std::exp(s_r(id)) * R_r(id);
        if (!std::isfinite(xi_r(id))) {
          stop("The adaptive scalar-r proposal variance became non-finite.");
        }
        xi_r(id) = std::max(1e-12, xi_r(id));
      }

    }

    // The experimental r-block and transported A-row proposals from the
    // enhanced sampler are intentionally excluded from this baseline build.

    //////////////
    // Update v //
    //////////////

    arma::vec pro_vv = safe_random_walk(vv, xi_vv, "stick weights");

    arma::vec pro_v = arma::vec(n_atoms, arma::fill::zeros);
    for (int k = 0; k < n_atoms; k++){
      pro_v(k) = inv_trans_par(pro_vv(k), 0.0, 1.0);
    }

    // Calculate probabilities
    arma::vec pro_probs = makeprobs(pro_v);

    // Calculate u
    arma::ivec pro_u = makeu(F_r, pro_probs);

    double MH_vv = 0.0;
    for (int atom = 0; atom < n_atoms; atom++) {
      MH_vv +=
        pro_vv(atom) - (alpha + 1.0) * log1pexp_stable(pro_vv(atom)) -
        vv(atom) + (alpha + 1.0) * log1pexp_stable(vv(atom));
    }

    for (int id = 0; id < nq; id++) {
      MH_vv += inverse_temperature * (
        y(id) * (theta(pro_u(id) - 1) - theta(u(id) - 1)) -
        E(id) * std::exp(X_beta(id)) *
          (std::exp(theta(pro_u(id) - 1)) - std::exp(theta(u(id) - 1)))
      );
    }

    MH_vv = finite_log_acceptance(MH_vv);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_vv)) {
      vv = pro_vv;
      v = pro_v;
      probs = pro_probs;
      u = pro_u;
      acceptv++;
    }

    if (adapt_proposals) {
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
    }

    /////////////////
    // Update taus //
    /////////////////

    taus = arma::randg<double>(arma::distr_param(0.5 * n_atoms + c,
                                                 1/(dot(theta,theta)/2 + d2)));

    ////////////////
    // Update rho //
    ////////////////

    arma::vec pro_rhorho = safe_random_walk(rhorho, xi_rhorho, "rho");

    // Apply inverse transformation
    arma::vec pro_rho = arma::vec(q, arma::fill::zeros);
    for (int d = 0; d < q; d++){
      pro_rho(d) = inv_trans_par(pro_rhorho(d), 0.0, 1.0);
    }

    std::vector<arma::mat> pro_Tau = Tau_new(pro_rho, n, q, W);
    std::vector<arma::mat> pro_Q = Q_new(pro_rho, n, q, W, pro_Tau);

    std::vector<arma::vec> rho_pro_Q_marginal_variances =
      compute_Q_marginal_variances(pro_Q);
    arma::vec rho_pro_marginal_Vr_diag = combine_marginal_Vr_diag(
      n, q, A, rho_pro_Q_marginal_variances
    );

    double pro_r_dens = ldmvnorm(n, q, b, A, pro_Q, pro_Tau);

    arma::vec rho_pro_F_r = arma::normcdf(
      r,
      mu_r,
      arma::sqrt(rho_pro_marginal_Vr_diag)
    );
    arma::ivec rho_pro_u = makeu(rho_pro_F_r, probs);
    int rho_u_changes = static_cast<int>(arma::accu(rho_pro_u != u));
    proposed_rho_u_changes += rho_u_changes;
    double rho_likelihood_difference =
      allocation_log_likelihood_difference(
        y, E, X_beta, theta, rho_pro_u, u
      );

    double MH_rhorho = pro_r_dens - r_dens +
      inverse_temperature * rho_likelihood_difference;
    for (int disease = 0; disease < q; disease++) {
      MH_rhorho +=
        pro_rhorho(disease) -
          2.0 * log1pexp_stable(pro_rhorho(disease)) -
        rhorho(disease) +
          2.0 * log1pexp_stable(rhorho(disease));
    }

    MH_rhorho = finite_log_acceptance(MH_rhorho);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_rhorho)) {
      rhorho = pro_rhorho;
      rho = pro_rho;
      Tau = pro_Tau;
      Q = pro_Q;
      Q_marginal_variances = rho_pro_Q_marginal_variances;
      marginal_Vr_diag = rho_pro_marginal_Vr_diag;
      r_dens = pro_r_dens;
      F_r = rho_pro_F_r;
      u = rho_pro_u;
      acceptrho++;
      accepted_rho_u_changes += rho_u_changes;
    }

    if (adapt_proposals) {
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
    }

    ////////////////
    // Update eta //
    ////////////////
if (cvrts == "adj" || cvrts == "meanadj") {
    arma::vec pro_etaeta = safe_random_walk(etaeta, xi_etaeta, "eta");

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

    std::vector<arma::mat> pro_W = W_new_from_eta(
      cvrts, q, Minc, Z1, Z2, Z3, pro_eta
    );

    pro_Tau = Tau_new(rho, n, q, pro_W);

    pro_Q = Q_new(rho, n, q, pro_W, pro_Tau);

    std::vector<arma::vec> eta_pro_Q_marginal_variances =
      compute_Q_marginal_variances(pro_Q);
    arma::vec eta_pro_marginal_Vr_diag = combine_marginal_Vr_diag(
      n, q, A, eta_pro_Q_marginal_variances
    );

    pro_r_dens = ldmvnorm(n, q, b, A, pro_Q, pro_Tau);

    arma::vec eta_pro_F_r = arma::normcdf(
      r,
      mu_r,
      arma::sqrt(eta_pro_marginal_Vr_diag)
    );
    arma::ivec eta_pro_u = makeu(eta_pro_F_r, probs);
    int eta_u_changes = static_cast<int>(arma::accu(eta_pro_u != u));
    proposed_eta_u_changes += eta_u_changes;
    double eta_likelihood_difference =
      allocation_log_likelihood_difference(
        y, E, X_beta, theta, eta_pro_u, u
      );

    double MH_etaeta = pro_r_dens - r_dens +
      inverse_temperature * eta_likelihood_difference;
    for (int coordinate = 0; coordinate < 3 * q; coordinate++) {
      MH_etaeta +=
        pro_etaeta(coordinate) -
          2.0 * log1pexp_stable(pro_etaeta(coordinate)) -
        etaeta(coordinate) +
          2.0 * log1pexp_stable(etaeta(coordinate));
    }

    MH_etaeta = finite_log_acceptance(MH_etaeta);

    // Metropolis-Hastings acceptance step
    if (arma::randu() < std::exp(MH_etaeta)) {
      etaeta = pro_etaeta;
      eta = pro_eta;
      W = pro_W;
      Tau = pro_Tau;
      Q = pro_Q;
      Q_marginal_variances = eta_pro_Q_marginal_variances;
      marginal_Vr_diag = eta_pro_marginal_Vr_diag;
      r_dens = pro_r_dens;
      F_r = eta_pro_F_r;
      u = eta_pro_u;
      accepteta++;
      accepted_eta_u_changes += eta_u_changes;
    }

    if (adapt_proposals) {
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

  }

    // The experimental joint rho--eta proposal is intentionally excluded
    // from this baseline build.

    //////////////
    // Update A //
    //////////////

    arma::vec pro_A_vec = safe_random_walk(A_vec, xi_A, "A");

    arma::mat pro_A(q, q, arma::fill::zeros);
    arma::mat inv_pro_A(q, q, arma::fill::zeros);
    arma::mat pro_b(n, q, arma::fill::zeros);
    arma::mat inv_pro_Sigma(q, q, arma::fill::zeros);
    arma::vec A_pro_marginal_Vr_diag(nq, arma::fill::ones);
    arma::vec A_pro_F_r = F_r;
    arma::ivec A_pro_u = u;
    double pro_lpA = -arma::datum::inf;
    double MH_A = -arma::datum::inf;
    int A_u_changes = 0;

    bool valid_A_proposal = pro_A_vec.is_finite();
    int idxproA = 0;
    for (int d = 0; d < q && valid_A_proposal; d++) {
      for (int h = 0; h <= d; h++) {
        if (h == d) {
          // Reject an unrepresentable Cholesky diagonal instead of allowing
          // exp() or the triangular inverse to throw during a long run.
          if (pro_A_vec(idxproA) < -700.0 ||
              pro_A_vec(idxproA) > 700.0) {
            valid_A_proposal = false;
            break;
          }
          pro_A(d, h) = std::exp(pro_A_vec(idxproA));
        } else {
          pro_A(d, h) = pro_A_vec(idxproA);
        }
        idxproA++;
      }
    }

    if (valid_A_proposal) {
      valid_A_proposal = arma::inv(inv_pro_A, arma::trimatl(pro_A));
    }
    if (valid_A_proposal && inv_pro_A.is_finite()) {
      pro_b = compute_b(n, q, inv_pro_A, r);
      pro_r_dens = ldmvnorm(n, q, pro_b, pro_A, Q, Tau);
      inv_pro_Sigma = inv_pro_A.t() * inv_pro_A;
      pro_lpA = log_density_A(pro_A, inv_pro_Sigma, R, nu);
      A_pro_marginal_Vr_diag = combine_marginal_Vr_diag(
        n, q, pro_A, Q_marginal_variances
      );
      A_pro_F_r = arma::normcdf(
        r,
        mu_r,
        arma::sqrt(A_pro_marginal_Vr_diag)
      );
      A_pro_u = makeu(A_pro_F_r, probs);
      A_u_changes = static_cast<int>(arma::accu(A_pro_u != u));
      proposed_A_u_changes += A_u_changes;
      double A_likelihood_difference =
        allocation_log_likelihood_difference(
          y, E, X_beta, theta, A_pro_u, u
        );
      MH_A = finite_log_acceptance(
        pro_r_dens - r_dens + pro_lpA - lpA +
          inverse_temperature * A_likelihood_difference
      );
    }

    if(arma::randu() < std::exp(MH_A)){
      A_vec = pro_A_vec;
      A = pro_A;
      inv_A = inv_pro_A;
      b = pro_b;
      marginal_Vr_diag = A_pro_marginal_Vr_diag;
      inv_Sigma = inv_pro_Sigma;
      lpA = pro_lpA;
      r_dens = pro_r_dens;
      F_r = A_pro_F_r;
      u = A_pro_u;
      acceptA++;
      accepted_A_u_changes += A_u_changes;
    }

    if (adapt_proposals) {
    // Update s_A
    s_A += iter07 * (std::exp(MH_A) - 0.234);

    // Broad numerical safeguards only.  The former lower bound at -4 was
    // active in the SEER pilot and prevented the adaptive random walk from
    // shrinking enough to approach its target acceptance rate.  On the
    // covariance scale exp(s_A), +/-12 remains finite while being far outside
    // the range expected under ordinary proposal adaptation.
    s_A = std::max(-12.0, std::min(12.0, s_A));

    // Calculate dA_2 and update m_A in a single loop
    arma::vec dA_2 = A_vec - m_A;
    m_A += iter07 * dA_2;

    // Update R_A and xi_A using vectorized operations
    R_A = (1 - iter07) * R_A + iter07 * (dA_2 * dA_2.t());
    xi_A = exp(s_A) * R_A;

    // Adjust diagonal of xi_A using vectorized operations
    xi_A.diag() += 1e-3 * xi_A.diag();
    }

    ////////////////////
    // record samples //
    ////////////////////

    if (iter >= burn && (iter - burn) % thin == 0) {
      assert_state_consistency(
        r,
        marginal_Vr_diag,
        F_r,
        u,
        probs,
        state_consistency_tolerance,
        chain_id,
        iter,
        maximum_F_gap_seen,
        maximum_u_mismatches_seen
      );
      assert_W_eta_consistency(
        W, cvrts, q, Minc, Z1, Z2, Z3, eta,
        chain_id, iter, maximum_W_eta_mismatches_seen
      );
      invariant_checks++;
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
      for (int disease = 0; disease < q; disease++) {
        sumW[disease] += W[disease];
        keepWcardinality(g, disease) = static_cast<int>(
          arma::accu(arma::trimatu(W[disease], 1))
        );
      }
      keepA.row(g) = A_vec.t();

      g++;
    }

  }

  assert_state_consistency(
    r,
    marginal_Vr_diag,
    F_r,
    u,
    probs,
    state_consistency_tolerance,
    chain_id,
    iterations,
    maximum_F_gap_seen,
    maximum_u_mismatches_seen
  );
  assert_W_eta_consistency(
    W, cvrts, q, Minc, Z1, Z2, Z3, eta,
    chain_id, iterations, maximum_W_eta_mismatches_seen
  );
  invariant_checks++;

  if (g != runs) {
    stop(
      "Chain %d retained %d draws, but runs requested %d.",
      chain_id,
      g,
      runs
    );
  }

  if (write_progress_log) {
    progress_log.close();
  }

  List keepWmean(q);
  for (int disease = 0; disease < q; disease++) {
    keepWmean[disease] = sumW[disease] / static_cast<double>(runs);
  }

  double final_log_likelihood = 0.0;
  for (int id = 0; id < nq; id++) {
    double linear_predictor = X_beta(id) + theta(u(id) - 1);
    final_log_likelihood +=
      y(id) * linear_predictor -
      E(id) * std::exp(linear_predictor);
  }

  int eta_proposed =
    (cvrts == "adj" || cvrts == "meanadj") ? iterations : 0;

  List acceptance;
  acceptance["beta_proposed"] = iterations;
  acceptance["beta_accepted"] = acceptbeta;
  acceptance["beta_rate"] =
    static_cast<double>(acceptbeta) / iterations;
  acceptance["theta_proposed"] = iterations;
  acceptance["theta_accepted"] = accepttheta;
  acceptance["theta_rate"] =
    static_cast<double>(accepttheta) / iterations;
  acceptance["r_proposed"] = iterations * nq;
  acceptance["r_accepted"] = arma::sum(acceptr);
  acceptance["r_rate"] =
    static_cast<double>(arma::sum(acceptr)) /
    static_cast<double>(iterations * nq);
  acceptance["r_accepted_by_coordinate"] = acceptr;
  acceptance["r_rate_by_coordinate"] =
    arma::conv_to<arma::vec>::from(acceptr) /
    static_cast<double>(iterations);
  acceptance["r_final_proposal_sd"] = arma::sqrt(xi_r);
  acceptance["V_proposed"] = iterations;
  acceptance["V_accepted"] = acceptv;
  acceptance["V_rate"] =
    static_cast<double>(acceptv) / iterations;
  acceptance["rho_proposed"] = iterations;
  acceptance["rho_accepted"] = acceptrho;
  acceptance["rho_rate"] =
    static_cast<double>(acceptrho) / iterations;
  acceptance["rho_mean_proposed_allocation_changes"] =
    proposed_rho_u_changes / static_cast<double>(iterations);
  acceptance["rho_mean_accepted_allocation_changes"] =
    acceptrho > 0 ?
      accepted_rho_u_changes / static_cast<double>(acceptrho) :
      NA_REAL;
  acceptance["eta_proposed"] = eta_proposed;
  acceptance["eta_accepted"] = accepteta;
  acceptance["eta_rate"] =
    eta_proposed > 0 ?
    static_cast<double>(accepteta) / eta_proposed :
    NA_REAL;
  acceptance["eta_mean_proposed_allocation_changes"] =
    eta_proposed > 0 ?
      proposed_eta_u_changes / static_cast<double>(eta_proposed) :
      NA_REAL;
  acceptance["eta_mean_accepted_allocation_changes"] =
    accepteta > 0 ?
      accepted_eta_u_changes / static_cast<double>(accepteta) :
      NA_REAL;
  acceptance["A_proposed"] = iterations;
  acceptance["A_accepted"] = acceptA;
  acceptance["A_rate"] =
    static_cast<double>(acceptA) / iterations;
  acceptance["A_mean_proposed_allocation_changes"] =
    proposed_A_u_changes / static_cast<double>(iterations);
  acceptance["A_mean_accepted_allocation_changes"] =
    acceptA > 0 ?
      accepted_A_u_changes / static_cast<double>(acceptA) :
      NA_REAL;

  List final_state = List::create(
    _["beta"] = beta,
    _["theta"] = theta,
    _["tau"] = taus,
    _["V"] = v,
    _["rho"] = rho,
    _["eta"] = eta,
    _["A"] = A_vec,
    _["r"] = r,
    _["F_r"] = F_r,
    _["u"] = u
  );

  List proposal_state;
  proposal_state["s_theta"] = s_theta;
  proposal_state["s_beta"] = s_beta;
  proposal_state["s_r"] = s_r;
  proposal_state["s_vv"] = s_vv;
  proposal_state["s_rhorho"] = s_rhorho;
  proposal_state["s_etaeta"] = s_etaeta;
  proposal_state["s_A"] = s_A;
  proposal_state["m_theta"] = m_theta;
  proposal_state["m_beta"] = m_beta;
  proposal_state["m_r"] = m_r;
  proposal_state["m_vv"] = m_vv;
  proposal_state["m_rhorho"] = m_rhorho;
  proposal_state["m_etaeta"] = m_etaeta;
  proposal_state["m_A"] = m_A;
  proposal_state["R_theta"] = R_theta;
  proposal_state["R_beta"] = R_beta;
  proposal_state["R_r"] = R_r;
  proposal_state["R_vv"] = R_vv;
  proposal_state["R_rhorho"] = R_rhorho;
  proposal_state["R_etaeta"] = R_etaeta;
  proposal_state["R_A"] = R_A;
  proposal_state["adaptation_iterations"] =
    adaptation_offset + (adapt_proposals ? iterations : 0);

  List state_consistency = List::create(
    _["checks"] = invariant_checks,
    _["tolerance"] = state_consistency_tolerance,
    _["maximum_F_gap_seen"] = maximum_F_gap_seen,
    _["maximum_u_mismatches_seen"] = maximum_u_mismatches_seen,
    _["maximum_W_eta_mismatches_seen"] = maximum_W_eta_mismatches_seen,
    _["all_checks_passed"] =
      maximum_F_gap_seen <= state_consistency_tolerance &&
      maximum_u_mismatches_seen == 0 &&
      maximum_W_eta_mismatches_seen == 0
  );

  // Build these longer lists by named assignment because older Rcpp releases
  // cap List::create() at 20 arguments.
  List settings;
  settings["chain_id"] = chain_id;
  settings["seed"] = chain_seed;
  settings["worker_pid"] = worker_pid;
  settings["initialization_regime"] = initialization_regime;
  settings["runs"] = runs;
  settings["burn"] = burn;
  settings["thin"] = thin;
  settings["iterations"] = iterations;
  settings["sampler"] = "support_consistent_tempered_segment";
  settings["adaptation"] = adapt_proposals ? "enabled" : "frozen";
  settings["adaptation_offset"] = adaptation_offset;
  settings["inverse_temperature"] = inverse_temperature;
  settings["state_consistency_tolerance"] = state_consistency_tolerance;
  settings["transition_blocks"] =
    "baseline blocks with Poisson-likelihood tempering only";

  List result;
  result["chain_id"] = chain_id;
  result["seed"] = chain_seed;
  result["settings"] = settings;
  result["initial_values"] = initial_values;
  result["acceptance"] = acceptance;
  result["state_consistency"] = state_consistency;
  result["final_state"] = final_state;
  result["proposal_state"] = proposal_state;
  result["log_likelihood"] = final_log_likelihood;
  result["beta"] = keepbeta;
  result["phi"] = keepphi;
  result["theta"] = keeptheta;
  result["u"] = keepu;
  result["rho"] = keeprho;
  result["V"] = keepv;
  result["r"] = keepr;
  result["F_r"] = keepFr;
  result["eta"] = keepeta;
  result["tau"] = keeptaus;
  result["W_mean"] = keepWmean;
  result["W_cardinality"] = keepWcardinality;
  result["A"] = keepA;

  return result;

}
