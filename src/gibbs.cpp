// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Sampling and exact-value counting retain the source order.

// sample an index in {0,1,...,K-1} with given (already normalized) probs
int sample_index(const NumericVector& p) {
  double u = R::runif(0.0, 1.0);
  double csum = 0.0;
  for (int i = 0; i < p.size(); ++i) {
    csum += p[i];
    if (u <= csum) return i;
  }
  return p.size() - 1; // numerical safety
}



// count unique values and their frequencies (exact equality, as in R code path)
void count_values(const NumericVector& v, NumericVector& uniques, IntegerVector& freqs) {
  // Count exactly equal locations in ascending order, as plyr::count does in R.
  int n = v.size();
  std::vector<double> tmp(n);
  for (int i = 0; i < n; ++i) tmp[i] = v[i];
  std::sort(tmp.begin(), tmp.end());
  std::vector<double> u;
  std::vector<int> f;
  for (int i = 0; i < n; ) {
    double val = tmp[i];
    int j = i + 1;
    while (j < n && tmp[j] == val) j++;
    u.push_back(val);
    f.push_back(j - i);
    i = j;
  }
  uniques = NumericVector(u.begin(), u.end());
  freqs  = IntegerVector(f.begin(), f.end());
}

// remove element i from vector (copy, like theta[-i] in R)
NumericVector drop_index(const NumericVector& v, int idx) {
  int n = v.size();
  NumericVector out(n - 1);
  for (int i = 0, k = 0; i < n; ++i) {
    if (i == idx) continue;
    out[k++] = v[i];
  }
  return out;
}

// [[Rcpp::export]]
List py_sampler_cpp(NumericVector standardized_observations,
                    int niter = 1000,
                    double concentration = 1.0,
                    double base_mean = 0.0,
                    double base_variance = 1.0,
                    double variance_shape = 1.0,
                    double variance_rate = 1.0,
                    double discount = 0.3
) {
  // Input is already standardized in R; no transformation is performed here.
  NumericVector observations = clone(standardized_observations);

  int n = observations.size();

  // storage
  NumericMatrix theta_store(niter, n);
  NumericVector sigma2_store(niter);
  IntegerVector k_store(niter);
  List pi_list(niter); // to mirror pi[[iter]] <- pi_vec (last i's probabilities each iter)

  // Initial distinct normal means and shared variance.
  NumericVector theta(n);
  for (int i = 0; i < n; ++i) theta[i] = R::rnorm(0.0, 1.0);
  double sigma2 = 1.0;


  // Allocation updates, variance update, then occupied-mean acceleration.
  for (int iter = 0; iter < niter; ++iter) {

    // update thetas
    NumericVector last_pi_vec; // to save after finishing i-loop
    for (int i = 0; i < n; ++i) {
      double Xi = observations[i];
      // theta_-i
      NumericVector theta_minus_i = drop_index(theta, i);

      // COUNT <- count(theta_-i)
      NumericVector theta_minus_i_star;
      IntegerVector n_minus_i_star;
      count_values(theta_minus_i, theta_minus_i_star, n_minus_i_star);
      int k_minus_i_star = theta_minus_i_star.size();

      // prepare probabilities
      NumericVector pi_vec_temp(k_minus_i_star + 1);

      // Gaussian conditional parameters and reference point for the marginal weight.
      double t0 = 0.0;
      double mui = (base_variance / (base_variance + sigma2)) * Xi + (sigma2 / (base_variance + sigma2)) * base_mean;
      double sigmai2 = 1.0 / (1.0 / base_variance + 1.0 / sigma2);

      // pi_0 (new cluster)
      // exp( log(concentration + k_-i^* * discount) - log(concentration + n - 1)
      //      + dnorm(Xi | t0, sqrt(sigma2), log=TRUE)
      //      + dnorm(t0 | base_mean, sqrt(base_variance), log=TRUE)
      //      - dnorm(t0 | mui, sqrt(sigmai2), log=TRUE) )
      double log_pi0 =
        std::log(concentration + k_minus_i_star * discount) - std::log(concentration + n - 1.0) +
        R::dnorm(Xi, t0, std::sqrt(sigma2), 1) +
        R::dnorm(t0, base_mean, std::sqrt(base_variance), 1) -
        R::dnorm(t0, mui, std::sqrt(sigmai2), 1);
      pi_vec_temp[0] = std::exp(log_pi0);

      // existing clusters j = 1..k_-i^*
      for (int j = 0; j < k_minus_i_star; ++j) {
        double theta_j = theta_minus_i_star[j];
        double log_p =
          std::log( (double)n_minus_i_star[j] - discount ) - std::log(concentration + n - 1.0) +
          R::dnorm(Xi, theta_j, std::sqrt(sigma2), 1);
        pi_vec_temp[j + 1] = std::exp(log_p);
      }

      // normalize
      double sumw = 0.0;
      for (int k = 0; k < pi_vec_temp.size(); ++k) sumw += pi_vec_temp[k];
      NumericVector pi_vec(pi_vec_temp.size());
      for (int k = 0; k < pi_vec.size(); ++k) pi_vec[k] = pi_vec_temp[k] / sumw;

      // sample index in {0,...,k_-i^*}
      int index = sample_index(pi_vec);

      if (index == 0) {
        // new: theta[i] ~ N(mui, sigmai2)
        double draw = R::rnorm(mui, std::sqrt(sigmai2));
        theta[i] = draw;
      } else {
        // old: copy existing cluster center
        theta[i] = theta_minus_i_star[index - 1];
      }

      last_pi_vec = pi_vec; // keep the last pi_vec from the loop (as in R: pi[[iter]] <- pi_vec)
    }

    // save last pi_vec for this iteration
    pi_list[iter] = last_pi_vec;

    // Store theta before occupied-mean acceleration.
    for (int i = 0; i < n; ++i) theta_store(iter, i) = theta[i];

    // update sigma2
    double a1 = variance_shape + n / 2.0;
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      double r = observations[i] - theta[i];
      ss += r * r;
    }
    double b1 = variance_rate + 0.5 * ss;

    // sigma2 <- 1 / rgamma(shape=a1, rate=b1)
    // R::rgamma takes scale; the R reference uses rate=b1, so scale=1/b1.
    double g = R::rgamma(a1, 1.0 / b1);
    sigma2 = 1.0 / g;

    sigma2_store[iter] = sigma2;

    // acceleration step
    // COUNT <- count(theta)
    NumericVector theta_star;
    IntegerVector n_star;
    count_values(theta, theta_star, n_star);
    int k_star = theta_star.size();

    for (int j = 0; j < k_star; ++j) {
      double thetaj_star = theta_star[j];
      // indices where theta == thetaj_star
      std::vector<int> idx;
      idx.reserve(n_star[j]);
      for (int i = 0; i < n; ++i) {
        if (theta[i] == thetaj_star) idx.push_back(i);
      }
      int nj = (int)idx.size();

      // Xbarj
      double sumx = 0.0;
      for (int t = 0; t < nj; ++t) sumx += observations[idx[t]];
      double Xbarj = sumx / nj;

      double sigmaj2_star = 1.0 / ( (nj / sigma2) + (1.0 / base_variance) );
      double muj_star = sigmaj2_star * ( (base_mean / base_variance) + (nj * Xbarj / sigma2) );

      double draw = R::rnorm(muj_star, std::sqrt(sigmaj2_star));
      for (int t = 0; t < nj; ++t) theta[idx[t]] = draw;
    }

    k_store[iter] = k_star;
  }

  return List::create(
    _["theta_store"] = theta_store,
    _["sigma2_store"] = sigma2_store,
    _["k_store"] = k_store,
    _["pi"] = pi_list,
    _["data_scaled"] = observations
  );
}

// Backward-compatible source arguments. nburn remains unused: all rows are saved.
// [[Rcpp::export]]
List dp_sampler_cpp(NumericVector galaxies, int niter = 1000, int nburn = 100,
                    double c = 1.0, double mu0 = 0.0, double sigma02 = 1.0,
                    double a = 1.0, double b = 1.0, double d = 0.3) {
  (void)nburn;
  return py_sampler_cpp(galaxies, niter, c, mu0, sigma02, a, b, d);
}
