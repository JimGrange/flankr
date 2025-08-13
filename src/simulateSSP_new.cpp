#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getSSP_new(NumericVector parms,
                         int trialType,
                         int nTrials,
                         double dt,
                         double var) {

  // fixed parameters
  const double sdRand = std::sqrt(dt * var);

  // unpack parameters
  const double A = parms[0], B = -parms[0];
  const double tEr = parms[1];
  const double p = parms[2];
  const double rd = parms[3];
  const double sda = parms[4];

  const double p_target = p;
  const double p_flanker = (trialType == 2) ? -p : p;

  const int m = 20000;
  NumericVector mu_target_vec(m), mu_flanker_vec(m);
  double* mu_target = REAL(mu_target_vec);
  double* mu_flanker = REAL(mu_flanker_vec);

  // precompute drift vectors
  for (int t = 0; t < m; ++t) {
    double sd_t = sda - (rd * t);
    if (sd_t <= 0.001) sd_t = 0.001;

    const double inv_sd = 1.0 / sd_t;

    const double a_target =
      R::pnorm5(0.5 * inv_sd, 0.0, 1.0, 1, 0) -
      R::pnorm5(-0.5 * inv_sd, 0.0, 1.0, 1, 0);

    const double a_flanker =
      R::pnorm5(10.0 * inv_sd, 0.0, 1.0, 1, 0) -
      R::pnorm5(0.5 * inv_sd, 0.0, 1.0, 1, 0);

    mu_target[t] = p_target * a_target * dt;
    mu_flanker[t] = 2.0 * p_flanker * a_flanker * dt;
  }

  NumericMatrix trialData(nTrials, 2);

  // simulate each trial
  for (int i = 0; i < nTrials; ++i) {
    double currEvidence = 0.0;
    int t = 0;

    while (currEvidence <= A && currEvidence >= B && t < m) {
      const double drift = mu_target[t] + mu_flanker[t];
      currEvidence += R::rnorm(drift, sdRand);
      ++t;
    }

    trialData(i, 0) = (t * dt) + tEr;
    trialData(i, 1) = (currEvidence >= A) ? 1.0 : 0.0;
  }

  return trialData;
}
