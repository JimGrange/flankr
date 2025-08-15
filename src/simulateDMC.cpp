#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getDMC(NumericVector parms,
                     int trialType,   // 1 = congruent, 2 = incongruent
                     int nTrials,
                     double dt) {

  Rcpp::RNGScope scope;              // respects set.seed in R

  // unpack parameters (order: b, mu_c, maxAmp, tau, alpha, ter, ster, sigma)
  const double bU    = parms[0];
  const double mu_c  = parms[1];
  const double maxAmp= parms[2];
  const double tau   = parms[3];
  const double alpha = parms[4];
  const double ter   = parms[5];
  const double ster  = parms[6];
  const double sigma = parms[7];

  const double bL   = -bU;
  const double rhs  = std::sqrt(dt) * sigma;

  // maximum number of steps; 2 s with dt = 0.1 ms = 20 000 steps
  const int maxIter = 20000;

  // precompute the automatic drift μ_a(t) (Eq. 5) and scale by dt
  NumericVector autoDrift(maxIter);
  double* autoPtr = REAL(autoDrift);
  const double a1   = alpha - 1.0;
  const double sign = (trialType == 1) ? 1.0 : -1.0;
  for (int j = 0; j < maxIter; ++j) {
    double t = (j + 1) * dt;
    double expTerm = std::exp(-t / tau);
    double powTerm = std::pow((t * std::exp(1.0)) / (a1 * tau), a1);
    double bracket = (a1 / t) - (1.0 / tau);
    // multiply by sign and dt so we add directly per step
    autoPtr[j] = sign * maxAmp * expTerm * powTerm * bracket * dt;
  }

  // controlled drift per step
  const double mu_c_dt = mu_c * dt;

  // output matrix (RT, accuracy)
  NumericMatrix trialData(nTrials, 2);

  for (int i = 0; i < nTrials; ++i) {

    // non-decision time: uniform on [ter - ster/2, ter + ster/2]
    const double ter_i = (ster > 1e-8)
    ? ter + (R::unif_rand() * ster) - (ster / 2.0)
      : ter;

    double x = 0.0;
    int tstep = 0;
    int resp  = -1;

    // simulate until bound crossing or maxIter
    while (x < bU && x > bL && tstep < maxIter) {
      // noise ~ N(0, 1) scaled by rhs
      double noise = R::norm_rand() * rhs;
      // add controlled drift, automatic drift, and noise
      x += mu_c_dt + autoPtr[tstep] + noise;
      ++tstep;
      if (x >= bU) { resp = 1; break; }
      if (x <= bL) { resp = 0; break; }
    }

    trialData(i, 0) = (tstep * dt) + ter_i; // RT
    trialData(i, 1) = resp;                 // response: 1 = upper, 0 = lower, -1 = none
  }

  return trialData;
}
