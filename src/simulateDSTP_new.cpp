#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getDSTP_new(NumericVector parms, int trialType, int nTrials,
                      double dt, double var) {

  // Fixed parameters
  double sdRand = sqrt(dt * var);

  // Free parameters
  double A = parms[0], B = -parms[0];
  double C = parms[1], D = -parms[1];
  double muTa = parms[2] * dt;
  double muFl = parms[3] * dt;
  double muSS = parms[4] * dt;
  double muRS2 = parms[5] * dt;
  double tEr = parms[6];

  // Pre-generate noise
  int noiseLen = 10000;
  NumericVector muNoise(noiseLen);
  if (trialType == 1) {
    muNoise = rnorm(noiseLen, muTa + muFl, sdRand);
  } else {
    muNoise = rnorm(noiseLen, muTa - muFl, sdRand);
  }

  NumericVector rsNoise = rnorm(noiseLen, muRS2, sdRand);
  NumericVector ssNoise = rnorm(noiseLen, muSS, sdRand);

  NumericMatrix trialData(nTrials, 2);

  srand(1);  // set RNG seed once

  for (int i = 0; i < nTrials; ++i) {

    bool stimSelected = false;
    int whichStim = 0;
    double currEvidenceResp = 0.0;
    double currEvidenceStim = 0.0;
    int j = 0;

    while (currEvidenceResp <= A && currEvidenceResp >= B) {
      int idx = rand() % noiseLen;

      if (!stimSelected) {
        currEvidenceResp += muNoise[idx];
      } else if (trialType == 2 && whichStim == 2) {
        currEvidenceResp -= rsNoise[idx];
      } else {
        currEvidenceResp += rsNoise[idx];
      }

      currEvidenceStim += ssNoise[idx];

      if (!stimSelected) {
        if (currEvidenceStim >= C) {
          whichStim = 1;
          stimSelected = true;
        } else if (currEvidenceStim <= D) {
          whichStim = 2;
          stimSelected = true;
        }
      }

      ++j;
    }

    trialData(i, 0) = (j * dt) + tEr;
    trialData(i, 1) = (currEvidenceResp >= A) ? 1.0 : 0.0;
  }

  return trialData;
}
