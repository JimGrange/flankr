#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getDSTP_new(NumericVector parms,
                          int trialType,
                          int nTrials,
                          double dt,
                          double var) {

  // initialise R's RNG (respects set.seed in R)
  Rcpp::RNGScope scope;

  double sdRand = sqrt(dt * var);
  double A = parms[0], B = -parms[0];
  double C = parms[1], D = -parms[1];
  double muTa = parms[2] * dt;
  double muFl = parms[3] * dt;
  double muSS = parms[4] * dt;
  double muRS2 = parms[5] * dt;
  double tEr = parms[6];

  // pre-compute drift for first stage based on trial type
  double muResp = (trialType == 1) ? (muTa + muFl) : (muTa - muFl);

  // pre-generate all noise vectors
  int noiseLen = 10000; // keep original length for behaviour parity
  NumericVector muNoise = rnorm(noiseLen, muResp, sdRand);
  NumericVector rsNoise = rnorm(noiseLen, muRS2, sdRand);
  NumericVector ssNoise = rnorm(noiseLen, muSS, sdRand);

  // use raw pointers for fast access
  double* muNoisePtr = REAL(muNoise);
  double* rsNoisePtr = REAL(rsNoise);
  double* ssNoisePtr = REAL(ssNoise);

  NumericMatrix trialData(nTrials, 2);

  // initialise random number generator
  // now using R's RNG (no srand/rand; CRAN-friendly and reproducible with set.seed)

  // loop over trials
  for (int i = 0; i < nTrials; ++i) {

    // reset everything
    double currEvidenceResp = 0.0;
    double currEvidenceStim = 0.0;
    int whichStim = 0;
    bool stimSelected = false;
    int j = 0;

    // diffusion process starts here
    while (currEvidenceResp <= A && currEvidenceResp >= B) {

      // select random noise
      // uniform index in [0, noiseLen)
      int idx = static_cast<int>(R::unif_rand() * noiseLen);

      // update the response drift rates based on
      // stimulus selection status
      if (!stimSelected) {
        currEvidenceResp += muNoisePtr[idx];
      } else {
        double delta = rsNoisePtr[idx];
        // apply multiplier for when flanker is selected on incongruent trials
        if (trialType == 2 && whichStim == 2) {
          delta = -delta;
        }
        currEvidenceResp += delta;
      }

      // update the stimulus selection drift rate
      currEvidenceStim += ssNoisePtr[idx];

      // check whether stimulus selection has occurred
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

    // collate trial level data
    trialData(i, 0) = (j * dt) + tEr;
    trialData(i, 1) = (currEvidenceResp >= A) ? 1.0 : 0.0;
  }

  return trialData;
}
