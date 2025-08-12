#include <Rcpp.h>
#include <Rcpp/Rmath.h>



using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getSSP(NumericVector parms, int trialType, int nTrials,
double dt, double var) {

  //trialType: 1 = congruent, 2 = incongruent

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //set fixed parameters
  double sdRand = sqrt(dt * var); //get standard deviation from variance

  //get the free parameters from the numeric vector parms
  double A = parms[0];
  double B = -parms[0];
  double tEr = parms[1];
  double p = parms[2];
  double rd = parms[3];
  double sda = parms[4];

  //declare random number variable for noise
  double noise = 0.0;

  //set empty matrix to store trial data
  int nRow = nTrials; //first set number of rows
  NumericMatrix trialData(nRow, 2); //set the matrix. 2 columns (RT & accuracy)

  //keep track of current evidence (response selection)
  double currEvidenceResp = 0.0;
  double t = 0.0; //to log number of diffusion steps
  double drift = 0.0; //the drift rate for the current moment


  //declare trial-specific parameters
  double p_target = 0.0;
  double p_flanker = 0.0;

  //get perceptual input for targets and flanker
  p_target = p;
  p_flanker = p;

  // flip the sign if current trial is incongruent
  if(trialType == 2){
    p_flanker = -p_flanker;
  }


  // pre-compute sd(t) vector, and from it the attentional areas
  // over target and flankers
  int m = 20000; // m needs to be larger than maximum expexted RT for the model
  NumericVector sd_t_vec(m);
  NumericVector a_target_vec(m);
  NumericVector a_flanker_vec(m);
  NumericVector mu_target_vec(m);
  NumericVector mu_flanker_vec(m);

  for(int t = 0; t < m; ++t) {
    sd_t_vec[t] = sda - (rd * t);
    if(sd_t_vec[t] <= 0.001) {
      sd_t_vec[t] = 0.001;
    }

    // attentional spotlights
    a_target_vec[t] = ::Rf_pnorm5(0.5, 0.0, sd_t_vec[t], 1, 0) - ::Rf_pnorm5(-0.5, 0.0, sd_t_vec[t], 1, 0);
    a_flanker_vec[t] = ::Rf_pnorm5(10.0, 0.0, sd_t_vec[t], 1, 0) - ::Rf_pnorm5(0.5, 0.0, sd_t_vec[t], 1, 0);

    // drift rate contributions from the target and the flankers
    mu_target_vec[t] = (p_target * a_target_vec[t]) * dt;
    mu_flanker_vec[t] = ((2 * p_flanker) * a_flanker_vec[t]) * dt;

  }


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  //start trial loop here
  for (int i=0; i<=nTrials - 1; i++){

    //keep track of the current evidence (response selection)
    currEvidenceResp = 0.0;

    //reset counter to log how many steps taken in diffusion process
    t = 0.0;

      //diffusion simulation starts here
      while((currEvidenceResp <= A) && (currEvidenceResp >= B)){

      //current drift rate
      drift = mu_flanker_vec[t] + mu_target_vec[t];

      //get random noise
      noise = ::Rf_rnorm(drift, sdRand);

      //update the response selection random walk
      currEvidenceResp = currEvidenceResp + noise;

      //update diffusion step number
      t++;

      } // while loop ends here

   trialData(i, 0) = (t * dt) + tEr;


    if(currEvidenceResp >= A){
      trialData(i, 1) = 1;
    }

      if(currEvidenceResp <= B){
        trialData(i, 1) = 0;
      }

  } // trial loop ends here
     return trialData;
}
