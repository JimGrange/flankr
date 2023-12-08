#include <Rcpp.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]

NumericMatrix getDSTP(NumericVector parms, int trialType, int nTrials,
   double dt, double var) {

  //trialType: 1 = congruent, 2 = incongruent

  //set fixed parameters
//  double dt = 0.001; //time per step of diffusion process
//  double var = 0.01;
  double sdRand = sqrt(dt * var);


  //get the free parameters from the numeric vector parms
  double A = parms[0];  // correct response selection boundary
  double B = -parms[0]; // error response selection boundary
  double C = parms[1];  // target stimulus selection boundary
  double D = -parms[1]; // flanker stimulus selection boundary
  double muTa = parms[2] * dt; // drift for target
  double muFl = parms[3] * dt; // drift for flanker
  double muSS = parms[4] * dt; // drift for stimulus selection
  double muRS2 = parms[5] * dt; // drift for second stage of response selection
  double tEr = parms[6];       // non-decision time


  //declare random number variables for noise
  NumericVector muNoise;
  if(trialType == 1){
    muNoise = rnorm(10000, (muTa + muFl), sdRand);
  } else {
    muNoise = rnorm(10000, (muTa + -muFl), sdRand);
  }

  NumericVector rsNoise = rnorm(10000, muRS2, sdRand);
  NumericVector ssNoise = rnorm(10000, muSS, sdRand);

  int randomIndex;
  srand (1);

  //set empty matrix to store trial data
  int nRow = nTrials; //first set number of rows
  NumericMatrix trialData(nRow, 2); //set the matrix. 2 columns (RT & accuracy)

  //initialise trial flags etc. for whether certain processes have finished
  //during diffusion
  bool stimSelected = false; //has stimulus been selected?

  int whichStim; //to flag whether target (1) or flanker (2) has been selected
  double currEvidenceResp = 0.0; //keep track of evidence (response selection)
  double currEvidenceStim = 0.0; //keep track of evidence (stimulus selection)
  double j = 0.0; //to log number of diffusion steps


  //##
  //start trial loop here

  for (int i=0; i<=nTrials - 1; i++){

    //reset stim selection to false
    stimSelected = false;

    //reset some trial parameters
    //which stim was selected by the stimulus drift (target [1] / flanker [2])?
    whichStim = 0;

    //keep track of the current evidence (response selection)
    currEvidenceResp = 0.0;

    //keep track of the current evidence (Stimulus selection)
    currEvidenceStim = 0.0;

    //reset counter to log how many steps taken in diffusion process
    j = 0.0;


    //##
    //diffusion simulation starts here
    while((currEvidenceResp <= A) && (currEvidenceResp >= B)){

      j++; //update diffusion step number

      randomIndex = rand() % muNoise.size();

      //set the correct drift rate
      if (stimSelected == false){
        currEvidenceResp = currEvidenceResp + muNoise[randomIndex];
        }


      if((stimSelected == true) && (trialType == 2)){
        if(whichStim == 1){
          currEvidenceResp = currEvidenceResp + rsNoise[randomIndex];
        }
        if(whichStim == 2){
          currEvidenceResp = currEvidenceResp + -rsNoise[randomIndex];
        }
      }

      if((stimSelected == true) && (trialType == 1)){
        currEvidenceResp = currEvidenceResp + rsNoise[randomIndex];
      }

        //update the stimulus selection drift rate
        currEvidenceStim = currEvidenceStim + ssNoise[randomIndex];

        //has stimulus selection finished??
        if(currEvidenceStim >= C){
          whichStim = 1;
          stimSelected = true;
        }

        if(currEvidenceStim <= D){
          whichStim = 2;
          stimSelected = true;
        }


    } // while loop ends here

    trialData(i, 0) = (j * dt) + tEr;


    if(currEvidenceResp >= A){
      trialData(i, 1) = 1;
    }
    else
    {
      trialData(i, 1) = 0;
    }

  } // trial loop ends here

  return trialData;

}



