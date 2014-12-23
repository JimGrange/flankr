#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getDSTP(NumericVector parms, int trialType, int nTrials,
   double dt, double var) {

  //trialType: 1 = congruent, 2 = incongruent

  //set fixed parameters
//  double dt = 0.001; //time per step of diffusion process
//  double var = 0.01;
  double sdRand = sqrt(dt * var);


  //get the free parameters from the numeric vector parms
  double A; //correct response selection boundary
    A = parms[0];
   double B; //error response selection boundary
     B = -parms[0];
  double C; //target stimulus selection boundary
    C = parms[1];
   double D; //flanker stimulus selection boundary
     D = -parms[1];
  double muTa; //drift for target
    muTa = parms[2] * dt;
  double muFl; //drift for flanker
    muFl = parms[3] * dt;
  double muSS; //drift for stimulus selection
     muSS = parms[4] * dt;
  double muRS2; //drift for second stage of response selection
    muRS2 = parms[5] * dt;
  double tEr; //non-decision time
    tEr = parms[6];

    //declare random number variable for noise
    NumericVector muNoise_cong;
    NumericVector muNoise_incong;
    NumericVector rsNoise;
    NumericVector ssNoise;
       muNoise_cong = rnorm(100000, (muTa + muFl), sdRand);
       muNoise_incong = rnorm(100000, (muTa + -muFl), sdRand);
       rsNoise = rnorm(100000, muRS2, sdRand);
       ssNoise = rnorm(100000, muSS, sdRand);

  int randomIndex;
      srand (1);

  //set empty matrix to store trial data
  int nRow = nTrials; //first set number of rows
  NumericMatrix trialData(nRow, 2); //set the matrix. 2 columns (RT & accuracy)

  //initialise trial flags etc. for whether certain processes have finished
  //during diffusion
  bool stimSelected; //has stimulus been selected?
    stimSelected = false;
  int whichStim; //to flag whether target (1) or flanker (2) has been selected
  double currEvidenceResp; //keep track of evidence (response selection)
    currEvidenceResp = 0.0;
  double currEvidenceStim; //keep track of evidence (stimulus selection)
    currEvidenceStim = 0.0;
  double j; //to log number of diffusion steps
    j = 0.0;


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


      randomIndex = rand() % muNoise_incong.size();

      j = j + 1.0; //update diffusion step number

      //set the correct drift rate
      if (stimSelected == false){
        if(trialType == 2){ //if trialType is incongruent
          currEvidenceResp = currEvidenceResp + muNoise_incong[randomIndex];
        }
        if(trialType == 1){ //if trialType is congruent
          currEvidenceResp = currEvidenceResp + muNoise_cong[randomIndex];
        }

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



