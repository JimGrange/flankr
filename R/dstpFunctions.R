###
# Functions to run the DSTP model itself. This includes functions to simulate
# data from the DSTP model, as well as to run the fitting routine.

#------------------------------------------------------------------------------
#' Obtain simulated response times and accuracy from the DSTP model
#'
#' \code{simulate.DSTP} generates synthetic data from the DSTP model in the
#' form of response time (RT) in seconds and accuracy for both congruent and
#' incongruent trials.
#'
#' This function can be employed by the user to generate synthetic data, but
#' its main purpose is to be used by the fitting procedure to generate model
#' predictions for a set of parameter values when trying to find the best-
#' fitting values.
#'
#' @param parms The set of parameters to use to simulate the data. Must be
#' contained in a vector in the order: \code{A}, \code{C}, \code{driftTarget},
#' \code{driftFlanker}, \code{diftStimSelection}, \code{driftRS2}, \code{ter}.
#' @param n How many trials to simulate per congruency condition.
#' @param var The variance of the diffusion process. By default this is set to
#' 0.01.
#' @param dt The diffusion scaling parameter (i.e., time steps). By default,
#' this is set to 0.001.
#' @param seed The value for the \code{set.seed} function to set random
#' generation state.
#'
#' @examples
#' # declare the parameters
#' parms <- c(0.070, 0.086, 0.045, 0.065, 0.368, 1.575, 0.225)
#'
#' # simulate the data
#' modelData <- simulateDSTP(parms, n = 1000)
#'
#' @return Returns a data frame with three columns: RT (response time) in
#' seconds, accuracy of the model's response (1 for correct, 0 for error), and
#' congruency condition.
#' @useDynLib flankr
#' @importFrom Rcpp sourceCpp
#' @export
simulateDSTP <- function(parms,  n, var = 0.01, dt = 1/1000, seed = 10){

  # Set random number seed, so same predictions occur every time. By default
  # this is set for the user.
  set.seed(seed, kind = NULL, normal.kind = NULL)

  # initialise empty matrix for simulation data with two columns
  # (RT & accuracy) and with rows = number of trials
  trialData <- matrix(0, nrow = n * 2, ncol = 3)
  colnames(trialData) <- c("RT", "Accuracy", "Congruency")
  trialData <- data.frame(trialData)

  # first generate congruent data by calling the C++ function
  trialData[1:n, 1:2] <- getDSTP(parms, trialType = 1, nTrials = n, dt, var)
  trialData[1:n, 3] <- "congruent"

  # now do incongruent data
  trialData[(n + 1):(n * 2), 1:2] <- getDSTP(parms, trialType = 2, nTrials = n,
                                             dt, var)
  trialData[(n + 1):(n * 2), 3] <- "incongruent"


  return(trialData);

}  # end of function
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#' fit the DSTP model to human data
#'
#' \code{fitDSTP} fits the DSTP model to human data.
#'
#' This function can be employed by the user to find best-fitting parameters of
#' the DSTP model to human data.
#'
#' @param data A data frame containing human data. See \code{?exampleData} for
#' data formatted correctly.
#'
#' @param conditionName A string (e.g., "present") which states which condition the
#' model is to fit currently.
#'
#' @param cdfs A vector of quantile values for cumulative distribution functions
#' to be estimated from the human data. The model will attempt to find the
#' best-fitting parameters that match this distributional data.
#'
#' @param cafs A vector of quantiles for conditional accuracy functions to be
#' estimated from the human data. The model will attempt to find the best-
#' fitting parameters that match this distributional data.
#'
#' @param nTrials An integer stating how many trials to simulate per iteration
#' of the fitting cycle for each congruency type.
#'
#' @export
fitDSTP <- function(data, conditionName, cdfs = c(.1, .3, .5, .7, .9),
                    cafs = c(.25, .50, .75),
                    parms = c(0.145, 0.08, 0.10, 0.07, 0.325, 1.30, 0.240),
                    maxParms = c(1, 1, 1, 1, 1, 2, 1), nTrials = 50000){


  # get the desired condition's data
  conditionData <- subset(data, data$condition == conditionName)

  # get all of the distribution & proportion information from human data.
  # This returns a list with all information in separate "cotainers" for ease
  # of access & generalisation to different CDF and CAF sizes.
  humanProportions <- getHumanProps(conditionData, cdfs, cafs)

  # perform the fit
  fit <- optim(parms, fn = fitFunctionDSTP, humanProportions = humanProportions,
               n = nTrials, maxParms = maxParms)


} # end of function
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Get the predicted proportions from the DSTP model

#'@export
predictionsDSTP <- function(parms, n, propsForModel, dt = 0.001, var = 0.01){

  # parms = parameters for the model run
  # n = number of trials per congruency condition
  # propsForModel = CDF & CAF distributional information

  # Run model to get congruent RTs
  set.seed(42)
  modelCon <- getDSTP(parms, trialType = 1, n = n, dt, var)
  modelConCDF <- getCDFProps(propsForModel$conCDF, modelCon)
  modelConCAF <- getCAFProps(propsForModel$conCAFsCutoff, modelCon)

  # Run model to get incontruent RTs
  set.seed(42)
  modelIncon <- getDSTP(parms, trialType = 2, n = n, dt, var)
  modelInconCDF <- getCDFProps(propsForModel$inconCDF, modelIncon)
  modelInconCAF <- getCAFProps(propsForModel$inconCAFsCutoff, modelIncon)

  modelProps <- list(modelConCDF = modelConCDF,
                     modelConCAF = modelConCAF,
                     modelInconCDF = modelInconCDF,
                     modelInconCAF = modelInconCAF)


  return(modelProps)

}
#------------------------------------------------------------------------------







