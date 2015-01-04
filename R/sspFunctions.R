###
# Functions to run the SSP model itself. This includes functions to simulate
# data from the SSP model, as well as to run the fitting routine.

#------------------------------------------------------------------------------
#' Obtain simulated response times and accuracy from the SSP model
#'
#' \code{simulateSSP} generates synthetic data from the DSTP model in the
#' form of response time (RT) in seconds and accuracy for both congruent and
#' incongruent trials.
#'
#' This function can be employed by the user to generate synthetic data, but
#' its main purpose is to be used by the fitting procedure to generate model
#' predictions for a set of parameter values when trying to find the best-
#' fitting values.
#'
#' @param parms The set of parameters to use to simulate the data. Must be
#' contained in a vector in the order: \code{A}, \code{ter},
#' \code{p}, \code{rd}, \code{sda}.
#' @param nTrials How many trials to simulate per congruency condition.
#' @param var The variance of the diffusion process. By default this is set to
#' 0.01.
#' @param dt The diffusion scaling parameter (i.e., time steps). By default,
#' this is set to 0.001.
#' @param seed The value for the \code{set.seed} function to set random
#' generation state.
#'
#' @examples
#' # declare the parameters
#' parms <- c(0.050, 0.300, 0.400, 0.040, 1.500)
#'
#' # simulate the data
#' modelData <- simulateSSP(parms, nTrials = 10000)
#'
#' @return Returns a data frame with three columns: rt (response time) in
#' seconds, accuracy of the model's response (1 for correct, 0 for error), and
#' congruency condition.
#' @useDynLib flankr
#' @importFrom Rcpp sourceCpp
#' @export
simulateSSP <- function(parms,  nTrials, var = 0.01, dt = 1/1000, seed = 42){

  # transfer nTrials to shorter name
  n <- nTrials

  # Set random number seed, so same predictions occur every time. By default
  # this is set for the user.
  set.seed(seed, kind = NULL, normal.kind = NULL)

  # initialise empty matrix for simulation data with two columns
  # (RT & accuracy) and with rows = number of trials
  trialData <- matrix(0, nrow = n * 2, ncol = 3)
  colnames(trialData) <- c("rt", "accuracy", "congruency")
  trialData <- data.frame(trialData)

  # first generate congruent data by calling the C++ function
  trialData[1:n, 1:2] <- getSSP(parms, trialType = 1, nTrials = n, dt, var)
  trialData[1:n, 3] <- "congruent"

  # now do incongruent data
  trialData[(n + 1):(n * 2), 1:2] <- getSSP(parms, trialType = 2, nTrials = n,
                                            dt, var)
  trialData[(n + 1):(n * 2), 3] <- "incongruent"


  return(trialData);

}  # end of function
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#' Fit the SSP model to human data
#'
#' \code{fitSSP} fits the SSP model to a single experimental condition of
#' human data (besides congruency, which it accounts for simutaneously).
#'
#' This function can be employed by the user to find the best-fitting
#' parameters of the SSP model to fit the human data of a single experimental
#' condition. The fitting procedure accounts for congruent and incongruent
#' trials simultaneously. The fit is obtained by a gradient-descent method
#' (using the Nelder-Mead method contained in R's \code{optim} function) and is
#' fit to the proportion of data contained in human CDF and CAF distributional
#' data.
#'
#' @param data A data frame containing human data. See \code{?exampleData} for
#' data formatted correctly.
#'
#' @param conditionName If there is an additional experimental manipulation
#' (i.e., other than target congruency) the model can only be fit to one at a
#' time. Tell the function which condition is currently being fit by passing
#' a string to the function (e.g., "present"). The function by default assumes
#' no additional condition (e.g., conditionName is set to NULL).
#'
#' @param cdfs A vector of quantile values for cumulative distribution functions
#' to be estimated from the human data. The model will attempt to find the
#' best-fitting parameters that match this distributional data.
#'
#' @param cafs A vector of quantiles for conditional accuracy functions to be
#' estimated from the human data. The model will attempt to find the best-
#' fitting parameters that match this distributional data.
#'
#' @param maxParms A vector containing upper limits on possible parameter
#' values.
#'
#' @param nTrials An integer stating how many trials to simulate per iteration
#' of the fitting cycle for each congruency type.
#'
#' @param multipleSubjects A boolean stating whether the fit is to multiple
#' subjects (multipleSubjects = TRUE) or to a single subject
#' (multipleSubjects = FALSE).
#'
#' @return \code{bestParameters} A vector of the best-fitting parameters found
#' by the current fit run.
#'
#' @return \code{chiSquared} The value of chi square obtained by the current
#' fit run.
#'
#' @return \code{bBIC} The value of the  Bayesian Information Criterion (BIC)
#' obtained by the current fit run. This is calculated using the BIC equation
#' for binned data, hence bBIC (binned BIC).
#'
#' @examples
#' # Load the example data the comes with the \code{flankr} package
#' data(exampleData)
#'
#' # Fit the model to the condition "present" in the example data set using
#' # the default settings in the model.
#'
#' fit <- fitSSP(data = exampleData, conditionName = "present")
#'
#' # Fit the model using different CDF and CAF values, and 100,000 trials per
#' # fit cycle
#' cdfs <- c(.2, .4, .6, .8)
#' cafs <- c(.2, .4, .6, .8)
#'
#' fit <- fitSSP(exampleData, conditionName = "present", cdfs = cdfs,
#'               cafs = cafs, nTrials = 100000)
#'
#'
#'@export
fitSSP<- function(data, conditionName = NULL,
                  parms = c(0.050, 0.300, 0.400, 0.050, 1.500),
                  cdfs = c(.1, .3, .5, .7, .9), cafs = c(.25, .50, .75),
                  maxParms = c(1, 1, 1, 1, 3), nTrials = 50000,
                  multipleSubjects = TRUE){


  # get the desired condition's data
  if(is.null(conditionName)){
    conditionData <- data
  } else{
    conditionData <- subset(data, data$condition == conditionName)
  }

  # get all of the distribution & proportion information from human data.
  # This returns a list with all information in separate "cotainers" for ease
  # of access & generalisation to different CDF and CAF sizes.
  if(multipleSubjects == TRUE){
    humanProportions <- getHumanProps(conditionData, cdfs, cafs)
  } else {
    humanProportions <- getHumanPropsSingle(conditionData, cdfs, cafs)
  }

  modelStart <- "Model Fit Running. Please Wait..."
  print(modelStart)

  # perform the fit
  fit <- optim(parms, fn = fitFunctionSSP, humanProportions = humanProportions,
               n = nTrials, maxParms = maxParms)

  # what are the best-fitting parameters?
  bestParameters <- round(fit$par, 3)

  # what is the fit statistic value?
  chiSquare <- fit$value

  # get the approximate BIC value
  bBIC <- bBIC(humanProportions, model = "SSP", parms = bestParameters)

  # put all results into a list, and return the list to the user
  modelFit <- list(bestParameters = bestParameters, chiSquare = chiSquare,
                   bBIC = bBIC)

  modelFinished <- "Model Fit Finished."
  print(modelFinished)

  return(modelFit)


} # end of function
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Get the predicted proportions from the SSP model

#'@export
predictionsSSP<- function(parms, n, propsForModel, dt = 0.001, var = 0.01){

  # parms = parameters for the model run
  # n = number of trials per congruency condition
  # propsForModel = CDF & CAF distributional information

  # Run model to get congruent RTs
  set.seed(42)
  modelCon <- getSSP(parms, trialType = 1, n = n, dt, var)
  modelConCDF <- getCDFProps(propsForModel$conCDF, modelCon)
  modelConCAF <- getCAFProps(propsForModel$conCAFsCutoff, modelCon)

  # Run model to get incontruent RTs
  set.seed(42)
  modelIncon <- getSSP(parms, trialType = 2, n = n, dt, var)
  modelInconCDF <- getCDFProps(propsForModel$inconCDF, modelIncon)
  modelInconCAF <- getCAFProps(propsForModel$inconCAFsCutoff, modelIncon)

  modelProps <- list(modelConCDF = modelConCDF,
                     modelConCAF = modelConCAF,
                     modelInconCDF = modelInconCDF,
                     modelInconCAF = modelInconCAF)


  return(modelProps)

}
#------------------------------------------------------------------------------







