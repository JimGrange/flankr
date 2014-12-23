#------------------------------------------------------------------------------
#' Obtain simulated response times and accuracy from the DSTP model
#'
#' \code{simulateDSTP} generates synthetic data from the DSTP model in the
#' form of response time (RT) in seconds and accuracy.
#'
#' This function can be employed by the user to generate synthetic data, but
#' its main purpose is to be used by the fitting procedure to generate model
#' predictions for a set of parameter values when trying to find the best-
#' fitting values.
#'
#' @param parms The set of parameters to use to simulate the data. Must be in
#' the order: \code{A}, \code{C}, \code{driftTarget}, \code{driftFlanker},
#' \code{diftStimSelection}, \code{driftRS2}, \code{ter}.
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
#' # simulate  data
#' getData <- simulateDSTP(parms, n = 1000)
#'
#' @return Returns a data frame with three columns: RT (response time) in
#' seconds, accuracy of the model's response (1 for correct, 0 for error), and
#' congruency condition.
#' @useDynLib flankerModels
#' @importFrom Rcpp sourceCpp
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#' @export
simulateDSTP <- function(parms,  n, var = 0.01, dt = 1/1000,
                         seed = 10){

  # Set random number seed, so same predictions occur every time. By default
  # this is set for the user.
  set.seed(seed, kind = NULL, normal.kind = NULL)

  # initialise empty matrix for simulation data with two columns
  # (RT & accuracy) and with rows = number of trials
  trialData <- matrix(0, nrow = n * 2, ncol = 3)
    colnames(trialData) <- c("RT", "Accuracy", "Congruency")
    trialData <- data.frame(trialData)

  # first generate congruent data by calling C++ function
  trialData[1:n, 1:2] <- getDSTP(parms, trialType = 1, n = n, dt, var)
    trialData[1:n, 3] <- "congruent"

  # now do incongruent data
  trialData[(n + 1):(n * 2), 1:2] <- getDSTP(parms, trialType = 2, n = n, dt, var)
    trialData[(n + 1):(n * 2), 3] <- "incongruent"


  return(trialData);

}  # end of function
#------------------------------------------------------------------------------
