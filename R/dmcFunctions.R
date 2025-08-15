
# simulate DMC ------------------------------------------------------------


#' Obtain simulated response times and accuracy from the DMC model
#'
#' \code{simulateDMC} generates synthetic data from the DMC model in the
#' form of response time (RT) in seconds and accuracy for both congruent and
#' incongruent trials.
#'
#' This function can be employed by the user to generate synthetic data, but
#' its main purpose is to be used by the fitting procedure to generate model
#' predictions for a set of parameter values when trying to find the best-
#' fitting values.d
#'
#' @param parms The set of parameters to use to simulate the data. Must be
#' contained in a vector in the order: \code{b}, \code{driftControlled},
#' \code{maxAmp}, \code{tau}, \code{alpha}, \code{ter}, \code{ster},
#' \code{sigma}.
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
#' \dontrun{
#' parms <- c(50, 0.69, 19.20, 118, 2.15,331,36,3.98)
#'
#' # simulate the data
#' modelData <- simulateDMC(parms, nTrials = 10000)
#'}
#'
#' @return Returns a data frame with three columns: rt (response time) in
#' seconds, accuracy of the model's response (1 for correct, 0 for error), and
#' congruency condition.
#' @useDynLib flankr
#' @importFrom Rcpp sourceCpp
#' @export
simulateDMC <- function(parms,
                        nTrials,
                        dt = 0.1,
                        seed = 42){

  # transfer nTrials to shorter name
  n <- nTrials

  # Set random number seed, so same predictions occur every time.
  set.seed(seed, kind = NULL, normal.kind = NULL)

  # initialise empty matrix for simulation data with two columns
  # (RT & accuracy) and with rows = number of trials
  trialData <- matrix(0, nrow = n * 2, ncol = 3)
  colnames(trialData) <- c("rt", "accuracy", "congruency")
  trialData <- data.frame(trialData)

  # first generate congruent data by calling the C++ function
  trialData[1:n, 1:2] <- getDMC(parms,
                                trialType = 1,
                                nTrials = n,
                                dt)
  trialData[1:n, 3] <- "congruent"

  # now do incongruent data
  trialData[(n + 1):(n * 2), 1:2] <- getDMC(parms,
                                            trialType = 2,
                                            nTrials = n,
                                            dt)

  trialData[(n + 1):(n * 2), 3] <- "incongruent"


  return(trialData);

}  # end of function
