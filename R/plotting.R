###
# Plotting functions

#------------------------------------------------------------------------------
# Plot the fit of the DSTP model

#'Plot the fit of the DSTP model to human data.
#'
#'\code{plotFitDSTP} will plot the fit of the model to human distributional
#'data.
#'
#'This function is passed the object obtained by the model fitting procedure,
#'as well as the human data and the condition that was fitted by the routine.
#'The function simulates 100,000 trials (by default) using the best-fitting
#'parameters found by the fit procedure. This synthetic data is then considered
#'as the model's best predictions. The function then provides a plot of the
#'model fit to cumulative distribution functions (CDFs) of correct response
#'time, and conditional accuracy functions (CAFs) to show fit to accuracy data.
#'The function also returns the data used to plot the fit so that the user can
#'use their own plotting methods.
#'
#' @param modelFit The object obtained by the model fit.
#'
#' @param data The data frame of human data.
#' @param conditionName The name of the condition that was fit. By default,
#' it is set to conditionName = NULL.
#' @param nTrials How many trials used to generate the model's best
#' predictions. This should be higher than that used to fit the model.
#' @param cdfs The cut-off points for the cumulative distribution functions.
#' @param cafs The cut-off points for the conditional accuracy functions.
#' @param multipleSubjects A boolean stating whether the fit is to multiple
#' subjects (multipleSubjects = TRUE) or to a single subject
#' (multipleSubjects = FALSE).
#'
#' @return \code{cdfs} The CDF values requested by the user.
#' @return \code{cafs} The CAF values requested by the user.
#' @return \code{humanConCDFs} The response time cut-off values for each CDF
#' bin for congruent human data.
#' @return \code{humanInconCDFs} The response time cut-off values for each CDF
#' bin for incongruent human data.
#' @return \code{humanConCAFsRT} The mean response times for each bin of the
#' CAF functions for congruent human data.
#' @return \code{humanInconCAFsRT} The mean response times for each bin of the
#' CAF functions for incongruent human data.
#' @return \code{humanConCAFsError} The percent accuracy for each bin of the
#' CAF functions for congruent human data.
#' @return \code{humanConCAFsError} The percent accuracy for each bin of the
#' CAF functions for congruent human data.
#' @return \code{modelConCDFs} The quantile cut-off points for the model
#' predictions for congruent data. A perfect fit would match the cdfs asked
#' for by the user (e.g., .1, .3, .5, .7, .9).
#' @return \code{modelInconCDFs} The quantile cut-off points for the model
#' predictions for incongruent data. A perfect fit would match the cdfs asked
#' for by the user (e.g., .1, .3, .5, .7, .9).
#' @return \code{modelConCAFs} The percentage accuracy predicted for each CAF
#' bin by the model for congruent data.
#' @return \code{modelInconCAFs} The percentage accuracy predicted for each CAF
#' bin by the model for incongruent data.
#'
#'@examples
#'# Assume that the model was just fit to the data contained in
#'# \code{exampleData} (condition "present") and saved to the variable called
#'# "fit", then we can obtain a plot of that fit by:
#'
#' plot <- plotFitDSTP(modelFit = fit, data = exampleData,
#'                     conditionName = "present")
#'
#'# We can also change the default CDF and CAF quantiles used, as well as the
#'# number of trials used to simulate the fitted data.
#'
#' plot <- plotFitDSTP(modelFit = fit, data = exampleData,
#'                     conditioName = "present", cdfs = c(.2, .4, .6, .8),
#'                     cafs = c(.2, .4, .6, .8), nTrials = 500000)
#'
#'

#'@export
plotFitDSTP <- function(modelFit, data, conditionName = NULL, nTrials = 50000,
                        cdfs = c(.1, .3, .5, .7, .9), cafs = c(.25, .50, .75),
                        multipleSubjects = TRUE){

  # Change the plotting window
  par(mfrow = c(1, 2))



  # get the desired condition's data-------------------------------------------
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




  # get model proportions by running the model---------------------------------

  # First, what were the best-fitting parameters?
  parms <- modelFit$bestParameters

  # simulate the DSTP model with these parameters
  modelData <- plotPredictionsDSTP(parms, nTrials,
                                   propsForModel = humanProportions)

  # Find model proportion predictions for
  # congruent and incongruent trials (CDF & CAFs)
  modelConCDF <- modelData$modelCongruentCDF
  modelConCAF <- modelData$modelCongruentCAF

  modelInconCDF <- modelData$modelIncongruentCDF
  modelInconCAF <- modelData$modelIncongruentCAF


  # Generate the return data to go back to the user ---------------------------
  # Pass the human proportions to an object with a shorter name to save typing
  x <- humanProportions
  returnData <- list(cdfs = cdfs,
                     cafs = cafs,
                     humanConCDFs = x$conCDFs,
                     humanInconCDFs = x$inconCDFs,
                     humanConCAFsRT = x$conCAFsRT,
                     humanInconCAFsRT = x$inconCAFsRT,
                     humanConCAFsError = x$conCAFsError,
                     humanInconCAFsError = x$inconCAFsError,
                     modelConCDFs = modelConCDF,
                     modelInconCDFs = modelInconCDF,
                     modelConCAFs = modelConCAF,
                     modelInconCAFs = modelInconCAF)


  # Plot the CDFs--------------------------------------------------------------

  # first, find the response time boundaries to ensure plots are in bounds
  minRT <- min(humanProportions$congruentCDFs, humanProportions$incongruentCDFs)
  maxRT <- max(humanProportions$congruentCDFs, humanProportions$incongruentCDFs)

  # Incongruent human first
  plot(x = humanProportions$incongruentCDFs, y = cdfs, xlab = "Response Time",
       ylab = "Cumulative Probability", pch = 1, ylim = c(0, 1),
       xlim = c(minRT, maxRT))
  lines(humanProportions$incongruentCDFs, modelInconCDF, type = "l",
        lty = 2)

  # Now congruent
  points(x = humanProportions$congruentCDFs, y = cdfs, pch = 19)
  lines(humanProportions$congruentCDFs, modelConCDF, type = "l", lty = 1)


  # Plot the CAFs--------------------------------------------------------------
  minRT <- min(humanProportions$congruentCAFsRT, humanProportions$incongruentCAFsRT)
  maxRT <- max(humanProportions$congruentCAFsRT, humanProportions$incongruentCAFsRT)

  # Incongruent data
  plot(x = humanProportions$incongruentCAFsRT, y = humanProportions$incongruentCAFsError,
       xlab = "Response Time", ylab = "Accuracy", pch = 1, ylim = c(0.5, 1),
       xlim = c(minRT, maxRT))
  lines(humanProportions$incongruentCAFsRT, modelInconCAF, type = "l", lty = 2)

  # Congruent data
  points(x = humanProportions$congruentCAFsRT, y = humanProportions$congruentCAFsError,
         pch = 19)
  lines(humanProportions$congruentCAFsRT, modelConCAF, type = "l", lty = 1)

  # Add legend
  legend(maxRT - 0.25, 0.65, c("Congruent","Incongruent"), cex=1, pch=c(19, 1),
         lty=1:2, bty="n");

  # Change the plotting window
  par(mfrow = c(1, 1))

  # Return the information used to plot the model so the user can use their own
  # software should they so please.
  return(returnData)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Plot the fit of the SSP model

#'Plot the fit of the SSP model to human data.
#'
#'\code{plotFitSSP} will plot the fit of the model to human distributional
#'data.
#'
#'This function is passed the object obtained by the model fitting procedure,
#'as well as the human data and the condition that was fitted by the routine.
#'The function simulates 100,000 trials (by default) using the best-fitting
#'parameters found by the fit procedure. This synthetic data is then considered
#'as the model's best predictions. The function then provides a plot of the
#'model fit to cumulative distribution functions (CDFs) of correct response
#'time, and conditional accuracy functions (CAFs) to show fit to accuracy data.
#'The function also returns the data used to plot the fit so that the user can
#'use their own plotting methods.
#'
#' @param modelFit The object obtained by the model fit.
#'
#' @param data The data frame of human data.
#' @param conditionName The name of the condition that was fit. By default,
#' it is set to conditionName = NULL.
#' @param nTrials How many trials used to generate the model's best
#' predictions. This should be higher than that used to fit the model.
#' @param cdfs The cut-off points for the cumulative distribution functions.
#' @param cafs The cut-off points for the conditional accuracy functions.
#' @param multipleSubjects A boolean stating whether the fit is to multiple
#' subjects (multipleSubjects = TRUE) or to a single subject
#' (multipleSubjects = FALSE).
#'
#' @return \code{cdfs} The CDF values requested by the user.
#' @return \code{cafs} The CAF values requested by the user.
#' @return \code{humanConCDFs} The response time cut-off values for each CDF
#' bin for congruent human data.
#' @return \code{humanInconCDFs} The response time cut-off values for each CDF
#' bin for incongruent human data.
#' @return \code{humanConCAFsRT} The mean response times for each bin of the
#' CAF functions for congruent human data.
#' @return \code{humanInconCAFsRT} The mean response times for each bin of the
#' CAF functions for incongruent human data.
#' @return \code{humanConCAFsError} The percent accuracy for each bin of the
#' CAF functions for congruent human data.
#' @return \code{humanConCAFsError} The percent accuracy for each bin of the
#' CAF functions for congruent human data.
#' @return \code{modelConCDFs} The quantile cut-off points for the model
#' predictions for congruent data. A perfect fit would match the cdfs asked
#' for by the user (e.g., .1, .3, .5, .7, .9).
#' @return \code{modelInconCDFs} The quantile cut-off points for the model
#' predictions for incongruent data. A perfect fit would match the cdfs asked
#' for by the user (e.g., .1, .3, .5, .7, .9).
#' @return \code{modelConCAFs} The percentage accuracy predicted for each CAF
#' bin by the model for congruent data.
#' @return \code{modelInconCAFs} The percentage accuracy predicted for each CAF
#' bin by the model for incongruent data.
#'
#'@examples
#'# Assume that the model was just fit to the data contained in
#'# \code{exampleData} (condition "present") and saved to the variable called
#'# "fit", then we can obtain a plot of that fit by:
#'
#' plot <- plotFitSSP(modelFit = fit, data = exampleData,
#'                    conditionName = "present")
#'
#'# We can also change the default CDF and CAF quantiles used, as well as the
#'# number of trials used to simulate the fitted data.
#'
#' plot <- plotFitSSP(modelFit = fit, data = exampleData,
#'                    conditioName = "present", cdfs = c(.2, .4, .6, .8),
#'                    cafs = c(.2, .4, .6, .8), nTrials = 500000)
#'
#'

#'@export
plotFitSSP <- function(modelFit, data, conditionName = NULL, nTrials = 50000,
                       cdfs = c(.1, .3, .5, .7, .9), cafs = c(.25, .50, .75),
                       multipleSubjects = TRUE){

  # Change the plotting window
  par(mfrow = c(1, 2))



  # get the desired condition's data-------------------------------------------
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




  # get model proportions by running the model---------------------------------

  # First, what were the best-fitting parameters?
  parms <- modelFit$bestParameters

  # simulate the DSTP model with these parameters
  modelData <- predictionsSSP(parms, nTrials,
                              propsForModel = humanProportions)

  # Find model proportion predictions for
  # congruent and incongruent trials (CDF & CAFs)
  modelConCDF <- proportionCDFs(modelData$modelConCDF)
  modelConCAF <- modelData$modelConCAF

  modelInconCDF <- proportionCDFs(modelData$modelInconCDF)
  modelInconCAF <- modelData$modelInconCAF


  # Generate the return data to go back to the user ---------------------------
  # Pass the human proportions to an object with a shorter name to save typing
  x <- humanProportions
  returnData <- list(cdfs = cdfs,
                     cafs = cafs,
                     humanConCDFs = x$conCDFs,
                     humanInconCDFs = x$inconCDFs,
                     humanConCAFsRT = x$conCAFsRT,
                     humanInconCAFsRT = x$inconCAFsRT,
                     humanConCAFsError = x$conCAFsError,
                     humanInconCAFsError = x$inconCAFsError,
                     modelConCDFs = modelConCDF,
                     modelInconCDFs = modelInconCDF,
                     modelConCAFs = modelConCAF,
                     modelInconCAFs = modelInconCAF)


  # Plot the CDFs--------------------------------------------------------------

  # first, find the response time boundaries to ensure plots are in bounds
  minRT <- min(humanProportions$conCDFs, humanProportions$inconCDFs)
  maxRT <- max(humanProportions$conCDFs, humanProportions$inconCDFs)

  # Incongruent human first
  plot(x = humanProportions$inconCDFs, y = cdfs, xlab = "Response Time",
       ylab = "Cumulative Probability", pch = 1, ylim = c(0, 1),
       xlim = c(minRT, maxRT))
  lines(humanProportions$inconCDFs, modelInconCDF, type = "l",
        lty = 2)

  # Now congruent
  points(x = humanProportions$conCDFs, y = cdfs, pch = 19)
  lines(humanProportions$conCDFs, modelConCDF, type = "l", lty = 1)


  # Plot the CAFs--------------------------------------------------------------
  minRT <- min(humanProportions$conCAFsRT, humanProportions$inconCAFsRT)
  maxRT <- max(humanProportions$conCAFsRT, humanProportions$inconCAFsRT)

  # Incongruent data
  plot(x = humanProportions$inconCAFsRT, y = humanProportions$inconCAFsError,
       xlab = "Response Time", ylab = "Accuracy", pch = 1, ylim = c(0.5, 1),
       xlim = c(minRT, maxRT))
  lines(humanProportions$inconCAFsRT, modelInconCAF, type = "l", lty = 2)

  # Congruent data
  points(x = humanProportions$conCAFsRT, y = humanProportions$conCAFsError,
         pch = 19)
  lines(humanProportions$conCAFsRT, modelConCAF, type = "l", lty = 1)

  # Add legend
  legend(maxRT - 0.25, 0.65, c("Congruent","Incongruent"), cex=1, pch=c(19, 1),
         lty=1:2, bty="n");

  # Change the plotting window
  par(mfrow = c(1, 1))

  # Return the information used to plot the model so the user can use their own
  # software should they so please.
  return(returnData)

}
#------------------------------------------------------------------------------
