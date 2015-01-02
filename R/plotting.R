###
# Plotting functions

#------------------------------------------------------------------------------
# Plot the fit of the DSTP model

#'@export
plotFitDSTP <- function(modelFit, data, conditionName, nTrials = 100000,
                          cdfs = c(.1, .3, .5, .7, .9), cafs = c(.25, .50, .75)){

  # Change the plotting window
  par(mfrow = c(1, 2))



  # get the desired condition's data-------------------------------------------
  conditionData <- subset(data, data$condition == conditionName)

  # get all of the distribution & proportion information from human data.
  # This returns a list with all information in separate "cotainers" for ease
  # of access & generalisation to different CDF and CAF sizes.
  humanProportions <- getHumanProps(conditionData, cdfs, cafs)




  # get model proportions by running the model---------------------------------

  # First, what were the best-fitting parameters?
  parms <- modelFit$bestParameters

  # simulate the DSTP model with these parameters
  modelData <- predictionsDSTP(parms, nTrials,
                               propsForModel = humanProportions)

  # Find model proportion predictions for
  # congruent and incongruent trials (CDF & CAFs)
  modelConCDF <- proportionCDFs(modelData$modelConCDF)
  modelConCAF <- modelData$modelConCAF

  modelInconCDF <- proportionCDFs(modelData$modelInconCDF)
  modelInconCAF <- modelData$modelInconCAF


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
  legend(maxRT - 0.25, 0.65, c("Congruent","Incongruent"), cex=1, pch=c(19, 1), lty=1:2, bty="n");

  # Change the plotting window
  par(mfrow = c(1, 1))

  ### FOR DEBUGGING
  return(modelData)

}
#------------------------------------------------------------------------------
