###
# Optimisation function for the different models

#------------------------------------------------------------------------------
# Fit function for the DSTP model
#'@export
fitFunctionDSTP <- function(humanProportions, parms, n, maxParms){


  # Get the model's predictions
  modelPrediction <- predictionsDSTP(parms, n,
                                     propsForModel = humanProportions)


  ## Change proportions in CAFs to log proportion of ERRORs. Currently they
  ## are logging proportion of CORRECT trials.
  humanProportions$conCAFProportions <- 1 -
    humanProportions$conCAFsProportions

  humanProportions$inconCAFProportions <- 1 -
    humanProportions$inconCAFsProportions

  modelPrediction$modelConCAF <- 1 - modelPrediction$modelConCAF
  modelPrediction$modelInconCAF <- 1 - modelPrediction$modelInconCAF



  # Put all human data into one vector, for ease of comparison with model's
  # prediction.
  humanProps <- c(humanProportions$conProportions,
                  humanProportions$inconProportions,
                  humanProportions$conCAFProportions,
                  humanProportions$inconCAFProportions)

  # Do the same for the model data.
  modelProps <- c(modelPrediction$modelConCDF,
                  modelPrediction$modelInconCDF,
                  modelPrediction$modelConCAF,
                  modelPrediction$modelInconCAF)

  # If any proportion is zero, change it to a very small number. This is
  # is because the fit statistic cannot handle zeros due to a division
  # by zero causing errors.
  humanProps[humanProps == 0] <- 0.00001
  modelProps[modelProps == 0] <- 0.00001

  # Do the Chi-squared test
  fitStatistic <- sum(100 * ((humanProps - modelProps) ^ 2) / modelProps)

  ##########################
  ## For Debugging
  print(fitStatistic)
  ##########################

  # If the parameters are below zero or are above maxParms, then return poor
  # fit
  if ((min(parms) < 0) | (min(maxParms - parms) < 0)){
    return(.Machine$double.xmax)
  } else {
    return(fitStatistic)
  }


}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Approximate BIC for binned data
#'@export
bBIC <- function(humanProportions, model, parms, n = 100000){

  # If the model selected is the DSTP model
  if(model == "DSTP"){

    # Get the model's predictions
    modelPrediction <- predictionsDSTP(parms, n,
                                       propsForModel = humanProportions)
  }

  ## Change proportions in CAFs to log proportion of ERRORs. Currently they
  ## are logging proportion of CORRECT trials.

  humanProportions$conCAFProportions <- 1 -
    humanProportions$conCAFsProportions

  humanProportions$inconCAFProportions <- 1 -
    humanProportions$inconCAFsProportions

  modelPrediction$modelConCAF <- 1 - modelPrediction$modelConCAF
  modelPrediction$modelInconCAF <- 1 - modelPrediction$modelInconCAF



  # Put all human data into one vector, for ease of comparison with model's
  # prediction.
  humanProps <- c(humanProportions$conProportions,
                  humanProportions$inconProportions,
                  humanProportions$conCAFProportions,
                  humanProportions$inconCAFProportions)

  # Do the same for the model data.
  modelProps <- c(modelPrediction$modelConCDF,
                  modelPrediction$modelInconCDF,
                  modelPrediction$modelConCAF,
                  modelPrediction$modelInconCAF)

  # sum the proportions of model versus human
  sumProps <- 100 * humanProps * log(modelProps)
  sumProps <- sum(sumProps)

  # set the required number of parameters
  if(model == "DSTP"){
    m <- 7
  } else {
    m <- 5
  }

  bic <- (-2 * sumProps) + m * log(100)

  return(bic)

}



#------------------------------------------------------------------------------


