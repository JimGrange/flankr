###
# Optimisation function for the different models

#------------------------------------------------------------------------------
# Fit function for the DSTP model
fitFunctionDSTP <- function(humanProportions, parms, n, maxParms){


  # Get the model's predictions
  modelPrediction <- predictionsDSTP(parms,
                                     n,
                                     propsForModel = humanProportions)


  # Put all human data into one vector, for ease of comparison with model's
  # prediction.
  humanProps <- c(humanProportions$congruentCDFProportions,
                  humanProportions$incongruentCDFProportions,
                  humanProportions$congruentCAFProportions,
                  humanProportions$incongruentCAFProportions)


  # Do the same for the model data.
  modelProps <- c(modelPrediction$modelCongruentCDF,
                  modelPrediction$modelIncongruentCDF,
                  modelPrediction$modelCongruentCAF,
                  modelPrediction$modelIncongruentCAF)


  # If any proportion is zero, change it to a very small number. This is
  # is because the fit statistic cannot handle zeros due to a division
  # by zero causing errors.
  humanProps[humanProps == 0] <- 0.0001
  modelProps[modelProps == 0] <- 0.0001

  # Calculate likelihood ratio chi-square statistic
  fitStatistic <- 2 * sum(250 * humanProps * log(humanProps / modelProps))

  if(fitStatistic == Inf){
    return(.Machine$double.xmax)
  }

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
# Fit function for the SSP model
fitFunctionSSP <- function(humanProportions, parms, n, maxParms){


  # Get the model's predictions
  modelPrediction <- predictionsSSP(parms, n,
                                    propsForModel = humanProportions)


  # Put all human data into one vector, for ease of comparison with model's
  # prediction.
  humanProps <- c(humanProportions$congruentCDFProportions,
                  humanProportions$incongruentCDFProportions,
                  humanProportions$congruentCAFProportions,
                  humanProportions$incongruentCAFProportions)


  # Do the same for the model data.
  modelProps <- c(modelPrediction$modelCongruentCDF,
                  modelPrediction$modelIncongruentCDF,
                  modelPrediction$modelCongruentCAF,
                  modelPrediction$modelIncongruentCAF)


  # If any proportion is zero, change it to a very small number. This is
  # is because the fit statistic cannot handle zeros due to a division
  # by zero causing errors.
  humanProps[humanProps == 0] <- 0.0001
  modelProps[modelProps == 0] <- 0.0001

  # Calculate likelihood ratio chi-square statistic
  fitStatistic <- 2 * sum(250 * humanProps * log(humanProps / modelProps))

  if(fitStatistic == Inf){
    return(.Machine$double.xmax)
  }

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
# optimisation for fixed parameter values for the DSTP model
# Code modified from the "optifix" function originally written by Barry
# Rowlingson:
# (http://geospaced.blogspot.co.uk/2011/10/optifix-optim-with-fixed-values.html)
optimFix_DSTP <- function(parms, fixed, humanProportions, n, maxParms){

  # which parameters are fixed/free?
  whichFixed <- parms[fixed]
  whichFree <- parms[!fixed]

  # define function to run the fit after parameters have been fixed
  fixedFit <- function(.parms, fixed, humanProportions, parms, n, maxParms){

    currParms <- rep(NA, sum(!fixed))
    currParms[!fixed] <- .parms
    currParms[fixed] <- whichFixed

    fitFunctionDSTP(humanProportions = humanProportions, parms = currParms,
                    n = n, maxParms = maxParms)

  }

  fit <- optim(whichFree, fixed = fixed, fn = fixedFit, method = "Nelder-Mead",
               humanProportions = humanProportions, n = n, maxParms = maxParms)

  # now populate the output
  fit$fullPars <- rep(NA, sum(!fixed))
  fit$fullPars[fixed] <- whichFixed
  fit$fullPars[!fixed] <- fit$par

  return(fit)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# optimisation for fixed parameter values for the SSP model
# Code modified from the "optifix" function originally written by Barry
# Rowlingson:
# (http://geospaced.blogspot.co.uk/2011/10/optifix-optim-with-fixed-values.html)
optimFix_SSP <- function(parms, fixed, humanProportions, n, maxParms){

  # which parameters are fixed/free?
  whichFixed <- parms[fixed]
  whichFree <- parms[!fixed]

  # define function to run the fit after parameters have been fixed
  fixedFit <- function(.parms, fixed, humanProportions, parms, n, maxParms){

    currParms <- rep(NA, sum(!fixed))
    currParms[!fixed] <- .parms
    currParms[fixed] <- whichFixed

    fitFunctionSSP(humanProportions = humanProportions, parms = currParms,
                   n = n, maxParms = maxParms)

  }

  fit <- optim(whichFree, fixed = fixed, fn = fixedFit, method = "Nelder-Mead",
               humanProportions = humanProportions, n = n, maxParms = maxParms)

  # now populate the output
  fit$fullPars <- rep(NA, sum(!fixed))
  fit$fullPars[fixed] <- whichFixed
  fit$fullPars[!fixed] <- fit$par

  return(fit)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# BIC for binned data
bBIC <- function(humanProportions, model, parms, nTrials){

  n = nTrials

  # If the model selected is the DSTP model
  if(model == "DSTP"){

    # Get the model's predictions
    modelPrediction <- predictionsDSTP(parms, n,
                                       propsForModel = humanProportions)
  } else {
    modelPrediction <- predictionsSSP(parms, n,
                                      propsForModel = humanProportions)
  }


  # Put all human data into one vector, for ease of comparison with model's
  # prediction.
  humanProps <- c(humanProportions$congruentCDFProportions,
                  humanProportions$incongruentCDFProportions,
                  humanProportions$congruentCAFProportions,
                  humanProportions$incongruentCAFProportions)


  # Do the same for the model data.
  modelProps <- c(modelPrediction$modelCongruentCDF,
                  modelPrediction$modelIncongruentCDF,
                  modelPrediction$modelCongruentCAF,
                  modelPrediction$modelIncongruentCAF)




  # If any proportion is zero, change it to a very small number. This is
  # is because the fit statistic cannot handle zeros due to a division
  # by zero causing errors.
  humanProps[humanProps == 0] <- 0.0001
  modelProps[modelProps == 0] <- 0.0001


  # sum the proportions of model versus human
  sumProps <- 250 * humanProps * log(modelProps)
  sumProps <- sum(sumProps)

  # set the required number of parameters
  if(model == "DSTP"){
    m <- 7
  } else {
    m <- 5
  }

  bic <- -2 * sumProps + m * log(250)

  return(bic)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# BIC for binned data with fixed model parameters
bBIC_fixed <- function(humanProportions, model, parms, fixed, nTrials){


  n = nTrials

  # If the model selected is the DSTP model
  if(model == "DSTP"){

    # Get the model's predictions
    modelPrediction <- predictionsDSTP(parms, n,
                                       propsForModel = humanProportions)
  } else {
    modelPrediction <- predictionsSSP(parms, n,
                                      propsForModel = humanProportions)
  }


  # Put all human data into one vector, for ease of comparison with model's
  # prediction.
  humanProps <- c(humanProportions$congruentCDFProportions,
                  humanProportions$incongruentCDFProportions,
                  humanProportions$congruentCAFProportions,
                  humanProportions$incongruentCAFProportions)


  # Do the same for the model data.
  modelProps <- c(modelPrediction$modelCongruentCDF,
                  modelPrediction$modelIncongruentCDF,
                  modelPrediction$modelCongruentCAF,
                  modelPrediction$modelIncongruentCAF)




  # If any proportion is zero, change it to a very small number. This is
  # is because the fit statistic cannot handle zeros due to a division
  # by zero causing errors.
  humanProps[humanProps == 0] <- 0.0001
  modelProps[modelProps == 0] <- 0.0001


  # sum the proportions of model versus human
  sumProps <- 250 * humanProps * log(modelProps)
  sumProps <- sum(sumProps)

  # set the required number of free parameters
  m <- length(fixed) - sum(fixed)

  bic <- -2 * sumProps + m * log(250)

  return(bic)
}
#------------------------------------------------------------------------------
