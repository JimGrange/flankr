###
# This file contains model preparation functions, such as calculating
# proportions in each distributional bin of human data, finding proportions
# predicted by a model etc.

#------------------------------------------------------------------------------
# A function to get human distributional bin proportions from data, given
# desired CDF and CAF qauntile values. These values will be used by the model
# during fitting.

#' @export
getHumanProps <- function(conditionData, cdfs, cafs){


  # split trials on congruency
  congruentData <- subset(conditionData,
                          conditionData$congruency == "congruent")
  incongruentData <- subset(conditionData,
                            conditionData$congruency == "incongruent")


  # CDFs-----------------------------------------------------------------------
  # get correct RT CDF proportions. They are a function of the desired CDF
  # quantile points, and the overall accuracy for the congruency condition
  congruentProps <- cdfBinsize(cdfs)
  incongruentProps <- cdfBinsize(cdfs)

  # scale by the proportion correct
  congruentCDFProps <- cdfProportions(congruentData, congruentProps)
  incongruentCDFProps <- cdfProportions(incongruentData, incongruentProps)

  # get the CDF cut-off values of response time
  congruentCDFs <- cdf(congruentData, quantiles = cdfs)
  incongruentCDFs <- cdf(incongruentData, quantiles = cdfs)

  # CAFs-----------------------------------------------------------------------
  # get the CAF proportions
  congruentCAFProps <- cafProportions(congruentData, quantiles = cafs)
  incongruentCAFProps <- cafProportions(incongruentData, quantiles = cafs)

  # get the CAF cut-off points
  congruentCAFsCutoff <- cdf(congruentData, quantiles = cafs,
                             correctTrials = 3)
  incongruentCAFsCutoff  <- cdf(incongruentData, quantiles = cafs,
                                correctTrials = 3)

  # now get the CAFs themselves
  congruentCAFs <- caf(congruentData, quantiles = cafs)
  incongruentCAFs <- caf(incongruentData, quantiles = cafs)


  #collate data and return-----------------------------------------------------
  humanProportions <- list(cdfQuantiles = cdfs,
                           cafQuantiles = cafs,
                           congruentCDFProportions = congruentCDFProps,
                           incongruentCDFProportions = incongruentCDFProps,
                           congruentCAFProportions = congruentCAFProps,
                           incongruentCAFProportions = incongruentCAFProps,
                           congruentCDFs = congruentCDFs,
                           incongruentCDFs = incongruentCDFs,
                           congruentCAFsCutoff = congruentCAFsCutoff,
                           incongruentCAFsCutoff =  incongruentCAFsCutoff,
                           congruentCAFsRT = congruentCAFs[1, ],
                           congruentCAFsError = congruentCAFs[2, ],
                           incongruentCAFsRT = incongruentCAFs[1, ],
                           incongruentCAFsError = incongruentCAFs[2, ])

  # return it
  return(humanProportions)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
########NEED TO MODIFY TO USE CORRECT PROPORTION FUNCTIONS



# A function to get SINGLE human distributional bin proportions from data,
# given desired CDF and CAF qauntile values. These values will be used by the
# model during fitting.

#' @export
getHumanPropsSingle <- function(conditionData, cdfs, cafs){


  # split trials on congruency
  congruentData <- subset(conditionData,
                          conditionData$congruency == "congruent")
  incongruentData <- subset(conditionData,
                            conditionData$congruency == "incongruent")

  # CDFs-----------------------------------------------------------------------
  # get correct RT CDF proportions. They are a function of the desired CDF
  # quantile points, and the overall accuracy for the congruency condition
  congruentProps <- cdfBinsize(cdfs)
  incongruentProps <- cdfBinsize(cdfs)

  # scale by the proportion correct
  congruentCDFProps <- cdfProportions(congruentData, congruentProps,
                                      multipleSubjects = FALSE)
  incongruentCDFProps <- cdfProportions(incongruentData, incongruentProps,
                                        multipleSubjects = FALSE)

  # get the CDF cut-off values of response time
  congruentCDFs <- cdf(congruentData, quantiles = cdfs,
                       multipleSubjects = FALSE)
  incongruentCDFs <- cdf(incongruentData, quantiles = cdfs,
                         multipleSubjects = FALSE)

  # CAFs-----------------------------------------------------------------------
  # get the CAF proportions
  congruentCAFProps <- cafProportions(congruentData, quantiles = cafs,
                                      multipleSubjects = FALSE)
  incongruentCAFProps <- cafProportions(incongruentData, quantiles = cafs,
                                        multipleSubjects = FALSE)

  # get the CAF cut-off points
  congruentCAFsCutoff <- cdf(congruentData, quantiles = cafs,
                             correctTrials = 3,multipleSubjects = FALSE)
  incongruentCAFsCutoff  <- cdf(incongruentData, quantiles = cafs,
                                correctTrials = 3, multipleSubjects = FALSE)

  # now get the CAFs themselves
  congruentCAFs <- caf(congruentData, quantiles = cafs,
                       multipleSubjects = FALSE)
  incongruentCAFs <- caf(incongruentData, quantiles = cafs,
                         multipleSubjects = FALSE)


  #collate data and return-----------------------------------------------------
  humanProportions <- list(cdfQuantiles = cdfs,
                           cafQuantiles = cafs,
                           congruentCDFProportions = congruentCDFProps,
                           incongruentCDFProportions = incongruentCDFProps,
                           congruentCAFProportions = congruentCAFProps,
                           incongruentCAFProportions = incongruentCAFProps,
                           congruentCDFs = congruentCDFs,
                           incongruentCDFs = incongruentCDFs,
                           congruentCAFsCutoff = congruentCAFsCutoff,
                           incongruentCAFsCutoff =  incongruentCAFsCutoff,
                           congruentCAFsRT = congruentCAFs[1, ],
                           congruentCAFsError = congruentCAFs[2, ],
                           incongruentCAFsRT = incongruentCAFs[1, ],
                           incongruentCAFsError = incongruentCAFs[2, ])

  # return it
  return(humanProportions)

}
#------------------------------------------------------------------------------







#------------------------------------------------------------------------------
# Get model's bin proportions for CDF of correct RT. This uses the human CDFs
# as cut-offs. It then calculates the proportion found in each bin, and
# a later function then compares this to the proportions in the human data.

#'@export
getCDFProps <- function(cdfs, modelData){

  # initialise vector to store proportions in
  props <- numeric(length(cdfs) + 1)

  # only analyse correct RT from the model
  data <- subset(modelData, modelData[, 2] == 1)



  # how many bins are there to work through?
  nBins <- length(props)

  # loop over each, find the bin boundaries, and calculate the proportion of RT
  # in each of the bins
  for(i in 1:nBins){

    # do the first one manually
    if(i == 1){
      # get the data in the current bin
      bin <- subset(data, data[, 1] <= cdfs[i])
      # find the proportion of data in this bin
      props[i] <- length(bin[, 1]) / nrow(modelData)
    }


    # do the middle ones automatically
    if(i > 1 & i <= nBins){
      bin <- subset(data, data[, 1] > cdfs[i - 1] & data[, 1] <= cdfs[i])
      props[i] <- length(bin[, 1]) / nrow(modelData)
    }


    # do the last one manually
    if(i == nBins){
      bin <- subset(data, data[, 1] > cdfs[i - 1])
      props[i] <- length(bin[, 1]) / nrow(modelData)
    }

  }


  # return the proportions
  return(props)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Get model's bin proportions for CAFs. This uses the CAF cutoff CDFs from the
# human data, and finds the proportion of correct responses in each bin
# predicted by the model.

#'@export
getCAFProps <- function(cdfs, modelData){

  # initialise vector to store proportions in
  props <- numeric(length(cdfs) + 1)

  # how many bins are there to calculate?
  nBins <- length(props)

  # loop over each CDF value, find the bin boundaries, source the data, and
  # find the proportion of correct responses in each bin.
  for(i in 1:nBins){


    # do the first bin manually
    if(i == 1){

      # get the current bin's data
      bin <- subset(modelData, modelData[, 1] <= cdfs[i])

      # if there are data points in this bin
      if(length(bin[, 2]) > 0) {

        # find the proportion of correct RTs
        props[i] <- (sum(bin[, 2] == 0)) / nrow(modelData)

      } else {
        props[i] <- 0
      }
    }


    # do middle ones automatically
    if(i > 1 & i <= nBins){

      bin <- subset(modelData, modelData[, 1] > cdfs[i - 1] &
                      modelData[, 1] <= cdfs[i])
      if(length(bin[, 2]) > 0) {
        props[i] <- (sum(bin[, 2] == 0)) / nrow(modelData)
      } else {
        props[i] <- 0
      }

    }


    # do the final bin manually
    if(i == nBins){
      bin <- subset(modelData, modelData[, 1] > cdfs[i - 1])
      if(length(bin[, 2]) > 0) {
        props[i] <- (sum(bin[, 2] == 0)) / nrow(modelData)
      } else {
        props[i] <- 9
      }
    }


  }

  # return the proportions
  return(props)

}
#------------------------------------------------------------------------------
