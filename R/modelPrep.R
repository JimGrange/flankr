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


  # get CDF proportions (they are a function of the requested CDF values)
  congruentProps <- cdf.proportions(cdfs)
  incongruentProps <- cdf.proportions(cdfs)

  # split trials on congruency
  congruentData <- subset(conditionData,
                          conditionData$congruency == "congruent")
  incongruentData <- subset(conditionData,
                            conditionData$congruency == "incongruent")

  # get the CDFs
  congruentCDFs <- cdf(congruentData, quantiles = cdfs)
  incongruentCDFs <- cdf(incongruentData, quantiles = cdfs)

  # get the CAF cut-off points
  congruentCAFsCutoff <- cdf(congruentData, quantiles = cafs,
                             correctTrials = 3)
  incongruentCAFsCutoff  <- cdf(incongruentData, quantiles = cafs,
                                correctTrials = 3)

  # now get the CAFs themselves
  congruentCAFs <- caf(congruentData, quantiles = cafs)
  incongruentCAFs <- caf(incongruentData, quantiles = cafs)


  #collate data and return-----------------------------------------------------

  # collate all the information for human proportions
  humanProps <- c(congruentProps, incongruentProps, congruentCAFs[2, ],
                  incongruentCAFs[2, ])

  # now get the human cut-off values for CAFs and CDFs
  humanCutoffs <- c(congruentCDFs, incongruentCDFs, congruentCAFsCutoff,
                    incongruentCAFsCutoff, congruentCAFs, incongruentCAFs)

  # collate ALL human proportion information
  humanProportions <- c(humanProps, humanCutoffs)


  # return it
  return(humanProportions)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Get model's bin proportions for CDF of correct RT. This uses the human CDFs
# as cut-offs. It then calculates the proportion found in each bin, and
# a later function then compares this to the proportions in the human data.

#'@export
cdfProportions <- function(cdfs, modelData){

  # initialise vector to store proportions in
  props <- numeric(length(cdfs) + 1)


}
