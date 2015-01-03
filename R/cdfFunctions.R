###
# functions for finding response time cumulative distribution functions (CDFs)


#------------------------------------------------------------------------------
#' Find cumulative distribution function (CDF) values for a single condition
#'
#' \code{cdf} takes a data frame for a single experimental condition and
#' returns a vector of requested CDF values.
#'
#' The function only deals with one experimental condition. There is another
#' function (\code{cdfAll}) which will return CDFs for all experimental
#' conditions. If there are more than one subject in the data frame being
#' passed to this function, the function first finds the CDF values for each
#' subject, and then takes the average for each quantile. This average is then
#' returned to the user.
#'
#' @param data A data frame containing the data to be passed to the function.
#' At the very least, the data frame must contain columns named "accuracy"
#' logging the accuracy (1 for correct, 0 for error) and "rt" containing the
#' response time data. If the user wishes to find the average CDFs across
#' multiple subjects, then another column must be included ("subject") with
#' numbers identifying unique subjects. See \code{?exampleData} for a data
#' frame formatted correctly.
#'
#' @param quantiles The quantile values to be found by the function. By
#' default, the function finds the .1, .3, .5, .7, and .9 CDF values.
#'
#' @param correctTrials If set to 1, the function will find the CDFs of
#' correct trials. Set to 2 to find the CDFs of error trials. Set to 3 to find
#' CDFs of ALL trials. Note, though, that CDFs of error trials may be less
#' accurate due to usually-low number of error trials.
#'
#' @param multipleSubjects Inform the function whether the data frame contains
#' data from multiple subjects. If set to TRUE, the function returns the
#' average CDF values across all subjects. If set to FALSE, the function
#' assumes all data being passed is just from one subject.
#'
#'
#' @examples
#' ### example of multiple subjects and default quantile values
#'
#' # only select the congruent data from the example data set
#' data <- subset(exampleData, exampleData$congruency == "congruent")
#'
#' # get the CDFs
#' getCDF <- cdf(data)
#'
#' ### example of single subject and different quantile values
#'
#' # only select subject 1 from the example data. Also, select only the
#' # "absent" condition and incongruent trials. This is an example when working
#' # with multiple conditions (besides target congruency).
#' data <- subset(exampleData, exampleData$subject == 1 &
#'     exampleData$condition == "absent" &
#'     exampleData$congruency == "incongruent")
#'
#' # set new quantile values
#' newQuantiles <- c(.1, .2, .3, .4, .5, .6, .7, .8,  .9)
#'
#' # get the CDFs
#' getCDF <- cdf(data, quantiles = newQuantiles, multipleSubjects = FALSE)
#'

#' @export
cdf <- function(data, quantiles = c(.1, .3, .5, .7, .9),
                correctTrials = 1, multipleSubjects = TRUE){

  # perform the simple operation of calculating CDFs if only one subject
  if(multipleSubjects == FALSE){


    # select whether the user wants correct trials or error trials (or all!)
    if(correctTrials == 1){
      tempData <- subset(data, data$accuracy == 1)
    }
    if(correctTrials == 2){
      tempData <- subset(data, data$accuracy == 0)
    }
    if(correctTrials == 3){
      tempData <- data
    }

    # calculate the CDFs
    cdfs <- as.numeric(quantile(tempData$rt, quantiles))

    # return them to the user
    return(cdfs)
  }


  # if there are multiple subjects, find average CDF across these subjects
  if(multipleSubjects == TRUE){

    # find the unique subject numbers
    subs <- unique(data$subject)

    # how many subjects are there?
    nSubs <- length(subs)

    # create a n*m matrix where rows (n) reflect quantile, and columns (m) are
    # subjects. At the end, return the average of each row (quantile)
    cdfData <- matrix(0, nrow = length(quantiles), ncol = nSubs)

    # loop over all subjects, find their CDFs, and place in cdfData matrix
    for(i in 1:nSubs){

      tempData <- subset(data, data$subject == subs[i])

        # select whether the user wants correct trials or error trials (or all!)
        if(correctTrials == 1){
          tempData <- subset(data, data$accuracy == 1)
        }
        if(correctTrials == 2){
          tempData <- subset(data, data$accuracy == 0)
        }
        if(correctTrials == 3){
          tempData <- data
        }


      # log the result
      cdfData[, i] <- quantile(tempData$rt, quantiles)

    }

    #calculate average CDFs across subjects
    averageCDF <- apply(cdfData, 1, mean)

  }

  # return them to the user
  return(averageCDF)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# given a set of quantiles for CDFs, return the proportion of data within each
# bin. For example, the CDFs c(.1, .3, .5, .7, .9) have proportions of
# c(.1, .2, .2, .2, .2, .1). This is required because the model will try to
# predict response times which match the proportions in the human data.

#' @export
cdfProportions <- function(cdfs){

  # get empty vector of the right length
  props <- numeric(length = (length(cdfs) + 1))

  # loop over all cdf values
  for(i in 1:length(cdfs)){

    # do the first one manually
    if(i == 1){
      props[i] <- cdfs[i] - 0
    }


    # do the intermediate bins automatically
    if(i > 1 & i <= length(cdfs)){
      props[i] <- cdfs[i] - cdfs[i - 1]
    }

    # do the final one manually
    if(i == length(cdfs)){
      props[i + 1] <- 1 - cdfs[i]
    }

  } # end of bin loop

  # return the proportions
  return(props)

} # end of function
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# The opposite of cdfProportions. Given a set of proportions, work out the CDFs
# For example, the proportions c(.1, .2, .2, .2, .2, .1) have CDFs of
# c(.1, .3, .5, .7, .9).
#'@export
proportionCDFs <- function(proportions){

  # initialise empty vector for cdfs
  cdfs <- numeric(length(proportions) - 1)

  for(i in 1:length(cdfs)){
    cdfs[i] <- sum(proportions[1:i])
  }
  return(cdfs)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
####THIS IS NOT CURRENTLY USED

#'@export
#get model proportions from human CDFs (the RTs, not proportions). Returns proportions

getModelCDFs <- function(modelData, cdfs){


  # only select the correct trials
  modelData <- subset(modelData, modelData$Accuracy == 1)

  # initiate empty vector to store model CDFs in
  props <- numeric(length(cdfs))

  # loop over each human CDF cutoff point, and find the proportion of model
  # data in each bin
  for(i in 1:length(cdfs)){
    x <- subset(modelData, modelData$RT <= cdfs[i])
    props[i] <- length(x$RT) / nrow(modelData)
  }

  return(props)

}
#------------------------------------------------------------------------------





