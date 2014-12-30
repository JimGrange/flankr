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
#' @param correctTrials If set to TRUE, the function will find the CDFs of
#' correct trials. Set to FALSE to find the CDFs of error trials. Note, though,
#' that CDFs of error trials may be less accurate due to usually-low number of
#' error trials.
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
                correctTrials = TRUE, multipleSubjects = TRUE){

  # perform the simple operation of calculating CDFs if only one subject
  if(multipleSubjects == FALSE){


    # select whether the user wants correct trials or error trials
    if(correctTrials == TRUE){
      tempData <- subset(data, data$accuracy == 1)
    } else {
      tempData <- subset(data, data$accuracy == 0)
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

        # select whether the user wants correct trials or error trials
        if(correctTrials == TRUE){
          tempData <- subset(tempData, data$accuracy == 1)
        } else {
          tempData <- subset(tempData, data$accuracy == 0)
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

