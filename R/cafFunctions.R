###
# functions for finding response time conditional accuracy functions (CAFs)


#------------------------------------------------------------------------------
#' Find conditional accuracy function (CAF) values for a single condition
#'
#' \code{caf} takes a data frame for a single experimental condition and
#' returns a vector of requested conditional accuracty function (CAF) values.
#'
#' The function only deals with one experimental condition. There is another
#' function (\code{cafAll}) which will return CAFs for all experimental
#' conditions. If there are more than one subject in the data frame being
#' passed to this function, the function first finds the CAF values for each
#' subject, and then takes the average for each quantile. This average is then
#' returned to the user.
#'
#' @param data A data frame containing the data to be passed to the function.
#' At the very least, the data frame must contain columns named "accuracy"
#' logging the accuracy (1 for correct, 0 for error) and "rt" containing the
#' response time data. If the user wishes to find the average CAFs across
#' multiple subjects, then another column must be included ("subject") with
#' numbers identifying unique subjects. See \code{?exampleData} for a data
#' frame formatted correctly.
#'
#' @param quantiles The quantile values to be found by the function. By
#' default, the function finds the accuracy for the .25, .5, and .75 quantiles.
#'
#' @param multipleSubjects Inform the function whether the data frame contains
#' data from multiple subjects. If set to TRUE, the function returns the
#' average CAF values across all subjects. If set to FALSE, the function
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
#' getCAF <- caf(data)
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
#' newQuantiles <- c(.2, .4, .6, .8)
#'
#' # get the CAFs
#' getCAF <- caf(data, quantiles = newQuantiles, multipleSubjects = FALSE)
#'

#' @export
caf <- function(data, quantiles = c(.25, .50, .75), multipleSubjects = TRUE){

  # single participant---------------------------------------------------------
  # perform simple CAF calculation if only one participant is being passed
  # to the function
  if(multipleSubjects == FALSE){

    # initialise empty vector to store data
    cafData <- numeric(length = (length(quantiles) * 2) + 2)

    # get the RT values for quantile cut-offs
    cdfs <- quantile(data$rt, quantiles)

    # calculate mean RT and proportion error for each bin----------------------
    for(i in 1:length(quantiles)){

      ## do the first one manually
      if(i == 1){
        # get the data
        temp <- subset(data, data$rt < cdfs[i])
        # log the RT
        cafData[i] <- mean(temp$rt)
        # log the accuracy
        cafData[i + length(quantiles) + 1] <- sum(temp$accuracy) /
          length(temp$accuracy)
      }

      ## do the rest of the slots automatically
      if(i > 1 & i <= length(quantiles)){

        # get the data
        temp <- subset(data, data$rt > cdfs[i - 1] & data$rt < cdfs[i])
        # log the RT
        cafData[i] <- mean(temp$rt)
        # log the accuracy
        cafData[i + length(quantiles) + 1] <- sum(temp$accuracy) /
          length(temp$accuracy)

      }

      ## do the last one manually, too
      if(i == length(quantiles)){
        # get the data
        temp <- subset(data, data$rt > cdfs[i])
        # log the RT
        cafData[i + 1] <- mean(temp$rt)
        # log the accuracy
        cafData[(i + 1) + length(quantiles) + 1] <- sum(temp$accuracy) /
          length(temp$accuracy)
      }

    } #end of loop over quantiles

    return(cafData)

  } #end of single-subject sub-function


  #multiple subjects-----------------------------------------------------------
  if(multipleSubjects == TRUE){

    # what are the unique subject numbers?
    subs <- unique(data$subject)

    # how many subjects are there?
    nSubs <- length(subs)

    #empty matrix to store CAF data in
    cafData <- matrix(0, nrow  = ((length(quantiles) * 2) + 2), ncol = nSubs)


    # loop over all subjects, get their CAFs, and store in cafData matrix
    for(j in 1:nSubs){

      # get the current subject's data
      subData <- subset(data, data$subject == subs[j])

      # get the current subject's CDF criteria
      subCDFs <- quantile(subData$rt, quantiles)


        # calculate mean RT and proportion error for each bin----------------------
        for(i in 1:length(quantiles)){

          ## do the first one manually
          if(i == 1){
            # get the data
            temp <- subset(subData, subData$rt < subCDFs[i])
            # log the RT
            cafData[i, j] <- mean(temp$rt)
            # log the accuracy
            cafData[i + length(quantiles) + 1, j] <- sum(temp$accuracy) /
              length(temp$accuracy)
          }

          ## do the rest of the slots automatically
          if(i > 1 & i <= length(quantiles)){

            # get the data
            temp <- subset(subData, subData$rt > subCDFs[i - 1] & subData$rt <
                             subCDFs[i])
            # log the RT
            cafData[i, j] <- mean(temp$rt)
            # log the accuracy
            cafData[i + length(quantiles) + 1, j] <- sum(temp$accuracy) /
              length(temp$accuracy)

          }

          ## do the last one manually, too
          if(i == length(quantiles)){
            # get the data
            temp <- subset(subData, subData$rt > subCDFs[i])
            # log the RT
            cafData[i + 1, j] <- mean(temp$rt)
            # log the accuracy
            cafData[(i + 1) + length(quantiles) + 1, j] <- sum(temp$accuracy) /
              length(temp$accuracy)
          }

        } # end of loop over quantiles


    } # end of loop over subjects

    # find the mean values
    cafData <- apply(cafData, 1, mean)

    # return the means
    return(cafData)

  } # end of multiple subjects loop


} # end of function
#------------------------------------------------------------------------------
