# ------------------------------------------------------------------------------
# function so the user can get their data from a dialog box rather than code
getData <- function(){

  data <- file.choose()
  data <- read.csv(data, header = TRUE)

  return(data)

}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# get matrix of random starting parameters for DSTP model around fixed means
getRandomParms <- function(startParms, varParms, maxParms, n){

  # initialise matrix
  finalParms <- matrix(0, n, length(startParms))


  # for each row
  for(i in 1:n){
    currParms <-  startParms + rnorm(length(startParms), 0, varParms)

    # make sure no parameter is negative and not above maximum allowed value
    while(((min(currParms) < 0) | (min(maxParms - currParms) < 0))){
      currParms <-  startParms + rnorm(length(startParms), 0, varParms)
    }

    # store the current vector
    finalParms[i, ] <- currParms
  }

  # change the first entry to match the starting parameters
  finalParms[1, ] <- startParms

  finalParms <- round(finalParms, 3)

  # return the matrix
  return(finalParms)

}
