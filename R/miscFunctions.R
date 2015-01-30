# ------------------------------------------------------------------------------
# function so the user can get their data from a dialog box rather than code
#' @export
getData <- function(){

  data <- file.choose()
  data <- read.csv(data, header = TRUE)

  return(data)

}
# ------------------------------------------------------------------------------
