#------------------------------------------------------------------------------
#' Example response time data set for multiple subjects.
#'
#' An example data set containing multiple participants' data for a response
#' time study involving two experimental conditions, as well as the congruency
#' manipulation. The data includes response time and accuracy.
#'
#' This data is taken from Grange, J.A. (in preparation). The effect of
#' accessory stimuli on choice response time in the flanker task.
#'
#' @format A data frame with 12524 rows and 5 variables:
#' \describe{
#'     \item{subject}{The subject identification number. Note that some
#'           participants were removed due to poor accuracy.}
#'     \item{condition}{The experimental condition (2 in this example). It
#'           logs the presence/absence of an auditory tone before stimulus
#'           onset.}
#'     \item{congruency}{The congruency of the flanker stimulus.}
#'     \item{accuracy}{Accuracy of the response; 1 = correct, 0 = error}
#'      \item{rt}{Response time, coded in seconds.}
#'
#'     }
"exampleData"
#------------------------------------------------------------------------------
