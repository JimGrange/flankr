% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cafFunctions.R
\name{caf}
\alias{caf}
\title{Find conditional accuracy function (CAF) values for a single condition}
\usage{
caf(data, quantiles = c(0.25, 0.5, 0.75), multipleSubjects = TRUE)
}
\arguments{
\item{data}{A data frame containing the data to be passed to the function.
At the very least, the data frame must contain columns named "accuracy"
logging the accuracy (1 for correct, 0 for error) and "rt" containing the
response time data. If the user wishes to find the average CAFs across
multiple subjects, then another column must be included ("subject") with
numbers identifying unique subjects. See \code{?exampleData} for a data
frame formatted correctly.}

\item{quantiles}{The quantile values to be found by the function. By
default, the function finds the accuracy for the .25, .5, and .75 quantiles.}

\item{multipleSubjects}{Inform the function whether the data frame contains
data from multiple subjects. If set to TRUE, the function returns the
average CAF values across all subjects. If set to FALSE, the function
assumes all data being passed is just from one subject.}
}
\description{
\code{caf} takes a data frame for a single experimental condition and
returns a vector of requested conditional accuracty function (CAF) values.
}
\details{
The function only deals with one experimental condition. There is another
function (\code{cafAll}) which will return CAFs for all experimental
conditions. If there are more than one subject in the data frame being
passed to this function, the function first finds the CAF values for each
subject, and then takes the average for each quantile. This average is then
returned to the user.
}
\examples{
### example of multiple subjects and default quantile values

# only select the congruent data from the example data set
data <- subset(exampleData, exampleData$congruency == "congruent")

# get the CDFs
getCAF <- caf(data)

### example of single subject and different quantile values

# only select subject 1 from the example data. Also, select only the
# "absent" condition and incongruent trials. This is an example when working
# with multiple conditions (besides target congruency).
data <- subset(exampleData, exampleData$subject == 1 &
    exampleData$condition == "absent" &
    exampleData$congruency == "incongruent")

# set new quantile values
newQuantiles <- c(.2, .4, .6, .8)

# get the CAFs
getCAF <- caf(data, quantiles = newQuantiles, multipleSubjects = FALSE)

}
