% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dstpFunctions.R
\name{simulateDSTP}
\alias{simulateDSTP}
\title{Obtain simulated response times and accuracy from the DSTP model}
\usage{
simulateDSTP(parms, nTrials, var = 0.01, dt = 1/1000, seed = NULL)
}
\arguments{
\item{parms}{The set of parameters to use to simulate the data. Must be
contained in a vector in the order: \code{A}, \code{C}, \code{driftTarget},
\code{driftFlanker}, \code{diftStimSelection}, \code{driftRS2}, \code{ter}.}

\item{nTrials}{How many trials to simulate per congruency condition.}

\item{var}{The variance of the diffusion process. By default this is set to
0.01.}

\item{dt}{The diffusion scaling parameter (i.e., time steps). By default,
this is set to 0.001.}

\item{seed}{The value for the \code{set.seed} function to set random
generation state.}
}
\value{
Returns a data frame with three columns: rt (response time) in
seconds, accuracy of the model's response (1 for correct, 0 for error), and
congruency condition.
}
\description{
\code{simulateDSTP} generates synthetic data from the DSTP model in the
form of response time (RT) in seconds and accuracy for both congruent and
incongruent trials.
}
\details{
This function can be employed by the user to generate synthetic data, but
its main purpose is to be used by the fitting procedure to generate model
predictions for a set of parameter values when trying to find the best-
fitting values.d
}
\examples{
# declare the parameters
parms <- c(0.070, 0.086, 0.045, 0.065, 0.368, 1.575, 0.225)

# simulate the data
modelData <- simulateDSTP(parms, nTrials = 10000)

}
