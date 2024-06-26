% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npdrLearner.R
\name{npdrLearner}
\alias{npdrLearner}
\title{npdrLearner}
\usage{
npdrLearner(
  train.outcome,
  train.data,
  test.outcome,
  test.data,
  nbd.method = "relieff",
  nbd.metric = "manhattan",
  msurf.sd.frac = 0.5,
  dopar.nn = FALSE,
  knn = 0
)
}
\arguments{
\item{train.outcome}{character name or length-m numeric outcome vector for
train data}

\item{train.data}{m x p matrix of m instances and p attributes of train data.
May also include outcome vector but then outcome should be a name.
Include attr names as colnames.}

\item{test.outcome}{character name or length-m numeric outcome vector of test data}

\item{test.data}{m x p matrix of m instances and p attributes for test data.
May also include outcome vector but then outcome should be a name.
Include attr names as colnames.}

\item{nbd.method}{neighborhood method: `multisurf` or `surf` (no k) or
`relieff` (specify k). Used by nearestNeighbors2().}

\item{nbd.metric}{used in npdrDistances2 for distance matrix between
instances, default: `manhattan` (numeric). Used by nearestNeighbors2().}

\item{msurf.sd.frac}{multiplier of the standard deviation from the mean
distances; subtracted from mean for SURF or multiSURF.
The multiSURF default is `msurf.sd.frac=0.5`: mean - sd/2.
Used by nearestNeighbors2().}

\item{dopar.nn}{whether or not neighborhood is computed in parallel,
default as FALSE.}

\item{knn}{number of constant nearest hits/misses for `relieff` (fixed-k).
Used by nearestNeighbors2().
The default knn=0 means use the expected SURF theoretical `k` with
`msurf.sd.frac` (0.5 by default)}
}
\value{
list: neighborhoods for each test instance, prediction for each
test instance, accuracy on test set
}
\description{
Uses npdr neighborhoods to learn a nearest neighbor classification or
regression model (latter not implemented but easy).
Finds the nearest neighbors of test instances to a training dataset.
Uses majority class of training neighbors for test prediction.
Allows adaptive Relief neighborhoods or specify k.
Regression would simply use the average value of the neighbor phenotypes.
Uses functions npdrDistances2 and nearestNeighbors2.
}
\examples{
test.results <- npdrLearner(
  train.outcome = "class", train.data = case.control.3sets$train,
  test.outcome = "class", test.data = case.control.3sets$validation,
  nbd.method = "relieff", nbd.metric = "manhattan",
  dopar.nn = FALSE, knn = 0
)
test.results$accuracy
}
