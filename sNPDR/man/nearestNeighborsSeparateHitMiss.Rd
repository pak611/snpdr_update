% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearestNeighbors.R
\name{nearestNeighborsSeparateHitMiss}
\alias{nearestNeighborsSeparateHitMiss}
\title{nearestNeighborsSeparateHitMiss}
\usage{
nearestNeighborsSeparateHitMiss(
  attr.mat,
  pheno.vec,
  nbd.method = "relieff",
  nbd.metric = "manhattan",
  sd.frac = 0.5,
  k = 0,
  neighbor.sampling = "none",
  att_to_remove = c(),
  fast.dist = FALSE,
  dopar.nn = FALSE
)
}
\arguments{
\item{attr.mat}{m x p matrix of m instances and p attributes}

\item{pheno.vec}{vector of class values for m instances}

\item{nbd.method}{neighborhood method \code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)}

\item{nbd.metric}{used in npdrDistances for distance matrix between instances, default: \code{"manhattan"} (numeric)}

\item{sd.frac}{multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2}

\item{k}{number of constant nearest hits/misses for \code{"relieff"} (fixed k).
The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.}

\item{neighbor.sampling}{"none" or \code{"unique"} if you want to return only unique neighbor pairs}

\item{att_to_remove}{attributes for removal (possible confounders) from the distance matrix calculation.}

\item{fast.dist}{whether or not distance is computed by faster algorithm in wordspace, default as F}

\item{dopar.nn}{whether or not neighborhood is computed in parallel, default as F}
}
\value{
Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
}
\description{
Find nearest neighbors of each instance using relief.method.
Treat the hit and miss distributions separately to circument potential hit bias.
ReliefF version makes hit/miss neighborhoods balanced. Surf and MultiSurf are still imbalanced.
Used for npdr (no hits or misses specified in neighbor function).
}
\examples{
# reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
# One can change the theoretical value by changing sd.frac (default 0.5)
neighbor.pairs.idx <- nearestNeighborsSeparateHitMiss(
  predictors.mat, case.control.3sets$train$class, # need attributes and pheno
  nbd.method = "relieff", nbd.metric = "manhattan",
  sd.frac = .5, k = 0
)
head(neighbor.pairs.idx)
}
