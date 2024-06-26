% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeData.R
\name{generate_survival_data}
\alias{generate_survival_data}
\title{Generate Survival Data}
\usage{
generate_survival_data(
  n = 1000,
  p = 5,
  max_observation_time = 100,
  intEff = FALSE
)
}
\arguments{
\item{n}{Number of observations. Default is 1000.}

\item{p}{Number of features. Default is 5.}

\item{max_observation_time}{Maximum observation time for censoring. Default is 100.}

\item{intEff}{Include interactions}
}
\value{
A dataframe containing the generated survival data, with columns for the features, survival times, censoring status, and the survival object.
}
\description{
This function generates synthetic survival data based on a combination of main effects,
interaction effects, and random noise. It also visualizes the generated Kaplan-Meier survival curves.
}
\examples{
\dontrun{
data <- generate_survival_data(n=1000, p=5, max_observation_time=100)
head(data)
}
}
