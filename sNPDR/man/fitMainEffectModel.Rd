% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regain.R
\name{fitMainEffectModel}
\alias{fitMainEffectModel}
\title{\code{fitMainEffectModel}
Get the main effect of a variable using generalized linear regression .}
\usage{
fitMainEffectModel(
  labelledDataFrame,
  variableName,
  depVarName,
  regressionFamily,
  numCovariates
)
}
\arguments{
\item{labelledDataFrame}{\code{data.frame} with variables in columns and samples in rows.}

\item{variableName}{\code{string} name of the variable to consider.}

\item{depVarName}{\code{string} name of the classLabel variable.}

\item{regressionFamily}{\code{string} glm regression family name.}

\item{numCovariates}{\code{numeric}  of included covariates.}
}
\value{
\code{data.frame} with variable, convergence status, beta coefficient,
p-value, standard error and standardized beta columns.
}
\description{
\code{fitMainEffectModel}
Get the main effect of a variable using generalized linear regression .
}
