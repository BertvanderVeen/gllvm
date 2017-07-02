% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.gllvm.R
\name{confint.gllvm}
\alias{confint.gllvm}
\title{Confidence intervals for model parameters}
\usage{
\method{confint}{gllvm}(object, level = 0.95, ...)
}
\arguments{
\item{object}{An object of class 'gllvm'.}

\item{level}{The confidence level. Scalar between 0 and 1.}

\item{...}{Not used.}
}
\description{
Computes confidence intervals for  parameters in a fitted GLLVM model.
}
\examples{
\dontrun{
library(mvabund) ## Load a dataset from the mvabund package
data(antTraits)
y <- as.matrix(antTraits$abund)
X <- as.matrix(antTraits$env[,4:5])
# Fit GLLVM model
fit <- gllvm(y = y, X = X, family = "negative.binomial")
# 95 \% confidence intervals
confint(fit, level = 0.95)
}
}
\author{
Jenni Niku <jenni.m.e.niku@jyu.fi>
}