% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markedPP_funcs.R
\name{sim_markedpp}
\alias{sim_markedpp}
\title{Simulate a pair of event time series.}
\usage{
sim_markedpp(lambda, W, nSamp, indep = FALSE, p, normal = TRUE, sigma,
  lambda2)
}
\arguments{
\item{lambda}{Intensity of first process (i.e. of Mark 1/A).}

\item{W}{Observation window \code{c(low, high)}.}

\item{nSamp}{Number of samples to draw, intentionally too large.}

\item{indep}{Logical indicating if the second process should depend on the
first; \code{default == FALSE}
\itemize{
  \item \code{TRUE}: homogeneous Poisson process of rate lambda2
  \item \code{FALSE}: draw \eqn{Bern(p)}, if 1 draw point from
        \eqn{N(t_1[i],sigma^2)}
}}

\item{p}{Bernoulli probability if \code{indep == FALSE}}

\item{normal}{Logical indicating distribution used to simulate dependent events;
\code{default == TRUE}
\itemize{
  \item \code{TRUE}: Gaussian(t_1[i], sigma^2)
  \item \code{FALSE}: Exponential(sigma)
}}

\item{sigma}{Standard deviation of normal if \code{indep == FALSE}}

\item{lambda2}{intensity of second process if  \code{indep == TRUE}}
}
\value{
Data.frame of \code{<t, m>}
}
\description{
Simulate a pair of event time series.
}
