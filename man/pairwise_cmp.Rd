% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_funcs.R
\name{pairwise_cmp}
\alias{pairwise_cmp}
\title{Wrapper to compute CMP for each combination of a given pair of user data}
\usage{
pairwise_cmp(data, n, W = c(0, 7), bidirectional = TRUE,
  samp = "empirical", sampSpace = NULL, rng = NULL, mark = 1)
}
\arguments{
\item{data}{List of data.frames for two pairs of event series (output of
\code{assocr::sessionize_data()}):
\itemize{
  \item \code{data}: sessionized data, \code{<id, m, sid, t>}
  \item \code{sessions}: summary of sessionized data,
        \code{<id, m, sid, n, t>}
}}

\item{n}{Number of samples to draw}

\item{W}{Vector containing observation window \code{c(low, high)}.}

\item{bidirectional}{Logical indicating if the inter-event times should be
bidirectional (versus only looking at events prior to event of interest)}

\item{samp}{Sampling technique to use. One of
\itemize{
  \item \code{empirical}: sessionized resampling from empirical distn of
    start times
  \item \code{gaussian}: sessionized resampling from Gaussian distn such
    that 99% of start times fall in \code{rng}
  \item \code{uniform}: sessionized resampling from Uniform distn over
    \code{rng}
}}

\item{sampSpace}{vector representation of empirical sample space for session
start times if \code{samp == "empirical"}}

\item{rng}{Data.frame of ranges for resampling; column names correspond to
ids and the 2 rows correspond to (lower, upper) limits}

\item{mark}{which mark to perform the seampling for;  \code{default == 1}}
}
\value{
Array of CMPs for each pairwise combination of event streams &
  score function
}
\description{
Wrapper to compute CMP for each combination of a given pair of user data
}
