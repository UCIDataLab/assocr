% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markedPP_funcs.R
\name{format_data}
\alias{format_data}
\title{Filter and format data for use with marked point process functions.}
\usage{
format_data(data, id.sum, dev = FALSE, low.bds)
}
\arguments{
\item{data}{Data.frame of event data,
\code{<id, dt, t, type, category, app_name, lbl_moshe, is_fb, t.center>}}

\item{id.sum}{Data.frame of id and counts by event type,
\code{<id, ev.dur, dur, tot, phone, comp, comp.fb, comp.nonfb, start, end>}}

\item{dev}{Logical indicating if marks correspond to device type
(1=phone, 2=comp); \code{default == FALSE}}

\item{low.bds}{Numeric array \code{c(n1, n2)} of criteria for minimum number
events of mark types 1 and 2, resptively.}
}
\value{
Data.frame of \code{<id, m, t>}
}
\description{
Filter and format data for use with marked point process functions.
}
