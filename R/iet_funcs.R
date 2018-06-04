###################################################################################
### This file contains functions to calculate statistics for the inter-event times
### of bivariate marked point processes.
###################################################################################

#' Return smallest non-zero element of a vector.
#'
#' @param x A vector of positive real numbers.
#' @return Minimum value of \code{x} such that \code{min(x) > 0}
#' @examples
#' nonzero.min(c(1,3,5,2,0))
nonzero.min <- function(x){
  return( min(x[x > 0]) )
}

###################################################################################
#' Calculate inter-event times for a pair of event series.
#'
#' \code{calc.inter_event_times} calculates inter-event times regardless of direction
#' (i.e. all positive even if closest event was earlier).
#'
#' @param x Data.frame of \code{<id, mark, time>} for pair of event series.
#' @return List of lists:
#'   \itemize{
#'     \item \code{AB} = set of iets for fixed mark 1 time to closest event of mark 2
#'     \item \cdoe{BA} = set of iets for fixed mark 2 time to closest event of mark 1
#'   }
calc.inter_event_times <- function(x, sum.stat=FALSE){
  x <- x[with(x, order(t)),]  # sort data by ascending time

  ### inter-event times for same marks
  xA <- x[which(x$m == 1),]

  xB <- x[which(x$m == 2),]

  ### inter-event times for different marks
  AB <- lapply(xA$t, function(e) nonzero.min(abs(xB$t - e)) )
  BA <- lapply(xB$t, function(e) nonzero.min(abs(xA$t - e)) )

  return( list(AB=AB, BA=BA) )
}

###################################################################################
#' Calculate inter-event times between one event and the next event forward in time.
#'
#' \code{calc.forward_inter_event_times} calculates inter-event times between one
#' event and the next event forward in time. (i.e. equiv to inter-event time in
#' Poisson process).
#'
#' @inheritParams calc.inter_event_times
#' @return List of lists:
#'   \itemize{
#'     \item \code{AB} = set of forward iets for fixed mark 1 time to closest event of mark 2
#'     \item \code{BA} = set of forward iets for fixed mark 2 time to closest event of mark 1
#'   }
calc.forward_inter_event_times <- function(x){
  x <- x[with(x, order(t)),]  # sort data by ascending time

  ### inter-event times for same marks
  xA <- x[which(x$m == 1),]
  xB <- x[which(x$m == 2),]

  ### inter-event times for different marks
  AB <- lapply(xA$t, function(e) nonzero.min(xB$t - e) )
  BA <- lapply(xB$t, function(e) nonzero.min(xA$t - e) )

  if(sum(is.infinite(as.numeric(AB))) > 0){ AB <- as.numeric(AB)[-which(is.infinite(as.numeric(AB)))] }
  if(sum(is.infinite(as.numeric(BA))) > 0){ BA <- as.numeric(BA)[-which(is.infinite(as.numeric(BA)))] }

  return( list(AB=AB, BA=BA) )
}

###################################################################################
#' Plot inter event time empirical cdf's for all pairwise combinations of event series.
#'
#' @param x Data.frame of inter-event times \code{<id.A, id.B, iet>}
#' @param main Title of the plot.
#' @param sub Subtitle of the plot.
#' @param samp.prob Probability of including diff source line to reduce overplottting; \code{deafult == 1}
#' @param xmax Max value for x-axis range; \code{default == max(x$iet)}
#' @param legend Logical indicating whether or not to include legend; \code{default == FALSE}
plot.iet <- function(x, main="", sub="", samp.prob=1, xmax=NULL, legend=FALSE){
  if(is.null(xmax)){xmax <- max(x$iet)}
  ids <- unique(x$id.fb)
  ind <- which(x$id.fb == ids[1] & x$id.nonfb == ids[1])
  plot(ecdf(x$iet[ind]), col=alpha(rgb(0,0,0), 0.4), xlim=c(0, xmax),
       main=main, xlab="Hours", xaxt="n", verticals=TRUE, do.points=FALSE)
  mtext(sub)
  axis(1, at=seq(0,xmax*24,1)/24, labels=seq(0,xmax*24,1))
  for(i in ids){
    for(j in ids){
      if(i == j & i == ids[1]){ next }  # already did this one
      if(i == j){
        c <- alpha(rgb(0,0,0), 0.4)
        lwd <- 1
      }else{
        if(runif(1) > samp.prob){ next }  # skip with prob samp.prob
        c <- alpha(rgb(1,0,0), 0.1)
        lwd <- 0.75
      }
      ind <- which(x$id.fb == i & x$id.nonfb == j)
      lines(ecdf(x$iet[ind]), col=c, lwd=lwd, verticals=TRUE, do.points=FALSE)
    }
  }
  if(legend){
    legend("bottomright", lty=1, col=c("black", "red"), lwd=1, legend=c("Same Source", "Diff Source"))
  }
}

###################################################################################
#' Plot densities for a given iet statistic and pair of event streams.
#'
#' @param within,between Data.frame of statistics returned from repeatedly running
#'   calc.inter_event_times and saving results with same/different user event streams.
#' @param which.stat Mean (\code{"mn"}) or median (\code{"md"}); \code{default == "mn"}
#' @param which.iet Which distribution to use (\code{"AB", "BA"}); \code{default == "BA"}
#' @param bw.type Bandwidth selection method for density() function; \code{default == "nrd"}
#' @param font.sz Size of font for axis; \code{default == 16}
plot.iet_dist <- function(within, between, which.stat="mn", which.iet="BA", bw.type="nrd",
                          font.sz=12){
  l <- 0
  u <- max( between[,paste(which.stat,".",which.iet,sep="")] )
  rng <- seq(l,u,0.0001)
  p.win <- approxfun(density(within[,paste(which.stat,".",which.iet,sep="")],
                             from=l, to=u, bw=bw.type))
  p.bt <- approxfun(density(between[,paste(which.stat,".",which.iet,sep="")],
                            from=l, to=u, bw=bw.type))
  plt.dat <- data.frame(x = rng*24*60, win = p.win(rng), bt = p.bt(rng))
  plt <- ggplot(plt.dat, aes(x=x, y=value)) +
               geom_line(aes(y=bt), linetype="solid", col="black") +
               geom_line(aes(y=win), linetype="dashed", col="black") +
               theme_bw() +
               ylab("Density") +
               xlab("Minutes") +
               ggtitle(which.iet) +
               theme(legend.position="none", axis.text=element_text(size=font.sz),
                     axis.title=element_text(size=font.sz))
  return(plt)
}

###################################################################################
#' Calculate score-based likelihood ratios.
#'
#' @param ids Vector of unique ids of the data considered.
#' @param within Data.frame of inter-event times for same-source pairs of event series,
#'   \code{<id, mn, md>}
#' @param between Data.frame of inter-event times for different-source pairs of event
#'   series, \code{<id.A, id.B, mn, md>}
#' @param bds.mn,bds.md Numeric array \code{c(lower, upper)} of limmits for KDE on
#'   mean/median IET.
#' @param match Logical indicating to compute scores of known matches or non-matches;
#'   \code{default == TRUE}
#' @param bw.type Bandwidth estimation method to use in density estimation
#'   one of \code{c("nrd", "nrd0", "ucv", "bcv", "SJ")}; \code{default == "nrd"}
#' @return Data.frame of score-based likelihood ratios:
#'   \itemize{
#'     \item \code{match == TRUE}: \code{<id, mn, md>}
#'     \item \code{match == FALSE}: \code{<id.A, id.B, mn, md>}
#'   }
calc.iet_SLR <- function(ids, within, between, bds.mn, bds.md, match = TRUE, bw.type = "nrd"){
  if(match == TRUE){  # compute SLR for known matches
    slr_mn <- slr_md <- rep(NA, length(ids))
    for(i in 1:length(ids)){
      ### "within" density not including student i
      ind.win <- which(within$id != ids[i])
      p.mn.win <- approxfun(density(within$mn, from=bds.mn[1], to=bds.mn[2], bw=bw.type))
      p.md.win <- approxfun(density(within$md, from=bds.md[1], to=bds.md[2], bw=bw.type))
      ### "between" density not including any stream from student i
      ind.bt <- which(between$id.fb != ids[i] & between$id.nonfb != ids[i])
      p.mn.bt <- approxfun(density(between$mn, from=bds.mn[1], to=bds.mn[2], bw=bw.type))
      p.md.bt <- approxfun(density(between$md, from=bds.md[1], to=bds.md[2], bw=bw.type))
      ### compute score-based likelihood ratio for each
      slr_mn[i] <- p.mn.win(within$mn[i]) / p.mn.bt(within$mn[i])
      slr_md[i] <- p.md.win(within$md[i]) / p.md.bt(within$md[i])
    }
    return( data.frame(id = ids, mn = slr_mn, md = slr_md) )
  } else{  # compute SLR for known non-matches
    slr_mn <- slr_md <- rep(NA, nrow(between))
    for(i in 1:nrow(between)){
      rm.ids <- between[i,1:2]  # get ids of current row
      ### calc densities for within
      ind.win <- which( !(within$id %in% rm.ids) )  # ids for within that arent current ids
      p.mn.win <- approxfun(density(within$mn[ind.win], from=bds.mn[1], to=bds.mn[2], bw=bw.type))
      p.md.win <- approxfun(density(within$md[ind.win], from=bds.md[1], to=bds.md[2], bw=bw.type))
      ### calc densities for between
      ind.bt <- which( !(between$id.fb %in% rm.ids) & !(between$id.nonfb %in% rm.ids) )  # same for between
      p.mn.bt <- approxfun(density(between$mn[ind.bt], from=bds.mn[1], to=bds.mn[2], bw=bw.type))
      p.md.bt <- approxfun(density(between$md[ind.bt], from=bds.md[1], to=bds.md[2], bw=bw.type))
      ### calc SLR
      slr_mn[i] <- p.mn.win(between$mn[i]) / p.mn.bt(between$mn[i])
      slr_md[i] <- p.md.win(between$md[i]) / p.md.bt(between$md[i])
    }
    return( data.frame(between[,1:2], mn = slr_mn, md = slr_md) )
  }
}
