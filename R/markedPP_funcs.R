###################################################################################
### This file contains functions to calculate statistics for bivariate marked point
### processes.
###################################################################################
#' Filter and format data for use with marked point process functions.
#'
#' @param data Data.frame of event data,
#'   \code{<id, dt, t, type, category, app_name, lbl_moshe, is_fb, t.center>}
#' @param dev Logical indicating if marks correspond to device type (1=phone, 2=comp);
#'   \code{default == FALSE}
#' @param id.sum Data.frame of id and counts by event type,
#'   \code{<id, ev.dur, dur, tot, phone, comp, comp.fb, comp.nonfb, start, end>}
#' @param low.bds Numeric array \code{c(n1, n2)} of criteria for minimum number events
#'   of mark types 1 and 2, resptively.
#' @return Data.frame of \code{<id, m, t>}
format_data <- function(data, id.sum, dev=FALSE, low.bds){
  if(dev == FALSE){  # marks corr to FB and non-FB events on computer
    attach(id.sum)
    incl <- unique( id.sum$id[ which( dur == 7
                                      & comp.fb >= low.bds[1]
                                      & comp.nonfb >= low.bds[2] ) ] )
    detach(id.sum)
    x <- data[ which(data$id %in% incl & data$type == "comp_web" ), c(1,8,9)]
    names(x) <- c("id","m","t")
    x$m <- as.character(x$m)
    ind <- which( x$m == 1 )
  } else{  # marks corr to device type
    incl <- unique( id.sum$id[ which( id.sum$dur == 7
                                      & id.sum$phone >= low.bds[1]
                                      & id.sum$comp >= low.bds[2] ) ] )
    x <- data[ which(data$id %in% incl), c(1,4,9) ]  # keep all events
    names(x) <- c("id","m","t")
    x$m <- as.character(x$m)
    ind <- which( x$m == "phone")
  }
  x$m[ind] <- 1
  x$m[-ind] <- 2
  x$m <- as.numeric(x$m)
  return(x)
}

###################################################################################
#' Group pairs of event series data into sessions.
#'
#' @param data Data.frame of event data, \code{<id, t, m>}
#' @param thres Number of inactive minutes to consider a new session start
#' @param timeScale Scale of time measurements of data; one of
#'   \code{c("s", "min", "hr", "day")}
#' @return List of data.frames:
#'   \itemize{
#'     \item \code{data}: sessionized data, \code{<id, m, sid, t>}
#'     \item \code{sessions}: summary of sessionized data, \code{<id, m, sid, n, t>}
#'   }
sessionize_data <- function(data, thres = 10, timeScale = "min"){
  # convert event times to minutes
  if (timeScale == "s"){
    data$t.min <- data$t / 60
  } else if (timeScale == "min") {
    data$t.min <- data$t
  } else if (timeScale == "hr"){
    data$t.min <- data$t * 60
  } else if (timeScale == "day"){
    data$t.min <- data$t * 24 * 60
  } else {
    stop("Invalid value for timeScale; enter one of (\"s\", \"min\", \"hr\", \"day\")")
  }

  ids <- unique(data$id)  # get unique ids from the data
  for (i in ids) {  # sessionize data for all users one at a time
    tmp.1 <- .find_sessions(data[data$id == i & data$m == 1, ])  # user i mark 1
    tmp.1$sessions$m <- rep(1, nrow(tmp.1$sessions))
    tmp.2 <- .find_sessions(data[data$id == i & data$m == 2, ])  # user i mark 2
    tmp.2$sessions$m <- rep(2, nrow(tmp.2$sessions))
    if (i == ids[1]) {  # if its the first user, make data.frames to export
      dat <- rbind(tmp.1$data, tmp.2$data)
      ses <- rbind(tmp.1$sessions, tmp.2$sessions)
    } else {  # data.frames already made, concatenate current user's data
      dat <- rbind(dat, rbind(tmp.1$data, tmp.2$data))
      ses <- rbind(ses, rbind(tmp.1$sessions, tmp.2$sessions))
    }
  }
  dat <- dat[order(dat$t), ]
  row.names(dat) <- 1:nrow(dat)
  ses <- ses[order(ses$t), ]
  row.names(ses) <- 1:nrow(ses)
  return ( list(data = dat[, c("id", "m", "sid", "t")],
                sessions = ses[, c("id", "m", "sid", "n", "t")]) )
}

###################################################################################
#' Group single event series (i.e., one user & one type of event) into sessions.
#'
#' @param data Data.frame of event data of same mark for one user, \code{<id, t, t.min>}
#' @param thres Number of inactive minutes to consider a new session start
#' @param timeScale Scale of time measurements of data; one of
#'   \code{c("s", "min", "hr", "day")}
#' @return List of data.frames:
#'   \itemize{
#'     \item \code{data}: sessionized data, \code{<t, id, t.min, diff, sid>}
#'     \item \code{sessions}: summary of sessions, \code{<id, t, t.min, sid, n>}
#'   }
.find_sessions <- function(data, thres = 10, timeScale = "min"){
  data$diff <- c(0, diff(data$t.min))  # get times between events
  row.names(data) <- 1:nrow(data)
  data.sessions <- data[data$diff == 0 | data$diff >= thres, c("id", "t", "t.min")]  # find session starts
  data.sessions$sid <- 1:nrow(data.sessions)  # number sessions
  data.ses.len <- diff(as.numeric(row.names(data.sessions)))
  data.ses.len <- c(data.ses.len, nrow(data) - sum(data.ses.len))
  data.sessions$n <- data.ses.len  # n events in each session
  data$sid <- rep(data.sessions$sid, data.ses.len)  # label sessions in raw data
  return ( list(data = data, sessions = data.sessions) )
}

###################################################################################
#' Simulate a pair of event time series.
#'
#' @param lambda Intensity of first process (i.e. of Mark 1/A).
#' @param W Observation window \code{c(low, high)}.
#' @param nSamp Number of samples to draw, intentionally too large.
#' @param indep Logical indicating if the second process should depend on the first;
#'   \code{default == FALSE}
#'   \itemize{
#'     \item \code{TRUE}: homogeneous Poisson process of rate lambda2
#'     \item \code{FALSE}: draw \eqn{Bern(p)}, if 1 draw point from \eqn{N(t_1[i],sigma^2)}
#'   }
#' @param p Bernoulli probability if \code{indep == FALSE}
#' @param normal Logical indicating distribution used to simulate dependent events;
#'   \code{default == TRUE}
#'   \itemize{
#'     \item \code{TRUE}: Gaussian(t_1[i], sigma^2)
#'     \item \code{FALSE}: Exponential(sigma)
#'   }
#' @param sigma Standard deviation of normal if \code{indep == FALSE}
#' @param lambda2 intensity of second process if  \code{indep == TRUE}
#' @return Data.frame of \code{<t, m>}
sim_markedpp <- function(lambda, W, nSamp, indep = FALSE, p, normal = TRUE, sigma, lambda2){
  x <- cumsum( rexp(n = nSamp, rate = lambda) )  # homogeneous poisson process
  x.w <- x[x <= W[2]]  # Poisson process restricted to window W
  n.x <- length(x.w)  # number of points in window W

  if(indep == FALSE){  # generate second point process dependent on first process
    d <- rbinom(nSamp, size=1, prob=p)  # wether or not to draw a point from second process
    y <- rep(0, nSamp)  # place holder for second pp
    for(i in 1:nSamp){
      if(d[i] == 1){
        if(normal == TRUE){
          y[i] <- rnorm(n = 1, mean = x[i], sd = sigma)
        }else{
          y[i] <- x[i] + rexp(n = 1, rate = 1/sigma)
        }
      }
      else{
        y[i] <- -1
      }
    }
    y.w <- y[y >= W[1] & y <= W[2]]  # second pp restricted to W
    y.w <- y.w[order(y.w)]  # ensure it is in increasing order
    n.y <- length(y.w)  # number of points in W
  } else{  # generate indep homogeneous Poisson process of rate lambda2
    y <- cumsum( rexp(n = nSamp, rate = lambda2) )  # homogeneous poisson process
    y.w <- y[y <= W[2]]  # Poisson process restricted to window W
    n.y <- length(y.w)  # number of points in window W
  }

  ### set up data frame
  dat <- data.frame( cbind(t=c(x.w, y.w), m=c(rep(1,n.x),rep(2,n.y))) )
  dat <- dat[order(dat$t),]
  dat$m <- as.factor(dat$m)
  rownames(dat) <- 1:nrow(dat)

  return(dat)
}

###################################################################################
#' Plot single pair of event series with cumulative counts.
#' @param x Data.frame of \code{<t, m>} single pair of event series.
#' @param W Observation window \code{c(low, high)}, equivalent to \code{xlim}.
#' @param c Numeric cex value for axis & legend magnification; \code{default == 1}
#' @param lbl Vector of labels for the legend corresponding to marks 1/A and 2/B,
#'   respectively; \code{default == c("Mark 1", "Mark 2")}
#' @param title Main title of plot; \code{default == ""}
#' @param xlab X-axis label; \code{default == "Time"}
plot_markedpp <- function(x, W, c=1, lbl=c("Mark 1", "Mark 2"), title="", xlab="Time"){
  ind <- which(x$m == 1)  # get index for type 1 points
  x.1 <- x$t[ind]
  n1 <- length(x.1)
  x.2 <- x$t[-ind]
  n2 <- length(x.2)
  plot(c(0,x.1), 0:n1, type = "s", xlim = W, ylim=c(0, max(n1,n2)),
       xlab = xlab, ylab = "Total Number of Events", cex.axis=c, cex.lab=c,
       main = title)
  lines(c(0,x.2), 0:n2, type = "s", col="red")
  points(x.1, rep(-1,n1), pch = 20, cex=0.45)
  points(x.2, rep(-5,n2), pch = 20, cex=0.45, col="red")
  legend("topleft", inset=c(.05,.05), bty="n", cex=c,
         legend=c(lbl[1],lbl[2]), fill=c("black","red"), border=c("black","red"))
}

###################################################################################
#' Plot multiple pairs of event series.
#'
#' @inheritParams plot_markedpp
#' @param x Data.frame of \code{<id, t, m>} multiple pairs of event series.
#' @param c Numeric cex value for point size; \code{default == 1}
#' @param cc Numeric cex value for axis & legend magnification; \code{default == 1}
#' @param leg Logical indicating whether or not to include legend; \code{default == TRUE}
#' @param leg.loc Where to locate legend.
plot_mult_markedpp <-function(x, W, c=1, leg=TRUE, lbl=c("Mark 1", "Mark 2"), leg.loc, cc=1){
  ids <- unique(x$id)
  plot(1,
       ylim = c(0.5, length(ids) + 0.5),
       xlim = W,
       type = "n",
       ylab = "User ID",
       xlab = "t",
       yaxt = "n",
       cex.lab = cc,
       cex.axis = cc )
  abline(h = 0.5,
         col = gray(0.8, alpha=0.9),
         lwd = 0.5)
  for(i in 1:length(ids)){
    abline(h = i + 0.5,
           col = gray(0.8, alpha=0.9),
           lwd = 0.5)
    ind <- x$id == ids[i] & x$m == 1
    points(x = x$t[ind],
           y = rep(i - 0.15, sum(ind)),
           pch = "|",
           cex = c,
           col = gray(0, alpha=0.4))
    ind <- x$id == ids[i] & x$m == 2
    points(x = x$t[ind],
           y = rep(i + 0.15, sum(ind)),
           pch = "|",
           cex = c,
           col = gray(0.6, alpha=0.4))
  }
  axis(2,
       at = 1:length(ids),
       labels = as.character(ids),
       las = 1,
       cex.axis=cc*0.75)
  if(leg == TRUE){
    legend(2.25, 11.75,
           legend = c("Mark 1", "Mark 2"),
           pch = "|",
           col = c(gray(0.6), gray(0)),
           horiz = TRUE,
           bty = "n",
           cex = cc,
           xpd = TRUE)
  }
}
