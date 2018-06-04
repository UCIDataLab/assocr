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
format.data <- function(data, id.sum, dev=FALSE, low.bds){
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
sim.markedPP <- function(lambda, W, nSamp, indep = FALSE, p, normal = TRUE, sigma, lambda2){
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
#' @param lbl Vector of labels for the legend corresponding to marks 1/A and 2/B, respectively;
#'   \code{default == c("Mark 1", "Mark 2")}
#' @param title Main title of plot; \code{default == ""}
#' @param xlab X-axis label; \code{default == "Time"}
plot.markedpp <- function(x, W, c=1, lbl=c("Mark 1", "Mark 2"), title="", xlab="Time"){
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
#' @inheritParams plot.markedpp
#' @param x Data.frame of \code{<id, t, m>} multiple pairs of event series.
#' @param c Numeric cex value for point size; \code{default == 1}
#' @param cc Numeric cex value for axis & legend magnification; \code{default == 1}
#' @param leg Logical indicating whether or not to include legend; \code{default == TRUE}
#' @param leg.loc Where to locate legend.
plot.mult.markedPP <-function(x, W, c=1, leg=TRUE, lbl=c("Mark 1", "Mark 2"), leg.loc, cc=1){
  ids <- unique(x$id)
  plot( 1, ylim=c(1,length(ids)), xlim=W, type="n", ylab="Student ID", xlab="Day", yaxt="n",
        cex.lab = cc, cex.axis = cc )
  for(i in 1:length(ids)){
    ind <- which(x$id==ids[i])
    abline(h=i, col=gray(0.8, alpha=0.2), lwd=0.5)
    ind <- which( x$id == ids[i] & x$m == 1 )
    points( x$t[ind], rep(i,length(ind)), pch="|", cex=c, col=gray(0, alpha=0.35))
    ind <- which( x$id == ids[i] & x$m == 2 )
    points( x$t[ind], rep(i-1/3,length(ind)), pch="|", cex=c, col=gray(0.6, alpha=0.4))
  }
  axis(2, at=1:length(ids), labels=as.character(ids), las=1, cex.axis=cc*0.75)
  if(leg == TRUE){
    legend(leg.loc[1], leg.loc[2],
           legend=lbl, pch="|", col=c(gray(0), gray(0.6)),
           horiz=TRUE, bty="n", cex=cc, xpd=TRUE)
  }
}

###################################################################################
#' Calculate first order characteristics from Illian et al.
#'
#' @inheritParams plot.markedpp
#' @return List containing the following:
#'   \itemize{
#'     \item \code{lambda}: vector of intensities of each process and both combined
#'     \item \code{p}: matrix of joint probs of nearest neighbor marks, \code{p_11, p_12, p_21, p_22}
#'     \item \code{p_dot}: vector of margin values (columns) of p, \code{c(p_dot1, p_dot2)}
#'     \item \code{markp}: vector of mark probabilities, \code{c(p_1,p_2)}
#'     \item \code{S}: coefficient of segregation
#'     \item \code{M}: vector of mingling index with \code{k=1} for observed process and theoretical
#'       value if indep, respectively, \code{c(M_k, M_ki)}
#'   }
calc.firstOrd <- function(x, W){
  n <- c(sum(x$m == 1), sum(x$m == 2))

  ### intensity
  lambda <- c(n, sum(n)) / W[2]

  ### compute nn probability table 5.2 p 313
  counts <- matrix(0, nrow=2, ncol=2)
  d <- as.matrix(dist(x$t))  # matrix of pairwise distances
  diag(d) <- NA
  nn <- matrix(apply(d, 2, which.min))
  for(i in 1:sum(n)){
    counts[x$m[i], x$m[nn[i]]] <- counts[x$m[i],x$m[nn[i]]] + 1
  }
  p <- counts / sum(n)  # rowSums() gives mark probs, colSums() gives p_dot
  p_dot <- colSums(p)
  markp <- rowSums(p)

  ### coefficient of segregation (p 314)
  S <- 1 - (p[1,2] + p[2,1]) / (markp[1]*p_dot[2] + markp[2]*p_dot[1])

  ### mingling index (estimate via 5.2.24)
  M_k <- (counts[1,2] + counts[2,1]) / sum(counts)
  M_ki <- 2 * markp[1] * markp[2]

  ### format & return
  names(lambda) <- c("lambda_1", "lambda_2", "lambda")
  return( list(lambda=lambda, p=p, p_dot=p_dot, markp=markp,
               S=S, M=c(M_k=M_k, M_ki=M_ki)) )
}

###################################################################################
#' Calculate score-based likelihood ratios.
#'
#' @inheritParams calc.iet_SLR
#' @param s.bds,m.bds Bounds on density estimation of segregation and mingling, respectively;
#'   \code{default == c(-1,1), default == c(0, 0.6)}
#'
#' @return Data.frame of score-based likelihood ratios:
#'   \itemize{
#'     \item \code{match == TRUE}: \code{<id, s, m>}
#'     \item \code{match == FALSE}: \code{<id.A, id.B, s, m>}
#'   }

calc.SLR <- function(ids, within, between, match = TRUE, bw.type = "nrd", s.bds = c(-1,1), m.bds = c(0,0.6)){
  if(match == TRUE){  # compute SLR for known matches
    slr_s <- rep(NA, length(ids))
    slr_m <- rep(NA, length(ids))
    for(i in 1:length(ids)){
      ### "within" density not including student i
      ind.win <- which(within$id != ids[i])
      p.s.win <- approxfun(density(within$s[ind.win], from=s.bds[1], to=s.bds[2], bw=bw.type))
      p.m.win <- approxfun(density(within$m_k[ind.win], from=m.bds[1], to=m.bds[2], bw=bw.type))
      ### "between" density not including any stream from student i
      ind.bt <- which(between$id.fb != ids[i] & between$id.nonfb != ids[i])
      p.s.bt <- approxfun(density(between$s[ind.bt], from=s.bds[1], to=s.bds[2], bw=bw.type))
      p.m.bt <- approxfun(density(between$m_k[ind.bt], from=m.bds[1], to=m.bds[2], bw=bw.type))
      ### compute score-based likelihood ratio for each
      slr_s[i] <- p.s.win(within$s[i]) / p.s.bt(within$s[i])
      slr_m[i] <- p.m.win(within$m_k[i]) / p.m.bt(within$m_k[i])
    }
    return( data.frame(id = ids, s = slr_s, m = slr_m) )
  } else{  # compute SLR for known non-matches
    slr_s <- rep(NA, nrow(between))
    slr_m <- rep(NA, nrow(between))
    for(i in 1:nrow(between)){
      ids <- between[i,1:2]  # get ids of current row
      ### calc densities for within
      ind.win <- which( !(within$id %in% ids) )  # ids for within that arent current ids
      p.s.win <- approxfun(density(within$s[ind.win], from=s.bds[1], to=s.bds[2], bw=bw.type))
      p.m.win <- approxfun(density(within$m_k[ind.win], from=m.bds[1], to=m.bds[2], bw=bw.type))
      ### calc densities for between
      ind.bt <- which( !(between$id.fb %in% ids) & !(between$id.nonfb %in% ids) )  # same for between
      p.s.bt <- approxfun(density(between$s[ind.bt], from=s.bds[1], to=s.bds[2], bw=bw.type))
      p.m.bt <- approxfun(density(between$m_k[ind.bt], from=m.bds[1], to=m.bds[2], bw=bw.type))
      ### calc SLR
      slr_s[i] <- p.s.win(between$s[i]) / p.s.bt(between$s[i])
      slr_m[i] <- p.m.win(between$m_k[i]) / p.m.bt(between$m_k[i])
    }
    return( data.frame(between[,1:2], s = slr_s, m = slr_m) )
  }
}
