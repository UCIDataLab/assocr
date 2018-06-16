################################################################################
### This file contains functions to calculate measures of assocaiation (score
### fucntions) for pairs of event series.
################################################################################

#' Return smallest non-zero element of a vector.
#'
#' @param x A vector of real numbers.
#' @return Minimum value of \code{x} such that \code{min(x) > 0}
#' @examples
#' nonzero_min(c(1,3,5,2,0))

nonzero_min <- function(x){
  return( min(x[x > 0]) )
}

################################################################################
#' Return the ceiling of a decimal number as a decimal with user-specified
#' precision.
#'
#' @param x Positive double.
#' @param level Number of decimal places to return; \code{default == 0}
#' @return Decimal ceiling of \code{x}.

ceiling_dec <- function(x, level = 0){
  # return(trunc(x*10^level + 0.5) / 10^level)
  return( round(x + 5*10^(-level-1), level) )
}

################################################################################
#' Count the number of leading zeros between decimal and first non-zero entry.
#'
#' @inheritParams ceiling_dec
#' @return Number of leading zeros.

count_zeros <- function(x){
  return( attr(regexpr("(?<=\\.)0+", format(x, scientific = FALSE), perl = TRUE),
               "match.length") )
}

################################################################################
#' Calculate the coefficient of segregation, mingling index, and inter-event
#' time summary statistics of a given pair of event series.
#'
#' @param x Data.frame of \code{<m, t>} for a single pair of event series.
#' @param W Vector containing observation window \code{c(low, high)}.
#' @param bidirectional Logical indicating if the inter-event times should be
#'   bidirectional (versus only looking at events prior to event of interest)
#' @return List containing
#'   \itemize{
#'     \item \code{iet.mn}: mean iet for fixed mark 2/B to mark 1/A
#'     \item \code{iet.md}: median iet for fixed mark 2/B to mark 1/A
#'     \item \code{s}: coefficient of segregation
#'     \item \code{m}: mingling index, k=1
#'     \item \code{nA}: number of events of type 1/A
#'     \item \code{nB}: number of events of type 2/B
#'   }

calc_score_funcs <- function(x, W = c(0, 7), bidirectional = TRUE){
  ### calculate iets
  x <- x[order(x$t), ]
  xA <- x[which(x$m == 1),]
  xB <- x[which(x$m == 2),]
  if(nrow(xB) == 0){
    return( list(iet.mn = NA,
                 iet.md = NA,
                 s = NA ,
                 m = NA,
                 nA = nrow(xA),
                 nB = 0)) }
  if(bidirectional == TRUE){
    BA <- sapply(xB$t, function(e) nonzero_min(abs(xA$t - e)) )
  }else{  # only compute iet to nearest event of type A before it
    BA <- sapply(xB$t, function(e) nonzero_min(e - xA$t))
    BA[is.infinite(BA)] <- NA  # replace obs w/o A events before B events w/ NAs
  }

  ### compute nn probability table 5.2 p 313 of Ilian et al
  n <- c(nrow(xA), nrow(xB))  # n of each type
  counts <- matrix(0, nrow=2, ncol=2)
  d <- as.matrix(dist(x$t))  # matrix of pairwise distances
  diag(d) <- NA
  nn <- matrix(apply(d, 2, which.min))
  for(i in 1:sum(n)){
    counts[x$m[i], x$m[nn[i]]] <- counts[x$m[i],x$m[nn[i]]] + 1
  }
  p <- counts / sum(n)  # joint probs
  p_dot <- colSums(p)  # marginalize reference point
  markp <- rowSums(p)  # mark probabiliies

  ### coefficient of segregation (p 314)
  s <- 1 - (p[1,2] + p[2,1]) / (markp[1]*p_dot[2] + markp[2]*p_dot[1])

  ### mingling index with k=1 (estimate via 5.2.24)
  m <- (counts[1,2] + counts[2,1]) / sum(counts)

  return( list(iet.mn = mean(BA, na.rm=TRUE),
               iet.md = median(BA, na.rm=TRUE),
               s = s,
               m = m,
               nA = nrow(xA),
               nB = nrow(xB)  ) )
}

################################################################################
#' Calculate score-based likelihood ratios & coincidental match probabilities
#' for simulated pairs of event series.
#'
#' @param same_src,diff_src Data.frame of score functions for same- and
#'   different-source pairs of event series, respectively, computed via
#'   \code{assocr::calc_score_funcs}; i.e., \code{<iet.mn, iet.md, s, m>}
#' @param bds.iet,bds.s,bds.m Vector of bounds \code{c(low, high)} for density
#'   estimation of mean/median inter-event times, segregation and mingling,
#'   respectively.
#' @param bw.type Bandwidth estimation method to use in density estimation,
#'   one of \code{c("nrd", "nrd0", "ucv", "bcv", "SJ")}
#' @return Data.frame of SLRs and CMPs for each score function where
#'   \code{indep == 1} if the pair was independent; i.e.,
#'   \code{<iet.mn, iet.md, s, m, indep, slr.iet.mn, slr.iet.md, slr.s, slr.m,
#'          cmp.iet.mn, cmp.iet.md, cmp.s, cmp.m>}

calc_slr_cmp <- function(same_src, diff_src, bds.iet, bds.s = c(-1,1),
                         bds.m = c(0,0.6), bw.type = "nrd"){
  rslt <- data.frame(matrix(NA,
                            nrow = nrow(same_src) + nrow(diff_src),
                            ncol = 9))
  names(rslt) <- c("indep", "slr.iet.mn", "slr.iet.md", "slr.s", "slr.m",
                   "cmp.iet.mn", "cmp.iet.md", "cmp.s", "cmp.m")
  rslt$indep <- c( rep(0, nrow(same_src)), rep(1, nrow(diff_src)) )
  rslt <- cbind(rbind(same_src, diff_src), rslt)

  for(i in 1:nrow(rslt)){
    tmp <-rslt[-i, ]  # remove current obs for calculation of density funcs
    ### same source densities
    same <- tmp[tmp$indep == 0, ]
    p.mn.win <- approxfun(density(same$iet.mn, from=bds.iet[1], to=bds.iet[2],
                                  bw=bw.type))
    p.md.win <- approxfun(density(same$iet.md, from=bds.iet[1], to=bds.iet[2],
                                  bw=bw.type))
    p.s.win <- approxfun(density(same$s, from=bds.s[1], to=bds.s[2], bw=bw.type))
    p.m.win <- approxfun(density(same$m, from=bds.m[1], to=bds.m[2], bw=bw.type))

    ### different source densities
    diff <- tmp[tmp$indep == 1, ]
    p.mn.bt <- approxfun(density(diff$iet.mn, from=bds.iet[1], to=bds.iet[2],
                                 bw=bw.type))
    p.md.bt <- approxfun(density(diff$iet.md, from=bds.iet[1], to=bds.iet[2],
                                 bw=bw.type))
    p.s.bt <- approxfun(density(diff$s, from=bds.s[1], to=bds.s[2], bw=bw.type))
    p.m.bt <- approxfun(density(diff$m, from=bds.m[1], to=bds.m[2], bw=bw.type))

    ### compute score-based likelihood ratio
    rslt$slr.iet.mn[i] <- p.mn.win(rslt$iet.mn[i]) / p.mn.bt(rslt$iet.mn[i])
    rslt$slr.iet.md[i] <- p.md.win(rslt$iet.md[i]) / p.md.bt(rslt$iet.md[i])
    rslt$slr.s[i] <- p.s.win(rslt$s[i]) / p.s.bt(rslt$s[i])
    rslt$slr.m[i] <- p.m.win(rslt$m[i]) / p.m.bt(rslt$m[i])

    ### compute coincidental match probability
    rslt$cmp.iet.mn[i] <- sum(diff$iet.mn < rslt$iet.mn[i]) / nrow(diff)
    rslt$cmp.iet.md[i] <- sum(diff$iet.md < rslt$iet.md[i]) / nrow(diff)
    rslt$cmp.s[i] <- sum(diff$s < rslt$s[i]) / nrow(diff)
    rslt$cmp.m[i] <- sum(diff$m > rslt$m[i]) / nrow(diff)
  }

  return( rslt )
}

################################################################################
#' Calculate the coincidental match probability for a single pairs of event
#' series.
#'
#' @inheritParams calc_score_funcs
#' @param data List of data.frames for one pair of event series (output of
#'   \code{assocr::sessionize_data()}):
#'   \itemize{
#'     \item \code{data}: sessionized data, \code{<id, m, sid, t>}
#'     \item \code{sessions}: summary of sessionized data,
#'           \code{<id, m, sid, n, t>}
#'   }
#' @param n Number of samples to draw
#' @param samp Sampling technique to use. One of
#'   \itemize{
#'     \item \code{empirical}: sessionized resampling from empirical distn of
#'       start times
#'     \item \code{periodic}: sessionized resampling from sinusoidal curve
#'       fitted to the data
#'   }
#' @param rng Vector of \code{c(low, high)} limits for sampling session start
#'   times if \code{samp == "periodic"}
#' @return Data.frame of CMPs for each score function;
#'   \code{<iet.mn, iet.md, s, m, cmp.iet.mn, cmp.iet.md, cmp.s, cmp.m>}
calc_cmp <- function(data, n, W = c(0,7), bidirectional = TRUE,
                     samp = "empirical", rng = NULL){
  sim <- data.frame(matrix(NA, nrow = n, ncol = 4))
  names(sim) <- c("iet.mn", "iet.md", "s", "m")
  dat.m2 <- data$data[data$data$m == 2, c("m", "t")]  # event series of mark 2
  n1 <- sum(data$data$m == 1)  # number of events of mark 1
  for (i in 1:n) {# simulate start times of mark 1 & calc score funcs for sim pair
    t <- .session_resampling(data, samp = samp, rng = rng)
    tmp <- rbind(dat.m2, cbind(t, m = rep(1, n1)))
    sim[i, ] <- as.numeric(calc_score_funcs(tmp, W, bidirectional)[1:4])
  }

  ### compute & return coincidental match probability
  obs <- calc_score_funcs(data$data, W, bidirectional)  # observed pair

  return( list(iet.mn = sum(sim$iet.mn < obs$iet.mn) / n,
               iet.md = sum(sim$iet.md < obs$iet.md) / n,
               s = sum(sim$s < obs$s) / n,
               m = sum(sim$m > obs$m) / n
               )
          )
}

################################################################################
#' Resample session start times.
#'
#' @inheritParams calc_cmp
#' @return Vector of resampled times for event series of mark 1.
.session_resampling <- function(data, samp = "empirical", rng = NULL){
  # sample new session start times for events of mark 1
  ind <- data$sessions$m == 1
  nSamp <- sum(ind)
  if (samp == "empirical") {
    stop("Haven't set this up yet!")
  } else if (samp == "periodic") {
    mu <- mean(rng)
    sigma <- (rng[2] - mu) / 3  # 99% of start times fall in rng
    t.ses <- rnorm(nSamp, mean = mu, sd = sigma)
  } else {
    stop("Invalid value for samp; enter one of (\"empirical\", \"periodic\")")
  }
  sim.t <- data$data$t[data$data$m == 1] -  # current times of mark 1 events
    with(data$sessions[ind, ], rep(t, n)) +  # subtract current session start times
    rep(t.ses, data$sessions$n[ind])  # add new start times in
  return (sim.t)
}

################################################################################
#' Wrapper to compute CMP for each combination of a given pair of user data
#'
#' @inheritParams calc_cmp
#' @param data List of data.frames for two pairs of event series (output of
#'   \code{assocr::sessionize_data()}):
#'   \itemize{
#'     \item \code{data}: sessionized data, \code{<id, m, sid, t>}
#'     \item \code{sessions}: summary of sessionized data,
#'           \code{<id, m, sid, n, t>}
#'   }
#' @return Array of CMPs for each pairwise combination of event streams &
#'   score function
pairwise_cmp <- function(data, n, W = c(0,7), bidirectional = TRUE,
                         samp = "empirical", rng = NULL){
  ids <- levels(data$data$id)

  r1 <- .filter_list(data = data, id = ids[1])
  r2 <- .filter_list(data = data, id = ids[2])
  r3 <- .filter_list(data = data, id = ids, m = c(1,2))
  r4 <- .filter_list(data = data, id = ids, m = c(2, 1))
  out <- rbind(calc_cmp(r1, n, W, bidirectional, samp, rng),
               calc_cmp(r2, n, W, bidirectional, samp, rng),
               calc_cmp(r3, n, W, bidirectional, samp, rng),
               calc_cmp(r4, n, W, bidirectional, samp, rng))
  rownames(out) <- c(ids[1],
                     ids[2],
                     paste0(ids[1], "m1_", ids[2], "m2"),
                     paste0(ids[2], "m1_", ids[1], "m2") )
  return ( out )
}

################################################################################
#' Filter list of data.frames to contain one pair of event series.
#'
#' @inheritParams pairwise_cmp
#' @param id Either of the following
#'   \itemize{
#'     \item String containing ID of user for same-source pair of event streams
#'     \item Vector of strings containing IDs for different-source pair of
#'       event streams
#'   }
#' @param m Array of c(mark of id 1, mark of id 2) used for different-source
#'   pairs of event streams if \code{id == c(id1, id2)}
.filter_list <- function(data, id, m = NULL){
  if (length(id) == 1) {  # streams are from same user
    data$data <- data$data[data$data$id == id, ]
    data$sessions <- data$sessions[data$sessions$id == id, ]
  } else {  # streams are from different users
    data$data <- data$data[(data$data$id == id[1] & data$data$m == m[1]) |
                           (data$data$id == id[2] & data$data$m == m[2]),]
    data$sessions <- data$sessions[
      (data$sessions$id == id[1] & data$sessions$m == m[1]) |
        (data$sessions$id == id[2] & data$sessions$m == m[2]), ]
  }
  return ( data )
}

################################################################################
# ### DESCRIPTION - Calculate quantiles from simulation
# ### INPUT - scores: data.frame of SLRs and cmps for each score function
# ###             where indep == 1 if the marked point process was independent
# ###             <iet.mn, iet.md, indep, slr.iet.mn, slr.iet.md, cmp.iet.mn,
# ###              cmp.iet.md>
# ###         which.score: (default == "slr.iet.mn") which score function to
# ###                      evaulate
# ###         p: relative frequency of B to A events
# ###         sigma: std dev of normal distn used to generate dependent events
# ###                (in sec)
# ### OUTPUT - data.frame of <p, sigma, indep, lower 2.5%, median, upper 97.5%>
#
# summarize_rslts <- function(scores, which.score = "slr.iet.mn", p, sigma){
#   col.ind <- c("indep", which.score)
#   rslt <- data.frame(matrix(NA, nrow=2, ncol=6))
#   names(rslt) <- c("p", "sigma", "indep", "l", "est", "u")
#   rslt$p <- p
#   rslt$sigma <- sigma
#   rslt$indep <- c(0,1)
#   rslt[1, 4:6] <- quantile(scores[scores$indep == 0, which.score],
#                            probs = c(0.025, 0.5, 0.975),
#                            na.rm = TRUE)
#   rslt[2, 4:6] <- quantile(scores[scores$indep == 1, which.score],
#                            probs = c(0.025, 0.5, 0.975),
#                            na.rm = TRUE)
#   return( rslt )
# }
