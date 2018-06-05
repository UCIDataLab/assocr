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
  xA <- x[which(x$m == 1),]
  xB <- x[which(x$m == 2),]
  if(nrow(xB) == 0){
    return( list(iet.mn = NA, iet.md = NA, s= NA , m = NA, nA = nrow(xA), nB = 0)) }
  if(bidirectional == TRUE){
    BA <- sapply(xB$t, function(e) nonzero_min(abs(xA$t - e)) )
  }else{  # only compute iet to nearest event of type A before it
    BA <- sapply(xB$t, function(e) nonzero_min(e - xA$t))
    BA[is.infinite(BA)] <- NA  # replace obs without A events before B events with NAs
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
#' Calculate score-based likelihood ratios & coincidental match probabilities.
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
  rslt <- data.frame(matrix(NA, nrow = nrow(same_src) + nrow(diff_src), ncol = 5))
  names(rslt) <- c("indep", "slr.iet.mn", "slr.iet.md", "slr.s", "slr.m",
                   "cmp.iet.mn", "cmp.iet.md", "cmp.s", "cmp.m")
  rslt$indep <- c( rep(0, nrow(same_src)), rep(1, nrow(diff_src)) )
  rslt <- cbind(rbind(same_src, diff_src), rslt)

  for(i in 1:nrow(rslt)){
    tmp <-rslt[-i, ]  # remove current obs for calculation of density funcs
    ### same source densities
    same <- tmp[tmp$indep == 0, ]
    p.mn.win <- approxfun(density(same$iet.mn, from=bds.iet[1], to=bds.iet[2], bw=bw.type))
    p.md.win <- approxfun(density(same$iet.md, from=bds.iet[1], to=bds.iet[2], bw=bw.type))
    p.s.win <- approxfun(density(same$s, from=bds.s[1], to=bds.s[2], bw=bw.type))
    p.m.win <- approxfun(density(same$m, from=bds.m[1], to=bds.m[2], bw=bw.type))

    ### different source densities
    diff <- tmp[tmp$indep == 1, ]
    p.mn.bt <- approxfun(density(diff$iet.mn, from=bds.iet[1], to=bds.iet[2], bw=bw.type))
    p.md.bt <- approxfun(density(diff$iet.md, from=bds.iet[1], to=bds.iet[2], bw=bw.type))
    p.s.bt <- approxfun(density(diff$s, from=bds.s[1], to=bds.s[2], bw=bw.type))
    p.m.bt <- approxfun(density(diff$m, from=bds.m[1], to=bds.m[2], bw=bw.type))

    ### compute score-based likelihood ratio
    rslt$slr.iet.mn[i] <- p.mn.win(rslt$iet.mn[i]) / p.mn.bt(rslt$iet.mn[i])
    rslt$slr.iet.md[i] <- p.md.win(rslt$iet.md[i]) / p.md.bt(rslt$iet.md[i])
    slr_s[i] <- p.s.win(rslt$s[i]) / p.s.bt(rslt$s[i])
    slr_m[i] <- p.m.win(rslt$m[i]) / p.m.bt(rslt$m[i])

    ### compute coincidental match probability
    rslt$cmp.iet.mn[i] <- sum(diff$iet.mn < rslt$iet.mn[i]) / nrow(diff)
    rslt$cmp.iet.md[i] <- sum(diff$iet.md < rslt$iet.md[i]) / nrow(diff)
    rslt$cmp.s[i] <- sum(diff$s < rslt$s[i]) / nrow(diff)
    rslt$cmp.m[i] <- sum(diff$m < rslt$m[i]) / nrow(diff)
  }

  return( rslt )
}

################################################################################
# ### DESCRIPTION - Calculate quantiles from simulation
# ### INPUT - scores: data.frame of SLRs and cmps for each score function
# ###             where indep == 1 if the marked point process was independent
# ###             <iet.mn, iet.md, indep, slr.iet.mn, slr.iet.md, cmp.iet.mn, cmp.iet.md>
# ###         which.score: (default == "slr.iet.mn") which score function to evaulate
# ###         p: relative frequency of B to A events
# ###         sigma: std dev of normal distn used to generate dependent events (in sec)
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
