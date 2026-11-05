#'------------------------------------------------------------------------------
#' @title Mantel test with Moran Spectral Randomization (MSR)
#'------------------------------------------------------------------------------

#' @description Spatially-constrained Mantel test. Replicates of m1 preserve its
#' spatial structure (power spectrum on the MEMs). This is a fast version of
#' adespatial::msr.mantelrtest (Crabot et al. 2019, Methods Ecol Evol,
#' 10:532-540) for the default "pair" method of adespatial::msr
#' (Wagner & Dray 2015, Methods Ecol Evol, 6:1169-1178). With the same seed, both
#' return the same replicates (see tests). Several matrices (m2) can be tested
#' against the same replicates of m1 in one pass.

#'------------------------------------------------------------------------------
#' Library

library(Rcpp)
library(RcppEigen)
library(RcppProgress)
library(ade4)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/mantel_msr_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' @param m1 Matrix of Euclidean distances to be randomized (e.g. sqrt(SAM))
#' @param m2 Matrix of distances kept fixed, or a named list of them
#' @param mem Orthobasis from adespatial::scores.listw(listw, MEM.autocor = "all")
#' or adespatial::dbmem(xy, MEM.autocor = "all")
#' @param nrepet Number of replicates
#' @param threads Number of threads
#' @param single Logic, use single precision for the products (faster)
#' @param block Size of the blocks of plots for the products
#' @param check_euclid Check that m1 is Euclidean (as adespatial::msr.mantelrtest).
#' "subset" checks three disjoint, evenly spaced subsets of up to 2000 plots (a
#' non-Euclidean subset implies a non-Euclidean matrix), "full" checks the whole
#' matrix (slow for large n), "none" skips the check
#' @param tol Stop if a replicate has squared distances below -tol, relative to
#' the mean squared distance to the centroid (rounding errors are far smaller)
#' @param progress Logic, provide progress bar

#'------------------------------------------------------------------------------
#' Random draws as in adespatial::msr.default (method = "pair", NCOL(x) > 1)

msr_draws <- function(m, nrepet) {

  n.pair <- m %/% 2
  mod <- m %% 2

  idx0 <- integer(nrepet)
  angles <- matrix(NA_real_, nrow = n.pair, ncol = nrepet)
  signs <- rep(1, nrepet)

  for(nr in 1:nrepet) {

    # adespatial:::cutbypair
    idx0[nr] <- sample(1:m, 1)

    # adespatial:::genbypair (one angle per pair)
    for(i in 1:n.pair) {
      angles[i, nr] <- runif(1, 0, 2*pi)
    }

    # adespatial:::genbysingleton (last row of pairs)
    if(mod > 0) {
      signs[nr] <- sample(c(-1, 1), 1)
    }
  }

  return(list(idx0 = idx0, angles = angles, signs = signs))

}

#'------------------------------------------------------------------------------
#' Function

mantel_msr <- function(m1,
                       m2,
                       mem,
                       nrepet = 999,
                       threads = 1,
                       single = TRUE,
                       block = 2048,
                       check_euclid = c("subset", "full", "none"),
                       tol = 1e-3,
                       progress = TRUE) {

  cl <- match.call()
  m1 <- as.matrix(m1)
  storage.mode(m1) <- "double"
  n <- nrow(m1)

  is_list <- is.list(m2) && !inherits(m2, "dist")
  if(!is_list) {
    m2 <- list(m2)
  }
  m2 <- lapply(m2, function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  })

  # Checks as in adespatial::msr.default and adespatial::msr.mantelrtest
  if(!identical(attr(mem, "weights"), rep(1/n, n))) {
    stop("Only implemented for uniform weights")
  }
  if(ncol(mem) != n - 1) {
    stop("orthobasis object should contain n-1 vectors")
  }
  check_euclid <- match.arg(check_euclid)
  if(check_euclid == "full") {
    subsets <- list(1:n)
  } else if(check_euclid == "subset") {
    subsets <- split(1:n, rep(1:3, length.out = n))
    subsets <- lapply(subsets, function(x) x[unique(round(seq(1, length(x), length.out = min(length(x), 2000))))])
  } else {
    subsets <- list()
  }
  for(sub in subsets) {
    if(!suppressWarnings(is.euclid(as.dist(m1[sub, sub])))) {
      stop("This test only works for Euclidean distances (m1 is not).")
    }
  }

  # MEMs with unit norm
  vectors <- as.matrix(mem) / sqrt(n)

  draws <- msr_draws(m = ncol(vectors), nrepet = nrepet)

  res <- mantel_msr_rcpp(d1 = m1,
                         d2 = m2,
                         mem = vectors,
                         idx0 = draws$idx0,
                         angles = draws$angles,
                         signs = draws$signs,
                         threads = threads,
                         block = block,
                         single = single,
                         tol = tol,
                         progress = progress)

  test <- lapply(seq_along(m2), function(k) {
    as.randtest(obs = res$obs[k],
                sim = res$sim[, k],
                output = "full",
                call = cl)
  })
  test <- lapply(test, function(x) {
    x$min_d2 <- res$min_d2
    x
  })
  names(test) <- names(m2)

  if(!is_list) {
    test <- test[[1]]
  }

  return(test)

}
