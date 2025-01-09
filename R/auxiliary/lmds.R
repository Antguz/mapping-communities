#'------------------------------------------------------------------------------
#' @title Landmark Multi-Dimensional Scaling
#'------------------------------------------------------------------------------

#' This are a set of modified functions from lmds R package to perform Landmark
#' Multi-Dimensional Scaling.

#'------------------------------------------------------------------------------
#' @Library

library(assertthat)
library(Rcpp)
library(RcppArmadillo)
library(irlba)
library(bit64)

#'------------------------------------------------------------------------------
#' @Source-code

sourceCpp("src/column_wise_euclidean.cpp")

#'------------------------------------------------------------------------------
#' @Function

select_landmarks <- function(x,
                             landmarks_vector,
                             threads = 25) {
  
  # naive -> just subsample the cell ids
  dist_2lm <- column_wise_euclidean(x[landmarks_vector, , drop = FALSE], x, n_threads = threads)
  dist_2lm <- zapsmall(dist_2lm)
  colnames(dist_2lm) <- colnames(x)[landmarks_vector]
  rownames(dist_2lm) <- colnames(x)[landmarks_vector]
  
  dist_2lm
  
}

cmdscale_landmarks <- function(dist_2lm, ndim = 3, rescale = TRUE) {
  
  assert_that(
    !is.null(attr(dist_2lm, "landmark_ix")) || (!is.null(rownames(dist_2lm)) && !is.null(colnames(dist_2lm))),
    is.numeric(ndim), length(ndim) == 1, ndim >= 1,
    is.logical(rescale), length(rescale) == 1, !is.na(rescale)
  )
  
  ix_lm <- attr(dist_2lm, "landmark_ix")
  if (is.null(ix_lm)) {
    ix_lm <- match(as.integer64(rownames(dist_2lm)), as.integer64(colnames(dist_2lm)))
  }
  
  # short hand notations
  x <- dist_2lm[, ix_lm, drop = FALSE]^2
  n <- as.integer(nrow(x))
  N <- as.integer(ncol(dist_2lm))
  
  # double center data
  mu_n <- rowMeans(x)
  mu <- mean(x)
  x_c <- sweep(x, 1, mu_n, "-")
  x_dc <- sweep(x_c, 2, mu_n, "-") + mu
  
  # classical MDS on landmarks
  if (ndim > 0.5 * min(nrow(x_dc), ncol(x_dc))) {
    e <- eigen(-x_dc / 2, symmetric = TRUE)
    e$values <- e$values[seq_len(ndim)]
    e$vectors <- e$vectors[, seq_len(ndim), drop = FALSE]
  
  } else {
    e <- irlba::partial_eigen(-x_dc / 2, symmetric = TRUE, n = ndim)
  }
  
  ev <- e$values
  evec <- e$vectors
  ndim1 <- sum(ev > 0)
  
  if (ndim1 < ndim) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", ndim1, ndim), domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
    ndim <- ndim1
  }
  
  # distance-based triangulation
  points_inv <- evec / rep(sqrt(ev), each = n)
  dimred <- -t(t(points_inv) %*% ((dist_2lm - rep(mu_n, each = N)) / 2))
  dimred <- (-t(dist_2lm - rep(mu_n, each = N)) / 2) %*% points_inv
  
  if (rescale) {
    dimred <- dynutils::scale_uniform(dimred)
  }
  
  colnames(dimred) <- paste0("comp_", seq_len(ndim))
  
  dimred
  
}