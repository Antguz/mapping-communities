#'------------------------------------------------------------------------------
#' @title Compute species approximation of the predicted values
#'------------------------------------------------------------------------------

#' @description It computes the approximation of the predicted values to species.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(bit64)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/approximation_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param communities_projections
#' @param species_projections

#'------------------------------------------------------------------------------
#' Function

approximation <- function(communities_projections,
                          species_projections) {
  
  #---------------------------------------------------------
  # Estimate distance
  distance <- approximation_rcpp(as.matrix(communities_projections[, 2:4]), 
                                 as.matrix(species_projections[, 2:4]))
  
  #Get the inverse
  for(i in 1:dim(distance)[3]) {
    
    subdistance <- distance[,,i]
    total <- rowSums(subdistance, na.rm = FALSE)
    subdistance <- (subdistance/total)*10000000
    subdistance <- round(subdistance, 0)
    colnames(subdistance) <- as.vector(as.matrix(species_projections[, 1]))
    distance[,,i] <- subdistance
    
  }
  
  return(distance)
  
}
