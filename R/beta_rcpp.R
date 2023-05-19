#'------------------------------------------------------------------------------
#' @title High performance beta diversity: a modification from BAT
#'------------------------------------------------------------------------------

#' @description It estimates beta diversity using phylogeny and abundance

#'------------------------------------------------------------------------------
#' Library

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)

#'------------------------------------------------------------------------------
#' Source code
sourceCpp("src/beta_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Source code

beta_rcpp <- function(comm, tree, threads = 1L, progress = TRUE) {
  
  treea <- xTree(tree)
  xtree_1 <- treea[[1]]
  xtree_2 <- treea[[2]]
    
  results <- beta_arma(comm = comm, 
                       xtree1 = xtree_1,
                       xtree2 = xtree_2,
                       threads = threads,
                       progress = progress)
    
  colnames(results) <- rownames(comm)
  rownames(results) <- rownames(comm)
  
  #as.dist(
  results <- list(Btotal = as.dist(results[,,1]),
                  Brepl = as.dist(results[,,2]),
                  Brich = as.dist(results[,,3]))
  
  return(results)
  
}

prep <- function(comm, xtree, abund = TRUE){
  len <- xtree[[1]] 							## length of each branch
  A <- xtree[[2]]									## matrix species X branches
  minBranch <- min(len[colSums(A)==1]) 	## minimum branch length of terminal branches
  if(is.data.frame(comm))
    comm = as.matrix(comm)
  BA <- comm%*%A 												## matrix samples X branches
  if (!abund)	BA = ifelse(BA >= 1, 1, 0)
  return (list(lenBranch = len, sampleBranch = BA, speciesBranch = A, minBranch = minBranch))
}

clean <- function(comm, tree = NA){
  if(is.vector(comm))
    comm <- matrix(comm, nrow = 1)
  comm <- as.matrix(comm)
  if (!missing(tree)){
    comm = reorderComm(comm, tree)
    tree <- xTree(tree)
  }
  return(list(comm, tree))
}

reorderComm <- function(comm, tree = NULL){
  
  if(is.vector(comm))
    comm = as.matrix(comm, nrow = 1)
  if(is(tree, "hclust"))
    tree = ape::as.phylo(tree)
  
  if(is(tree, "phylo")){
    if(!is.null(tree$tip.label) && !is.null(colnames(comm))){ ##if both tree and comm have species names match and reorder species (columns) in comm
      
      #if some species are missing from comm add 0s
      if(length(tree$tip.label) > ncol(comm) && all(colnames(comm) %in% tree$tip.label)){
        miss = tree$tip.label[which(!(tree$tip.label %in% colnames(comm)))]
        addComm = matrix(0, nrow = nrow(comm), ncol = length(miss))
        colnames(addComm) = miss
        comm = cbind(comm, addComm)
      }
      if(length(dim(comm)) == 2)
        comm <- comm[,match(tree$tip.label, colnames(comm)), drop = FALSE]
      else
        comm <- comm[,match(tree$tip.label, colnames(comm)),, drop = FALSE]
      if (any(tree$tip.label != colnames(comm)))
        warning("Species names of comm and tree do not match!")
    }
  } else if(is(tree, "dist")){
    if(!is.null(colnames(tree)) && !is.null(colnames(comm))){ ##if both tree and comm have species names match and reorder species (columns) in comm
      comm <- comm[,match(colnames(tree), colnames(comm))]
      if (any(colnames(tree) != colnames(comm)))
        warning("Species names of comm and distance do not match!")
    }
  }
  return(comm)
}

#Function from BAT
xTree <- function(tree) {
  if (class(tree) == "hclust"){
    nSpp <- nrow(as.data.frame(tree['order']))
    sppEdges <- matrix(0, nSpp, 2 * nSpp - 2)
    lenEdges <- vector("numeric", 2 * nSpp - 2)
    for(i in 1:(nSpp - 1)) {
      if(tree$merge[i, 1] < 0) {
        lenEdges[2 * i - 1] <- tree$height[order(tree$height)[i]]
        sppEdges[ - tree$merge[i, 1], 2 * i - 1] <- 1
      } else {
        lenEdges[2 * i - 1] <- tree$height[order(tree$height)[i]] - tree$height[order(tree$height)[tree$merge[i, 1]]]
        sppEdges[, 2 * i - 1] <- sppEdges[, 2 * tree$merge[i, 1] - 1] + sppEdges[ , 2 * tree$merge[i, 1]]
      }
      if(tree$merge[i, 2] < 0) {
        lenEdges[2 * i] <- tree$height[order(tree$height)[i]]
        sppEdges[ - tree$merge[i, 2], 2 * i] <- 1
      } else {
        lenEdges[2 * i] <- tree$height[order(tree$height)[i]] - tree$height[order(tree$height)[tree$merge[i, 2]]]
        sppEdges[, 2 * i] <- sppEdges[, 2 * tree$merge[i, 2] - 1] + sppEdges[, 2 *tree$merge[i, 2]]
      }
    }
    rownames(sppEdges) <- tree$labels
    list(lenEdges, sppEdges)
  } else if (class(tree) == "phylo"){
    lenEdges <- tree$edge.length
    nSpp <- length(tree$tip.label)
    nEdges <- length(tree$edge.length)
    root <- nSpp + 1
    sppEdges <- matrix(0, nSpp, nEdges)
    for(i in 1:nSpp){
      find = i                                    #start by finding the ith species
      repeat{
        row = which(tree$edge[,2] == find)        #locate in which row of the edge table is our species or edge to be found
        sppEdges[i, row] = 1
        find = tree$edge[row,1]                   #find next edge if any until reaching the root
        if(find == root) break                    #all edges of this species were found, go to next species
      }
    }
    rownames(sppEdges) <- tree$tip.label
    list(lenEdges, sppEdges)
  } else {
    cat("Unrecognized tree object!")
  }
}
