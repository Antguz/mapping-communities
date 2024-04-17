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

beta_rcpp <- function(comm, tree = NULL, taxonomic = TRUE, threads = 1L, progress = TRUE) {
  
  # Phylogenetic or Functional
  if(!taxonomic) {
    
    cleanData = clean(comm, tree)
    comm = cleanData[[1]]
    treea = cleanData[[2]]
    
    xtree_1 <- treea[[1]]
    xtree_2 <- treea[[2]]
  
  # Taxonomic
  } else {
    
    comm <- as.matrix(comm)
    xtree_1 <- 1
    xtree_2 <- as.matrix(2)
    
  }
  
  results <- beta_arma(comm = comm, 
                       xtree1 = xtree_1,
                       xtree2 = xtree_2,
                       taxonomic = taxonomic,
                       threads = threads,
                       progress = progress)
    
  colnames(results) <- rownames(comm)
  rownames(results) <- rownames(comm)
  
  # as.dist
  results <- list(Btotal = results[,,1],
                  Brepl = results[,,2],
                  Brich = results[,,3])
  
  return(results)
  
}

reorderComm <- function(comm, tree = NULL) {
  
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

clean <- function(comm, tree = NA) {
  if(is.vector(comm))
    comm <- matrix(comm, nrow = 1)
  comm <- as.matrix(comm)
  if (!missing(tree)){
    comm = reorderComm(comm, tree)
    tree <- xTree(tree)
  }
  return(list(comm, tree))
}

#Function from BAT
xTree <- function(tree) {
  if(is(tree, "hclust")){
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
    return(list(lenEdges, sppEdges))
    
  } else if(is(tree, "phylo")){
    
    edgeList = tree$edge # edge list
    edgeLength = tree$edge.length # edge length
    spp = tree$tip.label # species
    basal = min(edgeList[,1]) # Basal node
    mat = matrix(data = 0, nrow = length(spp), ncol = length(edgeLength)) # empty matrix
    for (i in 1:length(spp)) {
      x = i
      repeat {
        mat[i, which(edgeList[,2] == x)] = 1
        x = edgeList[which(edgeList[,2] == x), 1]
        if (x == basal){
          break
        }
      }
    }
    return(list(edgeLength, mat))
    
  } else {
    cat("Unrecognized tree object!")
  }
}

