# A function to normalize a tree from 0 to 1.
normalize_tree <- function(tree) {
  
  tree_copy <- tree
  
  # Get nodes
  node_depths <- node.depth.edgelength(tree_copy)
  
  # Max and min
  max_depth <- max(node_depths)
  normalized_depths <- node_depths / max_depth
  
  for (i in 1:nrow(tree_copy$edge)) {
    parent <- tree_copy$edge[i, 1]
    child <- tree_copy$edge[i, 2]
    tree_copy$edge.length[i] <- abs(normalized_depths[child] - normalized_depths[parent])
  }
  
  return(tree_copy)
  
}