is.tip <- function(tree, node) {
  return(node <= length(tree$tip.label))
}

mrca.phylo.or.unique.tip <- function(tree, node) {
  if (length(node) == 1) {
    return(node)
  } else {
    mrca <- mrca.phylo(tree, node)
    return(mrca)
  }
}

get.node.height <- function(tree, node){
  return(max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[node])
}

get.parents.vector <- function(tree, a.sample, interventions = NULL){
  if(is.null(interventions)){
    interventions <- rep(0, length(tree$tip.label) + tree$Nnode)
  }
  parents <- vector()
  for(host in 1:max(a.sample)){
    a.node <- which(a.sample==host)[1]
    parent.node <- Ancestors(tree, a.node, type="parent")
    move.up <- parent.node != 0
    if(move.up){
      move.up <- a.sample[a.node] == a.sample[parent.node]
    }
    while(move.up){
      a.node <- parent.node
      parent.node <- Ancestors(tree, a.node, type="parent")
      move.up <- parent.node != 0
      if(move.up){
        move.up <- a.sample[a.node] == a.sample[parent.node]
      }
    }
    if(interventions[a.node] == 0){
      if(a.node == getRoot(tree)){
        parents[host] <- 0
      } else {
        parents[host] <- a.sample[parent.node]
      }
    } else {
      current.host <- host
      for(intervening in 1:interventions[a.node]){
        current.host.count <<- current.host.count + 1
        parents[current.host] <- current.host.count
        current.host <- parents[current.host]
      }
      if(getRoot(tree) == a.node){
        parents[current.host] <- 0
      } else {
        parents[current.host] <- a.sample[parent.node]
      }
    }
  }
  return(parents)
}