node.count <- function(tree){
  tree$Nnode + length(tree$tip.label)
}

is.tip <- function(tree, node) {
  return(node <= length(tree$tip.label))
}

# mrca.phylo doesn't work on a single tip and we need that

mrca.phylo.or.unique.tip <- function(tree, node) {
  if (length(node) == 1) {
    node
  } else {
    mrca.phylo(tree, node)
  }
}

get.node.height <- function(tree, node){
  return(max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[node])
}

build.igraph <- function(tt.info, sample){

  tree <- tt.info$tree
  annotations <- sample$annotations
  interventions <- sample$hidden

  if(is.null(interventions)){
    interventions <- rep(0, length(tree$tip.label) + tree$Nnode)
  }
  parents <- vector()

  current.host.count <- max(annotations)

  for(host in 1:max(annotations)){
    a.node <- which(annotations==host)[1]
    parent.node <- Ancestors(tree, a.node, type="parent")
    move.up <- parent.node != 0
    if(move.up){
      move.up <- annotations[a.node] == annotations[parent.node]
    }
    while(move.up){
      a.node <- parent.node
      parent.node <- Ancestors(tree, a.node, type="parent")
      move.up <- parent.node != 0
      if(move.up){
        move.up <- annotations[a.node] == annotations[parent.node]
      }
    }
    if(interventions[a.node] == 0){
      if(a.node == getRoot(tree)){
        parents[host] <- 0
      } else {
        parents[host] <- annotations[parent.node]
      }
    } else {
      current.host <- host
      for(intervening in 1:interventions[a.node]){
        current.host.count <- current.host.count + 1
        parents[current.host] <- current.host.count
        current.host <- parents[current.host]
      }
      if(getRoot(tree) == a.node){
        parents[current.host] <- 0
      } else {
        parents[current.host] <- annotations[parent.node]
      }
    }
  }

  node.names <- tt.info$hosts
  parent.names <- sapply(parents, function(x) if(x==0) "root" else tt.info$hosts[x])

  print(node.names)
  print(parent.names)

  edgelist <- cbind(parent.names, node.names)

  return(graph_from_edgelist(edgelist))

}
