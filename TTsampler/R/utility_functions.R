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
    phangorn::mrca.phylo(tree, node)
  }
}

get.node.height <- function(tree, node){
  return(max(ape::node.depth.edgelength(tree)) - ape::node.depth.edgelength(tree)[node])
}

#' For a sample, produce the transmission tree as a \code{igraph} object
#'
#' @param generator A list of class \code{tt.generator} produced by \code{tt.generator}.
#' @param sample A list of class \code{tt.sample} produced by \code{sample.tt} or \code{sample.partial.tt}
#' @return An \code{igraph} object
#' @export build.igraph
#' @import igraph


build.igraph <- function(generator, sample){

  tree <- generator$tree
  annotations <- sample$annotations
  interventions <- sample$hidden

  if(is.null(interventions)){
    interventions <- rep(0, length(tree$tip.label) + tree$Nnode)
  }
  parents <- vector()

  current.host.count <- max(annotations)

  for(host in 1:max(annotations)){
    a.node <- which(annotations==host)[1]
    parent.node <- phangorn::Ancestors(tree, a.node, type="parent")
    move.up <- parent.node != 0
    if(move.up){
      move.up <- annotations[a.node] == annotations[parent.node]
    }
    while(move.up){
      a.node <- parent.node
      parent.node <- phangorn::Ancestors(tree, a.node, type="parent")
      move.up <- parent.node != 0
      if(move.up){
        move.up <- annotations[a.node] == annotations[parent.node]
      }
    }
    if(interventions[a.node] == 0){
      if(a.node == phangorn::getRoot(tree)){
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
      if(phangorn::getRoot(tree) == a.node){
        parents[current.host] <- 0
      } else {
        parents[current.host] <- annotations[parent.node]
      }
    }
  }

  node.names <- generator$hosts
  parent.names <- sapply(parents, function(x) if(x==0) "root" else generator$hosts[x])

  print(node.names)
  print(parent.names)

  edgelist <- cbind(parent.names, node.names)

  return(igraph::graph_from_edgelist(edgelist))

}

#' For a sample with no unsampled hosts, draw the annotated phylogeny using \code{ggtree}
#'
#' @param generator A list of class \code{tt.generator} produced by \code{tt.generator}.
#' @param sample A list of class \code{tt.sample} produced by \code{sample.tt} or \code{sample.partial.tt}
#' @return A \code{ggtree} object
#' @export draw.fully.sampled
#' @import ggtree
#' @import ggplot2

draw.fully.sampled <- function(generator, sample){

  if(generator$type=="unsampled"){
    stop("This generator uses unsampled hosts.")
  }

  tree <- generator$tree

  annots <- generator$hosts[sample$annotations]

  first.in.split <- sapply(1:node.count(tree), function(x) {
    if(phangorn::getRoot(tree) == x){
      return(T)
    }
    if(annots[x] == annots[phangorn::Ancestors(tree, x, type="parent")]){
      return(F)
    }
    return(T)
  })

  branch.colour.annots <- annots
  branch.colour.annots[which(first.in.split)] <- NA

  attr(tree, "annots") <- annots
  attr(tree, "branch.colour.annots") <- branch.colour.annots

  picture <- ggtree(tree, aes(col=branch.colour.annots), size=1.5) +
    geom_tiplab(aes(col=annots), hjust=-1) +
    geom_point(aes(col=annots), size=4) +
    scale_fill_hue(na.value = "grey") +
    scale_color_hue(na.value = "grey")

  picture
}

#' For a sample with unsampled hosts, draw the annotated phylogeny using \code{ggtree}
#'
#' @param generator A list of class \code{tt.generator} produced by \code{tt.generator}.
#' @param sample A list of class \code{tt.sample} produced by \code{sample.tt} or \code{sample.partial.tt}
#' @return A \code{ggtree} object
#' @export draw.incompletely.sampled
#' @import ggtree
#' @import ggplot2

draw.incompletely.sampled <- function(generator, sample){

  if(generator$type!="unsampled"){
    stop("This generator does not use unsampled hosts.")
  }


  node.annots <- generator$hosts[sample$annotations]
  branch.annots <- sample$hidden
  tree <- generator$tree

  first.in.split <- sapply(1:node.count(tree), function(x) {
    if(phangorn::getRoot(tree) == x){
      return(T)
    }
    if(node.annots[x] == node.annots[phangorn::Ancestors(tree, x, type="parent")]){
      return(F)
    }
    return(T)
  })

  branch.colour.annots <- node.annots
  branch.colour.annots[which(first.in.split)] <- NA

  branch.annots[which(!first.in.split)] <- NA

  adjustments <- rep(0.5, node.count(tree))
  adjustments[phangorn::getRoot(tree)] <- 1.5

  attr(tree, "node.annots") <- node.annots
  attr(tree, "branch.colour.annots") <- branch.colour.annots
  attr(tree, "branch.annots") <- branch.annots
  attr(tree, "adjustments") <- adjustments

  picture <- ggtree(tree, aes(col=branch.colour.annots), size=1.5) +
    geom_tiplab(aes(col=node.annots), hjust=-1) +
    geom_point(aes(col=node.annots), size=4) +
    scale_fill_hue(na.value = "grey") +
    scale_color_hue(na.value = "grey") +
    geom_text(aes(x=branch, label=branch.annots, hjust=adjustments), vjust=-0.5, col="black", na.rm=TRUE)

  return(picture)
}
