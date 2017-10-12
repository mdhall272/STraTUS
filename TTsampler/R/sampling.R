
#' Sample one or more transmission trees uniformly
#'
#' @param tt.info A list of class \code{tt.generator} produced by \code{tt.sampler}.
#' @param count How many transmission trees to sample.
#' @param unsampled The number of unsampled hosts in the transmission chain. A value >0 requires a \code{tt.info} list whose \code{type} is \code{unsampled}.
#' @param draw Use \code{ggtree} to draw a coloured phylogeny showing each transmission tree overlaid onto the phylogeny.
#' @param network Produce the transmission trees in \code{igraph} format.
#' @return A list, each of whose elements is a list of class \code{tt} with one or more of the following elements:
#' \itemize{
#' \item{\code{annotations}}{ Always present. A vector indicating which host (given by numbers corresponding to the ordering in \code{tt.info$hosts}) is assigned to each phylogeny node.}
#' \item{\code{hidden}}{ Present if \code{unsampled} is greater than 0. The number of "hidden" unsampled hosts (with no associated nodes) along each branch.}
#' \item{\code{picture}}{ Present if \code{draw} was specified; a \code{ggtree} object.}
#' \item{\code{igraph}}{ Present if \code{network} was specified; an \code{igraph} object.}
#' }


sample.tt <- function(tt.info, count = 1, unsampled = 0, draw = F, network = F){
  return(sample.partial.tt(tt.info, count, unsampled, getRoot(tt.info$tree),  NULL, F, draw, network))
}

#' Resample the subtree rooted at any tree node, keeping the annotations for the rest of the tree fixed
#'
#' @param tt.info A list of class \code{tt.generator} produced by \code{tt.sampler}.
#' @param count How many transmission trees to sample.
#' @param unsampled The number of unsampled hosts in the transmission chain. (The whole transmission chain, even if only part of the transmission tree is being resampled). A value >0 requires a \code{tt.info} list whose \code{type} is \code{unsampled}.
#' @param existing An existing list of class \code{tt}, representing a transmission tree to be modified. Usually these are produced by a \code{sample.tt} or \code{sample.partial.tt} call.
#' @param starting.node The root of the subtree to resample. If this is the root of the whole tree, then \code{existing} is irrelevent (but generally \code{sample.tt} should be used for this purpose).
#' @param check.integrity Whether to check if \code{existing} is indeed a valid transmission tree.
#' @param draw Use \code{ggtree} to draw a coloured phylogeny showing each transmission tree overload onto the phylogeny
#' @param network Produce the transmission trees in \code{igraph} format.
#' @return A list, each of whose elements is a list of class \code{tt} with one or more of the following elements:
#' \itemize{
#' \item{\code{annotations}}{ Always present. A vector indicating which host (given by numbers corresponding to the ordering in \code{tt.info$hosts}) is assigned to each phylogeny node.}
#' \item{\code{hidden}}{ Present if \code{unsampled} is greater than 0. The number of "hidden" unsampled hosts (with no associated nodes) along each branch.}
#' \item{\code{picture}}{ Present if \code{draw} was specified; a \code{ggtree} object.}
#' \item{\code{igraph}}{ Present if \code{network} was specified; an \code{igraph} object.}
#' }

sample.partial.tt <- function(tt.info, count = 1, unsampled = 0, starting.node = getRoot(tt.info$tree), existing=NULL, check.integrity = T, draw = F, network = F){
  if(!inherits(tt.info, "tt.generator")){
    stop("Input is not a list of class tt.generator")
  }

  if(unsampled > 0 & tt.info$type != "unsampled"){
    stop("This sampler will not generate trees with unsampled hosts.")
  }

  if(unsampled > 0){
    if ((ncol(tt.info$node.calculations[[1]]$v) - 1) < unsampled){
      stop("This sampler will not generate trees with this many unsampled hosts.")
    }
  }

  if(is.null(existing) & starting.node!=getRoot(tt.info$tree)){
    stop("An existing sample is required to resample the tree from any node other than the root.")
  }

  tree <- tt.info$tree

  existing.annot <- rep(0, node.count(tree))
  existing.hidden <- rep(0, node.count(tree))

  if(unsampled > 0 & !is.null(existing)){
    existing.annot <- existing$annotations
    existing.hidden <- existing$hidden
  } else {
    existing.annot <- existing$annotations
  }

  if(check.integrity){

    if(length(existing.annot)!=node.count(tree)){
      stop(paste0("Existing sample vector should have ",node.count(tree)," items."))
    }

    if(!is.numeric(existing.annot)){
      stop("Existing sample vector must be numerical")
    }

    if(max(existing.annot) > length(tt.info$hosts)){
      stop("Extra hosts amongst the annotations in the existing sample vector")
    }

    for(host in 1:(length(tt.info$hosts))){
      nodes <- which(existing.annot==host)

      ok <- T
      # a valid integrity check is whether the MRCA of any two nodes with the same annotation also has that annotation

      for(node.1 in nodes){
        for(node.2 in nodes){
          if(node.1!=node.2 & !(node.1 %in% Ancestors(tree, node.2)) & !(node.2 %in% Ancestors(tree, node.1))){
            mrca <- mrca.phylo(tree, c(node.1, node.2))
            if(existing.annot[mrca]!=host){
              ok <- F
              break
            }
          }
        }
        if(!ok){
          break
        }
      }

      if(!ok){
        break
      }
    }
    if(!ok){
      stop("Failed integrity check; specified existing transmission tree is not a valid partition ")
    }
  }

  if(tt.info$type=="basic"){
    results <- replicate(count, .basic.down.phase(tree, starting.node, existing.annot, tt.info$node.calculations))

    results <- lapply(seq_len(ncol(results)), function(x){
      item <- list()
      item$annotations <- results[,x]
      class(item) <- append(class(item), "tt")
      item
    })

    if(draw){
      results <- lapply(results, function(x){
        picture <- draw.fully.sampled(tt.info, x)
        x$picture <- picture
        x
      })
    }

    if(network){
      results <- lapply(results, function(x){
        x$igraph <- build.igraph(tt.info, x)
        x
      })
    }

    return(results)
  }
  if(tt.info$type=="multisampled"){
    if(is.null(tt.info$bridge)){
      stop("This tt.generator should have a bridge but it is not present.")
    }

    results <- replicate(count, .multiply.sampled.down.phase(tree, starting.node, existing.annot, tt.info$node.calculations, tt.info$bridge))

    results <- lapply(seq_len(ncol(results)), function(x){
      item <- list()
      item$annotations <- results[,x]
      class(item) <- append(class(item), "tt")
      item
    })

    if(draw){
      results <- lapply(results, function(x){
        x$picture <- draw.fully.sampled(tt.info, x)
        x
      })
    }

    if(network){
      results <- lapply(results, function(x){
        x$igraph <- build.igraph(tt.info, x)
        x
      })
    }

    return(results)
  }
  if(tt.info$type=="height.aware"){
    if(is.null(tt.info$height.limits)){
      stop("This tt.generator should have height.limits but they are not present.")
    }

    results <- replicate(count, .height.aware.down.phase(tree, starting.node, existing.annot, tt.info$node.calculations, tt.info$height.limits))

    results <- lapply(seq_len(ncol(results)), function(x){
      item <- list()
      item$annotations <- results[,x]
      class(item) <- append(class(item), "tt")
      item
    })

    if(draw){
      results <- lapply(results, function(x){
        x$picture <- draw.fully.sampled(tt.info, x)
        x
      })
    }

    if(network){
      results <- lapply(results, function(x){
        x$igraph <- build.igraph(tt.info, x)
        x
      })
    }

    return(results)
  }
  if(tt.info$type=="unsampled"){

    # uuuuurghhh...

    if(starting.node!=getRoot(tree)){

      unsampled.nos <- seq(length(tt.info$tree$tip.label) + 1, length(tt.info$tree$tip.label) + unsampled)

      subtree.nodes <- c(starting.node, unlist(Descendants(tree, starting.node, type="all")))
      other.nodes <- setdiff(1:node.count(tree), subtree.nodes)

      # We leave alone everything not involving the subtree rooted at starting.node.
      # Hidden hosts on the root branch of that ARE resampled (they have to be, because it may not end up an infection branch)
      # Existing unsampled hosts exist outside the subtree but can also exist inside it
      # Only visible hosts have numbers

      unsampled.nos.outside.subtree <- intersect(unsampled.nos, unique(existing.annot[other.nodes]))

      visible.existing.unsampled.hosts <- length(unsampled.nos.outside.subtree)
      hidden.existing.unsampled.hosts <- sum(existing.hidden[other.nodes])

      existing.unsampled.hosts <- visible.existing.unsampled.hosts + hidden.existing.unsampled.hosts

      remaining.unsampled.hosts <- unsampled - existing.unsampled.hosts

      if(remaining.unsampled.hosts < 0){
        stop(paste0("There are already more than ",unsampled, " unsampled hosts in the provided transmission tree"))
      }

      starting.current.host.count <- length(tree$tip.label) + visible.existing.unsampled.hosts

      # Unpleasant as this seems, it is probably easiest to renumber the existing hosts.

      existing.annot <- sapply(existing.annot, function(x){
        if(x %in% unsampled.nos.outside.subtree){
          length(tree$tip.label) + which(unsampled.nos.outside.subtree == x)
        } else {
          x
        }
      })
    } else {
      starting.current.host.count <- length(tree$tip.label)
      existing.unsampled.hosts <- 0
      subtree.nodes <- 1:node.count(tree)
      other.nodes <- vector()
      remaining.unsampled.hosts <- unsampled
    }

    root.forced <- F

    if(starting.node == getRoot(tree)){
      counts <- tt.info$node.calculations[[starting.node]]$p
      tip.hosts <- 1:length(tree$tip.label)
      tip.count <- length(tree$tip.label)

      visible.count.weights <- sapply(0:remaining.unsampled.hosts, function(x) counts[x+1]*choose(tip.count + remaining.unsampled.hosts - 1, tip.count + x - 1))
    } else {

      parent.host <- existing.annot[Ancestors(tree, starting.node, type="parent")]
      tip.hosts <- unlist(Descendants(tree, starting.node, type="tips"))
      tip.count <- length(tip.hosts)
      if(parent.host %in% tip.hosts){
        # Suppose that the parent of the starting node is assigned to a sampled host from within the subtree rooted at that node. Then our weights come from the
        # appropriate row of the v matrix, and there are the tip count of the subtree, minus 1, plus the column index, minus 1, infection branches in that subtree
        # which receive remaining.unsampled.hosts minus the column index

        root.forced <- T

        counts <- tt.info$node.calculations[[starting.node]]$v[parent.host,]
        visible.count.weights <- sapply(0:remaining.unsampled.hosts, function(x) counts[x+1]*choose(tip.count + remaining.unsampled.hosts - 2, tip.count + x - 2))

      } else {
        # If the partition element to which existing.node is not forced, then the parent is assigned either to an existing unsampled element or to the host on a
        # tip not from the subtree. The weights then come from pstar, and there are the tip count of the subtree, plus the column index, minus 1, infection branches
        # in the subtree which receive remaining.unsampled.hosts minus the column index

        counts <- tt.info$node.calculations[[starting.node]]$pstar
        visible.count.weights <- sapply(0:remaining.unsampled.hosts, function(x) counts[x+1]*choose(tip.count + remaining.unsampled.hosts - 1, tip.count + x - 1))
      }

    }

    annotations <- vector()
    hidden <- vector()

    results <- list()

    for(i in 1:count){
      out <- list()
      class(out) <- append(class(out), "tt")

      current.host.count <<- starting.current.host.count

      no.visible <- sample(0:remaining.unsampled.hosts, 1, prob=visible.count.weights)
      no.hidden <- remaining.unsampled.hosts - no.visible

      a.sample <- .unsampled.down.phase(tree, starting.node, existing.annot, tt.info$node.calculations, no.visible)

      out$annotations <- a.sample

      branch.us.position.choice <- vector()
      if(no.visible != remaining.unsampled.hosts){
        if(!root.forced){
          branch.us.position.options <- combinations(tip.count + no.visible, no.hidden, repeats.allowed = T )
        } else {
          branch.us.position.options <- combinations(tip.count + no.visible - 1, no.hidden, repeats.allowed = T )
        }
        branch.us.position.choice <- branch.us.position.options[sample(1:nrow(branch.us.position.options), 1),]
      }

      interventions <- existing.hidden

      # wipe clean the subtree

      interventions[subtree.nodes] <- 0

      if(starting.node == getRoot(tree)){
        need.new.ibs <- 1:(length(tree$tip.label) + no.visible)
      } else if(!root.forced){
        if(no.visible==0){
          need.new.ibs <- tip.hosts
        } else {
          need.new.ibs <- c(tip.hosts, starting.current.host.count + (1:no.visible))
        }
      } else {
        if(no.visible==0){
          need.new.ibs <- c(setdiff(tip.hosts, parent.host))
        } else {
          need.new.ibs <- c(setdiff(tip.hosts, parent.host), starting.current.host.count + (1:no.visible))
        }
      }

      interventions[which(a.sample %in% need.new.ibs)] <- 0

      for(host in 1:(length(tree$tip.label) + no.visible)){

        # find a node in this region
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
        # should now be at the root of this partition element
        interventions[a.node] <- sum(branch.us.position.choice == which(need.new.ibs==host))

      }
      out$hidden <- interventions

      if(draw){
        out$picture <- draw.incompletely.sampled(tt.info, out)
      }

      if(network){
        out$igraph <- build.igraph(tt.info, out)
      }

      results[[i]] <- out

    }

    return(results)
  }
  stop("Unknown tt.info$type")
}


.basic.down.phase <- function(tree, node, result.vector, info){
  if(node == getRoot(tree)){
    result.vector[node] <- sample(1:length(tree$tip.label), 1, prob=info[[node]]$v)
  } else {
    parent.choice <- result.vector[Ancestors(tree, node, type="parent")]

    desc.tips <- unlist(Descendants(tree, node))

    if(parent.choice %in% desc.tips){
      result.vector[node] <- parent.choice
    } else {
      host.weights <- info[[node]]$v
      host.weights[parent.choice] <- info[[node]]$pstar - info[[node]]$p

      result.vector[node] <- sample(1:length(tree$tip.label),1,prob=host.weights)
    }
  }
  for(child in Children(tree, node)){
    result.vector <- .basic.down.phase(tree, child, result.vector, info)
  }
  return(result.vector)
}

.unsampled.down.phase <- function(tree, node, result.vector, info, us.count){
  #  cat("At node",node,"with",us.count,"remaining unsampled regions\n")
  kids <- Children(tree, node)
  vmatrix <- info[[node]]$v
  if(node == getRoot(tree)){

    prob.weights <- vmatrix[,us.count+1]

    result <- sample(1:(length(tree$tip.label)+1),1, prob=prob.weights)

    if(result == length(tree$tip.label) + 1){
      result <- current.host.count + 1
      current.host.count <<- result
    }
    creep <- F

    result.vector[node] <- result

  } else {
    parent.choice <- result.vector[Ancestors(tree, node, type="parent")]

    desc.tips <- unlist(Descendants(tree, node))

    if(parent.choice %in% desc.tips){
      result <- parent.choice
      result.vector[node] <- result
      creep <- F
    } else {
      new.pat.weights <- vmatrix[,us.count+1]

      continuation.weight <- info[[node]]$pstar[us.count+1] - info[[node]]$p[us.count+1]

      result <- sample(1:(length(tree$tip.label)+2),1, prob=c(new.pat.weights, continuation.weight))

      creep <- F

      if(result == length(tree$tip.label) + 2){
        creep <- T
        #whatever the parent has
        result <- parent.choice
      } else if(result == length(tree$tip.label) + 1){
        #new unsampled host
        result <- current.host.count + 1
        current.host.count <<- result
      }

      result.vector[node] <- result
    }
  }
  if(!is.tip(tree, node)){
    if(us.count==0){
      result.vector <- .unsampled.down.phase (tree, kids[1], result.vector, info, 0)
      result.vector <- .unsampled.down.phase (tree, kids[2], result.vector, info, 0)
    } else {
      if(creep){
        # If this element (sampled or unsampled) isn't one associated with any descendant tip of this node
        i.weights <- sapply(0:us.count, function(i){
          j <- us.count - i
          return(info[[kids[1]]]$pstar[i+1]*info[[kids[2]]]$pstar[j+1])
        })

        if(sum(i.weights) != (info[[node]]$pstar[us.count + 1] - info[[node]]$p[us.count + 1])){
          cat("Encountered miscalculation 1.\n")
          quit("N")
        }

        chosen.i <- sample(0:us.count, 1, prob=i.weights)
        chosen.j <- us.count - chosen.i
      } else if(result > length(tree$tip.label)){
        # If this is a new unsampled node
        i.weights <- sapply(0:(us.count-1), function(i){
          j <- us.count - 1 - i
          return(info[[kids[1]]]$pstar[i+1]*info[[kids[2]]]$pstar[j+1])
        })

        if(sum(i.weights) != info[[node]]$pu[us.count + 1] ){
          cat("Encountered miscalculation 2.\n")
          quit("N")
        }

        chosen.i <- sample(0:(us.count-1), 1, prob=i.weights)
        chosen.j <- us.count - 1 - chosen.i
      } else {
        # Otherwise
        i.weights <- sapply(0:us.count, function(i){
          j <- us.count - i
          if(result %in% unlist(Descendants(tree, kids[1], type="tips"))){
            possibilities <- info[[kids[1]]]$v[result, i+1]*info[[kids[2]]]$pstar[j+1]
          } else {
            possibilities <- info[[kids[2]]]$v[result, j+1]*info[[kids[1]]]$pstar[i+1]
          }
          return(possibilities)
        })

        if(sum(i.weights) != info[[node]]$v[result, us.count + 1])  {
          cat("Encountered miscalculation 3.\n")
          quit("N")
        }

        chosen.i <- sample(0:us.count, 1, prob=i.weights)
        chosen.j <- us.count - chosen.i
      }
      # cat(chosen.i,"go to the left tree,",kids[1],"\n")
      # cat(chosen.j,"go to the right tree,",kids[2],"\n")
      result.vector <- .unsampled.down.phase(tree, kids[1], result.vector, info, chosen.i)
      result.vector <- .unsampled.down.phase(tree, kids[2], result.vector, info, chosen.j)
    }
  }

  return(result.vector)
}


.height.aware.down.phase <- function(tree, node, result.vector, info, height.limits){
  if(node == getRoot(tree)){
    result.vector[node] <- sample(1:length(tree$tip.label), 1, prob=info[[node]]$v)
  } else {
    parent.choice <- result.vector[Ancestors(tree, node, type="parent")]

    desc.tips <- unlist(Descendants(tree, node))

    if(parent.choice %in% desc.tips){
      result.vector[node] <- parent.choice
    } else {
      host.weights <- info[[node]]$v
      if(height.limits[parent.choice,1] <= get.node.height(tree, node)){
        host.weights[parent.choice] <- info[[node]]$pstar[parent.choice] - info[[node]]$p
      }
      result.vector[node] <- sample(1:length(tree$tip.label),1,prob=host.weights)
    }
  }
  for(child in Children(tree, node)){
    result.vector <- .height.aware.down.phase(tree, child, result.vector, info, height.limits)
  }
  return(result.vector)
}


.multiply.sampled.down.phase <- function(tree, node, result.vector, info, bridge){
  if(node == getRoot(tree)){
    result.vector[node] <- sample(1:length(unique(na.omit(bridge))), 1, prob=info[[node]]$v)
  } else {

    parent.choice <- result.vector[Ancestors(tree, node, type="parent")]
    desc.tips <- unlist(Descendants(tree, node))

    if(parent.choice %in% bridge[desc.tips]){
      result.vector[node] <- parent.choice
    } else {
      host.weights <- info[[node]]$v

      host.weights[parent.choice] <- info[[node]]$pstar - info[[node]]$p

      result.vector[node] <- sample(1:length(unique(na.omit(bridge))),1,prob=host.weights)
    }
  }
  for(child in Children(tree, node)){
    result.vector <- .multiply.sampled.down.phase(tree, child, result.vector, info, bridge)
  }
  return(result.vector)
}

draw.fully.sampled <- function(tt.info, sample){

  tree <- tt.info$tree

  annots <- tt.info$hosts[sample$annotations]

  first.in.split <- sapply(1:node.count(tree), function(x) {
    if(getRoot(tree) == x){
      return(T)
    }
    if(annots[x] == annots[Ancestors(tree, x, type="parent")]){
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

draw.incompletely.sampled <- function(tt.info, sample){
  node.annots <- tt.info$hosts[sample$annotations]
  branch.annots <- sample$hidden
  tree <- tt.info$tree

  first.in.split <- sapply(1:node.count(tree), function(x) {
    if(getRoot(tree) == x){
      return(T)
    }
    if(node.annots[x] == node.annots[Ancestors(tree, x, type="parent")]){
      return(F)
    }
    return(T)
  })

  branch.colour.annots <- node.annots
  branch.colour.annots[which(first.in.split)] <- NA

  branch.annots[which(!first.in.split)] <- NA

  adjustments <- rep(0.5, node.count(tree))
  adjustments[getRoot(tree)] <- 1.5

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
