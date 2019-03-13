
#' Sample one or more transmission trees uniformly
#'
#' @param generator A list of class \code{tt.generator} produced by \code{tt.generator}.
#' @param count How many transmission trees to sample.
#' @param unsampled The number of unsampled hosts in the transmission chain.
#' @param draw Use \code{ggtree} to draw a coloured phylogeny showing each transmission tree overlaid onto the phylogeny.
#' @param igraph Produce the transmission trees in \code{igraph} format.
#' @param verbose Verbose output
#' @return A list, each of whose elements is a list of class \code{tt} with one or more of the following elements:
#' \itemize{
#' \item{\code{annotations}}{ Always present. A vector indicating which host (given by numbers corresponding to the ordering in \code{generator$hosts}) is assigned to each phylogeny node.}
#' \item{\code{edgelist}}{ Always present. A \code{data.frame} giving the edge list; the first column are parents and the second children.}
#' \item{\code{hidden}}{ Present if \code{unsampled} is greater than 0. The number of "hidden" unsampled hosts (with no associated nodes) along each branch.}
#' \item{\code{picture}}{ Present if \code{draw} was TRUE; a \code{ggtree} object.}
#' \item{\code{igraph}}{ Present if \code{igraph} was TRUE; an \code{igraph} object.}
#' }
#' @export sample.tt
#' @import ggtree phangorn gmp RcppAlgos
#' @importFrom igraph graph_from_edgelist
#' @examples
#' # draw one sample from the uniform distribution
#' generator <- tt.generator(stratus.example.tree)
#' samples <- sample.tt(generator, 1, draw = TRUE)
#' samples[[1]]
#' # with unsampled.hosts
#' generator.us <- tt.generator(stratus.example.tree, max.unsampled = 2)
#' # note that you can ask for less unsampled hosts than the generator has (but not more)
#' samples.1us <- sample.tt(generator.us, 1, unsampled = 1,  draw = TRUE)
#' samples.1us[[1]]
#' # with multiply sampled hosts
#' generator.ms <- tt.generator(stratus.example.tree, tip.map = grouping.map)
#' samples.ms <- sample.tt(generator.ms, 1, draw = TRUE)


sample.tt <- function(generator, count = 1, unsampled = 0, draw = count==1, igraph = FALSE, verbose = FALSE){
  return(sample.partial.tt(generator, count, unsampled, phangorn::getRoot(generator$tree),  NULL, FALSE, draw, igraph, verbose))
}

#' Resample the subtree rooted at any tree node, keeping the annotations for the rest of the tree fixed
#'
#' @param generator A list of class \code{tt.generator} produced by \code{tt.generator}.
#' @param count How many transmission trees to sample.
#' @param unsampled The number of unsampled hosts in the transmission chain. (The whole transmission chain, even if only part of the transmission tree is being resampled). A value >0 requires a \code{generator} list whose \code{type} is \code{unsampled}.
#' @param existing An existing list of class \code{tt}, representing a transmission tree to be modified. Usually these are produced by a \code{sample.tt} or \code{sample.partial.tt} call.
#' @param starting.node The root of the subtree to resample. If this is the root of the whole tree, then \code{existing} is irrelevent (but generally \code{sample.tt} should be used for this purpose).
#' @param check.integrity Whether to check if \code{existing} is indeed a valid transmission tree.
#' @param draw Use \code{ggtree} to draw a coloured phylogeny showing each transmission tree overload onto the phylogeny
#' @param igraph Produce the transmission trees in \code{igraph} format.
#' @param verbose Verbose output
#' @return A list, each of whose elements is a list of class \code{tt} with one or more of the following elements:
#' \itemize{
#' \item{\code{annotations}}{ Always present. A vector indicating which host (given by numbers corresponding to the ordering in \code{generator$hosts}) is assigned to each phylogeny node.}
#' \item{\code{edgelist}}{ Always present. A \code{data.frame} giving the edge list; the first column are parents and the second children.}
#' \item{\code{hidden}}{ Present if \code{unsampled} is greater than 0. The number of "hidden" unsampled hosts (with no associated nodes) along each branch.}
#' \item{\code{picture}}{ Present if \code{draw} was TRUE; a \code{ggtree} object.}
#' \item{\code{igraph}}{ Present if \code{igraph} was TRUE; an \code{igraph} object.}
#' }
#' @export sample.partial.tt
#' @import ggtree phangorn 
#' @importFrom igraph graph_from_edgelist
#' @importFrom stats runif
#' @examples
#' # draw one sample from the uniform distribution
#' generator <- tt.generator(stratus.example.tree)
#' samples <- sample.tt(generator, 1, draw = TRUE)
#' original.tt <- samples[[1]]
#' # sample anew, from node 31 downwards
#' revised.tt <- sample.partial.tt(generator, 1, starting.node = 31, existing = original.tt, draw = TRUE)

sample.partial.tt <- function(generator, 
                              count = 1, 
                              unsampled = 0, 
                              starting.node = phangorn::getRoot(generator$tree), 
                              existing = NULL, 
                              check.integrity = TRUE, 
                              draw = count==1, 
                              igraph = FALSE, 
                              verbose = FALSE){
  
  
  if(!inherits(generator, "tt.generator")){
    stop("Input is not a list of class tt.generator")
  }
  
  if(generator$tt.count == 0){
    stop("There are no valid transmission trees produced by this generator.")
  }
  
  if(unsampled > 0){
    if ((ncol(generator$node.calculations[[1]]$v) - 1) < unsampled){
      stop("This sampler will not generate trees with this many unsampled hosts.")
    }
  }
  
  if(is.null(existing) & starting.node!=phangorn::getRoot(generator$tree)){
    stop("An existing sample is required to resample the tree from any node other than the root.")
  }
  
  tree <- generator$tree
  
  existing.annot <- rep(0, node.count(tree))
  existing.hidden <- rep(0, node.count(tree))
  
  if(!is.null(existing)){
    existing.annot <- existing$annotations
    if(unsampled > 0){
      existing.hidden <- existing$hidden
    }
  }
  
  if(check.integrity){
    
    if(length(existing.annot)!=node.count(tree)){
      stop(paste0("Existing sample vector should have ",node.count(tree)," items."))
    }
    
    if(!is.numeric(existing.annot)){
      stop("Existing sample vector must be numerical")
    }
    
    if(max(existing.annot) > length(generator$hosts)){
      stop("Extra hosts amongst the annotations in the existing sample vector")
    }
    
    for(host in 1:(length(generator$hosts))){
      nodes <- which(existing.annot==host)
      
      ok <- TRUE
      # a valid integrity check is whether the MRCA of any two nodes with the same annotation also has that annotation
      
      for(node.1 in nodes){
        for(node.2 in nodes){
          if(node.1!=node.2 & !(node.1 %in% phangorn::Ancestors(tree, node.2)) & !(node.2 %in% phangorn::Ancestors(tree, node.1))){
            mrca <- phangorn::mrca.phylo(tree, c(node.1, node.2))
            if(existing.annot[mrca]!=host){
              ok <- FALSE
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
  
  sampled.host.count <- length(unique(stats::na.omit(generator$bridge)))
  
  # double uuuuurghhh...
  
  if(starting.node!=phangorn::getRoot(tree)){
    
    subtree.nodes <- c(starting.node, unlist(phangorn::Descendants(tree, starting.node, type="all")))
    other.nodes <- setdiff(1:node.count(tree), subtree.nodes)
    
    unsampled.nos.outside.subtree <- vector()
    
    # these are all the unsampled hosts in the entire tree
    if(unsampled > 0){
      unsampled.nos <- seq(sampled.host.count + 1, sampled.host.count + unsampled)
      
      # We leave alone everything not involving the subtree rooted at starting.node.
      # Hidden hosts on the root branch of that ARE resampled (they have to be, because it may not end up an infection branch)
      # Existing unsampled hosts exist outside the subtree but can also exist inside it (if they creep down)
      # Only visible hosts have numbers
      
      unsampled.nos.outside.subtree <- intersect(unsampled.nos, unique(existing.annot[other.nodes]))
      
      visible.existing.unsampled.hosts <- length(unsampled.nos.outside.subtree)
      hidden.existing.unsampled.hosts <- sum(existing.hidden[other.nodes])
      
      existing.unsampled.hosts <- visible.existing.unsampled.hosts + hidden.existing.unsampled.hosts
      
      remaining.unsampled.hosts <- unsampled - existing.unsampled.hosts
      
      if(remaining.unsampled.hosts < 0){
        stop(paste0("There are already more than ", unsampled, " unsampled hosts in the provided transmission tree"))
      }
    } else {
      visible.existing.unsampled.hosts <- 0
      hidden.existing.unsampled.hosts <- 0
      existing.unsampled.hosts <- 0
      remaining.unsampled.hosts <- 0
    }
    
    # this is just for numbering
    
    starting.current.host.count <- sampled.host.count + visible.existing.unsampled.hosts
    
    # Unpleasant as this seems, it is probably easiest to renumber the existing hosts.
    
    existing.annot <- lapply(existing.annot, function(x){
      if(unsampled > 0 & x %in% unsampled.nos.outside.subtree){
        sampled.host.count + which(unsampled.nos.outside.subtree == x)
      } else {
        x
      }
    })
    
    existing.annot <- do.call(c, existing.annot)
    
  } else {
    starting.current.host.count <- sampled.host.count
    existing.unsampled.hosts <- 0
    subtree.nodes <- 1:node.count(tree)
    other.nodes <- vector()
    remaining.unsampled.hosts <- unsampled
  }
  
  root.forced <- FALSE
  
  if(starting.node == phangorn::getRoot(tree)){
    
    counts <- generator$node.calculations[[starting.node]]$p
    tip.hosts <- generator$bridge[1:length(tree$tip.label)]
    
    visible.count.weights <- lapply(0:remaining.unsampled.hosts, function(x) counts[x+1]*choose(sampled.host.count + remaining.unsampled.hosts - 1, sampled.host.count + x - 1))
    visible.count.weights <- do.call(c, visible.count.weights)
    
    subtree.sampled.host.count <- sampled.host.count
    
  } else {
    
    parent.host <- existing.annot[phangorn::Ancestors(tree, starting.node, type="parent")]
    tip.hosts <- generator$bridge[unlist(phangorn::Descendants(tree, starting.node, type="tips"))]
    
    subtree.sampled.host.count <- length(unique(tip.hosts))
    
    if(parent.host %in% tip.hosts){
      # Suppose that the parent of the starting node is assigned to a sampled host from within the subtree rooted at that node. 
      # Then our weights come from the appropriate row of the v matrix, and there are the sampled host count of the subtree, minus 
      # 1, plus the column index, minus 1, infection branches in that subtree which receive remaining.unsampled.hosts minus the 
      # column index.
      
      root.forced <- TRUE
      
      counts <- generator$node.calculations[[starting.node]]$v[parent.host,]
      visible.count.weights <- lapply(0:remaining.unsampled.hosts, function(x) counts[x+1]*choose(subtree.sampled.host.count + remaining.unsampled.hosts - 2, subtree.sampled.host.count + x - 2))
      visible.count.weights <- do.call(c, visible.count.weights)
      
    } else {
      # If the partition element to which existing.node belongs is not forced (so it is not a bridge node), then the parent is 
      # assigned either to an existing unsampled element or to the host on a tip not from the subtree. The weights then come 
      # from pstar, and there are the sampled host count of the subtree, plus the column index, minus 1, infection branches in 
      # the subtree which receive remaining.unsampled.hosts minus the column index.
      
      counts <- generator$node.calculations[[starting.node]]$pstar
      visible.count.weights <- lapply(0:remaining.unsampled.hosts, function(x) counts[parent.host,x+1]*choose(subtree.sampled.host.count + remaining.unsampled.hosts - 1, subtree.sampled.host.count + x - 1))
      visible.count.weights <- do.call(c, visible.count.weights)
    }
    
  }
  
  if(all(visible.count.weights==0)){
    stop("No valid transmission trees for this configuration (root must be unsampled given these height limits?)")
  }
  
  annotations <- vector()
  hidden <- vector()
  
  results <- list()
  
  for(i in 1:count){
    if(verbose) cat("Sample ",i,"\n", sep="")
    out <- list()
    class(out) <- append(class(out), "tt")
    
    current.host.count <- starting.current.host.count
    
    no.visible <- sample(0:remaining.unsampled.hosts, 1, prob = as.numeric(visible.count.weights))
    
    no.hidden <- remaining.unsampled.hosts - no.visible
    
    results <- .unified.down.phase(tree, starting.node, current.host.count, existing.annot, generator$node.calculations, no.visible, generator$height.limits, generator$bridge, verbose)
    
    a.sample <- results$result.vector
    
    out$annotations <- a.sample
    
    branch.us.position.choice <- vector()
    if(no.visible != remaining.unsampled.hosts){
      if(!root.forced){
        branch.us.position.choice <- tryCatch({comboSample(subtree.sampled.host.count + no.visible, no.hidden, repetition = TRUE, n=1)},
                                              error = function(e) stop("Integer too large; too many combinations"))
      } else {
        branch.us.position.choice <- tryCatch({comboSample(subtree.sampled.host.count + no.visible - 1, no.hidden, repetition = TRUE, n=1)},
                                              error = function(e) stop("Integer too large; too many combinations"))
      }
    }
    
    branch.us.position.choice <- c(branch.us.position.choice)
    
    interventions <- existing.hidden
    
    # wipe clean the subtree
    
    interventions[subtree.nodes] <- 0
    
    if(starting.node == phangorn::getRoot(tree)){
      need.new.ibs <- 1:(sampled.host.count + no.visible)
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
    
    for(host in need.new.ibs){
      
      # find a node in this region
      a.node <- which(a.sample==host)[1]
      parent.node <- phangorn::Ancestors(tree, a.node, type="parent")
      move.up <- parent.node != 0
      if(move.up){
        move.up <- a.sample[a.node] == a.sample[parent.node]
      }
      while(move.up){
        a.node <- parent.node
        parent.node <- phangorn::Ancestors(tree, a.node, type="parent")
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
      out$picture <- draw.incompletely.sampled(generator, out)
    }
    
    out$edgelist <- build.edgelist(generator, out)
    
    if(igraph){
      out$igraph <- graph_from_edgelist(out$edgelist)
    }
    
    results[[i]] <- out
    
  }
  return(results)
}



.unified.down.phase <- function(tree, node, current.host.count, result.vector, info, us.count, height.limits, bridge, verbose = FALSE){
  
  if(verbose) cat("At node",node,"with",us.count,"remaining unsampled regions\n")
  kids <- phangorn::Children(tree, node)
  vmatrix <- info[[node]]$v
  sampled.host.count <- length(unique(stats::na.omit(bridge)))
  
  if(node ==  phangorn::getRoot(tree)){
    
    prob.weights <- vmatrix[,us.count+1]
    
    if(all(prob.weights==0)){
      stop("No valid transmission trees for this configuration")
    }
    
    if(is.bigz(prob.weights)){
      prob.weights <- c(prob.weights)
      result <- bigz.sample(1:(sampled.host.count + 1), prob = prob.weights)
    } else {
      result <- sample(1:(sampled.host.count + 1), 1, prob = prob.weights)
    }
    
  
    if(result == (sampled.host.count + 1)){
      result <- current.host.count + 1
      current.host.count <- result
    }
    creep <- FALSE
    
    result.vector[node] <- result
    
  } else {
    parent.choice <- result.vector[phangorn::Ancestors(tree, node, type="parent")]
    
    desc.tips <- unlist(phangorn::Descendants(tree, node))
    
    if(parent.choice %in% bridge[desc.tips]){
      result <- parent.choice
      result.vector[node] <- result
      creep <- FALSE
    } else {
      
      parent.unsampled <- parent.choice > sampled.host.count
      
      new.pat.weights <- vmatrix[,us.count+1]
      
      if(!parent.unsampled){
        if(height.limits[parent.choice,1] <= get.node.height(tree, node)){
          continuation.weight <- info[[node]]$pstar[parent.choice, us.count+1] - info[[node]]$p[us.count+1]
        } else {
          # I _think_ this is redundant
          continuation.weight <- 0
        }
      } else {
        continuation.weight <- info[[node]]$pstar[sampled.host.count + 1, us.count+1] - info[[node]]$p[us.count+1]
      }
      
      if(all(c(new.pat.weights, continuation.weight)==0)){
        stop("No valid transmission trees for this configuration")
      }
      
      prob.weights <- c(new.pat.weights, continuation.weight)
      
      if(is.bigz(prob.weights)){
        prob.weights <- c(prob.weights)
        result <- bigz.sample(1:(sampled.host.count + 2), prob = prob.weights)
      } else {
        result <- sample(1:(sampled.host.count+2), 1, prob=prob.weights)
      }
      
      
      creep <- FALSE
      
      if(result == sampled.host.count + 2){
        creep <- TRUE
        #whatever the parent has
        result <- parent.choice
      } else if(result == sampled.host.count + 1){
        #new unsampled host
        result <- current.host.count + 1
        current.host.count <- result
      }
      
      result.vector[node] <- result
    }
    
  }
  
  if(verbose) cat("Assigned host is", result, "\n")
  
  if(!is.tip(tree, node)){
    if(us.count==0){
      for(child in kids){
        results <- .unified.down.phase(tree, child, current.host.count, result.vector, info, 0, height.limits, bridge, verbose)
        
        result.vector <- results$result.vector
        current.host.count <- results$current.host.count
      }
    } else {
      if(creep){
        # If this element (sampled or unsampled) isn't one associated with any descendant tip of this node
        
        what.comes.down <- result
        if(what.comes.down > sampled.host.count){
          what.comes.down <- sampled.host.count + 1
        }
        
        distribution.of.us <- divide.k.into.n(us.count, length(kids))
        
        column.weights <- lapply(1:ncol(distribution.of.us), function(i){
          temp <- lapply(1:nrow(distribution.of.us), function(j){
            info[[kids[j]]]$pstar[what.comes.down, (distribution.of.us[j,i]+1)]
          })
          
          temp <- do.call(c, temp)
          
          return(prod(temp))
        })
        column.weights <- do.call(c, column.weights)
        
        if(sum(column.weights) != (info[[node]]$pstar[what.comes.down, us.count + 1] - info[[node]]$p[us.count + 1])){
          warning(paste0("Encountered miscalculation 1; error is ", 
                         (abs(sum(column.weights) - (info[[node]]$pstar[what.comes.down, us.count + 1] - info[[node]]$p[us.count + 1])))/
                           (info[[node]]$pstar[what.comes.down, us.count + 1] - info[[node]]$p[us.count + 1]), ". This probably means that R has encountered an
                         integer too large for it to reliably handle. This is unlikely to significantly affect results if the error is small. To avoid it entirely 
                         at the cost of considerable slowdown, use the --bigz option to tt.generator.\n"))
        }
        
        chosen.col <- distribution.of.us[,sample(1:ncol(distribution.of.us), 1, prob=as.numeric(column.weights))]
        
      } else if(result > sampled.host.count){
        # If this is a new unsampled node
        distribution.of.us <- divide.k.into.n(us.count-1, length(kids))
        
        column.weights <- lapply(1:ncol(distribution.of.us), function(i){
          temp <- lapply(1:nrow(distribution.of.us), function(j){
            info[[kids[j]]]$pstar[sampled.host.count + 1, (distribution.of.us[j,i]+1)]
          })
          
          temp <- do.call(c, temp)
          
          return(prod(temp))
        })
        
        column.weights <- do.call(c, column.weights)
        
        if(sum(column.weights) != info[[node]]$pu[us.count + 1] ){
          warning(paste0("Encountered miscalculation 2; error is ", 
                         (abs(sum(column.weights) - info[[node]]$pu[us.count + 1]))/info[[node]]$pu[us.count + 1], ". This probably means that R has encountered an
                         integer too large for it to reliably handle. This is unlikely to significantly affect results if the error is small. To avoid it entirely 
                         at the cost of considerable slowdown, use the --bigz option to tt.generator.\n"))
        }
        
        chosen.col <- distribution.of.us[,sample(1:ncol(distribution.of.us), 1, prob=as.numeric(column.weights))]
        
      } else {
        # Otherwise
        
        distribution.of.us <- divide.k.into.n(us.count, length(kids))
        
        column.weights <- lapply(1:ncol(distribution.of.us), function(i){
          temp <- lapply(1:nrow(distribution.of.us), function(j){
            if(result %in% bridge[unlist(Descendants(tree, kids[j], type="tips"))]){
              info[[kids[j]]]$v[result, distribution.of.us[j,i]+1]
            } else {
              info[[kids[j]]]$pstar[result, distribution.of.us[j,i]+1]
            }
          })
          
          temp <- do.call(c, temp)
          
          return(prod(temp))
        })
        
        column.weights <- do.call(c, column.weights)
        
        if(sum(column.weights) != info[[node]]$v[result, us.count + 1])  {
          warning(paste0("Encountered miscalculation 3; error is ", 
                         (abs(sum(column.weights) - info[[node]]$v[result, us.count + 1]))/info[[node]]$v[result, us.count + 1], ". This probably means that R has encountered an
                         integer too large for it to reliably handle. This is unlikely to significantly affect results if the error is small. To avoid it entirely 
                         at the cost of considerable slowdown, use the --bigz option to tt.generator.\n"))
        }
        
        
        chosen.col <- distribution.of.us[,sample(1:ncol(distribution.of.us), 1, prob=as.numeric(column.weights))]
        
      }
      
      for(cn in 1:length(chosen.col)){
        if(verbose) cat(chosen.col[cn], "unsampled regions go to the subtree rooted at",kids[cn],"\n")
      }
      
      for(child.no in 1:length(kids)){
        results <- .unified.down.phase(tree, kids[child.no], current.host.count, result.vector, info, chosen.col[child.no], height.limits, bridge, verbose)
        
        result.vector <- results$result.vector
        current.host.count <- results$current.host.count
      }
    }
  }
  return(list(result.vector = result.vector, current.host.count = current.host.count))
}

