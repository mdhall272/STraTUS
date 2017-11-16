
#' Enumerate transmission trees for the given pathogen phylogeny, and provide a uniform sample generator
#'
#' This function produces a list of class \code{tt.generator} which can be used to randomly sample transmission trees for the input phylogeny, and contains information on the number of compatible transmission trees. Note that at present, no combinations of variations (unsampled hosts, multiply-sampled hosts, and time constraints) are supported.
#' @param tree A \code{phylo} object
#' @param max.unsampled The maximum number of unsampled hosts in the transmission chain. The default is 0.
#' @param max.infection.to.sampling The greatest time period (in tree branch length units) that can have elapsed between the infection of a host and a tip from that host appearing. The default is infinity, meaning that no such time limit exists.
#' @param max.sampling.to.noninfectious The greatest time period (in tree branch length units) that can have elapsed between a tip from a host appearing and that host becoming noninfectious. If this is 0, a host's infection ends at the time of its last tip. The default is infinity, meaning that no such time limit exists.
#' @param minimum.heights A vector of the same length as the set of sampled hosts (at present this is always the number of tips of the tree) dictating the minimum height at which nodes can be allocated to each host. The order is the same as the order of tips in \code{tree$tip.label}. If absent, no such restrictions will be placed. Each must be equal to or smaller than the height of the last tip from the corresponding host. This overrides the given value of \code{max.sampling.to.noninfectious}.
#' @param maximum.heights A vector of the same length as the set of sampled hosts (at present this is always the number of tips of the tree) dictating the maximum height at which nodes can be allocated to each host. The order is the same as the order of tips in \code{tree$tip.label}. If absent, no such restrictions will be placed. Each must be equal to or greater than the height of the last tip from the corresponding host. This overrides the given value of \code{max.infection.to.sampling}.
#' @param tip.map A vector of the same length as the tip set of the tree listing a string giving the host from which the corresponding sample was derived. If absent, each tip is assumed to come from a different host and the tip names are taken to be the host names.
#' @return A list of class \code{tt.info} with the following fields:
#' \itemize{
#' \item{\code{tree}}{ The input tree}
#' \item{\code{tt.count}}{The total number of possible transmission trees.}
#' \item{\code{hosts}}{ The vector of host names. The order of the elements of this vector is used in the output of \code{sample.tt}.}
#' \item{\code{height.limits}}{ Height constraints only. A matrix giving maximum and minimum node heights, in two columns. Rows are ordered by the order of hosts given in the \code{host} field.}
#' \item{\code{bridge}}{ Multiple sampling only. A vector with the same length as the node set of the tree, dictating which nodes have their annotation forced by the tip annotations. Entries are host numbers for nodes whose annotation must be that host, and NA for nodes which can take multiple hosts.}
#' \item{\code{type}}{ The variation used; \code{height.aware} if height constraints were specified, \code{multisampled} if a \code{tip.map} was given, \code{unsampled} if \code{max.unsampled} is greater than zero, \code{basic} otherwise.}
#' \item{\code{node.calculations}}{ A list with the same length as the number of nodes of the tree and whose entries are indexed in the same order. If \code{max.unsampled} is 0, each has the following fields (the terminology here comes from the Hall paper):
#' \itemize{
#' \item{\code{p}}{ The number of valid partitions of the subtree rooted at this node.}
#' \item{\code{pstar}}{ The number of valid partitions of the unrooted tree obtained by attaching a single extra tip to the root node of the subtree rooted at this node. Alternatively, if any height constraints are given, a vector of the same length as the set of hosts, giving the number of partitions of the unrooted tree if the extra partition element is subject to the same minimum (but not maximum) height constraint as each host in turn.}
#' \item{\code{v}}{ A list indexed by the set of hosts, whose entries are the number of valid partitions of the subtree rooted at this node where the root node is in the partition element from each host.}
#' }
#' Alternatively, if \code{max.unsampled} is greater than 0, the entries are:
#' \itemize{
#' \item{\code{p}}{ A vector of length 1 + \code{max.unsampled} giving the number of valid partitions of the subtree rooted at this node if there are between 0 and \code{max.unsampled} (in order) partition elements containing no tips.}
#' \item{\code{pstar}}{ A vector of length 1 + \code{max.unsampled} giving the number of valid partitions of the tree obtained from the subtree rooted at this node by adding an extra tip connected to the root node, if there are between 0 and \code{max.unsampled} (in order) partition elements containing no tips.}
#' \item{\code{ps}}{ As with \code{p}, except this counts only partitions that have the root node in a sampled component (one containing at least one tip).}
#' \item{\code{pu}}{ As with \code{p}, except this counts only partitions that have the have the root node in an unsampled component (one containing no tip).}
#' \item{\code{v}}{ A list indexed by the set of hosts and "unsampled", whose entries are, for each host and an unsampled host, a vector of length 1 + \code{max.unsampled} counting the number of partitions that have the root node in that host's component if there are between 0 and \code{max.unsampled} partition elements containing no tips.}
#' }
#' }
#' }
#' @export tt.generator
#' @import phangorn


tt.generator <- function(tree,
                       max.unsampled = 0,
                       max.infection.to.sampling = Inf,
                       max.sampling.to.noninfectious = Inf,
                       minimum.heights = NULL,
                       maximum.heights = NULL,
                       tip.map = NULL){

  if(!ape::is.binary(tree)){
    stop("Binary trees only!\n")
  }

  if(max.infection.to.sampling < Inf & !is.null(maximum.heights)){
    max.infection.to.sampling <- Inf
    warning(paste0("Using the specified maximum heights vector; ignoring the specified value of max.infection.to.sampling (",max.infection.to.sampling,")"))
  }
  
  if(max.sampling.to.noninfectious < Inf & !is.null(minimum.heights)){
    max.sampling.to.noninfectious <- Inf
    warning(paste0("Using the specified minimum heights vector; ignoring the specified value of max.sampling.to.noninfectious (",max.sampling.to.noninfectious,")"))
  }
  

  has.heights <- max.infection.to.sampling<Inf | max.sampling.to.noninfectious<Inf | !is.null(minimum.heights) | !is.null(maximum.heights)
  has.unsampled <- max.unsampled > 0
  has.multisampled <- !is.null(tip.map)

  if(sum(has.heights, has.unsampled, has.multisampled) > 1){
    stop("This combination of variations is not currently supported.")
  }

  height.limits <- NULL
  bridge <- NULL

  if(has.heights){

    hosts <- tree$tip.label

    if(max.sampling.to.noninfectious < Inf){
      minimum.heights <- sapply(1:length(tree$tip.label), function(x) {
        get.node.height(tree, x) - max.sampling.to.noninfectious
      })
    }
    
    if(max.infection.to.sampling < Inf){
      maximum.heights <- sapply(1:length(tree$tip.label), function(x) {
        get.node.height(tree, x) + max.infection.to.sampling
      })
    }

    if(is.null(minimum.heights)){
      minimum.heights <- rep(-Inf, length(tree$tip.label))
    }

    if(is.null(maximum.heights)){
      maximum.heights <- rep(Inf, length(tree$tip.label))
    }

    height.limits <- cbind(minimum.heights, maximum.heights)

    results <- .height.aware.up.phase(tree, phangorn::getRoot(tree), list(), height.limits)

    type <- "height.aware"

  } else if(has.multisampled){

    hosts <- unique(tip.map)

    bridge <- tryCatch(
      .build.bridge(tree, hosts, tip.map),
      error=function(cond){
        stop("No transmission trees (without superinfection) are compatible with this set of tip labels")
      }
    )

    results <- .multiply.sampled.up.phase(tree, phangorn::getRoot(tree), list(), bridge)

    type <- "multisampled"

  } else {

    if(max.unsampled == 0){
      hosts <- tree$tip.label
    } else{
      hosts <- c(tree$tip.label, paste("uh", 1:max.unsampled, sep=""))
    }

    results <- .up.phase(tree, phangorn::getRoot(tree), list(), max.unsampled)

    # Temporarily, convert to the expected format for zero unsampled hosts. This will go in the upcoming version where all variations are possible at once.

    if(max.unsampled == 0){
      results <- lapply(results, function(x){
        x$v <- x$v[1:length(tree$tip.label),1]
        x$pu <- NULL
        x$ps <- NULL
        x
      })
      type <- "basic"
    } else {
      type <- "unsampled"
    }
  }
  
  tt.count = sum(results[[phangorn::getRoot(tree)]]$p)

  out <- list(tree = tree, hosts = hosts, tt.count = tt.count, height.limits = height.limits, bridge = bridge, type=type, node.calculations = results)

  class(out) <- append(class(out), "tt.generator")

  return(out)
}

.up.phase <- function(tree, node, node.calculations, max.unsampled){
  if(is.tip(tree, node)){

    node.info <- list()

    # rows are hosts except the last row is "unsampled". Columns are the number of unsampled partition elements in the tree, shifted by 1 because zero is a thing

    v <- matrix(0, nrow = length(tree$tip.label) + 1, ncol = max.unsampled + 1)
    v[node, 1] <- 1

    node.info$v <- v
    node.info$p <- c(1, rep(0, max.unsampled))
    node.info$pstar <- c(1, rep(0, max.unsampled))
    if(ncol(v)>1){
      node.info$ps <- colSums(v[1:length(tree$tip.label),])
    } else {
      node.info$ps <- sum(v)
    }
    node.info$pu <- v[length(tree$tip.label) + 1,]

    node.calculations[[node]] <- node.info

  } else {
    for(child in phangorn::Children(tree, node)){
      node.calculations <- .up.phase(tree, child, node.calculations, max.unsampled)
    }

    node.info <- list()

    v <- matrix(0, nrow = length(tree$tip.label) + 1, ncol = max.unsampled + 1)

    # the row for the unsampled state

    v[length(tree$tip.label) + 1, ] <- sapply(0:max.unsampled, function(x){
      if(x==0){
        return(0)
      } else {
        # how many unsampled elements

        # how many that are not the root unsampled element
        remaining.us.elements <- x - 1
        out <- 0
        for(i in 0:remaining.us.elements){
          j = remaining.us.elements - i
          kids <- phangorn::Children(tree,node)

          term <- node.calculations[[kids[1]]]$pstar[i+1] * node.calculations[[kids[2]]]$pstar[j+1]
          out <- out + term
        }
        return(out)
      }
    } )

    node.info$pu <- v[length(tree$tip.label) + 1,]

    # the rest of the rows
    kids <- phangorn::Children(tree,node)

    for(host in 1:length(tree$tip.label)){
      temp <- sapply(0:max.unsampled, function(x){
        out <- 0
        for(i in 0:x){
          j = x - i
          term <- node.calculations[[kids[1]]]$v[host,i+1] * node.calculations[[kids[2]]]$pstar[j+1] +
            node.calculations[[kids[2]]]$v[host,i+1] * node.calculations[[kids[1]]]$pstar[j+1]
          out <- out + term
        }
        return(out)
      })
      v[host,] <- temp
    }

    if(ncol(v)>1){
      node.info$ps <- colSums(v[1:length(tree$tip.label),])
    } else {
      node.info$ps <- sum(v)
    }

    node.info$p <- node.info$pu + node.info$ps

    node.info$pstar <- sapply(0:max.unsampled, function(x) {
      out <- node.info$p[x+1]
      for(i in 0:x){
        j = x - i

        term <- node.calculations[[kids[1]]]$pstar[i+1] * node.calculations[[kids[2]]]$pstar[j+1]
        out <- out + term
      }
      return(out)
    })

    node.info$v <- v
    node.calculations[[node]] <- node.info
  }
  return(node.calculations)
}

.height.aware.up.phase <- function(tree, node, node.calculations, height.limits){
  if(is.tip(tree, node)){

    if(get.node.height(tree, node) < height.limits[node,1] | get.node.height(tree, node) > height.limits[node,2]){
      stop("Bad node height limits for node ",node, sep="")
    }

    node.info <- list()

    v <- rep(0, length(tree$tip.label))
    v[node] <- 1

    node.info$v <- v
    node.info$p <- 1
    node.info$pstar <- rep(1, length(tree$tip.label))
    node.calculations[[node]] <- node.info

  } else {
    for(child in phangorn::Children(tree, node)){
      node.calculations <- .height.aware.up.phase(tree, child, node.calculations, height.limits)
    }

    node.info <- list()

    v <- rep(0, length(tree$tip.label))
    pstarprod <- rep(1, length(tree$tip.label))

    kids <- phangorn::Children(tree,node)

    for(kid in kids){
      desc.tips <- unlist(phangorn::Descendants(tree, kid, type="tips"))
      v[desc.tips] <- node.calculations[[kid]]$v[desc.tips]

      for(kid2 in kids){
        if(kid != kid2){
          v[desc.tips] <- v[desc.tips]*node.calculations[[kid2]]$pstar[desc.tips]
        }
      }
      pstarprod <- pstarprod*node.calculations[[kid]]$pstar
    }

    height <- get.node.height(tree, node)

    # replace v for those out of range with 0

    v[which(height.limits[,1] > height | height.limits[,2] < height )] <- 0

    node.info$v <- v
    node.info$p <- sum(v)
    node.info$pstar <- rep(node.info$p, length(tree$tip.label))
    tips.lower <- which(height.limits[,1] <= get.node.height(tree,node))

    node.info$pstar[tips.lower] <- node.info$pstar[tips.lower] + pstarprod[tips.lower]
    node.calculations[[node]] <- node.info
  }
  return(node.calculations)
}

.multiply.sampled.up.phase <- function(tree, node, node.calculations, bridge){
  if(is.tip(tree, node)){

    node.info <- list()

    v <- rep(0, length(unique(stats::na.omit(bridge))))
    v[bridge[node]] <- 1

    node.info$v <- v
    node.info$p <- 1
    node.info$pstar <- 1
    node.calculations[[node]] <- node.info

  } else {
    for(child in phangorn::Children(tree, node)){
      node.calculations <- .multiply.sampled.up.phase(tree, child, node.calculations, bridge)
    }

    node.info <- list()

    v <- rep(0, length(unique(stats::na.omit(bridge))))
    pstarprod <- 1

    kids <- phangorn::Children(tree,node)

    desc.tips.1 <- unlist(phangorn::Descendants(tree, kids[1], type="tips"))
    desc.tips.2 <- unlist(phangorn::Descendants(tree, kids[2], type="tips"))
    if(length(intersect(bridge[desc.tips.1], bridge[desc.tips.2]))==1){
      # this is the situation where this is the top of a bridge

      only.valid.host <- intersect(bridge[desc.tips.1], bridge[desc.tips.2])

      v[only.valid.host] <- node.calculations[[kids[1]]]$v[only.valid.host] * node.calculations[[kids[2]]]$v[only.valid.host]

    } else {
      for(host.no in 1:max(stats::na.omit(bridge))){
        if(!is.na(bridge[node]) & bridge[node] != host.no){
          v[host.no] <- 0
        } else if(host.no %in% bridge[desc.tips.1]){
          v[host.no] <- node.calculations[[kids[1]]]$v[host.no] * node.calculations[[kids[2]]]$pstar
        } else if(host.no %in% bridge[desc.tips.2]) {
          v[host.no] <- node.calculations[[kids[2]]]$v[host.no] * node.calculations[[kids[1]]]$pstar
        }
      }
    }

    if(is.na(bridge[node])){
      pstar <- sum(v) + node.calculations[[kids[1]]]$pstar * node.calculations[[kids[2]]]$pstar
    } else {
      pstar <- sum(v)
    }

    node.info$v <- v
    node.info$p <- sum(v)
    node.info$pstar <- pstar
    node.calculations[[node]] <- node.info
  }
  return(node.calculations)
}

.build.bridge <- function(tree, hosts, tip.map){
  bridge <- rep(NA, node.count(tree))

  for(host in hosts){
    tips <- which(tip.map == host)
    mrca <- mrca.phylo.or.unique.tip(tree, tips)
    for(tip in tips){
      current.node <- tip
      repeat{
        bridge[current.node] <- which(hosts==host)
        if(current.node == mrca){
          break
        }
        current.node <- phangorn::Ancestors(tree, current.node, type="parent")
        if(!is.na(bridge[current.node])){
          if(bridge[current.node] ==  which(hosts==host)){
            break
          } else {
            stop("Overlapping bridges; no transmission trees are compatible with this tree")
          }
        }
      }
    }
  }
  return(bridge)
}
