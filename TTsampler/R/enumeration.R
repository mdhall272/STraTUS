
#' Enumerate transmission trees for the given pathogen phylogeny, and provide a uniform sample generator
#'
#' This function produces a list of class \code{tt.generator} which can be used to randomly sample transmission trees for the input phylogeny, and contains information on the number of compatible transmission trees.
#' @param tree A \code{phylo} object
#' @param max.unsampled The maximum number of unsampled hosts in the transmission chain. The default is 0.
#' @param max.infection.to.sampling The greatest time period (in tree branch length units) that can have elapsed between the infection of a host and a tip from that host appearing. The default is infinity, meaning that no such time limit exists.
#' @param max.sampling.to.noninfectious The greatest time period (in tree branch length units) that can have elapsed between a tip from a host appearing and that host becoming noninfectious. If this is 0, a host's infection ends at the time of its last tip. The default is infinity, meaning that no such time limit exists.
#' @param minimum.heights A vector of the same length as the set of sampled hosts (at present this is always the number of tips of the tree) dictating the minimum height at which nodes can be allocated to each host. The order is the same as the order of tips in \code{tree$tip.label}. If absent, no such restrictions will be placed. Each must be equal to or smaller than the height of the last tip from the corresponding host. This overrides the given value of \code{max.sampling.to.noninfectious}.
#' @param maximum.heights A vector of the same length as the set of sampled hosts (at present this is always the number of tips of the tree) dictating the maximum height at which nodes can be allocated to each host. The order is the same as the order of tips in \code{tree$tip.label}. If absent, no such restrictions will be placed. Each must be equal to or greater than the height of the last tip from the corresponding host. This overrides the given value of \code{max.infection.to.sampling}.
#' @param tip.map A vector of the same length as the tip set of the tree listing a string giving the host from which the corresponding sample was derived. If absent, each tip is assumed to come from a different host and the tip names are taken to be the host names.
#' @param bigz Use \code{bigz} from \code{gmp} for integers, recommended for large trees
#' @return A list of class \code{tt.info} with the following fields:
#' \itemize{
#' \item{\code{tree}}{ The input tree}
#' \item{\code{tt.count}}{The total number of possible transmission trees.}
#' \item{\code{hosts}}{ The vector of host names. The order of the elements of this vector is used in the output of \code{sample.tt}.}
#' \item{\code{height.limits}}{A matrix giving maximum and minimum node heights, in two columns. Rows are ordered by the order of hosts given in the \code{host} field.}
#' \item{\code{bridge}}{A vector with the same length as the node set of the tree, dictating which nodes have their annotation forced by the tip annotations. Entries are host numbers for nodes whose annotation must be that host, and NA for nodes which can take multiple hosts.}
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
#' @import phangorn gmp


tt.generator <- function(tree,
                         max.unsampled = 0,
                         max.infection.to.sampling = Inf,
                         max.sampling.to.noninfectious = Inf,
                         minimum.heights = NULL,
                         maximum.heights = NULL,
                         tip.map = 1:length(tree$tip.label),
                         bigz = F){
  
  if(max.infection.to.sampling < Inf & !is.null(maximum.heights)){
    max.infection.to.sampling <- Inf
    warning(paste0("Using the specified maximum heights vector; ignoring the specified value of max.infection.to.sampling (",max.infection.to.sampling,")"))
  }
  
  if(max.sampling.to.noninfectious < Inf & !is.null(minimum.heights)){
    max.sampling.to.noninfectious <- Inf
    warning(paste0("Using the specified minimum heights vector; ignoring the specified value of max.sampling.to.noninfectious (",max.sampling.to.noninfectious,")"))
  }
  
  host.nos <- unique(tip.map)
  
  has.heights <- max.infection.to.sampling<Inf | max.sampling.to.noninfectious<Inf | !is.null(minimum.heights) | !is.null(maximum.heights)
  
  if(has.heights){
    
    if(max.sampling.to.noninfectious < Inf){
      minimum.heights <- sapply(1:length(host.nos), function(x) {
        min(sapply(which(tip.map==x), function(y) get.node.height(tree, y) - max.sampling.to.noninfectious))
      })
    }
    
    if(max.infection.to.sampling < Inf){
      maximum.heights <- sapply(1:length(host.nos), function(x) {
        max(sapply(which(tip.map==x), function(y) get.node.height(tree, y) + max.infection.to.sampling))
      })
    }
    
    if(is.null(minimum.heights)){
      minimum.heights <- rep(-Inf, length(host.nos))
    }
    
    if(is.null(maximum.heights)){
      maximum.heights <- rep(Inf, length(host.nos))
    }
    
    height.limits <- cbind(minimum.heights, maximum.heights)
    
    # extra row for the unsampled state
    
    height.limits <- rbind(height.limits, c(-Inf, Inf))
    
  } else {
    height.limits <- cbind(rep(-Inf, length(host.nos)+1), rep(Inf, length(host.nos)+1))
  }
  
  bridge <- tryCatch(
    .build.bridge(tree, host.nos, tip.map),
    error=function(cond){
      stop("No transmission trees (without superinfection) are compatible with this set of tip labels")
    }
  )
  
  
  if(length(tree$tip.label)==length(host.nos)){
    hosts <- tree$tip.label
  } else {
    # you can do better than this _really_...
    
    hosts <- vector()
    
    for(i in 1:length(host.nos)){
      hosts[i] <- paste(tree$tip.label[which(tip.map==i)], collapse="-")
    }
    
  }
  
  if(max.unsampled > 0){
    hosts <- c(tree$tip.label, paste("uh", 1:max.unsampled, sep=""))
  }
  
  
  results <- .unified.up.phase(tree, phangorn::getRoot(tree), list(), max.unsampled, height.limits, bridge, bigz)
  
  tt.count = lapply(0:max.unsampled, function(x){
    
    sampled.hosts <- length(host.nos)
    visible.unsampled.hosts <- x
    
    results[[phangorn::getRoot(tree)]]$p[x+1] * choose(max.unsampled + sampled.hosts - 1, visible.unsampled.hosts + sampled.hosts - 1)
    
  })
  
  tt.count <- do.call(sum, tt.count)
  
  out <- list(tree = tree, hosts = hosts, tt.count = tt.count, height.limits = height.limits, bridge = bridge, node.calculations = results)
  
  class(out) <- append(class(out), "tt.generator")
  
  return(out)
}
# 
# .up.phase <- function(tree, node, node.calculations, max.unsampled){
#   if(is.tip(tree, node)){
#     
#     node.info <- list()
#     
#     # Rows are hosts except the last row is "unsampled". 
#     # Columns are the number of unsampled partition elements in the tree, shifted by 1 because zero is a thing
#     
#     v <- matrix(0, nrow = length(tree$tip.label) + 1, ncol = max.unsampled + 1)
#     v[node, 1] <- 1
#     
#     node.info$v <- v
#     node.info$p <- c(1, rep(0, max.unsampled))
#     node.info$pstar <- c(1, rep(0, max.unsampled))
#     if(ncol(v)>1){
#       node.info$ps <- colSums(v[1:length(tree$tip.label),])
#     } else {
#       node.info$ps <- sum(v)
#     }
#     node.info$pu <- v[length(tree$tip.label) + 1,]
#     
#     node.calculations[[node]] <- node.info
#     
#   } else {
#     for(child in phangorn::Children(tree, node)){
#       node.calculations <- .up.phase(tree, child, node.calculations, max.unsampled)
#     }
#     
#     node.info <- list()
#     
#     v <- matrix(0, nrow = length(tree$tip.label) + 1, ncol = max.unsampled + 1)
#     
#     # the row for the unsampled state
#     
#     v[length(tree$tip.label) + 1, ] <- sapply(0:max.unsampled, function(x){
#       if(x==0){
#         return(0)
#       } else {
#         # how many unsampled elements
#         
#         # how many that are not the root unsampled element
#         remaining.us.elements <- x - 1
#         out <- 0
#         for(i in 0:remaining.us.elements){
#           j = remaining.us.elements - i
#           kids <- phangorn::Children(tree,node)
#           
#           term <- node.calculations[[kids[1]]]$pstar[i+1] * node.calculations[[kids[2]]]$pstar[j+1]
#           out <- out + term
#         }
#         return(out)
#       }
#     } )
#     
#     node.info$pu <- v[length(tree$tip.label) + 1,]
#     
#     # the rest of the rows
#     kids <- phangorn::Children(tree,node)
#     
#     for(host in 1:length(tree$tip.label)){
#       temp <- sapply(0:max.unsampled, function(x){
#         out <- 0
#         for(i in 0:x){
#           j = x - i
#           term <- node.calculations[[kids[1]]]$v[host,i+1] * node.calculations[[kids[2]]]$pstar[j+1] +
#             node.calculations[[kids[2]]]$v[host,i+1] * node.calculations[[kids[1]]]$pstar[j+1]
#           out <- out + term
#         }
#         return(out)
#       })
#       v[host,] <- temp
#     }
#     
#     if(ncol(v)>1){
#       node.info$ps <- colSums(v[1:length(tree$tip.label),])
#     } else {
#       node.info$ps <- sum(v)
#     }
#     
#     node.info$p <- node.info$pu + node.info$ps
#     
#     node.info$pstar <- sapply(0:max.unsampled, function(x) {
#       out <- node.info$p[x+1]
#       for(i in 0:x){
#         j = x - i
#         
#         term <- node.calculations[[kids[1]]]$pstar[i+1] * node.calculations[[kids[2]]]$pstar[j+1]
#         out <- out + term
#       }
#       return(out)
#     })
#     
#     node.info$v <- v
#     node.calculations[[node]] <- node.info
#   }
#   return(node.calculations)
# }
# 
# .height.aware.up.phase <- function(tree, node, node.calculations, height.limits){
#   if(is.tip(tree, node)){
#     
#     if(get.node.height(tree, node) < height.limits[node,1] | get.node.height(tree, node) > height.limits[node,2]){
#       stop("Bad node height limits for node ",node, sep="")
#     }
#     
#     node.info <- list()
#     
#     v <- rep(0, length(tree$tip.label))
#     v[node] <- 1
#     
#     node.info$v <- v
#     node.info$p <- 1
#     node.info$pstar <- rep(1, length(tree$tip.label))
#     node.calculations[[node]] <- node.info
#     
#   } else {
#     for(child in phangorn::Children(tree, node)){
#       node.calculations <- .height.aware.up.phase(tree, child, node.calculations, height.limits)
#     }
#     
#     node.info <- list()
#     
#     v <- rep(0, length(tree$tip.label))
#     pstarprod <- rep(1, length(tree$tip.label))
#     
#     kids <- phangorn::Children(tree,node)
#     
#     for(kid in kids){
#       desc.tips <- unlist(phangorn::Descendants(tree, kid, type="tips"))
#       v[desc.tips] <- node.calculations[[kid]]$v[desc.tips]
#       
#       for(kid2 in kids){
#         if(kid != kid2){
#           v[desc.tips] <- v[desc.tips]*node.calculations[[kid2]]$pstar[desc.tips]
#         }
#       }
#       pstarprod <- pstarprod*node.calculations[[kid]]$pstar
#     }
#     
#     height <- get.node.height(tree, node)
#     
#     # replace v for those out of range with 0
#     
#     v[which(height.limits[,1] > height | height.limits[,2] < height )] <- 0
#     
#     node.info$v <- v
#     node.info$p <- sum(v)
#     node.info$pstar <- rep(node.info$p, length(tree$tip.label))
#     tips.lower <- which(height.limits[,1] <= get.node.height(tree,node))
#     
#     node.info$pstar[tips.lower] <- node.info$pstar[tips.lower] + pstarprod[tips.lower]
#     node.calculations[[node]] <- node.info
#   }
#   return(node.calculations)
# }
# 
# .multiply.sampled.up.phase <- function(tree, node, node.calculations, bridge){
#   if(is.tip(tree, node)){
#     
#     node.info <- list()
#     
#     v <- rep(0, length(unique(stats::na.omit(bridge))))
#     v[bridge[node]] <- 1
#     
#     node.info$v <- v
#     node.info$p <- 1
#     node.info$pstar <- 1
#     node.calculations[[node]] <- node.info
#     
#   } else {
#     for(child in phangorn::Children(tree, node)){
#       node.calculations <- .multiply.sampled.up.phase(tree, child, node.calculations, bridge)
#     }
#     
#     node.info <- list()
#     
#     v <- rep(0, length(unique(stats::na.omit(bridge))))
#     pstarprod <- 1
#     
#     kids <- phangorn::Children(tree,node)
#     
#     desc.tips.1 <- unlist(phangorn::Descendants(tree, kids[1], type="tips"))
#     desc.tips.2 <- unlist(phangorn::Descendants(tree, kids[2], type="tips"))
#     if(length(intersect(bridge[desc.tips.1], bridge[desc.tips.2]))==1){
#       # this is the situation where this is the top of a bridge
#       
#       only.valid.host <- intersect(bridge[desc.tips.1], bridge[desc.tips.2])
#       
#       v[only.valid.host] <- node.calculations[[kids[1]]]$v[only.valid.host] * node.calculations[[kids[2]]]$v[only.valid.host]
#       
#     } else {
#       for(host.no in 1:max(stats::na.omit(bridge))){
#         if(!is.na(bridge[node]) & bridge[node] != host.no){
#           v[host.no] <- 0
#         } else if(host.no %in% bridge[desc.tips.1]){
#           v[host.no] <- node.calculations[[kids[1]]]$v[host.no] * node.calculations[[kids[2]]]$pstar
#         } else if(host.no %in% bridge[desc.tips.2]) {
#           v[host.no] <- node.calculations[[kids[2]]]$v[host.no] * node.calculations[[kids[1]]]$pstar
#         }
#       }
#     }
#     
#     if(is.na(bridge[node])){
#       pstar <- sum(v) + node.calculations[[kids[1]]]$pstar * node.calculations[[kids[2]]]$pstar
#     } else {
#       pstar <- sum(v)
#     }
#     
#     node.info$v <- v
#     node.info$p <- sum(v)
#     node.info$pstar <- pstar
#     node.calculations[[node]] <- node.info
#   }
#   return(node.calculations)
# }

# height.limits must have an extra (-Inf, Inf) for the unsampled state

.unified.up.phase <- function(tree, 
                              node, 
                              node.calculations, 
                              max.unsampled = 0, 
                              height.limits = cbind(rep(-Inf, length(tree$tip.label)+1), rep(Inf, length(tree$tip.label)+1)),
                              bridge = c(1:(length(tree$tip.label)), rep(NA, tree$Nnode)),
                              bigz = F){
  
  nhosts <- length(unique(stats::na.omit(bridge)))
  
  if(is.tip(tree, node)){
    
    tiphost <- bridge[node]
    
    if(get.node.height(tree, node) < height.limits[tiphost,1] | get.node.height(tree, node) > height.limits[tiphost,2]){
      stop("Bad node height limits for node ",node, sep="")
    }
    
    node.info <- list()
    
    # Rows are hosts except the last row is "unsampled". 
    # Columns are the number of unsampled partition elements in the tree, shifted by 1 because zero is a thing
    
    v <- matrix(0, nrow = nhosts + 1, ncol = max.unsampled + 1)
    v[tiphost, 1] <- 1
    
    # pstar now has to be a matrix as it has to pay attention to both the unsampled host count and the time limits
    
    pstar <- matrix(0, nrow = nhosts + 1, ncol = max.unsampled + 1)
    pstar[,1] <- 1
    
    p <- c(1, rep(0, max.unsampled))
    
    if(bigz){
      v <- as.bigz(v)
      p <- as.bigz(p)
      pstar <- as.bigz(pstar)
      
    } 
    
    ps <- colSums.fixed(v[1:nhosts, , drop=F])
    
    node.info$v <- v
    node.info$p <- p
    node.info$pstar <- pstar
    node.info$ps <- ps
    node.info$pu <- v[nhosts + 1,]
    
    node.calculations[[node]] <- node.info
    
  } else {
    
    kids <- phangorn::Children(tree,node)
    
    for(child in kids){
      node.calculations <- .unified.up.phase(tree, child, node.calculations, max.unsampled, height.limits, bridge, bigz)
    }
    
    node.info <- list()
    
    v <- matrix(0, nrow = nhosts + 1, ncol = max.unsampled + 1)
    
    if(bigz){
      v <- as.bigz(v)
    }
    
    # the row for the unsampled state
    
    if(!is.na(bridge[node])){
      # never unsampled
      v[nhosts + 1, ] <- rep(0, max.unsampled + 1)
    } else {
      v[nhosts + 1, ] <- do.call(c, lapply(0:max.unsampled, function(x){
        if(x==0){
          return(0)
        } else {
          # How many unsampled elements?
          
          # One element is accounted for at the root
          remaining.us.elements <- x - 1
          
          out <- 0
          
          distribution.of.us <- divide.k.into.n(remaining.us.elements, length(kids))
          
          for(i in 1:ncol(distribution.of.us)){
            term <- 1
            for(j in 1:nrow(distribution.of.us)){
              term <- term * node.calculations[[kids[j]]]$pstar[nhosts + 1, distribution.of.us[j,i]+1]
            }
            out <- out + term
          }
          return(out)
        }
      } ))
    }
    
    node.info$pu <- v[nhosts + 1,]
    
    # the rest of the rows
    
    for(host in 1:nhosts){
      if(get.node.height(tree, node) < height.limits[host,1] | get.node.height(tree, node) > height.limits[host,2]){
        # this node is outside this host's height limits
        v[host,] <- 0
      } else if(!is.na(bridge[node]) & bridge[node]!=host) {
        # this node is in some other host's bridge
        v[host,] <- 0
      } else if(!(host %in% bridge[unlist(phangorn::Descendants(tree, node, type="tips"))] )){
        # this host is not a tip of this subtree. THIS NEEDS TO BE IN THE PAPER
        v[host,] <- 0
      } else {
        
        temp <- lapply(0:max.unsampled, function(x){
          out <- 0
          
          distribution.of.us <- divide.k.into.n(x, length(kids))
          
          for(i in 1:ncol(distribution.of.us)){
            term <- 1
            for(j in 1:nrow(distribution.of.us)){
              desc.tips <- unlist(phangorn::Descendants(tree, kids[j], type="tips"))
              
              if(host %in% bridge[desc.tips]){
                # if child j is an ancestor of a tip from host
                term <- term * node.calculations[[kids[j]]]$v[host, distribution.of.us[j,i]+1]
              } else {
                term <- term * node.calculations[[kids[j]]]$pstar[host, distribution.of.us[j,i]+1]
              }
            }
            out <- out + term
          }
          return(out)
        })
        
        temp <- do.call(c, temp)
        
        v[host,] <- temp
      }
    }
    
    node.info$ps <- colSums.fixed(v[1:nhosts,, drop=F])
    
    node.info$p <- node.info$pu + node.info$ps
    
    tips.lower <- height.limits[,1] <= get.node.height(tree,node)
    
    pstar <- lapply(1:(nhosts+1), function(x){
      ffs <- lapply(0:max.unsampled, function(y){
        
        out <- node.info$p[y+1]
        
        # That's all she wrote if this is a bridge node or it is outside the height limits. Otherwise:
        if(is.na(bridge[node]) & tips.lower[x]){
          
          distribution.of.us <- divide.k.into.n(y, length(kids))
          for(i in 1:ncol(distribution.of.us)){
            term <- 1
            for(j in 1:nrow(distribution.of.us)){
              term <- term * node.calculations[[kids[j]]]$pstar[x, distribution.of.us[j,i]+1]
            }
            out <- out + term
          }
          
        }
        return(out)
      })
      do.call(c, ffs)
    })
    
    pstar <- do.call(rbind, pstar)
    
    node.info$pstar <- pstar
    node.info$v <- v
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

