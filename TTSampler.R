#!/usr/bin/env Rscript

library(phangorn)
library(gtools)
library(argparse)

arg_parser = ArgumentParser(description="Randomly sample from the set of transmission trees compatible with the input phylogeny")
arg_parser$add_argument("-i", "--infectionEnds", action="store_true", default=F, help="If present, hosts are assumed to cease to be infected at the time at which their last tip appears. Cannot be used simultaneously with a timings file, unsampled hosts, or multiple samples per host.")
arg_parser$add_argument("-I", "--timingsFile", action="store", help="The argument expected here is a CSV file that dictates the earliest possible infection date and latest possible recovery date for each host. It should have a header and three columns - in order, tip label, label, earliest infection date as a tree height, and latest recovery date as a tree height. Missing tip labels will not be given time limits. NAs in either time column mean an infinite bound. Cannot be used simultaneously with the -i flag, unsampled hosts, or multiple samples per host.")
arg_parser$add_argument("-m", "--multipleSamplingFile", action="store", help="The argument expected here is a CSV file that maps each tip to a host which that tip was sampled from. A header and two columns are expected, in order, tip label and host ID. This is unnecessary if one tip was sampled per host. Cannot be used simultaneously with time limits on infections or unsampled hosts.")
arg_parser$add_argument("-u", "--unsampledHostCount", action="store", default=0, help="The number of unsampled hosts. Default is 0. An unsampled count of greater than 1 is currently incompatible with time limits for infectiousness. Cannot be used simultaneously with time limits on infections or multiple samples per host.")
arg_parser$add_argument("-s", "--sampleCount", action="store", default=1, help="How many transmission trees to sample. Default is 1.")
arg_parser$add_argument("-d", "--drawOutput", action="store_true", help="If present, draw every sample as an annotated tree PDF. Requires ggtree.")
arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names. Tree files are assumed to end in '.tree'. This must be the 'processed' tree produced by SplitTreesToSubtrees.R.")
arg_parser$add_argument("outputString", action="store", help="A string identifying output files")

args <- arg_parser$parse_args()

tree.file.name <- args$treeFileName
unsampled.count <- as.numeric(args$unsampledHostCount)
sample.count <- as.numeric(args$sampleCount)
output.string <- args$outputString

multiple.sampling <- !is.null(args$multipleSamplingFile)

if(args$infectionEnds & !is.null(args$timingsFile)){
  stop("You specified both infections ending at sampling and an infection times file. Please do one or the other.")
}

draw <- args$drawOutput
if(draw){
  library(ggtree)
}

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

up.phase <- function(tree, node, node.calculations, max.unsampled){
  if(is.tip(tree, node)){
    
    node.info <- list()
    
    # rows are patients except the last row is "unsampled". Columns are the number of unsampled partition elements in the tree, shifted by 1 because zero is a thing
    
    v <- matrix(0, nrow = tip.count + 1, ncol = max.unsampled + 1)
    v[node, 1] <- 1
    
    node.info$v <- v
    node.info$p <- c(1, rep(0, max.unsampled))
    node.info$pstar <- c(1, rep(0, max.unsampled))
    if(ncol(v)>1){
      node.info$ps <- colSums(v[1:tip.count,])
    } else {
      node.info$ps <- sum(v)
    }
    node.info$pu <- v[tip.count + 1,]
    
    node.calculations[[node]] <- node.info
    
  } else {
    for(child in Children(tree, node)){
      node.calculations <- up.phase(tree, child, node.calculations, max.unsampled)
    }
    
    node.info <- list()
    
    v <- matrix(0, nrow = tip.count + 1, ncol = max.unsampled + 1)
    
    # the row for the unsampled state
    
    v[tip.count + 1, ] <- sapply(0:max.unsampled, function(x){
      if(x==0){
        return(0)
      } else {
        # how many unsampled elements
        
        # how many that are not the root unsampled element
        remaining.us.elements <- x - 1
        out <- 0
        for(i in 0:remaining.us.elements){
          j = remaining.us.elements - i
          children <- Children(tree,node)
          
          term <- node.calculations[[children[1]]]$pstar[i+1] * node.calculations[[children[2]]]$pstar[j+1]
          out <- out + term
        }
        return(out)
      }
    } )
    
    node.info$pu <- v[tip.count + 1,]
    
    # the rest of the rows
    kids <- Children(tree,node)
    
    for(host in 1:tip.count){
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
      node.info$ps <- colSums(v[1:tip.count,])
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

down.phase <- function(tree, node, result.vector, info, us.count){
  #  cat("At node",node,"with",us.count,"remaining unsampled regions\n")
  kids <- Children(tree, node)
  vmatrix <- info[[node]]$v
  if(node == getRoot(tree)){
    
    prob.weights <- vmatrix[,us.count+1]
    
    result <- sample(1:(tip.count+1),1, prob=prob.weights)
    
    if(result == tip.count + 1){
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
      
      result <- sample(1:(tip.count+2),1, prob=c(new.pat.weights, continuation.weight))
      
      creep <- F
      
      if(result == tip.count + 2){
        creep <- T
        #whatever the parent has
        result <- parent.choice
      } else if(result == tip.count + 1){
        #new unsampled host
        result <- current.host.count + 1
        current.host.count <<- result
      } 
      
      result.vector[node] <- result
    }
  }
  if(!is.tip(tree, node)){
    if(us.count==0){
      result.vector <- down.phase(tree, kids[1], result.vector, info, 0)
      result.vector <- down.phase(tree, kids[2], result.vector, info, 0)
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
      } else if(result > tip.count){
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
      result.vector <- down.phase(tree, kids[1], result.vector, info, chosen.i)
      result.vector <- down.phase(tree, kids[2], result.vector, info, chosen.j)
    }
  }
  
  return(result.vector)
}

get.node.height <- function(tree, node){
  return(max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[node])
}

height.aware.up.phase <- function(tree, node, node.calculations, height.limits){
  if(is.tip(tree, node)){
    
    if(get.node.height(tree, node) < height.limits[node,1] | get.node.height(tree, node) > height.limits[node,2]){
      stop("Bad node height limits for node ",node, sep="")
    }
    
    node.info <- list()
    
    # rows are patients except the last row is "unsampled". Columns are the number of unsampled partition elements in the tree, shifted by 1 because zero is a thing
    
    v <- rep(0, tip.count)
    v[node] <- 1
    
    node.info$v <- v
    node.info$p <- 1
    node.info$pstar <- rep(1, tip.count)
    node.calculations[[node]] <- node.info
    
  } else {
    for(child in Children(tree, node)){
      node.calculations <- height.aware.up.phase(tree, child, node.calculations, height.limits)
    }
    
    node.info <- list()
    
    v <- rep(0, tip.count)
    pstarprod <- rep(1, tip.count)
    
    kids <- Children(tree,node)
    
    for(kid in kids){
      desc.tips <- unlist(Descendants(tree, kid, type="tips"))
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
    node.info$pstar <- rep(node.info$p, tip.count) 
    tips.lower <- which(height.limits[,1] <= get.node.height(tree,node))
    
    node.info$pstar[tips.lower] <- node.info$pstar[tips.lower] + pstarprod[tips.lower]
    node.calculations[[node]] <- node.info
  }
  return(node.calculations)
}

height.aware.down.phase <- function(tree, node, result.vector, info, height.limits){
  if(node == getRoot(tree)){
    result.vector[node] <- sample(1:tip.count, 1, prob=info[[node]]$v)
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
      result.vector[node] <- sample(1:tip.count,1,prob=host.weights)
    }
  }
  for(child in Children(tree, node)){
    result.vector <- height.aware.down.phase(tree, child, result.vector, info, height.limits)
  }
  return(result.vector)
}

multiply.sampled.up.phase <- function(tree, node, node.calculations, bridge){
  if(is.tip(tree, node)){
    
    node.info <- list()
    
    v <- rep(0, length(unique(na.omit(bridge))))
    v[bridge[node]] <- 1
    
    node.info$v <- v
    node.info$p <- 1
    node.info$pstar <- 1
    node.calculations[[node]] <- node.info
    
  } else {
    for(child in Children(tree, node)){
      node.calculations <- multiply.sampled.up.phase(tree, child, node.calculations, bridge)
    }
    
    node.info <- list()
    
    v <- rep(0, length(unique(na.omit(bridge))))
    pstarprod <- 1
    
    kids <- Children(tree,node)
    
    desc.tips.1 <- unlist(Descendants(tree, kids[1], type="tips"))
    desc.tips.2 <- unlist(Descendants(tree, kids[2], type="tips"))
    if (length(intersect(bridge[desc.tips.1], bridge[desc.tips.2]))==1){
      # this is the situation where this is the top of a bridge
      
      only.valid.host <- intersect(bridge[desc.tips.1], bridge[desc.tips.2])
      
      v[only.valid.host] <- node.calculations[[kids[1]]]$v[only.valid.host] * node.calculations[[kids[2]]]$v[only.valid.host]
      
    } else {
      for(host.no in 1:max(na.omit(bridge))){
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

multiply.sampled.down.phase <- function(tree, node, result.vector, info, bridge){
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
    result.vector <- multiply.sampled.down.phase(tree, child, result.vector, info, bridge)
  }
  return(result.vector)
}

# bl <- match(1:(length(tree$tip.label) + tree$Nnode), tree$edge[,2])

# ggtree(tree, aes(col=hosts[a.sample])) + geom_text2(aes(x=branch, label=round(tree$edge.length[bl],3)), vjust=1.25) + geom_tiplab()
# ggtree(tree) + geom_text2(aes(x=branch, label=round(tree$edge.length[bl],3)), vjust=1.25) + geom_tiplab()


# exhaustive.down.phase <- function(tree, node, result.matrix, info, us.count, columns){
#   cat("At node",node,"with",us.count,"remaining unsampled regions\n")
#   kids <- Children(tree, node)
#   vmatrix <- info[[node]]$v
#   if(node == getRoot(tree)){
# #    cat("At the root\n")
#     new.pat.weights <- vmatrix[,us.count+1]
#     
#     row.contents <- vector()
#     
#     for(weight.index in 1:length(new.pat.weights)){
#       row.contents <- c(row.contents, rep(weight.index, new.pat.weights[weight.index]))
#     }
#     new.index <- unique(current.host.counts) + 1
#     if(length(new.index)>1){
#       stop("Unsampled region counter is not working as intended")
#     }
#     current.host.counts[which(row.contents > tip.count)] <- new.index # add an extra -
#     
#     result.matrix[node,] <- row.contents
#   }
#   current.row <- result.matrix[node,columns ] 
#   
#   # otherwise you've already done this
#   for(parent.choice in unique(current.row)){
#     cat("Going down with parent ", parent.choice, "\n", sep="")
#     this.parent.columns <- columns[which(current.row==parent.choice)]
#     #    if(us.count == 0){
#     if(parent.choice %in% unlist(Descendants(tree, kids[1]))){
#       # creep on the right side        
#       left.repeats <- rep(parent.choice, info[[kids[1]]]$v[parent.choice,1])
#       
#       
#       right.repeats <- c(unlist(lapply(1:(tip.count+1), function(x) {
#         rep(x, info[[kids[2]]]$v[x,1])
#       })), rep(parent.choice,  info[[kids[2]]]$pstar[1] - info[[kids[2]]]$p[1] )
#       )
#       
#     } else if(parent.choice %in% unlist(Descendants(tree, kids[2]))){
#       # creep on the left side
#       right.repeats <- rep(parent.choice, info[[kids[2]]]$v[parent.choice,1])
#       
#       left.repeats <- c(unlist(lapply(1:(tip.count+1), function(x) {
#         rep(x, info[[kids[1]]]$v[x,1])
#       })), rep(parent.choice,  info[[kids[1]]]$pstar[1] - info[[kids[1]]]$p[1] )
#       )
#       
#     } else {
#       # creeping down 
#       left.repeats <- c(unlist(lapply(1:(tip.count+1), function(x) {
#         rep(x, info[[kids[1]]]$v[x,1])
#       })), rep(parent.choice,  info[[kids[1]]]$pstar[1] - info[[kids[1]]]$p[1] )
#       )
#       right.repeats <- c(unlist(lapply(1:(tip.count+1), function(x) {
#         rep(x, info[[kids[2]]]$v[x,1])
#       })), rep(parent.choice,  info[[kids[2]]]$pstar[1] - info[[kids[2]]]$p[1] )
#       )
#       
#     }
#     child.element.combinations <- expand.grid(left.repeats, right.repeats)
#     if(length(this.parent.columns) %% nrow(child.element.combinations)  != 0) {
#       stop("Something doesn't add up")
#     }
#     result.matrix[kids[1], this.parent.columns] <- child.element.combinations[,1]
#     result.matrix[kids[2], this.parent.columns] <- child.element.combinations[,2]
#     
#     for(kid in kids){
#       if(!is.tip(tree, kid)){
#         temp <- exhaustive.down.phase(tree, kid, result.matrix, info, us.count, this.parent.columns)
#         result.matrix <- temp
#       }
#     }
#     #    }
#   }
#   return(result.matrix)
# }


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

tree <- read.tree(tree.file.name)

if(args$infectionEnds){
  height.limits <- lapply(1:length(tree$tip.label), function(x) {
    c(get.node.height(tree, x), Inf)
  })
  height.limits <- do.call(rbind, height.limits)
  has.timings <- T
} else if(!is.null(args$timingsFile)){
  raw.timings <- read.csv(args$timingsFile, stringsAsFactors = F, header = T)
  height.limits <- lapply(1:length(tree$tip.label), function(x) {
    row.no <- which(raw.timings[,1] == tree$tip.label[x])
    if(length(row.no) > 1){
      stop(paste0("Tip labelled ",tree$tip.label[x]," occurs twice in the timings file."))
    } else if(length(row.no) == 1){
      if(as.numeric(raw.timings[row.no,2] <  raw.timings[row.no,3])){
        stop("Some timings in timings file are in the wrong order; they should be tree heights, with earliest infection time GREATER than latest recovery time")
      }
      row <- as.numeric(raw.timings[row.no,3:2])
      if(is.na(row[1])){
        row[1] <- -Inf
      }
      if(is.na(row[2])){
        row[2] <- Inf
      }
      
      return(row)
    } else {
      cat("Warning: tip ",tree$tip.label[x]," does not occur in timings file. No temporal limits will be applied to this host.\n", sep="")
      return(c(-Inf, Inf))
    }
  })
  height.limits <- do.call(rbind, height.limits)
  has.timings <- T
} else {
  has.timings <- F
}


if(multiple.sampling){
  raw.mi <- read.csv(args$multipleSamplingFile, stringsAsFactors = F, header = T)
  reordering <- match(tree$tip.label, raw.mi[,1])
  if(any(is.na(reordering))){
    stop("Not all tips appear in multiple infections file")
  }
  host.map <- raw.mi[reordering,2]
  host.labels <- unique(host.map)
  host.labels <- host.labels[order(host.labels)]
} 

# height.limits <- lapply(1:length(tree$tip.label), function(x){
#   h <- get.node.height(tree, x)
#   c(h, h+1)
# })
# height.limits <- do.call(rbind, height.limits)

if(length(which(c(multiple.sampling, has.timings, unsampled.count)>0))>1){
  stop("At present, only one of multiple sampling, incomplete sampling, and time limits can be specified")
}

if(!is.binary(tree)){
  stop("Binary trees only!\n")
}

tip.count <- length(tree$tip.label)

if(unsampled.count == 0 & !multiple.sampling){
  host.labels <- tree$tip.label
} else if(unsampled.count > 0){
  host.labels <- c(tree$tip.label, paste("uh", 1:unsampled.count, sep=""))
} 

annotations.df <- data.frame(node.no = seq(1,tip.count+tree$Nnode))

parents.df <- data.frame(child = host.labels)

if(unsampled.count==0 & !multiple.sampling){
  if(has.timings){
    
    node.calculations <- height.aware.up.phase(tree, getRoot(tree), list(), height.limits)
    
    if(node.calculations[[getRoot(tree)]]$p == 0){
      stop("Timings prevent any valid partitions of this tree")
    }
    
    for(i in seq(1:sample.count)){
      current.host.count <- tip.count
      a.sample <- height.aware.down.phase(tree, getRoot(tree), vector(), node.calculations, height.limits)
      annot.column.name <- paste0("Sample.",i,".annotation")
      annotations.df[[annot.column.name]] <- host.labels[a.sample]
      parents.vector <- get.parents.vector(tree, a.sample)
      par.column.name <- paste0("Sample.",i,".parent")
      
      parents.df[[par.column.name]] <- sapply(parents.vector, function(x) if(x!=0) return(host.labels[x]) else return("root") )
    }
  } else {
    node.calculations <- up.phase(tree, getRoot(tree), list(), 0)
    
    for(i in seq(1:sample.count)){
      current.host.count <- tip.count
      a.sample <- down.phase(tree, getRoot(tree), vector(), node.calculations, 0)
      annot.column.name <- paste0("Sample.",i,".annotation")
      annotations.df[[annot.column.name]] <- host.labels[a.sample]
      parents.vector <- get.parents.vector(tree, a.sample)
      par.column.name <- paste0("Sample.",i,".parent")
      
      parents.df[[par.column.name]] <- sapply(parents.vector, function(x) if(x!=0) return(host.labels[x]) else return("root") )
    }
    
  }
} else if(!multiple.sampling) {
  if(has.timings){
    stop("Combination of unsampled hosts and infectiousness time limits not currently supported.")
  } else {
    node.calculations <- up.phase(tree, getRoot(tree), list(), unsampled.count)
    
    counts <- node.calculations[[getRoot(tree)]]$p
    region.count.weights <- sapply(0:unsampled.count, function(x) counts[x+1]*choose(tip.count + x + unsampled.count - x - 1,unsampled.count - x))
    
    for(i in seq(1:sample.count)){
      
      no.regions <- sample(0:unsampled.count, 1, prob=region.count.weights)
      
      current.host.count <- tip.count
      
      a.sample <- down.phase(tree, getRoot(tree), vector(), node.calculations, no.regions)
      
      if(max(a.sample) != tip.count + no.regions){
        cat("Error - failed to insert enough unsampled hosts. Contact the author.\n")
        quit("N")
      }
      
      branch.us.position.choice <- vector()
      if(no.regions!=unsampled.count){
        branch.us.position.options <- combinations(tip.count + no.regions, unsampled.count - no.regions,  repeats.allowed = T )
        branch.us.position.choice <- branch.us.position.options[sample(1:nrow(branch.us.position.options), 1),]
      } 
      
      annot.column.name <- paste0("Sample.",i,".annotation")
      interv.column.name <- paste0("Sample.",i,".intervening")
      
      annotations.df[[annot.column.name]] <- host.labels[a.sample]
      
      interventions <- rep(0, tip.count + tree$Nnode)
      parents <- vector()
      for(host in 1:(tip.count + no.regions)){
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
        interventions[a.node] <- sum(branch.us.position.choice == host)
      }
      
      annotations.df[[interv.column.name]] <- interventions
      
      parents.vector <- get.parents.vector(tree, a.sample, interventions)
      par.column.name <- paste0("Sample.",i,".parent")
      parents.df[[par.column.name]] <- sapply(parents.vector, function(x) if(x!=0) return(host.labels[x]) else return("root") )
    }
  }
} else {
  
  # bridge-building
  bridge <- rep(NA, length(tree$tip.label) + tree$Nnode)
  for(host in host.labels){
    tips <- which(host.map == host)
    mrca <- mrca.phylo.or.unique.tip(tree, tips)
    for(tip in tips){
      current.node <- tip
      repeat{
        bridge[current.node] <- which(host.labels==host)
        if(current.node == mrca){
          break
        }
        current.node <- Ancestors(tree, current.node, type="parent")
        if(!is.na(bridge[current.node])){
          if(bridge[current.node] ==  which(host.labels==host)){
            break
          } else {
            stop("Overlapping bridges; no transmission trees are compatible")
          }
        }
      }
    }
  }
  node.calculations <- multiply.sampled.up.phase(tree, getRoot(tree), list(), bridge)
  for(i in seq(1:sample.count)){
    current.host.count <- tip.count
    a.sample <- multiply.sampled.down.phase(tree, getRoot(tree), vector(), node.calculations, bridge)
    annot.column.name <- paste0("Sample.",i,".annotation")
    annotations.df[[annot.column.name]] <- host.labels[a.sample]
    parents.vector <- get.parents.vector(tree, a.sample)
    par.column.name <- paste0("Sample.",i,".parent")
    parents.df[[par.column.name]] <- sapply(parents.vector, function(x) if(x!=0) return(host.labels[x]) else return("root") )
  }
  
}

write.csv(annotations.df, paste0(output.string,"_annotations.csv"), row.names = F, quote=F)
write.csv(parents.df, paste0(output.string,"_tt.csv"), row.names = F, quote=F)

if(draw){
  if(unsampled.count == 0){
    for(col in seq(2, ncol(annotations.df))){
      annots <- annotations.df[,col]
      first.in.split <- sapply(1:(tip.count + tree$Nnode), function(x) {
        if(getRoot(tree) == x){
          return(T)
        }
        if(annots[x] == annots[Ancestors(tree, x, type="parent")]){
          return(F)
        }
        return(T)
      })
      
      
      branch.colour.annots<- annots
      branch.colour.annots[which(first.in.split)] <- NA
      
      ggtree(tree, aes(col=branch.colour.annots), size=1.5) + 
        geom_tiplab(aes(col=annots), hjust=-1) + 
        geom_point(aes(col=annots), size=4) +
        scale_fill_hue(na.value = "grey") +
        scale_color_hue(na.value = "grey") 
      ggsave(paste0(output.string,"_sample_",(col-1),"_tree.pdf"))
    }
  } else {
    for(col in seq(2, ncol(annotations.df), by=2)){
      node.annots <- annotations.df[,col]
      branch.annots <- annotations.df[,col+1]
      
      first.in.split <- sapply(1:(tip.count + tree$Nnode), function(x) {
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
      
      adjustments <- rep(0.5, tip.count+tree$Nnode)
      adjustments[getRoot(tree)] <- 1.5
      
      ggtree(tree, aes(col=branch.colour.annots), size=1.5) + 
        geom_tiplab(aes(col=node.annots), hjust=-1) + 
        geom_point(aes(col=node.annots), size=4) +
        scale_fill_hue(na.value = "grey") +
        scale_color_hue(na.value = "grey") +
        geom_text(aes(x=branch, label=branch.annots, hjust=adjustments), vjust=-0.5, col="black")
      ggsave(paste0(output.string,"_sample_",col/2,"_tree.pdf"))
    }
  }
}
