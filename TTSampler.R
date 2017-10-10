#!/usr/bin/env Rscript

library(phangorn)
library(gtools)
library(argparse)

arg_parser = ArgumentParser(description="Randomly sample from the set of transmission trees compatible with the input phylogeny")
arg_parser$add_argument("-i", "--infectionEnds", action="store_true", default=F, help="If present, hosts are assumed to cease to be infected at the time at which their last tip appears. Cannot be used simultaneously with a timings file, unsampled hosts, or multiple samples per host.")
arg_parser$add_argument("-I", "--timingsFile", action="store", help="The argument expected here is a CSV file that dictates the earliest possible infection date and latest possible recovery date for each host. It should have a header and three columns - in order, tip label, earliest infection date as a tree height, and latest recovery date as a tree height. Missing tip labels will not be given time limits. NAs in either time column mean an infinite bound. Cannot be used simultaneously with the -i flag, unsampled hosts, or multiple samples per host.")
arg_parser$add_argument("-m", "--multipleSamplingFile", action="store", help="The argument expected here is a CSV file that maps each tip to a host which that tip was sampled from. A header and two columns are expected, in order, tip label and host ID. This is unnecessary if one tip was sampled per host. Cannot be used simultaneously with time limits on infections or unsampled hosts.")
arg_parser$add_argument("-u", "--unsampledHostCount", action="store", default=0, help="The number of unsampled hosts. Default is 0. Values >0 cannot be used simultaneously with time limits on infections or multiple samples per host.")
arg_parser$add_argument("-s", "--sampleCount", action="store", default=1, help="How many transmission trees to sample. Default is 1.")
arg_parser$add_argument("-d", "--drawOutput", action="store_true", help="If present, draw every sample as an annotated tree PDF. Requires ggtree.")
arg_parser$add_argument("treeFileName", action="store", help="File name for the input tree in Newick format.")
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
    region.count.weights <- sapply(0:unsampled.count, function(x) counts[x+1]*choose(tip.count + unsampled.count - 1, tip.count + x - 1))

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
