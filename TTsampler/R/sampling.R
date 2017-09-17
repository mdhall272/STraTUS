.down.phase <- function(tree, node, result.vector, info, us.count){
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


.height.aware.down.phase <- function(tree, node, result.vector, info, height.limits){
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
    result.vector <- multiply.sampled.down.phase(tree, child, result.vector, info, bridge)
  }
  return(result.vector)
}