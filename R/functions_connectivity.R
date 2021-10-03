# FUNCTION TO FIND UPSTREAM IDs -------------------------------------------------------------------------

# new column for master_table with concatenated next downstream ID
next_down <- function(tab){
  tb <- tab[,(ncol(tab)-1):ncol(tab)]
  colnames(tb) <- c('id','nd')
  d1 <- rep(0,nrow(tb))
  for(i in 1:nrow(tb)){
    if(tb$nd[i] != 0 & tb$nd[i] %in% tb$id){
      d1[i] <- unique(tb$nd[tb$id == tb$nd[i]])
    }else{
      d1[i] <- 0
    }
  }
  return(d1)
  
}

# function that takes in input the master table and selects the unique upstream IDs
upstream_sel <- function(t,id){
  return(
    unique(
      do.call('c',as.list(apply(t,1,function(x){
        sel <- which(x == id)
        if(length(sel) != 0) return(as.vector(x[1:sel]))
      })))
    )
  )
}

find_upstream_ids <- function(t,IDs){
  
  ### concatenate next_downstream basins ids
  c = 1
  while(sum(t[,ncol(t)]) != 0){
    t[,paste0('d',c)] <- next_down(t)
    c = c+1
  }
  
  l <- lapply(IDs,function(x) upstream_sel(t,x))
  names(l) <- IDs
  return(l)
}

# function that allocates the row to a different group based on hybas_id
# if no group is found, returns 0
assign_group <- function(hybas_id,group){
  a <- 0
  for(j in 1:length(group)){
    if(hybas_id %in% group[[j]]){
      a <- j
      break
    }
  }
  return(a) 
}
