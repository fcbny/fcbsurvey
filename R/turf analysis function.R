turf.analysis <- function(dataset,max.reach=0.99,seed=NULL,ignore=NULL,
                          max.k="n",plot.it=TRUE) {
  require(turfR)
  if (class(dataset) != "data.frame") {
  stop("data object is not a data frame")
  }
  dataset.names <- colnames(dataset)
  dataset.names[1] <- "id"
  dataset.names[2] <- "weight"
  colnames(dataset) <- dataset.names
  
  if (!(class(dataset$id) %in% c("numeric","integer"))) {
    stop("identifiers is not of type numeric or integer")
  }
  if (!(class(dataset$weight) %in% c("numeric","integer"))) {
    stop("weight variable is not of type numeric or integer")
  }
  
  response_types <- unique(sapply(dataset[,3:ncol(dataset)],class))
  if (TRUE %in% (c("character","logical") %in% response_types)) {
    stop("some or all of your response variables are not numeric or integer")
  }
  
  response_classes <- sapply(dataset[,3:ncol(dataset)],function(x) {
    length(unique(x))
  })
  
  if (FALSE %in% unique(response_classes == 2)) {
    stop("some or all of your response variables are not binary")
  }
  
  if (!is.null(ignore)) {
    if (ignore %in% names(dataset)) {
      dataset <- dataset[,-match(ignore,colnames(dataset))]
    } else {
      stop("ignore variable does not exist in data frame")
    }
  }
  
  if (!is.null(seed)) {
    if (seed %in% names(dataset)) {
      seed.values <- apply(dataset[,seed],1,function(x) {
        ifelse(sum(x) == length(x),1,0)
      })
      
      tmp.seeds <- paste0(seed,collapse="|")
      dataset <- dataset[,-match(seed,colnames(dataset))]
      dataset[,tmp.seeds] <- seed.values
      seed <- tmp.seeds
    } else {
      stop("seed variable does not exist in data frame")
    }
  }
  
  # Generate initial seed combination---
  # This can be seed = "default", where item with most reach is used
  # to seed this function. Otherwise, if seed != "default", then
  # the seed provided will be AND together, creating a new dummy
  # item that would seed the function
  initial_reach <- apply(dataset[,3:ncol(dataset)],2,sum)
  initial_reach <- data.frame(items=names(initial_reach),
                              reach=initial_reach / length(unique(dataset$id)))
  initial_reach <- initial_reach[order(initial_reach$reach,decreasing=TRUE),]
  
  # Generate elements for results---
  # items.reach is named vector, containing the "optimal" item bundle
  # with increasing reach if item i and i - 1, i - 2, etc. are all included
  # all.items.reach is a list, containing all possible combinations of each
  # item in the "optimal item bundle
  items.reach <- c(); all.items.reach <- c()
  items.reach[rownames(initial_reach)[1]] <- initial_reach$reach[1]
  all.items.reach <- list(item1=initial_reach)
  
  
  combo_id_items <- function(x,items) {
    r <- apply(x[,4:ncol(x)],1,function(x){items[x == 1]})
    r <- as.data.frame(t(r))
    r <- cbind(combo=x$combo,reach=x$rchX,frequency=x$frqX,r)
    
    return(r)
  }
  
  items <- colnames(dataset)[!grepl("id|weight",colnames(dataset))]
  if (max.k == "n") {
    max.k <- length(items)
  } else {
    if (max.k > length(items)) {
      stop("max.k is larger than number of available items")
    }
  }
  k <- 2
  repeat {
    combos <- turfR::turf.combos(max.k,k)
    
    # Subset relevant combinations---
    # Keep only combinations contain at least the items within current
    # "optimal" item bundle
    # e.g. If it was decided that the last 2 items in the optimal bundle
    #      is item A and item B, then combination with k=3 MUST be at least
    #      A u B u _, where _ is any other item NOT A or B
    combos.keep <- (
      apply(as.matrix(combos[[1]][,match(names(items.reach),items)]),1,
            function(x) {ifelse(sum(x) != length(items.reach),FALSE,TRUE)})
    )
    
    combos[[1]] <- combos[[1]][combos.keep,]
    r <- turfR::turf(data=dataset,n=(ncol(dataset) - 2),k=k,
                     combos=combos,depth=1)
    r <- combo_id_items(r$turf[[1]],items)
    
    best.r <- r[1,]
    
    all.items.reach[[paste0("item",k)]] <- (
      data.frame(items=apply(r[,4:ncol(r)],1,
                             function(x) {x[!(x %in% names(items.reach))]}),
                 reach=r$reach - max(items.reach)))
    items.k.combo <- as.vector(t(best.r[,4:ncol(best.r)]))
    
    items.reach[items.k.combo[!(items.k.combo %in% names(items.reach))]] <- (
      r$reach[1])
    
    if ((max(items.reach) >= max.reach) | (k >= max.k)) {break}  # repeat break
    
    k <- k + 1
  }

  result <- data.frame(items=names(items.reach),
                       initial_reach=items.reach - c(min(items.reach),diff(items.reach)),
                       incremental_reach=c(min(items.reach),diff(items.reach)),
                       total_reach=items.reach,
                       stringsAsFactors=F)
  rownames(result) <- 1:nrow(result)
  
  
  if (plot.it) {
    require(ggplot2);require(reshape2);require(scales)
    best_items <- result
    best_items$items <- factor(best_items$items,
                               levels=best_items$items,ordered=TRUE)
    best_items <- melt(best_items,id.vars="items")
    best_items <- best_items[best_items$variable != "total_reach",]
    
    p <- ggplot(best_items) + aes(x=factor(items),y=value,fill=variable) +
      geom_bar(stat="identity",color="black") +
      theme_bw() + scale_y_continuous(labels=percent) +
      labs(x="Total Reach",y="Items Selected",
           title="Recommended Items for Maximizing Reach") +
      theme(legend.position="top",
            plot.title=element_text(face="bold",size=12),
            legend.title=element_blank()) +
      scale_fill_manual(values=c("initial_reach"="darkgrey",
                                 "incremental_reach"="orange"),
                        labels=c("Initial Reach","Incremental Reach"))
      
    r <- list(best_items=result,all=all.items.reach,plot=p)
  } else {
    r <- list(best_items=result,all=all.items.reach)
  }
  
  return(r)
}