rm(list = ls())

library(ape)

#BigMEMOIR simulations

#simulates 4-state independent editing for a given period of time (1 day by default)
independent_sim_endpt <- function(subunits = 6,
                                   time = 30,
                                   start_time = 0,
                                   p,
                                   start_point = rep(0, 6)){
  results <- matrix(0, nrow = subunits, ncol = subunits + 1)
  times <- rep(0, subunits)
  
  results[,1] <- start_point
  times[1] <- start_time
  
  i <- 2
  
  #If everything edited, return
  if(times[i-1] < time & all(results[,i - 1] != 0)){
    return(list(results[,i-1], times[i-2], times[i-1], 0))
  }
  
  #Loop through until we exceed the given time or edit the barcode to completion
  while((times[i-1] < time) & !all(results[,i - 1] != 0)){
    
    #get the propensities for all possible transitions
    propensities <- p
    
    #make propensities 0 if the edit has already occurred at a position
    propensities[which(results[,i-1] != 0),] <- rep(0, length(propensities[1,]))
    
    #Store previously edited sites in the current row
    results[,i] <- results[,i-1]
    
    #Pick two random numbers
    rand_nums <- runif(2, 0, 1)
    
    #Find the next reaction time (exponential distribution)
    #The times vector stores cumulative times
    times[i] <- 1/sum(propensities) * log(1 / rand_nums[1]) + times[i-1]
    
    if(sum(propensities) > 0){
      
      #3 edit states used here
      m = 3
      
      #Figure out which site was edited
      site_probs <- rowSums(propensities) / sum(propensities)
      site_probs_upper <- sapply(1:length(site_probs), function(x) {
        return(sum(site_probs[1:x]))
      })
      site_probs_lower <- site_probs_upper - site_probs
      
      edit_site <- intersect(which(site_probs_lower < rand_nums[2]), 
                             which(site_probs_upper > rand_nums[2]))
      
      #Figure out what the edit outcome is
      norm_rand_2 <- (rand_nums[2] - site_probs_lower[edit_site]) / site_probs[edit_site]
      
      edit_state_upper <- sapply(1:m, function(x) {
        return(sum(propensities[edit_site, 1:x]) / sum(propensities[edit_site,]))
      })
      edit_state_lower <- edit_state_upper - (propensities[edit_site,] / sum(propensities[edit_site,]))
      
      edit_outcome <- intersect(which(edit_state_lower < norm_rand_2),
                                which(edit_state_upper > norm_rand_2))
      
      results[edit_site, i] <- edit_outcome
    }
    
    i <- i + 1
    
    #If everything edited, return
    if(times[i-1] < time & all(results[,i - 1] != 0)){
      return(list(results[,i-1], times[i-2], times[i-1], 0))
    }
  }
  
  #Returns the edited state and final editing time prior to the specified point
  #return(list(results, times))
  return(list(results[,i-2], times[i-2], times[i-1], sum(propensities)))
}

#simulates hypercascade editing over a given phylogeny
independent_cell_tree <- function(subunits = 6,
                                  p = list(),
                                  integrations = 60,
                                  tree_phylogeny,
                                  root_height){
  #For this implementation, build the tree based on a given phylogeny
  build_tree <- function(curr_node, next_node, prev_barcs){
    
    #node_depths <- node.depth.edgelength(tree_phylogeny)
    tree_edges <- tree_phylogeny$edge
    tree_edge_lengths <- tree_phylogeny$edge.length
    
    #Find the next time point
    if(curr_node == 0){
      next_time = 30 - root_height
    } else {
      next_time = tree_edge_lengths[which(tree_edges[,1] == curr_node & tree_edges[,2] == next_node)]
    }
    
    if(next_node > length(tree_phylogeny$tip.label)){
      
      #Get the edit state at the cell division time
      barc_states <- lapply(1:integrations, function (x) {
        barc_state <- getEdits(start_time = 0,
                               end_time = next_time,
                               start_point = prev_barcs[[x]])[[1]]
        return(barc_state)
      })
      
      #store the node name as the current index
      node_name <- next_node
      
      #find the left and right nodes
      left_next <- tree_edges[tree_edges[,1] == next_node,2][1]
      right_next <- tree_edges[tree_edges[,1] == next_node,2][2]
      
      #Determine the left and right branches
      LEFT <- build_tree(curr_node = next_node, next_node = left_next, prev_barcs = barc_states)
      RIGHT <- build_tree(curr_node = next_node, next_node = right_next, prev_barcs = barc_states)
      
      return()
    } else {
      #Get the edit state for the final branch
      barc_states <- lapply(1:integrations, function (x) {
        barc_state <- getEdits(start_time = 0,
                               end_time = next_time,
                               start_point = prev_barcs[[x]])[[1]]
        return(barc_state)
      })
      
      #store the barcode results in lists in the parent scope
      barcode_results[[next_node]] <<- barc_states
      
      return()
    }
  }
  
  #Helper function to determine the edited state of the internal node cell
  getEdits <- function(start_time, end_time, start_point){
    
    sim_res <- independent_sim_endpt(subunits = subunits,
                                      time = end_time, 
                                      start_time = start_time,
                                      p = p, 
                                      start_point = start_point)
    barcode <- sim_res[[1]]
    
    return(list(barcode))
  }
  
  #Initialize lists to store the barcodes for the nodes
  barcode_results <- list()
  
  #Simulate the tree
  simulated_editing <- build_tree(curr_node = 0, 
                                  next_node = length(curr_tree$tip.label) + 1,
                               prev_barcs = lapply(1:integrations, function(x) {
                                 return(rep(0, subunits))
                               }))
  return(barcode_results)
}

#Modify to desired working directory
setwd("D:/OneDrive_Backup_220312/OneDrive - California Institute of Technology/Lab_Notebook/Experiment_by_numbers/e165/")

#Read in tree data and log
yule_trees <- read.nexus(file = 'truth.trees')
tree_log <- read.table(file = 'truth.log', header = T)

#Load in the edit rate data - original data from e020
load('BigMEMOIR_edit_rates.RData')

#Loop over all trees (except burn in tree #1), simulate barcode editing, save results to .nex files
for(i in 2:101){
  #Pick a tree and get the root height
  curr_tree <- yule_trees[[i]]
  root_height <- tree_log$TreeHeight[i]
  
  #Traverse through the tree from root to tips
  root <- length(curr_tree$tip.label) + 1
  
  barcode_results <- independent_cell_tree(subunits = 6, 
                                 p = p_active*100*(3/30),
                                 integrations = 167, 
                                 tree_phylogeny = curr_tree, 
                                 root_height = root_height)
  
  barcode_results <- sapply(1:length(barcode_results), function(x) {
    barc <- unlist(barcode_results[[x]])
    barc <- paste(barc, collapse = "")
    barc <- gsub("3", "T", barc)
    barc <- gsub("2", "C", barc)
    barc <- gsub("1", "G", barc)
    barc <- gsub("0", "A", barc)
    return(barc)
  })
  
  names(barcode_results) <- curr_tree$tip.label
  
  barcode_results <- sapply(1:length(barcode_results), function(x) {
    curr_barc <- barcode_results[x]
    curr_barc <- strsplit(curr_barc, split = "")
    return(curr_barc)
  })
  write.nexus.data(barcode_results, file = paste("220731_BigMEMOIR_Sims_3day_rescale/Tree", i, "_Barcodes.nex", sep = ""))
}

