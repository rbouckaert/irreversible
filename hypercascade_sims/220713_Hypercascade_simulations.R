rm(list = ls())

library(ape)

#Hypercascade simulations

#simulates hypercascade editing for a given period of time (1 day by default)
hypercascade_sim_endpt <- function(subunits = 20,
                                   time = 1,
                                   start_time = 0,
                                   p = c(H01_for_nls[H01_for_nls$Layer == 1 &
                                                       H01_for_nls$Edit_Pattern == "000",]$Edit_Rate,
                                         H01_for_nls[H01_for_nls$Layer == 2 &
                                                       H01_for_nls$Edit_Pattern == "1001",]$Edit_Rate,
                                         H01_for_nls[H01_for_nls$Layer == 3 &
                                                       H01_for_nls$Edit_Pattern == "1011",]$Edit_Rate,
                                         H01_for_nls[H01_for_nls$Layer == 4 &
                                                       H01_for_nls$Edit_Pattern == "1111",]$Edit_Rate),
                                   start_point = rep(0, 20*4 - 6),
                                   edit_penalties = c(1 - coef(H01_model)[1],
                                                      1 - coef(H01_model)[2],
                                                      1 - coef(H01_model)[2],
                                                      1 - coef(H01_model)[3],
                                                      1 - coef(H01_model)[4])){
  results <- matrix(0, nrow = 4*subunits - 6, ncol = 4*subunits - 6 + 2)
  times <- rep(0, 4*subunits - 6 + 1)
  
  results[,1] <- start_point
  times[1] <- start_time
  
  i <- 2
  
  #If everything edited, return
  if(times[i-1] < time & all(results[,i - 1] != 0)){
    return(list(results[,i-1], times[i-2], times[i-1], 0))
  }
  
  #Loop through until we exceed the given time or edit the barcode to completion
  while(times[i-1] < time){
    
    if(i > 2 & (sum(results[,i-2]) >= (4*subunits - 6))){
      break
    }
    
    #Mismatch penalty matrix
    const_mat <- matrix(1, nrow = length(results[,1]), ncol = 5)
    
    #check the previous rows for constraint edits
    #Layer 1
    const_mat[1:(subunits-1),2] <- sapply(1:(subunits-1), function(x) !results[x + subunits,i-1])
    const_mat[1:(subunits-2),3] <- sapply(1:(subunits-2), function(x) !results[x + 2*subunits - 1,i-1])
    const_mat[1:(subunits-3),4] <- sapply(1:(subunits-3), function(x) !results[x + 3*subunits - 3,i-1])
    #Layer 2
    const_mat[(subunits+1):(2*subunits-1),1] <- sapply((subunits+1):(2*subunits-1), function(x) results[x - subunits,i-1])
    const_mat[(subunits+1):(2*subunits-2),2] <- sapply((subunits+1):(2*subunits-2), function(x) !results[x + subunits - 1,i-1])
    const_mat[(subunits+1):(2*subunits-3),3] <- sapply((subunits+1):(2*subunits-3), function(x) !results[x + 2*subunits - 3,i-1])
    const_mat[(subunits+1):(2*subunits-1),5] <- sapply((subunits+1):(2*subunits-1), function(x) results[x - subunits + 1,i-1])
    #Layer 3
    const_mat[(2*subunits):(3*subunits-3),1] <- sapply((2*subunits):(3*subunits-3), function(x) results[x - subunits + 1,i-1])
    const_mat[(2*subunits):(3*subunits-4),2] <- sapply((2*subunits):(3*subunits-4), function(x) !results[x + subunits - 2,i-1])
    const_mat[(2*subunits):(3*subunits-3),4] <- sapply((2*subunits):(3*subunits-3), function(x) results[x - 2*subunits + 2,i-1])
    const_mat[(2*subunits):(3*subunits-3),5] <- sapply((2*subunits):(3*subunits-3), function(x) results[x - subunits + 2,i-1])
    #Layer 4
    const_mat[(3*subunits-2):(4*subunits-6),1] <- sapply((3*subunits-2):(4*subunits-6), function(x) results[x - subunits + 2,i-1])
    const_mat[(3*subunits-2):(4*subunits-6),3] <- sapply((3*subunits-2):(4*subunits-6), function(x) results[x - 3*subunits + 4,i-1])
    const_mat[(3*subunits-2):(4*subunits-6),4] <- sapply((3*subunits-2):(4*subunits-6), function(x) results[x - 2*subunits + 4,i-1])
    const_mat[(3*subunits-2):(4*subunits-6),5] <- sapply((3*subunits-2):(4*subunits-6), function(x) results[x - subunits + 3,i-1])
    
    const_mat[grep(0,const_mat[,1]),1] <- edit_penalties[1]
    const_mat[grep(0,const_mat[,2]),2] <- edit_penalties[2]
    const_mat[grep(0,const_mat[,3]),3] <- edit_penalties[3]
    const_mat[grep(0,const_mat[,4]),4] <- edit_penalties[4]
    const_mat[grep(0,const_mat[,5]),5] <- edit_penalties[5]
    
    #Update the editing probability per site with the amount of suppression at the site (this is the propensities)
    propensities <- c(rep(p[1], (subunits)), 
                      rep(p[2], (subunits - 1)),
                      rep(p[3], (subunits - 2)),
                      rep(p[4], (subunits - 3))) * apply(const_mat, MARGIN = 1, prod)
    
    #make propensity 0 if the edit has already occurred
    propensities[grep(1, results[,i-1])] <- 0
    
    #Store previously edited sites in the current row
    results[,i] <- results[,i-1]
    
    #Pick two random numbers
    rand_nums <- runif(2, 0, 1)
    
    #Find the next reaction time (exponential distribution)
    #The times vector stores cumulative times
    
    times[i] <- 1/sum(propensities) * log(1 / rand_nums[1]) + times[i-1]
    
    #Figure out which site was edited
    if(sum(propensities) > 0){
      
      site_probs <- propensities / sum(propensities)
      
      if(rand_nums[2] < site_probs[1]){
        results[1, i] <- 1
      } else {
        for(j in 2:length(site_probs)){
          if(rand_nums[2] < sum(site_probs[1:j]) & rand_nums[2] >= sum(site_probs[1:(j - 1)])){
            results[j, i] <- 1
          }
        }
      }
    }
    
    i <- i + 1
    
    #If everything edited, return
    if(times[i-1] < time & all(results[,i - 1] != 0)){
      return(list(results[,i-1], times[i-2], times[i-1], 0))
    }
  }
  
  #Returns the edited state and final editing time prior to the specified point
  return(list(results[,i-2], times[i-2], times[i-1], sum(propensities)))
}
                                                         
#simulates hypercascade editing over a given phylogeny
hypercascade_cell_tree <- function(subunits = 20,
                                   p = c(H01_for_nls[H01_for_nls$Layer == 1 &
                                                       H01_for_nls$Edit_Pattern == "000",]$Edit_Rate,
                                         H01_for_nls[H01_for_nls$Layer == 2 &
                                                       H01_for_nls$Edit_Pattern == "1001",]$Edit_Rate,
                                         H01_for_nls[H01_for_nls$Layer == 3 &
                                                       H01_for_nls$Edit_Pattern == "1011",]$Edit_Rate,
                                         H01_for_nls[H01_for_nls$Layer == 4 &
                                                       H01_for_nls$Edit_Pattern == "1111",]$Edit_Rate),
                                   edit_penalties = c(1 - coef(H01_model)[1],
                                                      1 - coef(H01_model)[2],
                                                      1 - coef(H01_model)[2],
                                                      1 - coef(H01_model)[3],
                                                      1 - coef(H01_model)[4]),
                                   integrations = 20,
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
    
    sim_res <- hypercascade_sim_endpt(subunits = subunits,
                                      time = end_time, 
                                      start_time = start_time,
                                      p = p, 
                                      start_point = start_point,
                                      edit_penalties = edit_penalties)
    barcode <- sim_res[[1]]
    
    return(list(barcode))
  }
  
  #Initialize lists to store the barcodes for the nodes
  barcode_results <- list()
  
  #Simulate the tree
  simulated_editing <- build_tree(curr_node = 0, 
                                  next_node = length(curr_tree$tip.label) + 1,
                               prev_barcs = lapply(1:integrations, function(x) {
                                 return(rep(0, 4*subunits - 6))
                               }))
  return(barcode_results)
}

#Modify to desired working directory
setwd("D:/OneDrive_Backup_220312/OneDrive - California Institute of Technology/Lab_Notebook/Experiment_by_numbers/e165/")

#Read in tree data and log
yule_trees <- read.nexus(file = 'truth.trees')
tree_log <- read.table(file = 'truth.log', header = T)

#Load in the edit rate data
load('edit_rates_and_penalties.RData')

#Loop over all trees (except burn in tree #1), simulate barcode editing, save results to .nex files
for(i in 2:101){
  #Pick a tree and get the root height
  curr_tree <- yule_trees[[i]]
  root_height <- tree_log$TreeHeight[i]
  
  #Traverse through the tree from root to tips
  root <- length(curr_tree$tip.label) + 1
  
  #Set edit rates and penalties for testing
  base_editing_rates <- c(H01_iPSC_for_nls[H01_iPSC_for_nls$Layer == 1 &
                                             H01_iPSC_for_nls$Edit_Pattern == "000",]$Edit_Rate,
                          H01_iPSC_for_nls[H01_iPSC_for_nls$Layer == 2 &
                                             H01_iPSC_for_nls$Edit_Pattern == "1001",]$Edit_Rate,
                          H01_iPSC_for_nls[H01_iPSC_for_nls$Layer == 3 &
                                             H01_iPSC_for_nls$Edit_Pattern == "1011",]$Edit_Rate,
                          H01_iPSC_for_nls[H01_iPSC_for_nls$Layer == 4 &
                                             H01_iPSC_for_nls$Edit_Pattern == "1111",]$Edit_Rate)
  edit_penalties <- c(1 - coef(H01_iPSC_model)[1],
                      1 - coef(H01_iPSC_model)[2],
                      1 - coef(H01_iPSC_model)[2],
                      1 - coef(H01_iPSC_model)[3],
                      1 - coef(H01_iPSC_model)[4])
  
  barcode_results <- hypercascade_cell_tree(subunits = 20, 
                                 p = base_editing_rates, 
                                 edit_penalties = edit_penalties, 
                                 integrations = 1, 
                                 tree_phylogeny = curr_tree, 
                                 root_height = root_height)
  
  barcode_results <- sapply(1:length(barcode_results), function(x) {
    barc <- unlist(barcode_results[[x]])
    barc <- paste(barc, collapse = "")
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
  write.nexus.data(barcode_results, file = paste("Tree", i, "_Barcodes.nex", sep = ""))
}

