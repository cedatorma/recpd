###
#Additional Functions for Calculating RecPD-Derived Measures.
###

##nrecpd_calc - Calculate the normalized recpd of a feature distribution
#Divide the recombination-adjusted feature PD by its associated Faith's PD value
#(vertical inheritance).

nrecpd_calc <- function(recpd, tree, tip){
  #Input arguments:
  #recpd - recpd for a given feature.
  #tree - a species phylogenetic tree.
  #tip - a named matrix providing the presence/absence states of the feature
  #among the tips of the corresponding phylogenetic tree.

  #If the feature has prevalence == 1, assign nRecPD = 1 (vertical inheritance)
  if(sum(tip) == 1){
    nrecpd <- 1
  }
  #Otherwise, calculate Faith's PD:
  else{
    pd <- picante::pd(tip, tree, include.root = FALSE)$PD/sum(tree$edge.length)

    #If RecPD > Faith's PD, assign nRecPD to 1.
    if(recpd > pd){
      nrecpd <- 1
    }
    #Otherwise divide RecPD by Faith's:
    else{
      nrecpd <- recpd / pd
    }
  }

  return(nrecpd)
}


#Span:
#For a given trait distribution with prevalence (P):
#Recalculate the PD by summing all tree branch lengths between tips possessing a trait and the tree root;
#Then divide by the sum of branch lengths for the maximum spanning locus distribution with the same P.

#get_max_span() - Finds The Maximum Spanning Trait Distributions (Based on Their Total Sum of Branch Lengths) For A Given Phylogenetic Tree:
#*Note that span is not adjusted for potential recombination.

get_max_span <- function(tree){

  #tip_dist():
  #A function for finding the maximum spanning feature distributions (stored as the sum of tree branches at a given feature prevalence):


  tip_dist <- function(np){
    np_tmp <- rev(np)[-1]
    #Get the tip:
    t <- rev(np)[1]

    #Check which ancestral nodes haven't been included (n_tally != 2) in the previously found maximum spanning feature distributions:
    n <- np_tmp[which(!np_tmp %in% n_found[which(n_tally == 2)])]

    #Get the maximum distance between the tip and its ancestral node, keeping in mind to only include unique branches which haven't been previously included in previous maximum spanning distributions  (n_tally != 1). Keep in mind that if there are no nodes, as when a tip is the nearest-neighbour of a tip belonging to a previously identified maximum spanning feature distribution, then use its parent node:

    #If the tip is a descendant of a previously identified maximum spanning feature distribution, then only check the distance between it and its ancestral node, i.e. the first node in its node path which has been seen only once before (n_tally  == 1).

    a <- which(n %in% n_found[which(n_tally == 1)])[1]

    #If this node doesn't exist, then use the most proximal(nearest the root) node of the node path (not required after initialization):
    #if(length(a) == 0) a <- length(n)

    #If such an ancestral node exists, then get its distance to the tip:
    #Otherwise, use the most distal node in the node path:
    max_d <- max(dn[t, n[1:a]])

    #Return the distance of the tip to its deepest ancestral node, along with its nodepath:
    cbind(max_d = max_d,
          n = list(as.character(n[1:a])))
  }


  ###
  #Initialization of variables:
  ###
  #Get the root-tip nodepaths for the tree:
  np <- ape::nodepath(tree)

  #Matrix of node distances:
  dn <- ape::dist.nodes(tree)


  #An array for storing the nodes/branches already found in previous maximum spanning distributions:
  #Numeric array of nodes for searching:
  n_found <- (ape::Ntip(tree) + 1) : (ape::Ntip(tree) + ape::Nnode(tree))

  #The corresponding array which keeps a tally of nodes found in maximum spanning feature distributions:
  n_tally <- array(rep(0, length(n_found)), dimnames = list(n_found))


  #An array for storing the corresponding maximum spanning distances of feature distributions over the entire range of prevalence:
  max_span <- rep(0, ape::Ntip(tree))

  #Copy np to a a temporary variable which will be successively shortened as
  #prevalence increases:

  np_tmp <- np


  #Initilization (feature prevalence == 1):

  #Find the tip having the maximum distance to the root node:
  max_d <- sort(dn[-n_found, n_found[1]], decreasing = TRUE)[1]

  #Extract the tip/np index and its sum of branch-lengths to the root:
  i <- as.numeric(names(max_d))
  max_d <- as.numeric(max_d)

  #Get the nodepath:
  n <- as.character(rev(np[[i]])[-1])

  #Use its nodepath to update n_tally:
  n_tally[n] <- 1

  #Use its distance to update max_span:
  max_span[1] <- max_d

  #Remove the corresponding nodepath from np_tmp:
  np_tmp <- np_tmp[-i]


  #Now iterate through the remaining range of feature prevalence: 2 .. Ntip(tree)

  for(p in 2:ape::Ntip(tree)){

    #Run the tip lineage distance calculation:
    res <- sapply(np_tmp, tip_dist)


    #Identify the tip/index of the np list with maximum distance to its ancestral node:
    i <- which(unlist(res[1,]) == max(unlist(res[1,])))[1]

    #Get the maximum distance:
    max_d <- unlist(res[1,i])

    #Get the nodepath leading to it:
    n <- unlist(res[2,i])

    #Increment n_tally using the nodepath of the given tip (only if its distance is the maximum of all other remaining tips):
    n_tally[n] <- n_tally[n] + 1

    #Remove the selected nodepath from np:
    np_tmp <- np_tmp[-i]

    #Update max_span by adding the maximum tip-ancestral node distance to that of the previously stored value:
    max_span[p] <- max_span[p - 1] + max_d

    #Continue until all tips/prevalence have been examined.
  }

  return(as.list(array(max_span, dimnames = list(1:Ntip(tree)))))
}

#The old function:
get_max_span_old <- function(tree){

  #Start with a full tree (P = Ntip),
  #Remove a single tip and its accompanying ancestral branch, and recalculate the sum of branch lengths.
  #Select the tip which results in the lowest reduction in the sum of branch lengths to find the maximum spanning distribution for P = Ntip -1.
  #Use this tip distribution to repeat the step above until reaching 2 tips.

  np <- ape::nodepath(tree)
  dn <- ape::dist.nodes(tree)

  b_tot <- sum(tree$edge.length)

  n <- 1:ape::Ntip(tree)

  tmp <- np

  s_res <- list()
  for(p in (ape::Ntip(tree) - 1):2){
    #How much does the total tree branch lengths change when removing the giving tip:
    #Also have to keep track of the previous nodes which have been removed.
    #However, to initialize the span distribution, this could go in two ways:
    #Start from P-1 presence tips, find the single tip which has the smallest distance to the root (only considering distance to ancestral node ignores if the tip is present in a highly divergent clade), and remove it along with its nodepath; use this procedure to iteratively find the next tip to remove.
    #Or take all the nodepaths of every tip,

    #Get the resulting nodepaths remaining in the tree after removing a given tip:
    j <- lapply(1:length(tmp), function(x) unique(unlist(tmp[-x])))

    #For each removed tip,
    #use the remaining nodepaths to calculate the total tree branch lengths:
    d <- unlist(
      lapply(j,
             function(x){
               k <- which(tree$edge[,1] %in% x & tree$edge[,2] %in% x)
               sum(tree$edge.length[k])
             }
      )
    )

    #Select the tip to remove which reduces the total tree branch length the least:
    rm <- order(d, decreasing=TRUE)[1]

    tmp <- tmp[-rm]

    s_res[[as.character(p)]] <- list(max_span=d[rm], tips=sapply(tmp, function(x) rev(x)[1]))
  }

  #Return the maximum spanning distribution and sum of branch lengths by prevalence:
  return(s_res)
}


#span_calc() - Calculates The Span For A Given Trait Distribution Found On A Species Phylogeny Based On
#The Sum Of Branch Lengths Joining All the Tree Tips In Which The Trait is Found, Normalized By The
#Potentially Maximum Spanning Trait Distribution With The Same Prevalence (calculated by get_max_span()).
span_calc <- function(tree, t, max_span){
  #Get the prevalence of the given trait/locus:
  p <- sum(t)

  #Get the nodepaths of the tree:
  np <- ape::nodepath(tree)

  #Calculate the sum of branch lengths for the given locus distribution
  #Note, presence/absence matrix should be sorted by the tip label order in the tree:
  t <- which(t != 0)

  if(length(t) != 0 & p != 1){

    n <- unique(unlist(np[t]))

    i <- which(tree$edge[,1] %in% n & tree$edge[,2] %in% n)

    br_p <- sum(tree$edge.length[i])

    return(br_p/max_span[[as.character(p)]])
    #return(br_p/max_span[[as.character(p)]]$max_span)

  }
  else{
    if(p == 1){
      return(0)
    }
    else if(p == ape::Ntip(tree)){
      return(1)
    }
  }
}

####
#Clustering:
####

#cluster_calc() - Calculate The Branch-State Clustering Coefficient For A Trait Distribution:
#For Each Internal Node of the Phylogenetic Tree (Excluding Root) Identified As A Gain State by RecPD:
#Calculate the Proportion of Descendant Branches (or Tips) Which Are Also Annotated As Present State,
#Then Sum These Together and Normalize by The Maximum Clustering Value Expected For A Trait Distribution With
#The Same Prevalence.

cluster_calc <- function(tree, t, ns, bs){
  #ns - node_states predicted using a given ancestral state reconstruction method (includes split nodes).
  #bs - branch_states.

  clust_b <- 0
  clust_n <- 0

  clust_coeff_n <- 0
  clust_coeff_b <- 0

  #Descendant tips of nodes:
  desc_t <- phangorn::Descendants(tree, type='tips')

  #Descendant children of nodes:
  desc_c <- phangorn::Descendants(tree, type='children')

  #Extract Internal Nodes:
  n <- as.numeric(names(ns[-(1:(ape::Ntip(tree)))]))

  #First step for calculating clustering: identify potential internal nodes which delimit presence state clades (not gains),
  #And calculate the proportion of presence state tips descended from them.
  #Any internal nodes not matching this criterion are assigned 0.

  #Two ways of selecting presence state clades:
  #1)Based on the proportion of descendant presence state tips from each internal node:
  #2)If the descendant branches of the node are assigned gain/presence states based on ancestor -> descendant node state relationships.

  #return(desc_c)
  #return(lapply(desc_t, function(x) mean(t[as.character(x)])))

  if(length(n) != 0){

    #Check if the descendant branches of a given node are presence state:
    db <- sapply(n, function(i){
      sum(bs[which(tree$edge[,1] %in% i & tree$edge[,2] %in% desc_c[[i]])])
    })

    #Calculate the proportion of presence state tips descended from each node in the gain lineage?
    #clust_b <- sapply(1:length(n), function(i) ifelse(db[i] == 2, sum(t[desc_t[[n[i]]]])/length(desc_t[[n[i]]]), 0))

    #Otherwise, just tally whether or not the descendant branches are both present states (=1), or split (=0)
    clust_b <- sapply(1:length(n), function(i) ifelse(db[i] == 2, 1, 0))
    names(clust_b) <- n

    #Note for node based clustering, include split state tips?
    clust_n <- sapply(n, function(i) ifelse(ns[as.character(i)] != 0, sum(t[desc_t[[i]]])/length(desc_t[[i]]), 0))

  }
  else{
    j <- rep(0, length(n))
  }

  #Calculate the two kinds of global clustering (by branch state gain lineages, vs. node state gain lineages...)
  #coefficients for comparison, and normalize by the maximum clustering possible,
  #assuming that each internal node of a clade leads to two presence state lineages: (Prevalence - 1) / (# Nodes):

  #New version, sum the proportion of presence state tip descendants per internal node.
  clust_coeff_n <- sum(clust_n) / ape::Nnode(tree)

  if(sum(t) != 1) clust_coeff_b <- mean(clust_b) / ((sum(t)-1)/(ape::Nnode(tree)))


  #Old version:
  #Note, that only internal nodes with only present tip descendants (== 1) are included.
  #This will ignore instances of nodes where descendant lossess have occured!

  #clust_coeff_n <- (length(which(clust_n == 1))/length(clust_n)) / ((sum(t)-1)/(Nnode(tree)))
  #clust_coeff_b <- (length(which(clust_b == 1))/length(clust_b)) / ((sum(t)-1)/(Nnode(tree)))

  ##There is also another thing to consider, which is partitioning, or the number of distinct clades the trait appears in.
  #For one needs to identify the earliest ancestral nodes (descendants from the root) where the trait is first gained, and the proportion of presence state tips descended from them.

  return(list(clust_coeff_n=clust_coeff_n,
              clust_coeff_b=clust_coeff_b,
              clust_n=clust_n,
              clust_b=clust_b))
}



#longevity_calc() - Calculate Longevity and Lability of Traits Based on Trait Lineages Identified by Ancestral State Reconstruction:
longevity_calc <- function(tree, ns){
  ts <- ns$tip_state
  bs <- ns$branch_state
  ns <- ns$node_state

  np <- ape::nodepath(tree)


  #Also add the root node if it has a presence branch joining it:
  r <- which(tree$edge[,1] == ape::Ntip(tree) + 1)
  r_br <- which(bs[r] == 1)

  if(length(r_br) != 0 & !as.character(ape::Ntip(tree) + 1) %in% names(ns)){
    ns <- append(1, ns)
    names(ns)[1] <- ape::Ntip(tree) + 1
  }

  n <- as.numeric(names(ns))

  #Get potential gain and loss nodes:
  gain_n <- which(ns == 1)
  loss_n <- which(ns == 0)

  #return(list(gain_n, loss_n))

  #Initialize output variables:
  longevity <- 0
  lability <- 0
  lability1 <- 0
  lability2 <- 0
  lability3 <- 0
  gain_times <- NULL
  loss_times <- NULL

  #Get the ancestral gain nodes for each presence tip:
  if(length(gain_n) != 0){
    ##Longevity:
    #For each presence tip, find the ancestral nodepath leading to their ancestral gain tip,
    #i.e. identify where the branch state changes from presence to absence, or if it leads directly to the root.

    #Then calculate the sum of branchlength distances  and normalize by the maximu of all tip to root node branch
    #lengths, to get the trait longevity distribution.
    #Like RecPD, except multiple branches can be counted, and they are taken for each individual locus,
    #to later come up with a mean, or distribution...

    dist <- ape::dist.nodes(tree)[1:ape::Ntip(tree), -(1:ape::Ntip(tree))]
    #dist <- dist.nodes(tree)

    #longevity:

    #For presence tips, get sum of consecutive presence state branches.
    #Define a function to do this (should only be done once, outside of this function, as a helper):
    path_state <- function(x){
      x <- rev(x)
      i <- sapply(1:(length(x)-1), function(i){
        which(tree$edge[,2] %in% x[i] & tree$edge[,1] %in% x[i+1])
      })

      b <- bs[i]
      j <- which(b[1:(length(b)-1)] != b[2:length(b)])[1]

      if(is.na(j)) j <- length(b)

      return(sum(tree$edge.length[i[1:j]])/max(dist[,1]))
      #return(list(j, bs, r_tr$edge.length[i], sum(r_tr$edge.length[i[1:j]])/max(dn) ))
    }

    longevity <- array(unlist(lapply(np[which(ts == 1)], path_state)), dimnames=list(which(ts == 1)))

    #The old longevity calculation with node states: will not work in the cases where loss nodes occur
    #on gain lineages.
    # longevity <- sapply(np[which(ts == 1)], function(x){
    #   i <- ns[as.character(rev(x))]
    #   i <- which(i == 1)[1]
    #   p <- x[1:i]
    #   #t <- as.character(rev(x)[1])
    #   #n <- as.character(rev(x)[i])
    #   e <- which(tree$edge[,1] %in% p & tree$edge[,2] %in% p)
    #   bl <- sum(tree$edge.length[e])
    #
    #   #return(list(tip=x, time=bl/max(dist[,1]), nodes=c(x[i], x[1]), edge=p))
    #   return(bl/max(dist[,1]))
    #   })




    #Gain times:
    #When do gain nodes appear, relative to the root... normalized by the maximum root-to-tip distance.

    gain_times <- sapply(n[gain_n], function(x){
      g <- ape::nodepath(tree, x, ape::Ntip(tree) + 1)
      #t <- as.character(rev(x)[1])
      #n <- as.character(rev(x)[i])
      e <- which(tree$edge[,1] %in% g & tree$edge[,2] %in% g)
      bl <- sum(tree$edge.length[e])
      return(bl/max(dist[,1]))
    })

    #return(longevity)

    #Lability:
    #Calculate lability as the normalized number of state changes occurring in the tree (excluding the root):
    #the sum of internal node gain and loss events, divided by the number of internal tree nodes where events could potentially occur;

    #If the root appears as a node-state change, remove it?
    lability1 <- ifelse(as.character(ape::Ntip(tree) + 1) %in% names(ns), length(ns) - 1, length(ns)) / (ape::Nnode(tree) - 1)


    #Or divided by the total number of internal nodes of each trait lineage where events could occur:
    ##To do this, need to identify the gain lineages (from gain_n) and loss lineages (from loss_n), and get the total number of unique
    ##descendant internal nodes from each...

    #Get gain lineage subtree nodes:
    s_n <- c(n[gain_n], n[loss_n])

    #If no loss lineages are identified, remove them:
    s_n <- s_n[which(!is.na(s_n))]

    #Get descendant nodes of each node state change lineage.
    #Note some gain/loss nodes may be descended from other ancestral gain/loss nodes, remove redundant nodes.
    st_n <- unlist(phangorn::Descendants(tree, s_n, type='all'))
    st_n <- unique(st_n[which(st_n >= ape::Ntip(tree) + 2)]) # +2 - remove root?

    #Distinguish the ancestral gain and loss lineage nodes from the descendant nodes of each gain or loss lineages,
    #those will be the actual node state change nodes:

    lability2 <- length(which(s_n %in% st_n))/length(st_n)


    #Or, just take the number of internal gain and loss event nodes and
    #divide by the number of gain and loss tips (excluding absences), descended from them.
    #Note, lability == 1 when all gain/loss states result directly from a gain/loss event, i.e. later occuring events,
    #while lower levels of lability indicates gains/losses occuring earlier in the tree.
    #and lability == 0 when all gain nodes are directly descended from an
    #ancestral gain event (have to subtract each category by 1 if only a single gain node exists).

    if(length(gain_n) == 1){
      lability3 <- ((length(which(ns == 1)) - 1) + length(which(ns == 0)))/ length(which(ts != -1))
    }
    else{
      lability3 <- ((length(which(ns == 1))) + length(which(ns == 0)))/ length(which(ts != -1))
    }


    if(length(loss_n) != 0){

      #If any internal loss nodes are found, also add them to the longevity array.
      #Get the sum of branch lengths between ancestral gain -> loss node lineages.
      #and normalize by max tip to root distance.

      #Get the branch state distances from loss nodes to ancestral gains:
      #i <- as.numeric(names(ns[[2]]$node_state)[which(ns[[2]]$node_state == 0)])

      l <- lapply(n[loss_n], function(x) ape::nodepath(phy=tree, from=ape::Ntip(tree) + 1, to=x))

      longevity <- append(longevity, array(unlist(lapply(l, path_state)), dimnames=list(n[loss_n])))


      # long_internal <- sapply(n[loss_n], function(x){
      #   i <- nodepath(tree, x, Ntip(tree)+1)
      #   g <- which(ns[as.character(i)] != 0)[1]
      #   #l <- which(ns[as.character(i)] == 0)[1]
      #   e <- which(tree$edge[,1] %in% i[1:g] & tree$edge[,2] %in% i[1:g])
      #   bl <- sum(tree$edge.length[e])
      #   return(bl/max(dist[,1]))
      #   #return(list(node=x, time=bl/max(dist[,1]), nodes=c(i[g], x), edge=i[1:g]))
      # })
      #
      # longevity <- append(longevity, long_internal)

      #Loss times:
      #Similar to gain times, but is the distance from the internal loss node relative to the root:
      loss_times <- sapply(n[loss_n], function(x){
        l <- ape::nodepath(tree, x, ape::Ntip(tree) + 1)
        #t <- as.character(rev(x)[1])
        #n <- as.character(rev(x)[i])
        e <- which(tree$edge[,1] %in% l & tree$edge[,2] %in% l)
        bl <- sum(tree$edge.length[e])
        return(bl/max(dist[,1]))
      })
    }
  }

  return(list(longevity=longevity,
              lability1=lability1,
              lability2=lability2,
              lability3=lability3,
              gain_times=gain_times,
              loss_times=loss_times))

}
