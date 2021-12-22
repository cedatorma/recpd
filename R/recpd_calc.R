#recpd_calc() - main wrapper function for performing RecPD calculations and
#branch state identification, and other cladistic and evolutionary distance
#based metrics.

#' Calculate the RecPD of a feature and other associated measures.
#'
#' See \href{https://www.biorxiv.org/content/10.1101/2021.10.01.462747v1}{Bundalovic-Torma, Guttman; (2021).}
#'
#' @param tree A phylogenetic tree, in phylo format.
#' @param tab A feature(s) presence/absence array, matrix, or data.frame:
#'   columns should be named according to the tip labels of the tree.
#' @param option Ancestral state reconstruction method to apply, either 'nn',
#'   'mpr', or 'ace'.
#' @param params Additional paramaters to pass to ancestral reconstruction
#'   methods.
#' @param calc logical. If TRUE, calculate additional measures: span,
#'   clustering, longevity, and lability.
#' @param prune logical. If TRUE and option == 'nn', will provide a conservative
#'   reconstruction of feature lineages.
#'
#' @return A named list including metrics and ancestral state annotations for
#'   tree nodes and branches for each feature in tab.
#' @export
#'
#' @examples
#' library(ape)
#' library(dplyr)
#'
#' ##Test case for a randomly generated tree (10 tips) and 10 randomly generated
#' ##feature presence/absence (1, 0) data.frame.
#'
#' #Generate a random phylogenetic tree:
#' n_tip <- 10
#' tree <- rtree(n_tip)
#'
#' #Generate 10 randomly distributed features (rows = features, columns =
#' #phylogenetic tree tips):
#' tab <- replicate(10,
#'                  sapply(10, function(x) sample(c(0,1), x, replace=TRUE)),
#'                  simplify='matrix') %>%
#'   t() %>%
#'   data.frame(check.names=FALSE) %>%
#'   rename_with(function(x) tree$tip.label[as.numeric(x)])
#'
#' #Note: columns should be labelled by tip.label of the tree. This also means that
#' #all of the tips in the tree should have an associated feature presence absence
#' #state. If tips lack this information, the tree should be parsed (see
#' #ape::keep.tip()).
#'
#'
#' #Run recpd_calc() without calculating additional measures (span, cluster,
#' #longevity, lability):
#' recpd_calc(tree, tab, calc=FALSE)
#'
#' #Run recpd_calc() with additional measures calculated:
#' recpd_calc(tree, tab, calc=TRUE)
#'
#' @author Cedoljub Bundalovic-Torma
#'
#' @references Insert reference with Rdpack.
#'
recpd_calc <- function(tree,
                       tab,
                       option = c('nn', 'mpr', 'ace'),
                       params = list(nn_select=c('min', 'median'), model_type=c('ER', 'ARD'), kappa=1, lik_cut=0.99),
                       calc = FALSE,
                       prune = TRUE){

  #Input arguments:

  #tree - the species tree tab - a presence absence matrix in dataframe format:
  #rows are a given trait/locus distribution, columns are the corresponding
  #species names matching the tip labels of the tree #Note: all species in the
  #presence/absence matrix (t) should be in the tree


  #method - the method for assigning internal node states - 3 options:
  #1) 'nn' : nearest-neighbours approach (by default)
  #2) 'mpr' : Most Parsimonious Reconstruction (using the MPR() function from the ape package)
  #3) 'ace' : Maximum Likelihood Ancestral Character Estimation (using the ace() function
  ##from the ape package).

  #Note that the ACE state ancestral likelihood measurements (marginal likelihoods)
  #are referred to as empirical bayesian posterior probabilities (see:
  #https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/ancestral-states-2/)

  #In the future: 4) a bayesian approach? anc.Bayes() from the phytools package?
  #This only works for continuous features.


  #Additional parameters are supplied using the params list:
  #If option = 'nn':
  #nn_select = choice of evolutionary distance cutoff to select nearest-neighbour nodes - 'min' for minimum pairwise distance, 'median' for median pairwise distance.
  ##If option = 'ace':
  #model_type -  can choose between rates equal ('ER') or all rates different ('ARD') models.
  #kappa - exponent for transforming tree branch-lengths.

  #calc - TRUE/FALSE if you want to perform topology or evolutionary metric calculations.



  #If calc == TRUE, get the maximum spanning distributions here, to prevent
  #rerunning multiple times for multiple trait distributions - majority of time,
  #is spent on this step, especially for larger species trees.
  if(calc == TRUE) max_span <- get_max_span(tree)


  #Output is a list, for now...:
  res <- list()

  #An alternative: results is an object (of class 'recpd' - ?) which stores a
  #data.frame of RecPD and associated measures (measures) calculated for each feature, with
  #the following attributes:

  #1) anc_new - list of arrays storing finalized trait lineage phylogenetic tree
  #tip/node/branch state mappings.

  #2) anc_old - list of arrays storing preliminary trait phylogenetic tree
  #tip/node/branch state mappings.

  #3) tree - the user-provided phylogenetic tree phylo object.

  #Initialize the measures data.frame:
  measures <- data.frame()

  #Initialize anc_old and anc_new lists:
  anc_new <- list()
  anc_old <- list()



  #If tab supplied is a named vector (nrow = NULL), convert it to a matrix format:
  if(is.null(nrow(tab))){
    tab <- matrix(tab, nrow=1, ncol=length(tab), dimnames=list(make.names(1), names(tab)))
  }


  #Iterate through the rows/features of tab,
  #perform ancestral state reconstruction and identify feature lineages:
  for(i in 1:nrow(tab)){
    #Make sure tips in the presence absence matrix are ordered according to the
    #tip ordering in the phylogenetic tree: Names of columns in t must match the
    #tree! Have to have some error checking to make sure that names match...
    #Resulting array should have strain names....
    t <- tab[i,tree$tip.label]

    #Initialize an array to store the states of terminal tip and internal nodes
    node_state <- NULL

    #Assign states based on nearest neighbour descendant tips of a given node

    #Option 1: Nearest neighbour descendants:
    if(option[1] == 'nn'){
      #Only generate the species-tree interal-node nearest-neighbour tip list once!
      if(i == 1){
        #Check if any input parameters have been changed:
        nn_select <- 'min'
        if(!is.null(params$nn_select[1])) nn_select <- params$nn_select[1]

        node_sp <- nn_nodes(tree, option=nn_select)
      }
      node_state <- node_state_nn(tree, t, node_sp)

    }
    #Option 2: Most Parsimonious Reconstruction (from ape package)
    else if(option[1] == 'mpr'){
      node_state <- node_state_mpr(tree, t)

    }
    #Option 3: Maximum Likelihood Ancestral Character Estimation *from ape package)
    else if(option[1] == 'ace'){

      #Update input parameters for ace, if provided
      if(i == 1){
        #defaul parameters:
        model_type <- 'ER'
        kappa <- 1
        lik_cut <- 0.99

        #Check if any of the model paramaters have been changed:
        if(!is.null(params$model_type[1])) model_type <- params$model_type[1]
        if(!is.null(params$kappa)) kappa <- params$kappa
        if(!is.null(params$lik_cut)) lik_cut <- params$lik_cut
      }

      #ace will only work if there is at least one presence or absence state on the tips:
      if(sum(t) < length(t)){

        node_state <- node_state_ace(tree, t, model_type, kappa, lik_cut)
      }
      else{
        node_state <- array(rep(1, ape::Ntip(tree) + ape::Nnode(tree)),
                            dimnames=list(1:(ape::Ntip(tree) + ape::Nnode(tree))))
      }
    }

    #Assign branch states based on ancestral node state reconstruction:
    branch_state <- branch_state_assign(tree, node_state)

    #Merge descendant and ancestral nsode states together to identify putative
    #gain and loss nodes: also output all of the updated state annotations
    #(gain, loss, absence) for tips, branches, and internal nodes, for both
    #plotting and additional calculations.
    node_state_merge <- node_ident5(tree, node_state, branch_state)

    #Update node, tip, and branch state annotations, pruning excessive internal
    #gain lineage branches resulting from 'nn':
    if(option[1] == 'nn' & prune == TRUE) node_state_merge <- prune_state(tree, node_state_merge)


    #Calculate RecPD and nRecPD:
    recpd <- sum(tree$edge.length*branch_state)/sum(tree$edge.length)

    nrecpd <- nrecpd_calc(recpd, tree, t)


    ##Perform other calculations?:

    span <- NA
    cluster <- NA
    evo <- NA
    longevity <- NA
    lability <- NA

    if(calc){

      #Topology:
      ##Span:
      span <- span_calc(tree, t, max_span)

      ##Clustering (has to take the entire set of annotated internal node states):
      cluster <- cluster_calc(tree, t, node_state, branch_state)[[2]]

      #Evolutionary:
      ##Longevity, lability, gain and loss times.
      evo <- longevity_calc(tree, node_state_merge)
      longevity <- stats::median(evo$longevity)
      lability <- evo$lability3

    }

    #Old:
    # res[[rownames(tab)[i]]] <- list('metrics' = list(recpd=recpd,
    #                                                  span=span,
    #                                                  cluster=cluster[[2]],
    #                                                  #Note, longevity is taken
    #                                                  #as the median of trait
    #                                                  #lineage branches, leading
    #                                                  #from gain nodes to both
    #                                                  #presence state tips and
    #                                                  #descendant loss nodes (if
    #                                                  #present).
    #                                                  longevity=stats::median(evo$longevity),
    #                                                  lability=evo$lability3),
    # 'node_state' = node_state_merge,
    # 'ns_old' = list(tip_state=node_state[(1:ape::Ntip(tree))],
    #                 branch_state=branch_state,
    #                 node_state=node_state[-(1:ape::Ntip(tree))]))


    #Update measures, anc_old, and anc_new:
    measures <- rbind.data.frame(measures,
                                 data.frame(feature = rownames(tab)[i],
                                            prevalence = sum(tab[i,]),
                                            recpd = signif(recpd, 3),
                                            nrecpd = signif(nrecpd, 3),
                                            span = signif(span, 3),
                                            cluster = signif(cluster, 3),
                                            longevity = signif(longevity, 3),
                                            lability = signif(lability, 3))
                                 )


    anc_old[[rownames(tab)[i]]] <- list(tip_state=node_state[(1:ape::Ntip(tree))],
                                        branch_state=branch_state,
                                        node_state=node_state[-(1:ape::Ntip(tree))])

    anc_new[[rownames(tab)[i]]] <- node_state_merge
  }

  #Remove uncalculated measures from the measures data.frame:
  if(calc == FALSE) measures <- measures[, 1:4]

  #Assign anc_old, anc_new, and tree as attributes of the measures data.frame:
  measures <- structure(measures,
                        anc_old = anc_old,
                        anc_new = anc_new,
                        tree = tree)

  return(measures)
}


####
##The following functions are helpers called by recpd_calc().
####

##nn_nodes() - Identify Nearest Neighbour Tips Descended From Each Internal Node
#of the Phylogenetic Tree. The Preliminary Step For Nearest-Neighbours Ancestral
#State Reconstruction (required by the node_state_nn()).

#For each internal node of the tree, identifying the subset of tips of each
#child-node descendant, then ompare the pairwise phylogenetic distances between
#each tip subset, and select those tips which are: 1) Equal to the minimum
#pairwise distance, or 2) <= median pairwise distance.

nn_nodes <- function(tree, option=c('min', 'median')){
  #Get the pairwise tip distance matrix
  #dist <- cophenetic.phylo(tree)

  #Get the distances between tips and nodes:
  dist <- ape::dist.nodes(tree)

  #Get internal node indices
  nodes <- ape::Ntip(tree) + seq(1, ape::Nnode(tree))
  node_sp <- list()

  for(n in nodes){
    #Get descendant children nodes for each internal node
    descend <- phangorn::Descendants(tree, n, type='children')

    i <- which(descend <= ape::Ntip(tree))

    #Check if descendant nodes are tips:
    #If nodes include 2 tips:
    if(length(i) > 1){

      t1 <- descend[1]
      t2 <- descend[2]

      close_t <- descend

    }
    #If only one descendant is a tip, get the descendant tips of the other child
    #internal node
    else if(length(i) == 1){
      #If nodes include 1 tip:
      t1 <- descend[i]

      #Get the descendant tips of the other child descendant node:
      t2 <- unlist(phangorn::Descendants(tree, descend[-i], type='tips'))

      #Get their pairwise evolutionary distances:
      d_sub <- dist[t1, t2]

      #Find the closest tip to the child descendant node:
      if(option[1] == 'min'){
        close_t <- c(t1, as.numeric(names(d_sub)[which(d_sub == min(d_sub))]))
      }
      #Or, find the pairs of tips <= the 25%-tile (?median) pairwise distance:
      else{
        close_t <- c(t1, as.numeric(names(d_sub)[which(d_sub <= stats::median(d_sub))]))
      }

    }
    #If both descendant nodes are internal-nodes, get the descendant tips for
    #each:
    else{
      #If nodes are internal, find the nearest neighbour tips between the
      #descendant nodes:
      t1 <- unlist(phangorn::Descendants(tree, descend[1], type='tips'))
      t2 <- unlist(phangorn::Descendants(tree, descend[2], type='tips'))

      #Get their parwise evolutionary distances:
      d_sub <- dist[t1, t2]

      #Take the closest related tips:
      if(1){

        #Find the closest tip pair to the child descendant node:
        if(option[1] == 'min'){
          m_ind <- which(d_sub == min(d_sub), arr.ind=TRUE)
          #close_t <- as.numeric(c(rownames(d_sub)[m_ind[1,1]], colnames(d_sub)[m_ind[1,2]]))
        }
        #Or, find the pairs of tips <= the 25%-tile (?median) pairwise distance:
        else{
          m_ind <- which(d_sub < stats::median(d_sub), arr.ind=TRUE)
          #close_t <- as.numeric(c(rownames(d_sub)[unique(m_ind[,1])], colnames(d_sub)[unique(m_ind[,2])]))

        }
        close_t <- as.numeric(c(rownames(d_sub)[unique(m_ind[,1])], colnames(d_sub)[unique(m_ind[,2])]))
      }
      else{
        #Or choose the closest related pairs of strains between descendant nodes.
        m_ind <- apply(d_sub, 1, function(x) which(x == min(x)))
        close_t <- as.numeric(c(rownames(d_sub)), unique(as.numeric(colnames(d_sub)[m_ind])))
      }
    }

    #Order the pair of tips based on increasing distance to the given node:
    #Sort tips by increasing distance to their parent node:
    d <- dist[n, close_t]

    #Store closest related descendant species of a node:
    node_sp[[as.character(n)]] <- tree$tip.label[close_t[order(d)]] #list(close_t, tree$tip.label[close_t[order(d)]], sort(d))
  }

  return(node_sp)
}


####
#Ancestral Reconstruction Functions of Internal Node Trait States, Used For
#Branch Selection to Calculate RecPD and Associated Metrics:
# node_state_nn() - The Default / Original Implemenatation Based on Nearest-Neighbours:
####

node_state_nn <- function(tree, t, node_sp){

  #Define a function to assign states to internal nodes:
  assign_state <- function(j, node_sp, presabs){
    close_sp <- as.character(unlist(node_sp[j]))

    #Assign split states based on the states of nearest neighbour descendant
    #nodes:

    #Version 1: look at the closest related pair of strains descended
    #from each node, ensuring that they come from different child nodes. If they
    #both have the same state, then assign the parent node to that state,
    #otherwise annotate the parent node as a split (potential gain/loss event).

    #Version 2: only look at the first, closest descendant of the parent node.

    s <- sum(presabs[1,close_sp])

    state <- 0
    #If all tips are present state:
    if(s == length(close_sp)){
      state <- 1
    }
    else if(s != 0){
      state <- 2
    }

    return(state)
  }

  #Presence / absence state of locus/trait
  presabs <- ifelse(t != 0, 1, 0)

  #return(presabs)

  #An array to store the states of internal nodes and tips
  node_state <- NULL

  n <- 1:(ape::Ntip(tree) + ape::Nnode(tree))

  for(i in n){
    if(i <= ape::Ntip(tree)){
      node_state <- append(node_state, presabs[1,tree$tip.label[i]])
    }
    else{
      node_state <- append(node_state, assign_state(as.character(i), node_sp, presabs))
    }
  }
  names(node_state) <- n
  return(node_state)

}



#node_state_mpr() - Ancestral Node State Reconstruction Using Most Parsimonious
#Reconstruction: *requires the MPR() function from the 'ape' package.
node_state_mpr <- function(tree, t){
  #Presence / absence state of locus/trait
  c <- ifelse(t[tree$tip.label] != 0, 1, 0)
  names(c) <- 1:ape::Ntip(tree)

  #Add the internal node labels to the tree, replacing


  #Use the Most Parsimonious Reconstruction (MPR) function in the ape package to
  #assign ancestral states for a given T3SE to all nodes in the tree. The MPR
  #algorithm will produce a matrix of nodes (in their order of appearance in the
  #tree?) and their corresponding lower and upper value character state
  #prediction for each node. Therefore, in the case of a binary discrete
  #character, e.g. presence = 1, absence = 0:

  #- Nodes with [1,1] indicate that the ancestral state of the is predicted to be present.
  #- Nodes with [0,0] indicate the ancestral state of the trait is predicted to be absent.
  #- Nodes with [0,1] indicate an ambiguous ancestral state.

  #Note: the tree must be unrooted.

  #The states assigned according to the appearance of parent nodes in the edge
  #matrix of the unrooted tree: Note, have to supply an unrooted tree but also
  #designate a tip to use as a root. Given that midpoint rooted and chronograms
  #(usually) do not have a tip to designate as a root, so the resulting
  #node-state mappings have to be remapped to the original tree to be visualized
  #properly. To make this process easier by preserving the original clade
  #structuring of the tree, employ the following trick:

  #Replace node labels (usually bootstrap supports) with node IDs,
  #So output from mpr() is intelligible.
  tree$node.label <- (ape::Ntip(tree)+1):(ape::Ntip(tree) + ape::Nnode(tree))

  #Add an artificial outgroup, attaching it to the root node of the original
  #tree, then root the new tree on it, which will preserve the :
  r_tree <- phytools::bind.tip(tree, 'outgroup', edge.length=1)
  r_tree <- ape::root(r_tree, 'outgroup')

  r <- 'outgroup'

  #Make all branch lengths the same length.
  r_tree$edge.length <- rep(1, length(r_tree$edge.length))

  #Add the outgroup to the original tip presence/absence state array:
  t2 <- c(t, `outgroup`=0)

  #Now perform MPR:
  parsim <- ape::MPR(ifelse(t2[r_tree$tip.label] != 0, 1, 0), ape::unroot(r_tree), r)

  #Assign node labels to the node parsimony state array according to their
  #appearance in the rooted tree (note: in some imported trees, the original
  #node labels assigned by MPR are the node bootstrap support values, but
  #generally they will be in the order of the node numbers of the tree).

  #Convert the table into an array, assigning states based on correspondence of
  #lower and upper state predictions.
  node_state <- NULL
  node_state <- apply(parsim, 1, function(x) if(x[1] == x[2]){return(x[1])}else{return(2)})

  #Assign internal node IDs to the node_state array
  names(node_state) <- rownames(parsim)

  #Also append the presence/absence states of the tips:
  node_state <- append(c, node_state)

  #Check if there are any presence state tips without ancestral nodes
  #annotated: will occur for trait distributions of low prevalence,
  #or if the trait appears in a clade where states are largely absent.
  #If there are, assign as a split (= 2).

  a <- phangorn::Ancestors(tree, which(c == 1), type='parent')

  i <- which(node_state[as.character(a)] == 0)

  if(length(i) != 0) node_state[as.character(a)[i]] <- 2

  return(node_state)
}



#node_state_ace() -  Ancestral Node State Reconstruction Using ML (Ancestral
#Character Estimation): *Requires the ace() function from the 'ape' package.
node_state_ace <- function(tree, t, model_type='ER', kappa=1, lik_cut=0.9){
  #lik_cut - a cutoff for selecting states based on the scaled state likelihood
  #of each node.

  #This will effect the corresponding node state assignments, if it is too high,
  #all nodes will be splits.

  #Presence / absence state of locus/trait
  c <- ifelse(t[tree$tip.label] != 0, 1, 0)
  names(c) <- 1:ape::Ntip(tree)

  #ace will throw an error if there are 0 or negative length branches.
  #Convert all of these edges to a minimum edge length. Let's say 1e-10
  if(length(which(tree$edge.length == 0)) != 0){
    tree$edge.length[which(tree$edge.length == 0)] <- 1e-9
  }

  #Results will be quite different if branch lengths are all set to 1, and RecPD
  #appears to be close to the other approaches
  #tree$edge.length <- rep(1,length(tree$edge.length))

  #Equal Rates assumed for state transitions?
  if(model_type == 'ER'){
    ace_res <- ape::ace(ifelse(t[tree$tip.label] != 0, 1, 0),
                        tree,
                        type='discrete',
                        method="ML",
                        model='ER',
                        marginal=TRUE,
                        kappa=kappa)
  }
  else if(model_type == 'ARD'){
    #Unequal Rates?

    #The 'model' argument allows you to specify the number of
    #different rate categories (integer) to describe transitions between states.
    #Rates are estimated using ML using the best fit-model.

    #Note that given the size of the tree and locus/state tip distribution, some
    #rates (e.g. gains) might not be estimated. Changing the scaling of the
    #branch-lengths (kappa) will result in a multiple rates being estimated,
    #although this means that the estimated rates have to be rescaled. This
    #rescaling also changes the resulting likelihood estimates as well.... The
    #intuition is that by raising the branch-lengths of a tree by a certain
    #exponent to approximate speciational change leads to better ancestral
    #character estimates (This applies also to the ER model).

    ace_res <- ape::ace(ifelse(t[tree$tip.label] != 0, 1, 0),
                        tree,
                        type='discrete',
                        method="ML",
                        model='ARD',
                        marginal=TRUE,
                        kappa=kappa)
  }


  #Ancestral state scaled likelihoods:
  anc <- ace_res$lik.anc

  #An array for storing node states
  node_state <- NULL

  #Assign presence, absence, or split state to nodes if their corresponding
  #state likelihoods are >= lik_cut, <= 1-lik_cut, or not, respectively.
  for(i in 1:nrow(anc)){
    if(anc[i,1] >= lik_cut){
      node_state <- append(node_state, 0)
    }
    else if(anc[i,1] <= 1-lik_cut){
      node_state <- append(node_state, 1)
    }
    else{
      node_state <- append(node_state, 2)
    }
  }

  #Nodes in sequence:
  n_ord <- seq(ape::Ntip(tree) + 1, ape::Ntip(tree) + ape::Nnode(tree))

  names(node_state) <- n_ord

  #Add the tip node states
  node_state <- append(c, node_state)

  return(node_state)
}





#node_ident5() - consolidation of previously generated phylogenetic tree
#ancestral state node state assignments.

#Assigns split state nodes to gain or loss events, if they correspond to
#the following ancestor -> descendant branch state changes:

#E.g. : branch state ancestor -> branch state descendant = node state:
#absent -> present = gain node,
#present -> absent = loss node.

#Also merges any loss or gain nodes of idential state which occurr within the
#same phylogenetic tree lineages.

node_ident5 <- function(tree, node_state, branch_state){
  #Get the tip states:
  tip_state <- node_state[1:ape::Ntip(tree)]

  #Get the internal node states:
  node_state <- node_state[-(1:ape::Ntip(tree))]

  #Get the edge matrix of the tree:
  edge <- tree$edge


  #Identify state change nodes based on changes of ancestral and descendant
  #branches of each node:

  #Go through the edge matrix, for each parent node, get its ancestral edge, and
  #check if the branch states change, if they do then keep that node as a state
  #change.

  #Define an array for storing updated node states:
  ns_merge <- array(rep(NA, ape::Nnode(tree)),
                    dimnames=list((ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree))))


  for(i in 1:nrow(edge)){
    #Get parent and child nodes for each edge:
    p <- edge[i,1]
    c <- edge[i,2]
    #Exclude the root:
    if(p != ape::Ntip(tree) + 1){
      #Get the ancestral node:
      a <- phangorn::Ancestors(tree, p, 'parent')

      #Find the ancestral edge:
      j <- which(edge[,1] == a & edge[,2] == p)

      #Compare branch states of the ancestral and decendant branches of the parent node:

      #If the state changes, then assign the node state to:

      #Node = Gain if ancestral branch is absent (0), and descendent branch is
      #present (1) = gain lingeage branches

      #Node = Loss if ancestral branch is present (1), and descendent branch is
      #absent (0) = loss lineage branches

      if(branch_state[i] == 1 & branch_state[j] == 0) ns_merge[as.character(p)] <- 1
      if(branch_state[i] == 0 & branch_state[j] == 1) ns_merge[as.character(p)] <- 0
    }
  }

  #Keep only the nodes where states change:
  rm <- which(is.na(ns_merge))
  ns_merge <- ns_merge[-rm]

  #Also check if the root node is assigned to an absence state. Then check
  #descendant nodepaths from the root which are all absent and assign to -1
  #state, i.e. absent lineages.

  r <- which(edge[,1] == ape::Ntip(tree) + 1)

  if(branch_state[r[1]] == 0 | branch_state[r[2]] == 0){
    #Get all nodepaths:
    np <- ape::nodepath(tree)

    #Traverse the nodepaths from the root, and identify any nodes on absent
    #lineages, assign them to -1:
    bs_abs <- list()

    for(i in 1:length(np)){
      x <- np[[i]]
      #Get the edge indices for consecutive nodes in each nodepath:
      j <- which(edge[,1] %in% x[-length(x)] & edge[,2] %in% x[-1])

      bs <- branch_state[j]
      #Find the consecutive absent state edge indices,
      #i.e. where does the branch state change to a gain?
      k <- which(bs[-length(bs)] != bs[-1])

      if(bs[1] == 0){
        if(length(k) != 0){
          bs_abs[[i]] <- j[1:k[1]]
        }
        else{
          bs_abs[[i]] <- j
        }
      }

      #Also assign any tips to absent states if they appear on a absent lineage
      #Descended directly from the root:
      if(bs[length(j)] == 0 & length(k) == 0){
        tip_state[i] <- -1
      }
    }

    #Assign the identified branches to absent state:
    branch_state[unique(unlist(bs_abs))] <- -1
  }

  #Add the root node if at least one descendant branch has been assigned as presence:
  r_d <- which(edge[,1] == ape::Ntip(tree) + 1)

  if(length(which(branch_state[r_d] == 1)) != 0){
    ns_merge <- append(array(1, dimnames=list(ape::Ntip(tree) + 1)), ns_merge)
  }

  return(list(tip_state=tip_state,
              branch_state=branch_state,
              node_state=ns_merge)
  )
}





#branch_state_assign() - branch state annotation function based on the ancestral
#state reconstruction of neighbouring internal nodes of the species phylogenetic
#tree.

branch_state_assign <- function(tree, node_state){
  branch_state <- NULL

  for(j in 1:nrow(tree$edge)){
    e <- tree$edge[j,]
    e1 <- e[1]
    e2 <- e[2]

    s1 <- node_state[e1]
    s2 <- node_state[e2]

    #If the branch ancestor node is the root, only include it if its node state
    #is present (=1)?

    #Call branches present if ancestor and child nodes are both assigned as
    #presence, split, or a combination of presence/split:

    #Call branches absent if ancestor and child nodes are both assigned as
    #absence, split->absence, or absence->split.

    #Only keep branches leading to split nodes if its parent is a gain, and vice
    #versa?

    if(s1 != 2 & s1 == s2){
      branch_state <- append(branch_state, as.numeric(s1))
    }
    # else if(s1 != 2 && s2 == 2){
    #   branch_state <- append(branch_state, 1)
    # }
    else if(s1 != 0 && s2 != 0){
      branch_state <- append(branch_state, 1)
    }
    else{
      branch_state <- append(branch_state, 0)
    }
  }

  return(branch_state)
}




#prune_state() - prune excessive internal gain lineage branches resulting from
#nearest-neighbours annotation:

#1) Get the gain nodes, and identify their children descendant branch states:
#2) Remove them if the descendant branches are split between loss and gain
#lineages;
#3) Stop at the last split state node before a gain (both children descendant
#branches present state), or if one of the children nodes is a presence-state
#tip;
#4) Also, reannotate any internal branches which are descended from an absence
#state ancestral node.

prune_state <- function(tree, states){

  #Input: the recpd_res[[i]]$node_state list

  #Extract the previously defined node, branch, and tip states:
  ns <- states$node_state
  bs <- states$branch_state
  ts <- states$tip_state

  edge <- tree$edge

  #Define variables for updated node, branch, and tip states:
  ns_new <- ns
  bs_new <- bs
  ts_new <- ts


  #Get the gain nodes:
  gain_n <- as.numeric(names(ns)[ns == 1])

  #Exclude any gain nodes which have only a pair of tip descendants:
  nt_d <- unlist(lapply(phangorn::Children(tree, gain_n),
                        function(x) length(which(x <= ape::Ntip(tree)))))

  #Condition for checking single gain nodes:
  if(length(gain_n) == 1 & length(nt_d) == 2) gain_n = NULL else gain_n <- gain_n[which(nt_d != 2)]

  #Need a condition for pruning lineages descended from the root...
  #For now, just remove them:
  #if((Ntip(tree) + 1) %in% gain_n) gain_n <- gain_n[-which(gain_n == Ntip(tree) + 1)]


  #An array for storing any modified/pruned gain nodes:
  rm_nodes <- NULL


  if(length(gain_n) != 0){
    for(n in gain_n){
      while(1){
        #Get the children nodes and branch states of the present gain node:
        c <- phangorn::Children(tree, n)
        c_branch <- which(edge[,1] %in% n & edge[,2] %in% c)

        #Find which child node exists on the gain lineage based on its ancestral
        #branch state:
        branch_pres <- which(bs[c_branch] == 1)

        #Get the ancestral node and its branch state leading to the current
        #node:
        a <- phangorn::Ancestors(tree, n, 'parent')
        a_branch <- which(edge[,1] %in% a & edge[,2] %in% n)


        #If the current node is a split (one child node is connected with a
        #present, the other absent), then remove the given node, and reassign
        #its branch state based on the anestral branch-state leading to it:
        if(length(branch_pres) == 1){
          #Remove the current node:
          rm_n <- which(names(ns_new) == as.character(n))
          ns_new <- ns_new[-rm_n]

          #Store the initially annotated gain nodes which will be removed to
          #later update their descendant branches:
          if(n %in% gain_n) rm_nodes <- append(rm_nodes, n)

          #Reassign ancestral branch states:
          #If the current node is the root, then use the node state of the other
          #child node:
          if(n == ape::Ntip(tree) + 1){
            bs_new[c_branch[branch_pres]] <- bs_new[c_branch[-branch_pres]]

          }
          #Otherwise use the branch-state leading to its ancestral node:
          else{
            bs_new[c_branch[branch_pres]] <- bs_new[a_branch]

          }


          #Update the current node to the present state child:
          n <- edge[c_branch[branch_pres], 2]
        }
        #If the children descendant nodes are both on the present state
        #lineages, then assign the ancestral node as a gain and its branch as
        #present (if not root node):
        else{
          if(n != (ape::Ntip(tree) + 1)){
            ns_new[as.character(a)] <- 1
            bs_new[a_branch] <- 1
            ns_new <- ns_new[order(as.numeric(names(ns_new)))]
          }
          break
        }
      }
    }


    #Now that excess present state branches have been trimmed, check if any
    #previously inferred loss lineages should be changed to absent state.

    #Define a function to update branch state annotations on the tree by
    #checking all loss state branches and changing them to absent if they are
    #directly descended from an absent state node:
    state_change <- function(n){
      c <- phangorn::Children(tree, n)

      if(length(c) != 0){
        e <- which(edge[,1] %in% n & edge[,2] %in% c)

        e <- e[which(bs_new[e] != 1)]

        if(length(e) != 0){
          #Update branch states
          bs_new[e] <<- -1

          #Also update tip states if descendants are tips:
          t <- which(edge[e,2] <= ape::Ntip(tree))

          #Note, have convert the edge-tip indicies to their appropriate tip numbers;
          if(length(t) != 0) ts_new[as.character(edge[e[t],2])] <<- -1

          for(n in c){
            state_change(n)
          }
          #print(c)
        }
      }
    }

    #First check if the root state is absent:
    r_state <- bs_new[which(edge[,1] == ape::Ntip(tree) + 1)]

    #Reassign the branches for the pruned internal gain-state nodes previously
    #identified to the corresponding state of their ancestral branches:
    if(length(which(r_state == -1)) != 0){
      for(n in rm_nodes){
        state_change(n)
      }
    }
  }

  return(list(tip_state = ts_new,
              node_state = ns_new,
              branch_state = bs_new)
  )
}
