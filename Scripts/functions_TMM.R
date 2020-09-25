
####################################Functions used for our study##################################

# Function to simulate gene trees for our study 
Gene_trees <- function(nsim, ntips){
  
  # function simualtes trees under the constant rate birth-death process
  trees <- sim.bd.taxa(n=ntips, 1, numbsim = nsim,  lambda=0.4, mu=0.1, complete=FALSE)
  return(trees)
}

##Function created to take into account heterogeneous substitution rates of the speciation and duplication events
## Function to create substitution rate tree from time tree
simulate.rates_heterogeneous <-function (tree, rate, noise) 
{
  #rate <- params$rate
  #noise <- params$noise
  rate <- rate[colnames(tree$mapped.edge)]
  data.matrix <- get.tree.data.matrix(tree)
  branch.rates <- rep(0, times = nrow(tree$edge))
  branch.rates <- ifelse(tree$mapped.edge[,"D"]>0, 
                         abs(branch.rates + rate_vector[["D"]]+ rnorm(length(tree$edge.length), mean = 0, sd = noise)),
                         abs(branch.rates + rate_vector[["S"]]+ rnorm(length(tree$edge.length),mean = 0, sd = noise)))
  data.matrix[, 5] <- branch.rates
  data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
  tree$edge.length <- data.matrix[, 6]
  tree$mapped.edge <- NULL
  tree$maps <- NULL
  res <- list(tree, data.matrix)
  names(res) <- c("phylogram", "tree.data.matrix")
  class(res) <- "ratesim"
  return(res)
}

## Function used for heterogeneous rates of trait (in our case tissue specificity measure tau with bound values of 0 and 1) evolution for two types of events
#sim.rates_bound <- function (tree, sig2, anc=0, nsim = 1, internal = F, plot = F) {
 sim.rates_bound <- function (tree, sig2, anc=0, nsim = 1, internal = F, plot = F) {
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (is.null(tree$mapped.edge)) {
    message("tree does not contain a mapped discrete character history, using fastBM")
    X <- fastBM(tree, sig2 = sig2[1], anc = anc, nsim = nsim, 
                internal = internal)
  }
  else {
    if (is.null(names(sig2))) {
      message("names absent from sig2: assuming same order as $mapped.edge")
      if (length(sig2) == ncol(tree$mapped.edge)) 
        names(sig2) <- colnames(tree$mapped.edge)
      else stop("the number of elements in sig2 should match the number of rows in mapped.edge")
    }
    sig2 <- sig2[colnames(tree$mapped.edge)]
    edge.length <- rep(0, nrow(tree$edge))
    for (i in 1:ncol(tree$mapped.edge)) edge.length <- edge.length + 
      sig2[i] * tree$mapped.edge[, i]
    names(edge.length) <- NULL
    tree_new <- list(Nnode = tree$Nnode, edge = tree$edge, tip.label = tree$tip.label, 
                     edge.length = edge.length)
    class(tree_new) <- "phylo"
    if (plot) 
      plot(tree_new)
    X <- fastBM(tree_new, anc = anc, nsim = nsim, internal = internal, bounds=c(0,1))
  }
  X
}

## internal function does BM simulation
## written by Liam J. Revell 2011, 2013
simBM<-function(tree,a,mu,sig2,bounds,internal,nsim){
  if(bounds[2]<bounds[1]){
    warning("bounds[2] must be > bounds[1]. Simulating without bounds.")
    bounds<-c(-Inf,Inf)
  }
  if(bounds[1]==-Inf&&bounds[2]==Inf) no.bounds=TRUE
  else no.bounds=FALSE
  if(a<bounds[1]||a>bounds[2]){
    warning("a must be bounds[1]<a<bounds[2]. Setting a to midpoint of bounds.")
    a<-bounds[1]+(bounds[2]-bounds[1])/2
  }
  if(sig2<0){
    warning("sig2 must be > 0.  Setting sig2 to 1.0.")
    sig2=1.0
  }
  # function for reflection off bounds
  reflect<-function(yy,bounds){
    while(yy<bounds[1]||yy>bounds[2]){
      if(yy<bounds[1]) yy<-2*bounds[1]-yy
      if(yy>bounds[2]) yy<-2*bounds[2]-yy
    }
    return(yy)
  }
  # how many species?
  n<-length(tree$tip)
  # first simulate changes along each branch
  x<-matrix(data=rnorm(n=length(tree$edge.length)*nsim,mean=rep(mu*tree$edge.length,nsim),sd=rep(sqrt(sig2*tree$edge.length),nsim)),length(tree$edge.length),nsim)
  # now add them up
  y<-array(0,dim=c(nrow(tree$edge),ncol(tree$edge),nsim))
  for(i in 1:nrow(x)){
    if(tree$edge[i,1]==(n+1))
      y[i,1,]<-a
    else
      y[i,1,]<-y[match(tree$edge[i,1],tree$edge[,2]),2,]
    
    y[i,2,]<-y[i,1,]+x[i,]
    if(!no.bounds) y[i,2,]<-apply(as.matrix(y[i,2,]),1,function(yy) reflect(yy,bounds))
  }
  rm(x); x<-matrix(data=rbind(y[1,1,],as.matrix(y[,2,])),length(tree$edge.length)+1,nsim)
  rownames(x)<-c(n+1,tree$edge[,2])
  x<-as.matrix(x[as.character(1:(n+tree$Nnode)),])
  rownames(x)[1:n]<-tree$tip.label
  # return simulated data
  if(internal==TRUE)
    return(x[1:nrow(x),]) # include internal nodes
  else
    return(x[1:length(tree$tip.label),]) # tip nodes only
}
## This function aims to assign proportion of duplication events and also paints duplication branches of time tree symmetrically and asymmetrically
## We chose duplication and speciation events randomly for each tree
Initial_timetree_pre <- function(tree1,p_dup)
{

  ##Initializing the list
  tree_painted_symmetric=vector('list',length(tree1))
  tree_painted_asymmetric=vector('list',length(tree1))
  dup_nodes_paint_tree= vector('list',length(tree1))
  spe_nodes_paint_tree= vector('list',length(tree1))
  dup_edges_symmetric= vector('list',length(tree1))
  dup_edges_asymmetric= vector('list',length(tree1))

  for(i in 1:length(tree1))
  {
    print(i)
    tree<-tree1[[i]]
    events <- rep("speciation", times=tree$Nnode)
    dup_nodes <- sample(x= c(1:length(events)), size=round(tree$Nnode*p_dup), replace=F)
    events [dup_nodes ] <- "duplication"
    names(events) <- c((ntips+1) : (ntips+tree$Nnode))

    # Painting both the the edges or branches that are duplication
    dup_nodes_paint <- as.numeric(names( events [ which(events =="duplication") ]) )
    dup_nodes_paint_tree[[i]] <- dup_nodes_paint
    spe_nodes_paint <- as.numeric(names( events [ which(events =="speciation") ]) )
    spe_nodes_paint_tree[[i]] <- spe_nodes_paint
    tree_edges<- unique(tree$edge[which(tree$edge[,1] %in% dup_nodes_paint), c(1,2)])
    tree_edges<-tree_edges[order(tree_edges[,1],tree_edges[,2]),]
    dup_edges <- unique(tree$edge[ which(tree$edge[,1] %in% dup_nodes_paint), 2])
    dup_edges_symmetric[[i]]<- dup_edges

    ##Painting initial time tree with certain proportions of duplication nodes
    tree_painted_symmetric[[i]]<- paintBranches (tree, edge=dup_edges, "D", anc.state="S")

    ##Randomly choosing one of the two duplicate branches for painting
    dup_edge_paint<- vector()
    for(x in 1:length(dup_nodes_paint))
    {
      dup_edge_selected <- sample(tree$edge[ which(tree$edge[,1] %in% dup_nodes_paint[x]), 2], 1, replace = F)
      dup_edge_paint<-append(dup_edge_paint,dup_edge_selected)
    }
    dup_edges_asymmetric[[i]]<- dup_edge_paint
    tree_painted_asymmetric[[i]] <- paintBranches (tree, edge=dup_edge_paint, "D", anc.state="S")
  }
  return(list(tree_painted_symmetric, tree_painted_asymmetric, dup_nodes_paint_tree, spe_nodes_paint_tree, dup_edges_symmetric, dup_edges_asymmetric))
}




## This function aims to assign proportion of duplication events and also paints duplication branches of time tree symmetrically and asymmetrically
## However we need to fix some speciation nodes which can further be used as focal speciation nodes for time calibration 
## Since our tree has 100 tips, so we have 99 internal nodes
## Therefore we fix 50% of 0.2 speciation event i.e. 10 nodes which will be present as focal speciation nodes in all simulated trees
Initial_timetree <- function(tree1,p_dup)
{
  
  ##Initializing the list
  tree_painted_symmetric <-vector('list',length(tree1)) 
  tree_painted_asymmetric<-vector('list',length(tree1))
  dup_nodes_paint_tree<- vector('list',length(tree1))
  spe_nodes_paint_tree <-vector('list',length(tree1))
  dup_edges_symmetric<- vector('list',length(tree1))
  dup_edges_asymmetric<- vector('list',length(tree1))
  
  for(i in 1:length(tree1))
  {
    print(i)
    tree<-tree1[[i]]
    events <- rep("speciation", times=tree$Nnode)
   
    ## All internal nodes
    all_nodes<- c((ntips+1) : (ntips+tree$Nnode))
    names(events)<-all_nodes
    #fixed_focal_spe_nodes<- sample(x = all_nodes,size = 10, replace = F ) ##Once selected using i=1, then used the same speciation time point for calibrating all the trees
    fixed_focal_spe_nodes<-c(102, 141, 153, 111, 129, 175, 194, 164, 150, 143)
    #events [fixed_focal_spe_nodes ] <- "speciation"
    nodes_excluding_focal_nodes<-as.numeric(all_nodes[!(all_nodes%in%fixed_focal_spe_nodes)])
    dup_nodes <- sample(x= nodes_excluding_focal_nodes, size=round(tree$Nnode*p_dup), replace=F) 
    events [which(names(events)%in% dup_nodes)] <- "duplication"
    
    # Painting both the the edges or branches that are duplication 
    dup_nodes_paint <- as.numeric(names( events [ which(events =="duplication") ]) )
    dup_nodes_paint_tree[[i]] <- dup_nodes_paint
    spe_nodes_paint <- as.numeric(names( events [ which(events =="speciation") ]) )
    spe_nodes_paint_tree[[i]] <- spe_nodes_paint
    tree_edges<- unique(tree$edge[which(tree$edge[,1] %in% dup_nodes_paint), c(1,2)])
    tree_edges<-tree_edges[order(tree_edges[,1],tree_edges[,2]),]
    dup_edges <- unique(tree$edge[ which(tree$edge[,1] %in% dup_nodes_paint), 2])
    dup_edges_symmetric[[i]]<- dup_edges
    
    ##Painting initial time tree with certain proportions of duplication nodes
    tree_painted_symmetric[[i]]<- paintBranches (tree, edge=dup_edges, "D", anc.state="S")
    
    ##Randomly choosing one of the two duplicate branches for painting 
    dup_edge_paint<- vector()
    for(x in 1:length(dup_nodes_paint))
    {
      dup_edge_selected <- sample(tree$edge[ which(tree$edge[,1] %in% dup_nodes_paint[x]), 2], 1, replace = F)
      dup_edge_paint<-append(dup_edge_paint,dup_edge_selected)
    }
    dup_edges_asymmetric[[i]]<- dup_edge_paint
    tree_painted_asymmetric[[i]] <- paintBranches (tree, edge=dup_edge_paint, "D", anc.state="S")
  }
  return(list(tree_painted_symmetric, tree_painted_asymmetric, dup_nodes_paint_tree, spe_nodes_paint_tree, dup_edges_symmetric, dup_edges_asymmetric))
}


## Function to provide different substitution rates for the two different events base on painted time tree 
##Providing different substitution rates for the two different events (speciation and duplication) base on painted time tree
Substitution_tree <- function(tree,rate,noise)
{
  tree_name_sym<-tree[[1]]
  tree_name_asym<-tree[[2]]
  tree_name_dup_edge_sym<-tree[[5]]
  tree_name_dup_edge_asym<-tree[[6]]
  tree_sim_painted_sym=vector('list',length(tree_name_sym))
  tree_sim_painted_asym=vector('list',length(tree_name_asym))
  
  for(i in 1:length(tree[[1]]))
  {
    print(i)
    Timetree_sym<- tree_name_sym[[i]]
    Timetree_asym<- tree_name_asym[[i]]
    dup_edges_sym<- tree_name_dup_edge_sym[[i]]
    dup_edges_asym <- tree_name_dup_edge_asym[[i]]
    
    if( "phylo" %in% class( Timetree_sym))
    {
      tree_rate_sim<-simulate.rates_heterogeneous(Timetree_sym, rate=rate_vector, noise)
      tree_sim <- tree_rate_sim[[1]]
      
      ##Painting substitution rate tree
      tree_sim_painted_sym[[i]] <- paintBranches(tree_sim, edge=dup_edges_sym, "D", anc.state="S") 
    }
    if( "phylo" %in% class( Timetree_asym))
    {
      tree_rate_sim_new<-simulate.rates_heterogeneous(Timetree_asym, rate=rate_vector, noise)
      tree_sim_new <- tree_rate_sim_new[[1]]
      tree_sim_painted_asym[[i]] <- paintBranches (tree_sim_new, edge=dup_edges_asym, "D", anc.state="S") 
    }
  }
  return(list(tree_sim_painted_sym, tree_sim_painted_asym))
}

## Calibrating the substitution rate trees
## We used speciation nodes as the focal speciation node and all the trees are calibrated based on the focal speciation nodes
Calibration_tree <- function(Timetree,Substree) ##here the tree is substitution rate tree
{
  sym_calibrated<-0
  asym_calibrated<-0
  calibrated_painted_sym<-vector('list',length(Timetree[[1]]))
  calibrated_painted_asym<-vector('list',length(Timetree[[2]]))
  

  for(i in 1:length(Substree[[1]]))
  {
    print(i)
    Timetree_sym<- Timetree[[1]][[i]]
    Timetree_asym<- Timetree[[2]][[i]]
    spe_nodes<- Timetree[[4]][[i]]
    dup_edges_sym<- Timetree[[5]][[i]]
    dup_edges_asym <- Timetree[[6]][[i]]
    SubsTree_sym <- Substree[[1]][[i]]
    SubsTree_asym <- Substree[[2]][[i]]
    
    ##getting age of all nodes
    #age_all=(round(branching.times(Timetree_sym),5)) 
    age_all<-branching.times(Timetree_sym) 
    
    
    ##Calibrating time of the substitution tree to obtain psuedo-time tree (Using chronos function to calibrate times for focal speciation nodes)
    ##But if the tree starts with duplication, then chronos fails to calibrate the root node; it often changes the scale drastically
    root_age<-max(age_all)
    root_node<-as.integer(names(age_all[which(age_all%in%root_age)]))
    fixed_focal_spe_nodes<-c(102, 141, 153, 111, 129, 175, 194, 164, 150, 143)
    spe_nodes_new<- NULL
    #spe_nodes_new<-unique(append(root_node,spe_nodes)) #we plan to incorporate the root age (although in real life this root age we never know!!!) for our tree to maintain the scale
    spe_nodes_new<-sort(fixed_focal_spe_nodes)
    speciation_age <- as.numeric(age_all[(names(age_all) %in% spe_nodes_new)])
    #speciation_age <- age_all[(names(age_all) %in% spe_nodes)]
    
    focal_calibration_times<-data.frame(node=spe_nodes_new,age=speciation_age)
    calibration_matrix <- data.frame(node=spe_nodes_new,age.min=speciation_age, age.max=speciation_age, soft.bound=NA)
    
    class(SubsTree_sym) <-"phylo"
    SubsTree_sym$mapped.edge <- NULL
    SubsTree_sym$maps <- NULL
    class(SubsTree_asym) <-"phylo"
    SubsTree_asym$mapped.edge <- NULL
    SubsTree_asym$maps <- NULL
    ctree_sym <- try(ape::chronos(SubsTree_sym,calibration = calibration_matrix,  model="correlated" ))
    ctree_asym <- try(ape::chronos(SubsTree_asym,calibration = calibration_matrix,  model="correlated" ))
    
    ##Trees those are not time calibrated can not be used further painting "duplication" and "speciation" branch.
    ##To avoid error due to non calibrated tree we did the following
    if( "phylo" %in% class(ctree_sym)) {class( ctree_sym ) <- "phylo"}
    else{ctree_sym<-NA}
    if( "phylo" %in% class(ctree_asym)) {class( ctree_asym ) <- "phylo"}
    else{ctree_asym<-NA}
    #print (calibrated_tree)
    if(!is.na(ctree_sym ))
    {
      if( "phylo" %in% class( ctree_sym) ) 
      {
        calibrated_painted_sym[[i]]<- paintBranches (ctree_sym, edge=dup_edges_sym, "D", anc.state="S")
        sym_calibrated<-sym_calibrated+1
        print(sym_calibrated) 
      }
    }
    if(!is.na(ctree_asym ))
    {
      if( "phylo" %in% class( ctree_asym) ) 
      {
        calibrated_painted_asym[[i]]<- paintBranches (ctree_asym, edge=dup_edges_asym, "D", anc.state="S")
        asym_calibrated<-asym_calibrated+1
        print(asym_calibrated) 
      }
    }
  }
  return(list(calibrated_painted_sym, calibrated_painted_asym))
}
## Calibrating the substitution rate trees
## We used speciation nodes as the focal speciation node and all the trees are calibrated based on the focal speciation nodes
Calibration_tree2 <- function(Timetree,Substree) ##here the tree is substitution rate tree
{
  sym_calibrated<-0
  asym_calibrated<-0
  calibrated_painted_sym<-vector('list',length(Timetree[[1]]))
  calibrated_painted_asym<-vector('list',length(Timetree[[2]]))
  
  
  for(i in 1:length(Substree[[1]]))
  {
    print(i)
    Timetree_sym<- Timetree[[1]][[i]]
    Timetree_asym<- Timetree[[2]][[i]]
    spe_nodes<- Timetree[[4]][[i]]
    dup_edges_sym<- Timetree[[5]][[i]]
    dup_edges_asym <- Timetree[[6]][[i]]
    SubsTree_sym <- Substree[[1]][[i]]
    SubsTree_asym <- Substree[[2]][[i]]
    
    ##getting age of all nodes
    #age_all=(round(branching.times(Timetree_sym),5)) 
    age_all<-branching.times(Timetree_sym) 
    
    
    ##Calibrating time of the substitution tree to obtain psuedo-time tree (Using chronos function to calibrate times for focal speciation nodes)
    ##But if the tree starts with duplication, then chronos fails to calibrate the root node; it often changes the scale drastically
    root_age<-max(age_all)
    root_node<-as.integer(names(age_all[which(age_all%in%root_age)]))
    
    ## Using root node as well for calibration
    #fixed_focal_spe_nodes<- sample(x = spe_nodes,size = 10, replace = F ) ##Once selected using i=1, then used the same speciation time point for calibrating all the tree
    #fixed_focal_spe_nodes<- spe_nodes
    #fixed_focal_spe_nodes<-c(102, 141, 153, 111, 129, 175, 194, 164, 150, 143)
    spe_nodes_new<- NULL
    spe_nodes_new<-unique(append(root_node,spe_nodes)) #we plan to incorporate the root age (although in real life this root age we never know!!!) for our tree to maintain the scale
    #spe_nodes_new<-sort(fixed_focal_spe_nodes)
    speciation_age <- as.numeric(age_all[(names(age_all) %in% spe_nodes_new)])
    #speciation_age <- age_all[(names(age_all) %in% spe_nodes)]
    
    focal_calibration_times<-data.frame(node=spe_nodes_new,age=speciation_age)
    calibration_matrix <- data.frame(node=spe_nodes_new,age.min=speciation_age, age.max=speciation_age, soft.bound=NA)
    
    class(SubsTree_sym) <-"phylo"
    SubsTree_sym$mapped.edge <- NULL
    SubsTree_sym$maps <- NULL
    class(SubsTree_asym) <-"phylo"
    SubsTree_asym$mapped.edge <- NULL
    SubsTree_asym$maps <- NULL
    ctree_sym <- try(ape::chronos(SubsTree_sym,calibration = calibration_matrix,  model="correlated" ))
    ctree_asym <- try(ape::chronos(SubsTree_asym,calibration = calibration_matrix,  model="correlated" ))
    
    ##Trees those are not time calibrated can not be used further painting "duplication" and "speciation" branch.
    ##To avoid error due to non calibrated tree we did the following
    if( "phylo" %in% class(ctree_sym)) {class( ctree_sym ) <- "phylo"}
    else{ctree_sym<-NA}
    if( "phylo" %in% class(ctree_asym)) {class( ctree_asym ) <- "phylo"}
    else{ctree_asym<-NA}
    #print (calibrated_tree)
    if(!is.na(ctree_sym ))
    {
      if( "phylo" %in% class( ctree_sym) ) 
      {
        calibrated_painted_sym[[i]]<- paintBranches (ctree_sym, edge=dup_edges_sym, "D", anc.state="S")
        sym_calibrated<-sym_calibrated+1
        print(sym_calibrated) 
      }
    }
    if(!is.na(ctree_asym ))
    {
      if( "phylo" %in% class( ctree_asym) ) 
      {
        calibrated_painted_asym[[i]]<- paintBranches (ctree_asym, edge=dup_edges_asym, "D", anc.state="S")
        asym_calibrated<-asym_calibrated+1
        print(asym_calibrated) 
      }
    }
  }
  return(list(calibrated_painted_sym, calibrated_painted_asym))
}

##This function is created to work for progressive jump model 
##replace i_node_x by i_node for some data
##This function is created to work for progressive jump model
recursive.speciation <- function(tree,new_tree,i_node,tip_label_new,trait_jumped,new_trait_i,sig2,events,ntips)
{
  i_node_des<-as.vector(tree$edge[ which(tree$edge[,1] == names(new_trait_i)) ,2]) 
  for(z in 1:length(i_node_des))
  {
    i_node_n = i_node_des[z]
    
    # taking the branch length of the branch between the nodes and assign to a new single-branch fake tree, for the BM calculation
    new_tree$edge.length<-NULL
    new_tree$edge.length <- tree$edge.length[which(tree$edge[,1] == i_node & tree$edge[,2]==i_node_n)]
    # print(new_tree$edge.length) 
    
    # starting trait from the simulated at the beginning, simulating with it's new ancestral trait value 
    new_trait <- trait_jumped[which(names(trait_jumped)==i_node)]
    
    # Estimate the new trait based on a BM process, the branch length of the randomly selected branch, and the jumped ancestral state i_trait
    new_trait_s <- simBM(new_tree, a=new_trait,mu=0, sig2=sig2, bounds=c(0,1),internal=F, nsim=1)
    names(new_trait_s) <- i_node_n
    
    # This new trait gets assigned to the node i_node_n 
    trait_jumped[ which(names(trait_jumped) == i_node_n)]<- new_trait_s
    ##To change the tip value at the terminal branches i.e. in tips
    if(i_node_n <= ntips)
    {
      i_node_n1 <-tip_label_new[[i_node_n]]
      trait_jumped[ which(names(trait_jumped) == i_node_n1)]<- new_trait_s
    }
    
    else if((i_node_n > ntips) & (events[which(names(events)%in%i_node_n)]=="duplication"))
    {
      next
      #return (trait_jumped)
    }
    
    else if((i_node_n > ntips) & (events[which(names(events)%in%i_node_n)]=="speciation"))
    {
      trait_jumped<-recursive.speciation(tree,i_node_n,tip_label_new,trait_jumped,new_trait_s,sig2,events,ntips)
      #return(trait_jumped)
    }
  }
  return (trait_jumped) 
}


##This function is created to work for progressive jump model where trait jump happens some of the duplication nodes and some of the speciation nodes
##replace i_node_x by i_node for some data
##This function is created to work for progressive jump model
recursive_mixed <- function(tree,new_tree,i_node,tip_label_new,trait_jumped,new_trait_i,sig2,events,ntips,dups_no_jump,spe_no_jump)
{
  #options(expressions=100)
  i_node_des<-as.vector(tree$edge[ which(tree$edge[,1] == names(new_trait_i)) ,2]) 
  for(z in 1:length(i_node_des))
  {
    i_node_n = i_node_des[z]
    
    # taking the branch length of the branch between the nodes and assign to a new single-branch fake tree, for the BM calculation
    new_tree$edge.length<-NULL
    new_tree$edge.length <- tree$edge.length[which(tree$edge[,1] == i_node & tree$edge[,2]==i_node_n)]
    
    # starting trait from the simulated at the beginning, simulating with it's new ancestral trait value 
    new_trait <- trait_jumped[which(names(trait_jumped)==i_node)]
    
    # Estimate the new trait based on a BM process, the branch length of the randomly selected branch, and the jumped ancestral state i_trait
    new_trait_s <- simBM(new_tree, a=new_trait, sig2=sig2, bounds=c(0,1), nsim=1, mu=0, internal=F)
    names(new_trait_s) <- i_node_n
    
    # This new trait gets assigned to the node i_node_n 
    trait_jumped[ which(names(trait_jumped) == i_node_n)]<- new_trait_s
    
    ##To change the tip value at the terminal branches i.e. in tips
    if(i_node_n <= ntips)
    {
      i_node_n1 <-tip_label_new[[i_node_n]]
      trait_jumped[ which(names(trait_jumped) == i_node_n1)]<- new_trait_s
    }
    else if((i_node_n > ntips) & (events[which(names(events)%in%i_node_n)]=="duplication") & (i_node_n%in%dups_no_jump))
    {
      trait_jumped<-recursive_mixed(tree,new_tree,i_node_n,tip_label_new,trait_jumped,new_trait_s,sig2,events,ntips,dups_no_jump,spe_no_jump)
    }
    
    else if((i_node_n > ntips) & (events[which(names(events)%in%i_node_n)]=="duplication"))
    {
      next
      #return (trait_jumped)
    }
    
    else if((i_node_n > ntips) & (events[which(names(events)%in%i_node_n)]=="speciation") & (i_node_n%in%spe_no_jump))
    {
      #return(recursive.speciation(tree,i_node_n,tip_label_new,trait_jumped,new_trait_s,sig2,events,ntips,dups_no_jump,spe_no_jump))
      trait_jumped<-recursive_mixed(tree,new_tree,i_node_n,tip_label_new,trait_jumped,new_trait_s,sig2,events,ntips,dups_no_jump,spe_no_jump)
    }
    else if((i_node_n > ntips) & (events[which(names(events)%in%i_node_n)]=="speciation"))
    {
      next
    }
  }
  return (trait_jumped) 
}


## original paintBranches function of phytools modified by Tina
paintBranches.tina<-function(tree,edge,state,anc.state="1",t){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(is.null(tree$maps)) maps<-lapply(tree$edge.length,function(x) c(setNames(x,anc.state),setNames(x,state)))
  else maps<-tree$maps
  ii<-sapply(edge,function(x,y) which(y==x),y=tree$edge[,2])
  for(i in 1:length(ii)) maps[[ii[i]]]<-c(setNames(t,state),setNames((tree$edge.length[[ii[i]]]-t),anc.state)) 
  ## build mapped.edge matrix
  s<-vector()
  for(i in 1:nrow(tree$edge)) s<-c(s,names(maps[[i]]))
  s<-unique(s)
  mapped.edge<-matrix(0,length(tree$edge.length),length(s),dimnames=list(edge=apply(tree$edge,1,function(x) paste(x,collapse=",")),state=s))
  for(i in 1:length(maps)) for(j in 1:length(maps[[i]])) mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
  
  ## add attributes to the tree
  tree$mapped.edge<-mapped.edge
  tree$maps<-maps
  class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
  tree
}

## Function written to perform recursive painting for a given time period
recursive.painting<-function(ptree,ntip,parentnode,pedgelength,input.time,all.event)
{
  nodes.descendant<-as.vector(ptree$edge[which(ptree$edge[,1]== parentnode),2])
  n_desc_node<-NULL
  edgelength.new<-NULL
  
  for(j in 1:length(nodes.descendant))
  {
    edgelength.new<-pedgelength
    n_desc_node<-nodes.descendant[j]
    desc.edgelength<-NULL
    desc.edgelength<-ptree$edge.length[which(ptree$edge[,1]==parentnode & ptree$edge[,2]==n_desc_node)]
    edgelength.new<-edgelength.new + desc.edgelength
    
    if((edgelength.new <= input.time) & (n_desc_node > ntip))
    {
      ptree<-paintBranches(ptree,n_desc_node,state = "D",anc.state = "S")
    }
    if((edgelength.new > input.time) & (n_desc_node > ntip))
    {
      ## If the combined length is greater than time, need to paint a segment 
      ##that makes up the paren length to time with "D" state and rest of the segment with "S" state
      newtime<-input.time - pedgelength
      ptree<-paintBranches.tina(ptree,n_desc_node,state = "D",anc.state = "S",newtime)
    }
    
    ##Conditions to test
    if(n_desc_node <= ntip) 
    {
      ## If the combined length equals to time, need to paint the whole descendant branch with "D"state
      if(edgelength.new <= input.time)
      {
        ptree<-paintBranches(ptree,n_desc_node,state = "D",anc.state = "S")
      }
      if(edgelength.new > input.time)
      {
        ## If the combined length is greater than time, need to paint a segment 
        ##that makes up the paren length to time with "D" state and rest of the segment with "S" state
        ntime<-input.time - pedgelength
        
        ptree<-paintBranches.tina(ptree,n_desc_node,state = "D",anc.state = "S",ntime)
      }
    }   
    else if((n_desc_node > ntip) & (all.event[which(names(all.event)%in% n_desc_node)]=="duplication"))
    {
      next
    } 
    else if((n_desc_node > ntip) & (all.event[which(names(all.event)%in% n_desc_node)]=="speciation"))
    {
      
      parent.edgelength<-edgelength.new
      new_parent.node<-n_desc_node
      if(parent.edgelength < input.time)
      {
        ptree<-recursive.painting(ptree,ntip,new_parent.node,parent.edgelength,input.time,all.event)
      }     
    }
    else 
    { 
      stop("should never reach")
    }
  }
  return (ptree) 
}  

## Function to compute Wilcoxon one-tailed test of data1 and data2 using one sided test
one_tailed_wilcox<-function (data1, data2)
{
  wilcox_oc_one_tailed <- wilcox.test(data1,data2,alternative="greater")$p.value
  
  #star <-stars.pval(wilcox_oc_one_tailed)
  #star0<- stars.pval(0)
  if (wilcox_oc_one_tailed == 0){label_p = paste0("P < 2.2e-16")}
  if ((wilcox_oc_one_tailed != 0))
  {
    wilcox_oc_one_tailed <- format(wilcox_oc_one_tailed, digits= 3, scientific = TRUE)
    label_p =paste0("P = ",wilcox_oc_one_tailed)
  }
  return(label_p)
} 

## Function to compute Wilcoxon two-tailed test of data1 and data2 using two sided test
two_tailed_wilcox<-function (data1, data2)
{
  wilcox_oc_two_tailed <- wilcox.test(data1,data2,alternative="two.sided")$p.value
  
  #star <-stars.pval(wilcox_oc_two_tailed)
  #star0<- stars.pval(0)
  if (wilcox_oc_two_tailed == 0){label_p = paste0("P < 2.2e-16")}
  if ((wilcox_oc_two_tailed != 0))
  {
    wilcox_oc_two_tailed <- format(wilcox_oc_two_tailed, digits= 3, scientific = TRUE)
    label_p =paste0("P = ",wilcox_oc_two_tailed)
  }
  return(label_p)
} 

## Correlation anaysis function
## Function to compute correlation in a time interval for pairwise data
corr_analysis_pw<-function(test_OC_pairwise)
{
  test_OC_pairwise$Time.MRCA.new<-round(test_OC_pairwise$Time.MRCA,3)
  test_OC_pairwise<-test_OC_pairwise[order(test_OC_pairwise$Event, test_OC_pairwise$Time.MRCA.new),] ##Sorting based on events and on MRCA age
  
  xs <- data.frame(R=NA, Mya=NA, Event=NA)
  xs_all <- data.frame(R=NA, Mya=NA, Event=NA)

  ##dataframe for interval data
  ##Dividing our result as per time frame
  xs_interval<-data.frame(R=NA, Mya=NA, all1.Event=NA,n=NA)

  ##Initialization of variables and data frames for further use
  all_corr<- data.frame(R=NA, P=NA, Mya=NA, Event=NA,n=NA)
  all_dupl<- data.frame(R=NA, P=NA, Mya=NA, Event=NA,n=NA) ##To draw ablines
  all_spec<- data.frame(R=NA, P=NA, Mya=NA, Event=NA,n=NA) ##To draw ablines
  all1<-NULL
  dup_c<-NULL
  spe_c<-NULL
  New_dup<-NULL
  New_spe<-NULL
  
  all<- test_OC_pairwise[order(test_OC_pairwise$Time.MRCA.new),]
  
  ## Making a time interval bin of 0.05 My for analysis
  ##Dividing our result as per time frames
  time_interval<-data.frame(start.range=seq(0,20, by=0.051),
                            # c(0,0.1,1,1.5,2,2.5,3,3.5), 
                            end.range=seq(0.05,20.05,by=0.051))
  #c(0.1,0.21,1.5,2,2.5,3,3.5,4),
  #range=c("0-<0.5","0.5-<1","1-<1.5","1.5-<2","2-<2.5","2.5-<3","3-<3.5","3.5-<4"))
  for(t in 1:nrow(time_interval))
  {
    print(t)
    all1<-NULL
    all1_result<-NULL
    all1_result1<-NULL
    Dup<-NULL
    Spe<-NULL
    dup_c<-NULL
    spe_c<-NULL
    New_dup<-NULL
    New_spe<-NULL
    all1<-subset(all,all$Time.MRCA.new>=time_interval[t,1] & all$Time.MRCA.new<time_interval[t,2])
    
    all1_result <-by(all1, all1$Event,function(X){X<-cor(X['Tau_tip1'],X['Tau_tip2'], method = "pearson")}) ##but we can not get P value against each correlation by this method
    all1_result1 <- data.frame(R=as.vector(all1_result),Mya=median(c(time_interval[t,1],time_interval[t,2])),Event=as.vector(dimnames(all1_result)),n=as.vector(tabulate(all1$Event)))
    xs_interval<-rbind(xs_interval,all1_result1)
    xs_interval<-xs_interval[complete.cases(xs_interval),]
    #xs_interval<-xs_interval[-1,]
    xs_interval <- subset(xs_interval,xs_interval$n>=20) ##Considering atleast 20 points to get the correlation
    
    ##Another way to get P values against each correlation
    Dup<-(all1[all1$Event=="duplication",]) ##Identifying duplication events
    if(nrow(Dup)>=20){
      dup_c <-cor.test(Dup$Tau_tip1,Dup$Tau_tip2)
      New_dup<-data.frame(R=as.numeric(dup_c$estimate),P= as.numeric(dup_c$p.value),Mya=median(c(time_interval[t,1],time_interval[t,2])),Event=as.character(unique(Dup$Event)),n=tabulate(Dup$Event))
      all_corr<-rbind(all_corr,New_dup)
      all_dupl<-rbind(all_dupl,New_dup)}
    
    if(nrow(Dup)<20){
      New_dup<-data.frame(R=NA,P= NA,Mya=median(c(time_interval[t,1],time_interval[t,2])),Event="duplication",n=0)
      all_corr<-rbind(all_corr,New_dup)
      all_dupl<-rbind(all_dupl,New_dup)}
    
    Spe<-(all1[all1$Event=="speciation",]) ##Identifying speciation events
    
    if(nrow(Spe)>=20){
      spe_c <-cor.test(Spe$Tau_tip1,Spe$Tau_tip2)
      New_Spe<-data.frame(R=as.numeric(spe_c$estimate),P= as.numeric(spe_c$p.value),Mya=median(c(time_interval[t,1],time_interval[t,2])),Event=as.character(unique(Spe$Event)),n=tabulate(Spe$Event)[-1])
      all_corr<-rbind(all_corr,New_Spe)
      all_spec<-rbind(all_spec, New_Spe)}
    
    if(nrow(Spe)<20){
      New_Spe<-data.frame(R=NA,P= NA,Mya=median(c(time_interval[t,1],time_interval[t,2])),Event="speciation",n=0)
      all_corr<-rbind(all_corr,New_Spe)
      all_spec<-rbind(all_spec, New_Spe)}
  }
  ## Collecting correlation data for both events 
  all_corr<-all_corr[complete.cases(all_corr),]##New dataframe with R and P value against each time point
  all_corr<-unique(all_corr)
  
  return(all_corr) 
}


## Selecting best fit model for fitting curve in pairwise analysis
model_selection<-function(data)
{
  model1<-lm(R~poly(Mya,1,raw=TRUE), data=data) ##linear model
  model1.summary<-summary(model1)
  model2<-lm(R~poly(Mya,2,raw=TRUE), data=data) ##polynomial model with degree 2 
  model2.summary<-summary(model2)
  model3<-lm(R~poly(Mya,3,raw=TRUE), data=data) ##polynomial model with degree 2 
  model3.summary<-summary(model3)
  
  model.selection<-anova(model1,model2,model3)
  model.summary<-data.frame(model=NA,R.square=NA)
  
  capture.output(cat("\t\n"), file="Genomicus_empirical.txt") ## After first time, make it comment line so that all results get append to it
  capture.output(cat("\n\nModel selection for Empirical data \n\n"), append=TRUE, file="Genomicus_empirical.txt") 
  capture.output(cat(" model1 summary:",model1.summary$adj.r.squared,"\n model2 summary:",model2.summary$adj.r.squared,"\n model3 summary:",model3.summary$adj.r.squared,"\n"),append=TRUE, file="Genomicus_empirical.txt")
  #cat(" model1 adjusted R2:",model1.summary$adj.r.squared,"\n model2 adjusted R2:",model2.summary$adj.r.squared,"\n model3 adjusted R2:",model3.summary$adj.r.squared,"\n")
  model.summary<-rbind(model.summary,data.frame(model="linear",R.square=model1.summary$adj.r.squared))
  model.summary<-rbind(model.summary,data.frame(model="quadratic",R.square=model2.summary$adj.r.squared))
  model.summary<-rbind(model.summary,data.frame(model="cubic",R.square=model3.summary$adj.r.squared))
  model.summary<-model.summary[-1,]
  
  ##Printing best fit model
  if(model.selection$`Pr(>F)`[2] < 0.05 ) 
  {
    if(model.selection$`Pr(>F)`[3] < 0.05)
    {
      ##To avoid over-fitting problem
      if((model1.summary$adj.r.squared < model2.summary$adj.r.squared) & (model2.summary$adj.r.squared < model2.summary$adj.r.squared))
      {
        delta1_diff= model2.summary$adj.r.squared - model1.summary$adj.r.squared
        delta2_diff= model3.summary$adj.r.squared - model2.summary$adj.r.squared
        cat("\ndelta1_diff:",delta1_diff,"\ndelta2_diff:",delta2_diff)
        if(delta2_diff <= 0.1)
        {
          capture.output(cat("\nWe consider model2 i.e. quadratic model as the best fit model although polynomial model with degree 3 is fitting better\n"),append=TRUE, file="Genomicus_empirical.txt")
        }
        if(delta2_diff > 0.1)
        {
          capture.output(cat("\nWe consider model3 i.e. cubic model as the best fit model"),append=TRUE, file="Genomicus_empirical.txt")
        }
        
      }
      if((model1.summary$adj.r.squared < model2.summary$adj.r.squared) & (model2.summary$adj.r.squared >= model2.summary$adj.r.squared))
      {
        capture.output(cat("\nmodel2 i.e. quadratic model is the best fit model\n"),append=TRUE, file="Genomicus_empirical.txt")
      }
    } 
    if(model.selection$`Pr(>F)`[3] >= 0.05)
    {
      capture.output(cat("\nmodel2 i.e. quadratic model is the best fit model\n"),append=TRUE, file="Genomicus_empirical.txt")
    }
    
  }
  if(model.selection$`Pr(>F)`[2] >= 0.05 )
  {
    capture.output(cat("\nmodel1 i.e. linear model is the best fit model\n"),append=TRUE, file="Genomicus_empirical.txt")
  }
  return(model.summary) 
}  

## ANCOVA test for linear model fit
ANCOVA_pairwise <- function (data)
{
  capture.output(cat("\n\n ANCOVA result \n\n"), append=TRUE, file="Genomicus_empirical.txt") 
  
  ## Separating speciation and duplication evenets from the data
  all_spec<-data[data$Event=="speciation",]
  all_dupl<- data[data$Event=="duplication",]
  
  AP1<-NULL
  AP1<- lm(all_dupl$R~all_dupl$Mya)
  AO1<-NULL
  AO1<- lm(all_spec$R~all_spec$Mya)
  
  
  capture.output(cat("ANCOVA result for our Null model \n\n"), append=TRUE, file="Genomicus_empirical.txt") 
  
  ## Linear regression model1 with interaction term
  reg_imp1<-NULL
  reg_imp2<-NULL
  reg_imp1 <-aov(R~Mya+Event+Mya*Event, data = data)
  interaction_effect<-summary(reg_imp1)
  interaction_P<-interaction_effect[[1]][["Pr(>F)"]][3]
  event_effect_P<- interaction_effect[[1]][["Pr(>F)"]][2]
  covariate_effect_P<- interaction_effect[[1]][["Pr(>F)"]][1]
  
  capture.output(cat("ANCOVA result with interaction term \n"), append=TRUE, file="Genomicus_empirical.txt") 
  capture.output(summary(reg_imp1), append=TRUE, file="Genomicus_empirical.txt") 
  
  ## Linear regression model2 without interaction term
  reg_imp2 <-aov(R~Mya+Event, data = data)
  without_interaction_effect<-summary(reg_imp2)
  event_effect_P1<- without_interaction_effect[[1]][["Pr(>F)"]][2]
  covariate_effect_P1<- without_interaction_effect[[1]][["Pr(>F)"]][1]
  
  capture.output(cat("\nANCOVA result without interaction term \n"), append=TRUE, file="Genomicus_empirical.txt") 
  capture.output(summary(reg_imp2), append=TRUE, file="Genomicus_empirical.txt") 
  
  ## testing which model is better between thw two
  test<-anova(reg_imp1,reg_imp2)
  capture.output(cat("\nComparing two ANCOVA models\n"), append=TRUE, file="Genomicus_empirical.txt")
  capture.output(test, append=TRUE, file="Genomicus_empirical.txt") 
  
  pValue_comp <- test[6][2,1]
  capture.output(cat("\n"), append=TRUE, file="Genomicus_empirical.txt") 
  
  
  if(pValue_comp >= 0.05) 
  {
    capture.output(cat("\nHence,the most parsimonious model is model 2\n"), append=TRUE, file="Genomicus_empirical.txt")
    if(event_effect_P1 >= 0.05) 
    {
      star <-stars.pval(event_effect_P1)
      event_effect_P1 <- format(event_effect_P1, digits= 3, scientific = TRUE)
      label_p1 =paste0("ANCOVA, P = ",event_effect_P1, star)
      capture.output(cat("\nThis means there is a single regression line: no significant difference\n\n"), append=TRUE, file="Genomicus_empirical.txt")
    }
    if((event_effect_P1 < 0.05))
    {
      capture.output(cat("\nThis means slope is similar for both events\n"), append=TRUE, file="Genomicus_empirical.txt")
      star <-stars.pval(event_effect_P1)
      star0<- stars.pval(0)
      if (event_effect_P1 == 0){label_p1 = paste0("ANCOVA, P < 2.2e-16",star0)}
      if ((event_effect_P1 != 0))
      {
        event_effect_P1 <- format(event_effect_P1, digits= 3, scientific = TRUE)
        label_p1 =paste0("ANCOVA, P = ",event_effect_P1, star)
      }
      
      ##Now to check whether there is a significant difference in the intercepts
      ##We will fit linear regressions separately for duplication and speciation in this case
      dupli<-NULL
      speci<-NULL
      dupli<-summary(AP1)
      dupli_intercept <-dupli[[4]][1]
      speci<-summary(AO1)
      speci_intercept <-speci[[4]][1]
      if(dupli_intercept > speci_intercept)
      {
        capture.output(cat("\nThis means paralogs have a higher intercept than orthologs\n\n"), append=TRUE, file="Genomicus_empirical.txt")
      }
      if(dupli_intercept < speci_intercept)
      {
        capture.output(cat("\nThis means orthologs have a higher intercept than paralogs\n\n"), append=TRUE, file="Genomicus_empirical.txt")
      }
    }
    
  }
  if(pValue_comp < 0.05) 
  {
    capture.output(cat("\nHence,the most parsimonious model is model 1\n Hence,the interaction term significantly affect the model\n\n"), append=TRUE, file="Genomicus_empirical.txt")
    star <-stars.pval(event_effect_P)
    star0<- stars.pval(0)
    if (event_effect_P == 0){label_p1 = paste0("ANCOVA, P < 2.2e-16",star0)}
    if ((event_effect_P != 0))
    {
      event_effect_P <- format(event_effect_P, digits= 3, scientific = TRUE)
      label_p1 =paste0("ANCOVA, P = ",event_effect_P, star)
    }
    
    dupli<-NULL
    speci<-NULL
    dupli<-summary(AP1)
    dupli_intercept <-dupli[[4]][1]
    speci<-summary(AO1)
    speci_intercept <-speci[[4]][1]
    if(dupli_intercept > speci_intercept)
    {
      capture.output(cat("\nThis means paralogs have a higher intercept than orthologs\n\n"), append=TRUE, file="Genomicus_empirical.txt")
    }
    if(dupli_intercept < speci_intercept)
    {
      capture.output(cat("\nThis means orthologs have a higher intercept than paralogs\n\n"), append=TRUE, file="Genomicus_empirical.txt")
    }
  }
  return(label_p1)
}

##  For nonlinear quadratic fit model
Nonlinear_pairwise <- function (data,file)
{
  capture.output(cat("\n\n Nonlinear result \n\n"), append=TRUE, file=file) 
  
  
  ## Separating speciation and duplication evenets from the data
  all_spec<-data[data$Event=="speciation",]
  all_dupl<- data[data$Event=="duplication",]
  
  ##Model1 with interaction term
  ##Model2 without interaction term
  model_gnm1<-NULL
  model_gnm2<-NULL
  m2_spec<- NULL
  m2_dupl<- NULL
  
  model_gnm1<-gnm(R~ Mya + I(Mya^2) + Event + Mya*Event, data = data) ##not using Mult(Mya,Event)
  interaction_effect<-summary(model_gnm1)
  interaction_P<-interaction_effect$coefficients[5,4]
  event_effect_P<- interaction_effect$coefficients[4,4]
  covariate_effect_P<-interaction_effect$coefficients[3,4]
  
  model_gnm2<-gnm(R~Mya + I(Mya^2) + Event, data = data)
  without_interaction_effect<-summary(model_gnm2)
  event_effect_P1<-without_interaction_effect$coefficients[4,4]
  covariate_effect_P1<-without_interaction_effect$coefficients[3,4]
  
  ## Model comparison
  test<-anova(model_gnm1,model_gnm2, test = "F")
  pValue_comp <- test$`Pr(>F)`[2]
  
  ##Considering nonlinear regression separately for speciation and duplication
  m2_spec<-lm(R~poly(Mya,2,raw=TRUE), data=all_spec)
  m2_dupl<-lm(R~poly(Mya,2,raw=TRUE), data=all_dupl)
  
  capture.output(cat("\nnonlinear_regression result with interaction term \n\n"), append=TRUE, file=file) 
  capture.output(summary(model_gnm1), append=TRUE, file=file)
  capture.output(cat("\nnonlinear_regression result without interaction term \n"), append=TRUE, file=file) 
  capture.output(summary(model_gnm2), append=TRUE, file=file)
  capture.output(cat("\nComparing two nonlinear_regression models\n"), append=TRUE, file=file)
  capture.output(test, append=TRUE, file=file) 
  
  
  if(pValue_comp >= 0.05) 
  {
    capture.output(cat("\nHence,the most parsimonious model is model 2\n"), append=TRUE, file=file)
    if(event_effect_P1 >= 0.05) 
    {
      label_p1 = str_c( "p = ", signif(event_effect_P1,4))
      capture.output(cat("\nThis means there is a single regression line: no significant difference\n\n"), append=TRUE, file=file)
    }
    if((event_effect_P1 < 0.05) & (event_effect_P < 0.05) | ((event_effect_P1 < 0.05) & (event_effect_P > 0.05)))
    {
      capture.output(cat("\nThis means slope is similar for both events\n"), append=TRUE, file=file)
      label_p1 = str_c( "p = ", signif(event_effect_P1,4))
      
      ##Now to check whether there is a significant difference in the intercepts
      ##We will fit linear regressions separately for duplication and speciation in this case
      dupli<-NULL
      speci<-NULL
      dupli<-summary(m2_dupl)
      dupli_intercept <-dupli$coefficients[1]
      speci<-summary(m2_spec)
      speci_intercept <-speci$coefficients[1]
      if(dupli_intercept > speci_intercept)
      {
        capture.output(cat("\nThis means paralogs have a higher intercept than orthologs\n\n"), append=TRUE, file=file)
      }
      if(dupli_intercept < speci_intercept)
      {
        capture.output(cat("\nThis means orthologs have a higher intercept than paralogs\n\n"), append=TRUE, file=file)
      }
    }
    
  }
  if(pValue_comp < 0.05) 
  {
    capture.output(cat("\nHence,the most parsimonious model is model 1\n Hence,the interaction term significantly affect the model\n\n"), append=TRUE, file=file)
    label_p1 = str_c( "p = ", signif(event_effect_P,4))
    dupli<-NULL
    speci<-NULL
    dupli<-summary(m2_dupl)
    dupli_intercept <-dupli$coefficients[1]
    speci<-summary(m2_spec)
    speci_intercept <-speci$coefficients[1]
    if(dupli_intercept > speci_intercept)
    {
      capture.output(cat("\nThis means paralogs have a higher intercept than orthologs\n\n"), append=TRUE, file=file)
    }
    if(dupli_intercept < speci_intercept)
    {
      capture.output(cat("\nThis means orthologs have a higher intercept than paralogs\n\n"), append=TRUE, file=file)
    }
  }
  return(label_p1)
}

##  For nonlinear cubic  model fit 
Nonlinear_pairwise2 <- function (data,file) 
{ 
  capture.output(cat("\n\n Nonlinear result \n\n"), append=TRUE, file=file)  
   
   
  ## Separating speciation and duplication evenets from the data 
  all_spec<-data[data$Event=="speciation",] 
  all_dupl<- data[data$Event=="duplication",] 
   
  ##Model1 with interaction term 
  ##Model2 without interaction term 
  model_gnm1<-NULL 
  model_gnm2<-NULL 
  m2_spec<- NULL 
  m2_dupl<- NULL 
   
  model_gnm1<-gnm(R~ Mya + I(Mya^2) + I(Mya^3) + Event + Mya*Event, data = data) ##not using Mult(Mya,Event) 
  interaction_effect<-summary(model_gnm1) 
  interaction_P<-interaction_effect$coefficients[6,4] 
  event_effect_P<- interaction_effect$coefficients[5,4] 
  covariate_effect_P<-interaction_effect$coefficients[4,4] 
   
  model_gnm2<-gnm(R~Mya + I(Mya^2) + I(Mya^3)+ Event, data = data) 
  without_interaction_effect<-summary(model_gnm2) 
  event_effect_P1<-without_interaction_effect$coefficients[5,4] 
  covariate_effect_P1<-without_interaction_effect$coefficients[4,4] 
   
  ## Model comparison 
  test<-anova(model_gnm1,model_gnm2, test = "F") 
  pValue_comp <- test$`Pr(>F)`[2] 
   
  ##Considering nonlinear regression separately for speciation and duplication 
  m2_spec<-lm(R~poly(Mya,3,raw=TRUE), data=all_spec) 
  m2_dupl<-lm(R~poly(Mya,3,raw=TRUE), data=all_dupl) 
   
  capture.output(cat("\nnonlinear_regression result with interaction term \n\n"), append=TRUE, file=file)  
  capture.output(summary(model_gnm1), append=TRUE, file=file) 
  capture.output(cat("\nnonlinear_regression result without interaction term \n"), append=TRUE, file=file)  
  capture.output(summary(model_gnm2), append=TRUE, file=file) 
  capture.output(cat("\nComparing two nonlinear_regression models\n"), append=TRUE, file=file) 
  capture.output(test, append=TRUE, file=file)  
   
  if(pValue_comp >= 0.05)  
  { 
    capture.output(cat("\nHence,the most parsimonious model is model 2\n"), append=TRUE, file=file) 
    if(event_effect_P1 >= 0.05)  
    { 
      label_p1 = str_c( "p = ", signif(event_effect_P1,4)) 
      capture.output(cat("\nThis means there is a single regression line: no significant difference\n\n"), append=TRUE, file=file) 
    } 
    if((event_effect_P1 < 0.05) & (event_effect_P < 0.05) | ((event_effect_P1 < 0.05) & (event_effect_P > 0.05))) 
    { 
      capture.output(cat("\nThis means slope is similar for both events\n"), append=TRUE, file=file) 
      label_p1 = str_c( "p = ", signif(event_effect_P1,4)) 
       
      ##Now to check whether there is a significant difference in the intercepts 
      ##We will fit linear regressions separately for duplication and speciation in this case 
      dupli<-NULL 
      speci<-NULL 
      dupli<-summary(m2_dupl) 
      dupli_intercept <-dupli$coefficients[1] 
      speci<-summary(m2_spec) 
      speci_intercept <-speci$coefficients[1] 
      if(dupli_intercept > speci_intercept) 
      { 
        capture.output(cat("\nThis means paralogs have a higher intercept than orthologs\n\n"), append=TRUE, file=file) 
      } 
      if(dupli_intercept < speci_intercept) 
      { 
        capture.output(cat("\nThis means orthologs have a higher intercept than paralogs\n\n"), append=TRUE, file=file) 
      } 
    } 
     
  } 
  if(pValue_comp < 0.05)  
  { 
    capture.output(cat("\nHence,the most parsimonious model is model 1\n Hence,the interaction term significantly affect the model\n\n"), append=TRUE, file=file) 
    label_p1 = str_c( "p = ", signif(event_effect_P,4)) 
    dupli<-NULL 
    speci<-NULL 
    dupli<-summary(m2_dupl) 
    dupli_intercept <-dupli$coefficients[1] 
    speci<-summary(m2_spec) 
    speci_intercept <-speci$coefficients[1] 
    if(dupli_intercept > speci_intercept) 
    { 
      capture.output(cat("\nThis means paralogs have a higher intercept than orthologs\n\n"), append=TRUE, file=file) 
    } 
    if(dupli_intercept < speci_intercept) 
    { 
      capture.output(cat("\nThis means orthologs have a higher intercept than paralogs\n\n"), append=TRUE, file=file) 
    } 
  } 
  return(label_p1) 
}


## This function is used to plot phylogenetic tree

plot_phylogeny<-function(tree)
{
  ## Setting colors
  colSimmap <- c("#F8766D", "#00BFC4")
  names(colSimmap) <- c("S", "D")
  p1<-plotSimmap(tree,col=colSimmap,mar=c(5,1,1,2))+
  axisPhylo()
  return(p1)
}
