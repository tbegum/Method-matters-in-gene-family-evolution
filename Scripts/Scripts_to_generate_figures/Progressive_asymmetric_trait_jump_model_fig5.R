
############## Simulating jump in trait asymmetrically ####################

## Due to the instant trait jump, the values at tips also changes

##Loading required libraries
library(geiger)
library(TreeSim)
library(phytools)

ntips<-10 ### Number of tips 
nsim<-1### Number of simulations
p_dup<-0.2   ### Proportion of duplications
sig2<-0.02 ##Trait evolutionary rate


##Setting working directory to use functions
setwd("/Users/admin/Desktop/")
source( "functions_TMM.R" )
set.seed(1234)

# function simualtes trees under the constant rate birth-death process
tree<-sim.bd.taxa(n=ntips, nsim, lambda=0.4, mu=0.1, complete=FALSE)[[1]]

##Using a single tree of our simulated dataset
#tree<-tree1[[1]]

# simulating Tau with bounds between 0,1 and rate = sig2
trait <- fastBM(tree, sig2=sig2, bounds=c(0,1), internal=T)
trait_names_tips<- grep("^t.", names(trait), value = T)##Just to get the trait values at the tips 
tip_names<-trait_names_tips ## Initial tipnames
trait_tips<-trait[tip_names] ##Tau valus at tips

##Creating a hash table to get the tiplabels 
tip_label_new <- list()
for(x in 1:length(tip_names))
{
  tip_label_new[[x]] <-tip_names[x]
}  
#tree$tip.label <- round(trait_tips,2)

# randomizing types of events in nodes
events <- rep("speciation", times=tree$Nnode)
dup_nodes <- sample(x= c(1:length(events)), size=round(tree$Nnode*p_dup), replace=F) 

## duplication can not be assigned to root, otherwise in OC the mapping doesn't work 
#if (length( which(dup_nodes == 1)) == 1) {
# dup_nodes <- dup_nodes[-which(dup_nodes == 1)]
#}

events [dup_nodes ] <- "duplication"
names(events) <- c((ntips+1) : (ntips+tree$Nnode))

dups<- as.numeric(names(events)[ which(events == "duplication")])
spe <-as.numeric(names(events)[ which(events == "speciation")])
size_dup<-round(length(dups)*0.5)
size_spe<-round(length(spe)*0.2)
dup_edges <- unique(tree$edge[ which(tree$edge[,1] %in% dups), 2])
tree_painted <- paintBranches (tree, edge=dup_edges, "D", anc.state="S")

##Selecting duplication and speciation nodes for jump
dups_jump<-sample(dups,size_dup,replace = F)
spe_jump<-sample(spe,size_spe,replace = F)
nodes_jump<-NULL
nodes_jump<-append(dups_jump, spe_jump)

##Identifying duplication nodes where jump does not take place
##For these nodes trait simulation (with new ancestral trait value) will be followed as per speciation nodes after jump 
dups_no_jump<-dups[!(dups%in%dups_jump)]
spe_no_jump <-spe[!(spe%in%spe_jump)]

#dup_edges <- unique(tree$edge[ which(tree$edge[,1] %in% dups), 2])
#tree_painted <- paintBranches (tree, edge=dup_edges, "D", anc.state="S")
bt<- branching.times(tree)
bt_all<-sort(bt, decreasing = T)
bt_jumps <- sort(bt [ which(names(bt) %in% nodes_jump) ],decreasing=T)
node_paint<-as.integer(names(bt_jumps))

## loop to re-estimate each state after the jump and loop into the new trait
trait_jumped <- trait 
#trait_jumped_tips<-trait_tips

tip_edge_nodes<- unique(tree$edge[which(tree$edge[,2] <= ntips), 1])
tip_edges_all<- unique(tree$edge[which(tree$edge[,2] <= ntips), c(1,2)])
#tip_trait_value <- tips(tree,tip_edges_all[x,2]))

new_tree <- drop.tip(tree, tip=tree$tip.label[-1])
new_tree$tip.label <- c("test")


##Starting from the oldest duplication/speciation node for trait jump
for (t in 1:length(bt_jumps))
{ 
  #print(t)
  # go over all nodes by time (oldest-> newest)
  i_bt <- bt_jumps[t] 
  
  # selecting the nodes that descend on the i duplication node, and sampling randomly only one of those to be affected by the jump, as assymetrical
  i_node <-sample(tree$edge[ which(tree$edge[,1] == names(i_bt)) ,2], 1)
  #print(i_node)
  
  # taking the branch length of the branch between the nodes and assign to a new single-branch fake tree, for the BM calculation
  new_tree$edge.length <- tree$edge.length[which(tree$edge[,1] == names(i_bt) & tree$edge[,2]==i_node)]
  
  # starting trait from the simulated at the beginning, but assuming the jump model, it changes instantaneously so, is the starting trait is drawn from random new value. 
  i_trait <- runif(1,min=0,max=1)
  
  # Estimate the new trait based on a BM process, the branch length of the randomly selected branch, and the jumped ancestral state i_trait
  new_trait_i <- simBM(new_tree, a=i_trait, sig2=sig2, bounds=c(0,1), nsim=1, mu=0, internal=F)
  names(new_trait_i) <- i_node
  
  
  ##To check whether this duplication has happend at the terminal node (near tips) and to change the tip value at the terminal branches i.e. in tips
  if(i_node <= ntips)
  {
    i_node1 <-tip_label_new[[i_node]]
    # This new trait value gets assigned to the node i_node 
    trait_jumped[ which(names(trait_jumped) == i_node1)]<- new_trait_i
  }
  
  ##If the descendant node is a speciation node; then we need to reestimate the trait value as now the ancestral trait value has changed 
  else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="speciation") & (i_node%in%spe_no_jump))
  {
    trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i
    trait_jumped<-recursive_mixed(tree,new_tree,i_node,tip_label_new,trait_jumped,new_trait_i,sig2,events,ntips,dups_no_jump,spe_no_jump)
   
  }
  
  ##If the descendant node is a speciation node for jump; we go to next loop
  else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="speciation") & (i_node%in%spe_jump))
  {
    trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i###check
  }
  
  ##If the descendant node is a duplication node where jump does not take place; then we need to reestimate the trait value as now the ancestral trait value has changed 
  else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="duplication") & (i_node%in%dups_no_jump))
  {
    trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i
    trait_jumped<-recursive_mixed(tree,new_tree,i_node,tip_label_new,trait_jumped,new_trait_i,sig2,events,ntips,dups_no_jump,spe_no_jump)
  }
  
  ##If the descendant node is a duplication node for jump; we go to next loop
  else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="duplication") & (i_node%in%dups_jump))
  {
    trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i###check
  }
  
}
trait_jumped_tips<- grep("^t.", names(trait_jumped), value = T)
jumped_tips_trait<-trait_jumped[trait_jumped_tips] ##Tau valus at tips

# To compare the results

#png("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Plots/jump_model_example.png",width=500,height=300)
#pdf("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Plots/jump_assym_dup_example.pdf", width=4, height=6)

par(mfrow=c(1,3),mar=c(5,4,2,3))
colSimmap<-c("#F8766D","#00BFC4")
names(colSimmap)<-c("S","D")

##Plotting Time time and substitution tree
#plotTree.wBars(tree_painted,trait_oc,method="plotSimmap",col= "cornflowerblue",mar=c(5,1,5,2))
# legend(0,28,legend=c("Trait evolution rate",str_c("Duplication= ",sig2*dup_adjust),str_c("Speciation= ",sig2)), pch=19, col=c("white","black","red"),bty="n")
plotSimmap(tree_painted,col=colSimmap,ftype="reg",fsize=0.9,mar = c(5,4,2,3)) 
axisPhylo(font=2, cex=1) 
title(main = "A. Example of time tree", cex.main= 1.2)
mtext("Time in My", font = 2, side = 1, line = 3,cex = 0.7)
#legend(2,21,legend=c("Trait (Tau)"), text.font=2, ncol=1, pch=19, col=c("white"), cex=1, bty='n' )
#legend(0.2,5,legend=c("Speciation","Duplication"), ncol=1, pch=19, col=c("#F8766D","#00BFC4"), cex=0.9,text.font = 2, bty='n')
legend(0.4,10,legend=c("Speciation","Duplication"), ncol=1, pch=16, col=c("#F8766D","#00BFC4"), cex=0.8,text.font = 1, bty='n')
nodelabels(node=nodes_jump,frame = "none",col = "purple", adj=c(-0.05,1.2),font=2, cex=1.2 )
nodelabels(node=nodes_jump,frame="none", pch=16, cex= 1.3, col="blue")

phenogram(tree_painted, trait, main="B. Phenogram without jump", mar = c(5,4,2,3), cex.main= 1.3,colors = colSimmap,spread.labels=T,spread.cost=c(0,0),ylim=c(0.0,1),ftype ="b",font.lab=2, ylab=expression(bold(paste("Phenotypic trait, ", tau))),xlab=expression(bold("Time in My")),cex.lab=2)
#fancyTree(tree=tree,type="phenogram95",x=trait_tips,main="Simulation without trait jump",spread.labels=TRUE, spread.cost=c(1,0),ylim=c(0.0,1),ylab=expression(bold("Phenotype")),xlab=expression(bold("My")))
nodelabels(node=nodes_jump,frame = "none",col = "purple", adj=c(-0.05,1.2),font=2, cex=1.2 )
nodelabels(node=nodes_jump,frame="none", pch=16, cex= 1.3, col="blue")

phenogram(tree_painted, trait_jumped,main="C. Phenogram after jump", mar = c(5,4,2,3), cex.main= 1.3,colors = colSimmap,spread.labels=T,spread.cost=c(0,0),ylim=c(0.0,1),ftype ="b",font.lab=2,cex.axis=1.2, ylab=expression(bold(paste("Phenotypic trait, ", tau))),xlab=expression(bold("Time in My")))
#fancyTree(tree=tree,type="phenogram95",x=jumped_tips_trait,main="Simulation after trait jump",spread.labels=TRUE, spread.cost=c(1,0),ylim=c(0.0,1),ylab=expression(bold("Phenotype")),xlab=expression(bold("My")))
nodelabels(node=nodes_jump,frame = "none",col = "purple", adj=c(-0.05,1.2),font=2, cex=1.2 )
nodelabels(node=nodes_jump,frame="none", pch=16, cex= 1.3, col="blue")
#dev.off()


