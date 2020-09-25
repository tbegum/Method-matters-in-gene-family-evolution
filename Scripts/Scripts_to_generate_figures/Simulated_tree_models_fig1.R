## This codes is written to generate simulated tree models to generate figure 1

## Required libraries
library(NELSI)
library(geiger)
library(TreeSim)
library(ape)
library(phytools)
library(ggplot2)

set.seed(1234) 

##Setting working directory to use functions
setwd("/Users/admin/Desktop/")

## The functions written for simulation study
source("functions_TMM.R")

## Parameters setting
sig2 <- 0.02 ##Trait evolutionary rate
rate<- 0.02 ###Sequence evolutionary rate
fast<- 8 ###Adjustment of seq evol rate acceleration following duplication
dup_adjust<-2### Adjustment of trait evol rate acceleration following duplication
p_dup <- 0.2 ### Proportion of duplications
ntips <- 20 ### Number of tips 
noise <-0 ## noise parameter
sigma_vector_oc <- c(sig2, sig2*dup_adjust) ##Trait evolutionary rate vector of events
names(sigma_vector_oc) <- c("S", "D")
rate_vector <- c(rate, rate*fast) ##Sequence evolutionary rate vector of events
names(rate_vector) <- c("S", "D")

## Simulating a tree with 10 tips
test_tree <- sim.bd.taxa(n=ntips, 1, lambda=0.4, mu=0.1, complete=FALSE)[[1]]
events <- rep("speciation", times=test_tree$Nnode)

## All internal nodes
all_nodes<- c((ntips+1) : (ntips+test_tree$Nnode))
names(events)<-all_nodes

## Fixing by randomly chosen focal speciation nodes & annotating events
fixed_focal_spe_nodes<-sample(x = all_nodes,size = 4, replace = F )
nodes_excluding_focal_nodes<-as.numeric(all_nodes[!(all_nodes%in%fixed_focal_spe_nodes)])
dup_nodes <- sample(x= nodes_excluding_focal_nodes, size=round(test_tree$Nnode*p_dup), replace=F) 
events[which(names(events)%in% dup_nodes)] <- "duplication"

##painting Tree asymmetrically
##Randomly choosing one of the two duplicate branches for painting 
dup_edge_paint<- vector()
for(x in 1:length(dup_nodes))
{
  dup_edge_selected <- sample(test_tree$edge[ which(test_tree$edge[,1] %in% dup_nodes[x]), 2], 1, replace = F)
  dup_edge_paint<-append(dup_edge_paint,dup_edge_selected)
}
painted_tree<-paintBranches(test_tree,dup_edge_paint,state = "D", anc.state = "S")

##Simulating Tau
simulated_Tau<-round(sim.rates_bound(painted_tree, sig2=sigma_vector_oc),1)
painted_tree$tip.label<- simulated_Tau

par(mfrow=c(1,3))

##Plotting time_tree
plot_phylogeny(painted_tree)+
  mtext("Time in My", side=1, line=3,font = 2, adj=0.5, col=NA, cex= 0.8)+
  title(main = "A. Original time trees",cex.main= 1.5)+
  legend(0.1,12,legend=c("Speciation","Duplication"), ncol=1, pch=19, col=c("#F8766D","#00BFC4"), cex=1, text.font=2, bty='n')

## Substitution tree
tree_rate_sim<-simulate.rates_heterogeneous(painted_tree, rate_vector, noise)
tree_sim <- tree_rate_sim[[1]]
tree_sim_painted <- paintBranches(tree_sim ,dup_edge_paint,state = "D", anc.state = "S")

##Plotting substitution rate tree
plot_phylogeny(tree_sim_painted)+
  mtext("Substitution rates", side=1, line=3,font = 2, adj=0.5, col=NA, cex= 0.8)+
  title(main = "B. Substitution trees",cex.main= 1.5)

##Time calibrating substitution rate tree

## Age of all nodes
age_all<-branching.times(painted_tree) 

##Calibrating time of the substitution tree to obtain psuedo-time tree (Using chronos function to calibrate times for focal speciation nodes)

spe_nodes_new<- NULL
spe_nodes_new<-sort(fixed_focal_spe_nodes)
speciation_age <- as.numeric(age_all[(names(age_all) %in% spe_nodes_new)])

## generation of calibration matrix
focal_calibration_times<-data.frame(node=spe_nodes_new,age=speciation_age)
calibration_matrix <- data.frame(node=spe_nodes_new,age.min=speciation_age, age.max=speciation_age, soft.bound=NA)

## Pesudo time calibrated tree
class(tree_sim_painted) <-"phylo"
tree_sim_painted$mapped.edge <- NULL
tree_sim_painted$maps <- NULL
calibrated_timetree_new <- try(ape::chronos(tree_sim_painted,calibration = calibration_matrix,  model="correlated" ))

##Painting tree and plotting it
Calibrated_tree_painted <- paintBranches (calibrated_timetree_new,edge=dup_edge_paint, "D", anc.state="S")
plot_phylogeny(Calibrated_tree_painted)+
  mtext("Time in My", side=1, line=3,font = 2, adj=0.5, col=NA,cex= 0.8)+
  title(main = "C. Calibrated time trees",cex.main= 1.5)

dev.off()
