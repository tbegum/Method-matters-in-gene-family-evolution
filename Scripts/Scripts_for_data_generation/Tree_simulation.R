
############################Calling functions to generate original simulated time trees, substitution rate trees, and pseudo time calibrated trees for this study #############################
## This study uses 3 different proportions of duplication events, i.e., 0.2, 0.5, and 0.8
##The value of the variable p_dup needs to be changed for different proportions of duplication events

## Libraries need to be loaded
library(NELSI)  ##NELSI simulates rates of molecular evolution along phylogenetic tree
library(geiger)
library(TreeSim)
library(ape)
library(phytools)
library(stringr)


##Initializing variables
nsim<-10000 ##Number of simulations
ntips<-100 ##Number of tips
# Proportion of duplications = p_dup
sig2<-0.02 ##Trait evolutionary rate
rate<-0.02 ##Trait evolutionary rate
noise <- 0.00 ## Random noise


##Setting working directory to use functions
setwd("/Users/admin/Desktop/")
source( "functions_TMM.R" )
set.seed(1234)


##Simulated trees, and stored the tree sets for further use
## These simulated trees are original time trees
tree1<-Gene_trees(nsim, ntips)
save.image("Timetree_10000sim_100tip.rda")
#load("Timetree_10000sim_100tip.rda")


################## Painted, and stored original time trees (symmetrically and asymmetrically) ###################

######## Following the description of 6 elements of each set of 10000 trees
## 1st element <- list of tree_painted_symmetric
## 2nd element <- list of tree_painted_asymmetric
## 3rd element <- list of dup_nodes_paint_tree (example: Timetree_0.2_painted1[[3]][[1]] provides the duplication nodes or the 3rd element of the first painted tree)
## 4th element <- list of spe_nodes_paint_tree
## 5th element <- list of dup_edges_symmetric
## 6th element <- list of dup_edges_asymmetric


## Annotations of events in timetree according to proportion of duplication events
Timetree_painted1<-Initial_timetree_pre(tree1,p_dup)
assign(paste("Timetree_",p_dup,"_painted1",sep=""), Timetree_painted1)

#save.image(paste("Timetree_painted_pdup",p_dup,"_sig",sig2,"_two_ways_latest.rda",sep="")) 
#load("Timetree_painted_pdup0.2_sig0.02_two_ways_latest.rda") # for pdup=0.2
#load("Timetree_painted_pdup0.5_sig0.02_two_ways_latest.rda") # for pdup=0.5
#load("Timetree_painted_pdup0.8_sig0.02_two_ways_latest.rda") # for pdup=0.8


#################### Generating substitution rates tree ######################
##Dataframe with different combination of trait and seq evol rate under different conditions
dt<-data.frame(sig2=rep(sig2,6), ##Trait evolutionary rate 
               rate=rep(rate,6), ###Sequence evolutionary rate 
               dup_adjust=c(1,1,5,8,8,2), ### Adjustment of trait evol rate acceleration
               fast=c(1,5,1,2,8,8), ###Adjustment of seq evol rate acceleration
               condition=c("Scene1: Null hypothesis","Scene2: Fixed trait OC(Tau)", "Scene3: OC_Fixed rate (K)", "Scene4: OC (Tau > K)", "Scene5: OC (Tau = K)", "Scene6: OC (Tau < K)"))

##Substitution rate trees based on six different criterion for each proportions of duplication events
for(x in 1:nrow(dt))
{
  condition<-dt$condition[x]
  rate_vector<-c(dt$rate[x],(dt$rate[x]*dt$fast[x]))
  names(rate_vector) <- c("S", "D")
  if(x == 1)
  {
    null<-NULL
    #Substree_0.2_0.02painted_null<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    null<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    assign(paste("Substree_",p_dup,"_",sig2,"painted_null",sep=""), null)
    rm(null)
  }
  if(x == 2)
  {
    #Substree_0.2_0.02painted_ftrait_OC<- Substitution_tree(Timetree_painted1, rate_vector, noise)
    ftrait<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    assign(paste("Substree_",p_dup,"_",sig2,"painted_ftrait_OC",sep=""), ftrait)
    rm(ftrait)
  }
  if(x == 3)
  {
    #Substree_0.2_0.02painted_frate_OC<-Substitution_tree(Timetree_painted1, rate_vector, noise)
    frate<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    assign(paste("Substree_",p_dup,"_",sig2,"painted_frate_OC",sep=""), frate)
    rm(frate)
  }
  if(x == 4)
  {
    #Substree_0.2_0.02painted_OC1<-Substitution_tree(Timetree_painted1, rate_vector, noise)
    OC1<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    assign(paste("Substree_",p_dup,"_",sig2,"painted_OC1",sep=""), OC1)
    rm(OC1)
  }
  if(x == 5)
  {
    #Substree_0.2_0.02painted_OC2<-Substitution_tree(Timetree_painted1, rate_vector, noise)
    OC2<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    assign(paste("Substree_",p_dup,"_",sig2,"painted_OC2",sep=""), OC2)
    rm(OC2)
  }
  if(x == 6)
  {
    #Substree_0.2_0.02painted_OC3<-Substitution_tree(Timetree_painted1, rate_vector, noise)
    OC3<- Substitution_tree(Timetree_painted1, rate_vector, noise) 
    assign(paste("Substree_",p_dup,"_",sig2,"painted_OC3",sep=""), OC3)
    rm(OC3)
  }
}

##Clearing memory
rm(dt)
#save.image(paste("Substitution_tree_painted_pdup",p_dup,"_sig",sig2,"_two_ways_latest.rda",sep=""))

###To check the plots 
par(mfrow=c(1,3))
plotSimmap(Timetree_0.2_painted1[[1]][[100]],mar=c(5,2,2,2))
axisPhylo()
plotSimmap(Substree_0.2_0.02painted_frate_OC[[1]][[100]],mar=c(5,2,2,2))
axisPhylo()
plotSimmap(Pseudotree_0.2_0.02painted_null_new[[1]][[100]],mar=c(5,2,2,2))
axisPhylo()
remove(tree_sim_painted_asym)
remove(tree_sim_painted_sym)

#save.image("Substitution_tree_painted_pdup0.2_sig0.02_two_ways_latest.rda") ##for pdup=0.2
#save.image("Substitution_tree_painted_pdup0.5_sig0.02_two_ways_latest.rda") ##for pdup=0.5
#save.image("Substitution_tree_painted_pdup0.8_sig0.02_two_ways_latest.rda") ##for pdup=0.8

#################### Pseudo time calibrated substitution rates tree ######################
## We need to load the corresponding substitution rates tree set according to required proportions of duplications (pdup)
## Then we need to run the time calibration step separately for each criteria to save, and to further use 
########For pdup=0.2
rm(Timetree_painted1)
Pseudotree_0.2_0.02painted_null_new<-Calibration_tree2(Timetree_0.2_painted1,Substree_0.2_0.02painted_null)
save.image("Calibration_tree_painted_pdup0.2_sig0.02_NULL_new.rda")
Pseudotree_0.2_0.02painted_ftrait_OC_new<-Calibration_tree2(Timetree_0.2_painted1,Substree_0.2_0.02painted_ftrait_OC)
save.image("Calibration_tree_painted_pdup0.2_sig0.02_ftrait_new.rda")
Pseudotree_0.2_0.02painted_frate_OC_new<-Calibration_tree2(Timetree_0.2_painted1,Substree_0.2_0.02painted_frate_OC)
save.image("Calibration_tree_painted_pdup0.2_sig0.02_frate_new.rda")
Pseudotree_0.2_0.02painted_OC1_new<-Calibration_tree2(Timetree_0.2_painted1,Substree_0.2_0.02painted_OC1)
save.image("Calibration_tree_painted_pdup0.2_sig0.02_OC1_new.rda")
Pseudotree_0.2_0.02painted_OC2_new<-Calibration_tree2(Timetree_0.2_painted1,Substree_0.2_0.02painted_OC2)
save.image("Calibration_tree_painted_pdup0.2_sig0.02_OC2_new.rda")
Pseudotree_0.2_0.02painted_OC3_new<-Calibration_tree2(Timetree_0.2_painted1,Substree_0.2_0.02painted_OC3)
save.image("Calibration_tree_painted_pdup0.2_sig0.02_OC3_new.rda")
remove(calibrated_painted_sym)
remove(calibrated_painted_asym)
save.image("Calibration_tree_painted_0.2_0.02two_ways_new.rda") ## Contains pseudo time calibrated trees for 6 different criteria for pdup = 0.2 


########For pdup=0.5
Psuedotree_0.5_0.02painted_null_new<-Calibration_tree2(Timetree_0.5_painted1,Substree_0.5_0.02painted_null)
save.image("Calibration_tree_painted_pdup0.5_sig0.02_NULL_new.rda")
Psuedotree_0.5_0.02painted_ftrait_OC_new<-Calibration_tree2(Timetree_0.5_painted1,Substree_0.5_0.02painted_ftrait_OC)
save.image("Calibration_tree_painted_pdup0.5_sig0.02_ftrait_new.rda")
Psuedotree_0.5_0.02painted_frate_OC_new<-Calibration_tree2(Timetree_0.5_painted1,Substree_0.5_0.02painted_frate_OC)
save.image("Calibration_tree_painted_pdup0.5_sig0.02_frate_new.rda")
Psuedotree_0.5_0.02painted_OC1_new<-Calibration_tree2(Timetree_0.5_painted1,Substree_0.5_0.02painted_OC1)
save.image("Calibration_tree_painted_pdup0.5_sig0.02_OC1_new.rda")
Psuedotree_0.5_0.02painted_OC2_new<-Calibration_tree2(Timetree_0.5_painted1,Substree_0.5_0.02painted_OC2)
save.image("Calibration_tree_painted_pdup0.5_sig0.02_OC2_new.rda")
Psuedotree_0.5_0.02painted_OC3_new<-Calibration_tree2(Timetree_0.5_painted1,Substree_0.5_0.02painted_OC3)
save.image("Calibration_tree_painted_pdup0.5_sig0.02_OC3_new.rda")
remove(calibrated_painted_sym)
remove(calibrated_painted_asym)
save.image("Calibration_tree_painted_0.5_0.02two_ways_new.rda") ## Contains pseudo time calibrated trees for 6 different criteria for pdup =0.5 


########For pdup=0.8
Psuedotree_0.8_0.02painted_null_new<-Calibration_tree2(Timetree_0.8_painted1,Substree_0.8_0.02painted_null)
save.image("Calibration_tree_painted_pdup0.8_sig0.02_NULL_new.rda")
Psuedotree_0.8_0.02painted_ftrait_OC_new<-Calibration_tree2(Timetree_0.8_painted1,Substree_0.8_0.02painted_ftrait_OC)
save.image("Calibration_tree_painted_pdup0.8_sig0.02_ftrait_new.rda")
Psuedotree_0.8_0.02painted_frate_OC_new<-Calibration_tree2(Timetree_0.8_painted1,Substree_0.8_0.02painted_frate_OC)
save.image("Calibration_tree_painted_pdup0.8_sig0.02_frate_new.rda")
Psuedotree_0.8_0.02painted_OC1_new<-Calibration_tree2(Timetree_0.8_painted1,Substree_0.8_0.02painted_OC1)
save.image("Calibration_tree_painted_pdup0.8_sig0.02_OC1_new.rda")
Psuedotree_0.8_0.02painted_OC2_new<-Calibration_tree2(Timetree_0.8_painted1,Substree_0.8_0.02painted_OC2)
save.image("Calibration_tree_painted_pdup0.8_sig0.02_OC2_new.rda")
Psuedotree_0.8_0.02painted_OC3_new<-Calibration_tree2(Timetree_0.8_painted1,Substree_0.8_0.02painted_OC3)
save.image("Calibration_tree_painted_pdup0.8_sig0.02_OC3_new.rda")
remove(calibrated_painted_sym)
remove(calibrated_painted_asym)
save.image("Calibration_tree_painted_0.8_0.02two_ways_new.rda") ## Contains pseudo time calibrated trees for 6 different criteria for pdup =0.8 


