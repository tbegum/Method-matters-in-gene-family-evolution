
## This program is written for phylogenetic independent contrasts (PICs) analyses to test the ortholog conjecture in 6 hypothetical scenario
##Loading required libraries
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtools) ## for stars.pval


## Loading previously stored data
load("Calibration_tree_painted_0.5_0.02two_ways_new.rda")
source("functions_TMM.R")

## To reproduce result, we set the seed number
set.seed(1234)


##Parameters setting 
pdup<-0.5 ### Proportion of duplications


## Function for asymmetric acceleration of sequence and trait evolutionary rates following duplication
asymmetric_function_PIC <- function (sig2,ntips,dup_adjust,Timetree,Substree,Pseudotree) { 
  
  sigma_vector_oc <- c(sig2, sig2*dup_adjust) 
  names(sigma_vector_oc) <- c("S", "D")
  results <- array(dim=c(ntips-1, 12, length(tree1)))
  
  ## Running loop for all 10000 simulated trees
  for(i in 1:length(Timetree[[1]]))
  {
    calibrated_painted <-Pseudotree[[2]][[i]] #Asymmetric model
    
    ## Excluding trees failed calibration step
    if(!is.null(calibrated_painted)) 
    {
      print (i)
      tree<-tree1[[i]]
      Timetree_asym<- Timetree[[2]][[i]]
      tree_painted <- Timetree_asym
      dup_nodes<- as.numeric(Timetree[[3]][[i]])
      spe_nodes<- as.numeric(Timetree[[4]][[i]])
      dup_edges_asym<- as.numeric(Timetree[[6]][[i]])
      tree_sim_painted <- Substree[[2]][[i]] ## Painted asymmetric substitution tree 
      bt<- branching.times(tree)
      events <- ifelse(names(bt)%in%dup_nodes,"duplication", "speciation") ##Identifying the two events as previously stored
      names(events) <- names(bt)
      
      # Simulating trait with values between 0,1 and rate = sig2 for all over the tree
      trait_oc<-sim.rates_bound(tree_painted, sig2=sigma_vector_oc)
      
      ## Estimating the Phylogenetic Independent Contrasts (pics) for all three types of trees   
      pic_Ttree <- pic( trait_oc, tree_painted, var.contrasts=TRUE )  #For original time tree,
      pic_Stree <- pic( trait_oc, tree_sim_painted, var.contrasts=TRUE ) # For substitution tree
      pic_pseudoTtree <- pic( trait_oc,calibrated_painted, var.contrasts = TRUE) # For pseudo calibrated time tree
      
      ## Logging result
      temp_results <- matrix(ncol=12, nrow=tree_sim_painted$Nnode) ##First 4 columns for Time tree data 
      temp_results[,1] <- c((ntips+1) : (ntips+tree_painted$Nnode))
      temp_results[,2] <- events
      temp_results[,3] <- pic_Ttree[,1]
      temp_results[,4] <- branching.times(tree_painted)
      temp_results[,5] <- c((ntips+1) : (ntips+tree_sim_painted$Nnode)) #Second 4 columns for Substitution tree data
      temp_results[,6] <- events
      temp_results[,7] <- pic_Stree[,1]
      temp_results[,8] <- branching.times(tree_sim_painted)
      temp_results[,9] <- c((ntips+1) : (ntips+calibrated_painted$Nnode)) ## Third 4 columns for Pseudo time tree data
      temp_results[,10] <- events
      temp_results[,11] <- pic_pseudoTtree[,1]
      temp_results[,12] <- branching.times(calibrated_painted)
      
      colnames(temp_results) <- c("node_Ttree", "event_Ttree", "PIC_Ttree", "b_time_Ttree", "node_Stree", "event_Stree", "PIC_Stree", "b_time_Stree", "node_pseudotree", "event_pseudotree", "PIC_pseudotree", "b_time_Psuedotree")
      results[,,i] <- temp_results
    }
  }  
  return(results)
}

print("calling function")

#Calling the function
test_OC_null_0.5 <- asymmetric_function_PIC(0.02, 100, 1, Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new) 
test_OC_ftrait_0.5 <- asymmetric_function_PIC(0.02, 100, 1, Timetree_0.5_painted1, Substree_0.5_0.02painted_ftrait_OC, Pseudotree_0.5_0.02painted_ftrait_OC_new) 
test_OC_frate_0.5 <- asymmetric_function_PIC(0.02, 100, 5, Timetree_0.5_painted1, Substree_0.5_0.02painted_frate_OC, Pseudotree_0.5_0.02painted_frate_OC_new) 
test_OC_OC1_0.5 <- asymmetric_function_PIC(0.02, 100, 8, Timetree_0.5_painted1, Substree_0.5_0.02painted_OC1, Pseudotree_0.5_0.02painted_OC1_new) 
test_OC_OC2_0.5 <- asymmetric_function_PIC(0.02, 100, 8, Timetree_0.5_painted1, Substree_0.5_0.02painted_OC2, Pseudotree_0.5_0.02painted_OC2_new) 
test_OC_OC3_0.5 <- asymmetric_function_PIC(0.02, 100, 2, Timetree_0.5_painted1, Substree_0.5_0.02painted_OC3, Pseudotree_0.5_0.02painted_OC3_new) 

## Saving data for further analyses and plotting
save.image("Asym_PIC_pdup0.5_sig0.02_all.rda")

## Cleaning memory
rm(dt)
rm(tree1)
rm(Timetree_painted1)
rm(Timetree_0.5_painted1)
rm(Substree_0.5_0.02painted_null)
rm(Substree_0.5_0.02painted_ftrait_OC)
rm(Substree_0.5_0.02painted_frate_OC)
rm(Substree_0.5_0.02painted_OC1)
rm(Substree_0.5_0.02painted_OC2)
rm(Substree_0.5_0.02painted_OC3)
rm(Pseudotree_0.5_0.02painted_null_new)
rm(Pseudotree_0.5_0.02painted_ftrait_OC_new)
rm(Pseudotree_0.5_0.02painted_frate_OC_new)
rm(Pseudotree_0.5_0.02painted_OC1_new)
rm(Pseudotree_0.5_0.02painted_OC2_new)
rm(Pseudotree_0.5_0.02painted_OC3_new)

## Saving data for further analyses and plotting
save.image("Asym_PIC_pdup0.5_sig0.02_plot_data1.rda")

############ Analyses further ###############
name<-c("null","ftrait","frate","OC1","OC2","OC3")
dff<-list(test_OC_null_0.5,test_OC_ftrait_0.5,test_OC_frate_0.5,test_OC_OC1_0.5,test_OC_OC2_0.5,test_OC_OC3_0.5)

for(x in 1:length(dff))
{
  test_OC<-NULL
  test_OC<-dff[[x]] ## Every time change the dataframe for analyses until plotting

  ##Initializing vectors 
  pic_speciation_Ttree <- vector()
  pic_duplication_Ttree <- vector()
  pic_speciation_Stree <- vector()
  pic_duplication_Stree <- vector()
  pic_speciation_pseudotree <- vector()
  pic_duplication_pseudotree <- vector()

  ##Collecting the PIC values for speciation and duplication events separately for real time tree, substitution tree and for pseudo time tree
  ##To test the OC, we also are interested in magnitude rather than direction. So, we took the absolute value of the contrast 

  print("timetree")
  ## For original time tree 
  for (a in 1:dim(test_OC)[3]) { 
  temp_Ttree <- abs(as.numeric(test_OC[ which(test_OC[,2,a]=="duplication") , 3 , a]))
  pic_duplication_Ttree <- append (pic_duplication_Ttree, temp_Ttree)
  }

  for (a in 1:dim(test_OC)[3]) { 
  temp_Ttree <- abs(as.numeric(test_OC [ which(test_OC[,2,a]=="speciation") , 3 , a]))
  pic_speciation_Ttree <- append (pic_speciation_Ttree, temp_Ttree)
  }
   
  print("substree")
  ##For substitution tree
  for (a in 1:dim(test_OC)[3]) { 
  temp_Stree <- abs(as.numeric(test_OC [ which(test_OC[,6,a]=="duplication") , 7 , a]))
  pic_duplication_Stree <- append (pic_duplication_Stree, temp_Stree)
  }

  for (a in 1:dim(test_OC)[3]) { 
  temp_Stree <- abs(as.numeric(test_OC [ which(test_OC[,6,a]=="speciation") , 7 , a]))
  pic_speciation_Stree <- append (pic_speciation_Stree, temp_Stree)
  }
 
 
  print("pseutree")
  ##for pseudoTime Tree
  for (a in 1:dim(test_OC)[3]) { 
  temp_PseudoTtree <- abs(as.numeric(test_OC [ which(test_OC[,10,a]=="duplication") , 11 , a]))
  pic_duplication_pseudotree <- append (pic_duplication_pseudotree, temp_PseudoTtree)
  }

  for (a in 1:dim(test_OC)[3]) { 
  temp_PseudoTtree <- abs(as.numeric(test_OC [ which(test_OC[,10,a]=="speciation") , 11 , a]))
  pic_speciation_pseudotree <- append (pic_speciation_pseudotree, temp_PseudoTtree)
  }
  
  print("now_comp_diff")

  label_p_ot<-NULL
  label_p_st<-NULL
  label_p_pseut<-NULL
  cutoff<- as.numeric(0.05/84)
  star<-stars.pval(0.05)
  
  ##Two sided Wilcoxon rank sum tests  
  label_p_ot1<-wilcox.test(pic_duplication_Ttree, pic_speciation_Ttree)$p.value
  if (label_p_ot1 == 0){label_p_ot <- paste0("P < 2.2e-16",star)}
  if ((label_p_ot1 != 0))
  {
    label_p_ot1 <- format(label_p_ot1, digits= 3, scientific = TRUE)
    if(label_p_ot1 < cutoff){label_p_ot<-paste0("P = ",label_p_ot1,star)}
    if(label_p_ot1 >= cutoff){label_p_ot<-paste0("P = ",label_p_ot1)}
  }
  #label_p_ot<-two_tailed_wilcox(pic_duplication_Ttree, pic_speciation_Ttree)
  assign(paste("label_p_ot",p_dup,"_",name[x],sep=""), label_p_ot)
  
  label_p_st1<-wilcox.test(pic_duplication_Stree, pic_speciation_Stree)$p.value
  if (label_p_st1 == 0){label_p_st <- paste0("P < 2.2e-16",star)}
  if ((label_p_st1 != 0))
  {
    label_p_st1 <- format(label_p_st1, digits= 3, scientific = TRUE)
    if(label_p_st1 < cutoff){label_p_st<-paste0("P = ",label_p_st1,star)}
    if(label_p_st1 >= cutoff){label_p_st<-paste0("P = ",label_p_st1)}
  }
  assign(paste("label_p_st",p_dup,"_",name[x],sep=""), label_p_st)

  label_p_pseut1<-wilcox.test(pic_duplication_pseudotree,pic_speciation_pseudotree)$p.value
  if (label_p_pseut1 == 0){label_p_pseut <- paste0("P < 2.2e-16",star)}
  if ((label_p_pseut1 != 0))
  { 
    label_p_pseut1 <- format(label_p_pseut1, digits= 3, scientific = TRUE)
    if(label_p_pseut1 < cutoff){label_p_pseut<-paste0("P = ",label_p_pseut1,star)}
    if(label_p_pseut1 >= cutoff){label_p_pseut<-paste0("P = ",label_p_pseut1)}
  }
  assign(paste("label_p_pseut",p_dup,"_",name[x],sep=""), label_p_pseut)

  ## Generating dataframe for plotting
  simulation.df<-data.frame(pic=NA, Event=NA,label=NA)
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_speciation_Ttree, Event="Speciation", label="Original time trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_duplication_Ttree, Event="Duplication", label="Original time trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_speciation_Stree, Event="Speciation", label="Substitution trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_duplication_Stree, Event="Duplication", label="Substitution trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_speciation_pseudotree, Event="Speciation", label="Calibrated time trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_duplication_pseudotree, Event="Duplication", label="Calibrated time trees"))
  simulation.df<-simulation.df[-1,]
  simulation.df$label <- factor(simulation.df$label, levels=c("Original time trees","Substitution trees","Calibrated time trees"))
  simulation.df$Event <- factor( simulation.df$Event, levels=c( "Speciation", "Duplication" ) )

  ## Plotting figure with the simulation data
  dodge <- position_dodge(width = 0.51)

  print("now_plotting")
  plot<-NULL
  plot<-ggplot(simulation.df,aes( x=label, y=pic, fill=Event)) + 
  guides(colour = guide_legend( override.aes = list( shape = 16 ))) +
  geom_boxplot( width=0.5,outlier.colour=NA, position = dodge, notch = T) +
  xlab( NULL ) +
  ylab(expression(bold(paste("PICs of ", tau)))) +
  ylim(0, 2.2) +
  theme_classic()+
  theme(legend.title=element_blank(),legend.position=c(0.9,0.9)) +
  theme(axis.text=element_text(size=10,face="bold")) +
  theme(legend.text=element_text(size=10,face="bold")) +
  annotate("text", x = 1.0, y = 0.5, label= label_p_ot, fontface = 4,size=4) +
  annotate("text", x = 2.0, y = 2.0, label= label_p_st, fontface = 4,size=4) +
  annotate("text", x = 3.0, y = 0.5, label= label_p_pseut, fontface = 4,size=4)
  assign(paste("PIC_plot_",p_dup,"_",name[x],sep=""), plot)

}

## Saving data for backup
save.image("Asym_PIC_pdup0.5_sig0.02_plot_data2.rda")

##Cleaning memory
rm(test_OC)
rm(rate_vector)
rm(simulation.df)
rm(label_p_ot)
rm(label_p_st)
rm(label_p_pseut)
rm(temp_Ttree)
rm(temp_Stree)
rm(temp_PseudoTtree)
rm(pic_duplication_Ttree)
rm(pic_speciation_Ttree)
rm(pic_duplication_Stree)
rm(pic_speciation_Stree)
rm(pic_duplication_pseudotree)
rm(pic_speciation_pseudotree)
rm(plot)
rm(dff)
rm(name)
rm(test_OC_null_0.5)
rm(test_OC_ftrait_0.5)
rm(test_OC_frate_0.5)
rm(test_OC_OC1_0.5)
rm(test_OC_OC2_0.5)
rm(test_OC_OC3_0.5)

## Saving data for further analyses and plotting
save.image("Asym_PIC_pdup0.5_sig0.02_plot_data_latest.rda")

