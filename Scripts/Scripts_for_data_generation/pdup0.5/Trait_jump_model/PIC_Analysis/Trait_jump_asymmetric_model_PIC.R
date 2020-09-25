

############## Simulating jump in trait asymetrically by PIC method ####################
### This model considers asymmetric trait jumps after 30% of the randomly chosen duplication, and 30% of the speciation nodes under a uniform null scenario
### This model considers asymmetric trait jumps after 50% of the randomly chosen duplication, and 20% of the speciation nodes under an ortholog conjecture scenario
### Hence, we find progressive changes in the trait values of descendant speciation and duplication nodes, where trait jump takes place 

##Loading required libraries
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(ggplot2)
library(gtools)

## Loading previously stored data
load("Calibration_tree_painted_pdup0.5_sig0.02_NULL_new.rda")

## Functions to use
source("functions_TMM.R")
set.seed(1234)## For reproducibility


## Parameters setting 
noise <- 0.00 ## Random noise
pdup<-0.5 ###Proportion of duplication events
pdup_jump<-NULL ## Proportion of jumps for duplication nodes
pspe_jump<-NULL ## Proportion of jumps for speciation nodes
#sig2<-Trait evolutionary rates
#ntips<- Number of tips
#Timetree<-Original simulated time tree
#Substree<-Substitution rates tree
#Pseudotree<- Time calibrated substitution rates tree

## Function for asymmetric trait jump model
asymmetric_jump_mixed <- function (sig2,ntips,Timetree,Substree,Pseudotree,pdup_jump,pspe_jump) { 
  
  results <- array(dim=c(ntips-1, 12, length(tree1)))
  
  for(i in 1:length(Timetree[[1]]))
  {
    calibrated_painted <-Pseudotree[[1]][[i]]
    
    if(!is.null(calibrated_painted)) 
    {
      print (i)
      tree<-tree1[[i]]
      Timetree_sym<- Timetree[[1]][[i]]
      tree_painted <- Timetree_sym
      dup_nodes<- as.numeric(Timetree[[3]][[i]])
      spe_nodes<- as.numeric(Timetree[[4]][[i]])
      dup_edges_sym<- as.numeric(Timetree[[5]][[i]])
      tree_sim_painted <- Substree[[1]][[i]] ## as the rates remain same for speciation as well as for duplication
      bt<- branching.times(tree)
      events <- ifelse(names(bt)%in%dup_nodes,"duplication", "speciation") ##Identifying the two events
      names(events) <- names(bt)
      
      # Simulating trait with bounds between 0,1 and rate = sig2 foreach tree
      trait <- fastBM(tree, sig2=sig2, bounds=c(0,1), internal=T)
      trait_names_tips<- grep("^t.", names(trait), value = T)##Just to get the trait values at the tips 
      tip_names<-trait_names_tips ## Initial tipnames
      trait_tips<-trait[tip_names] ##Trait values at tips
      
      ##Creating a hash table to match integer tiplabels (5) to alphanumeric tiplabels (like t25)
      tip_label_new <- list()
      for(x in 1:length(tip_names))
      {
        tip_label_new[[x]] <-tip_names[x]
      }
      
      ##Sampling nodes where trait jump takes place
      size_dup<-round(length(dup_nodes)*pdup_jump)
      size_spe<-round(length(spe_nodes)*pspe_jump)
      
      ##Randomly selecting duplication and speciation nodes for trait jump
      dups_jump<-sample(dup_nodes,size_dup,replace = F)
      spe_jump<-sample(spe_nodes,size_spe,replace = F)
      
      nodes_jump<-NULL
      nodes_jump<-append(dups_jump, spe_jump)## Nodes where jump in trait takes place
      
      ##Identifying duplication and speciation nodes where jump in trait do not take place
      dups_no_jump<-dup_nodes[!(dup_nodes%in%dups_jump)]
      spe_no_jump <-spe_nodes[!(spe_nodes%in%spe_jump)]
      
      ##To arrange the branching times of the nodes in descending order to make the jump model work
      bt<- branching.times(tree)
      bt_jumps <- sort(bt [ which(names(bt) %in% nodes_jump) ],decreasing=T)
      node_paint<-as.integer(names(bt_jumps))
      
      ## Loop to re-estimate each state after the jump and loop into the new trait
      trait_jumped <- trait 
      new_tree <- drop.tip(tree, tip=tree$tip.label[-1])
      new_tree$tip.label <- c("test")
      
      ##Collecting new trait values after jump
      for (t in 1:length(bt_jumps))
      { 
        #print(t)
        # go over all nodes by time (oldest-> newest)
        i_bt <- bt_jumps[t] 
        
        # Selecting the nodes that descend on the i duplication node, and sampling randomly only one of those to be affected by the jump, for asymmetric jump model
        i_node <-sample(tree$edge[ which(tree$edge[,1] == names(i_bt)) ,2], 1)
        #print(i_node)
        new_tree$edge.length <- tree$edge.length[which(tree$edge[,1] == names(i_bt) & tree$edge[,2]==i_node)]
        
        #Starting from the simulated trait at the beginning, but assuming the jump model, it changes instantaneously 
        #so, the trait value between 0 to 1 after jump is drawn at random. 
        i_trait <- runif(1,min=0,max=1)
        
        # Estimating the new trait based on a Brownian process, the branch length of the randomly selected branch, and the jumped ancestral state i_trait
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
        
        ##If the descendant node is a speciation node, where trait jump takes place; we go to next loop
        else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="speciation") & (i_node%in%spe_jump))
        {
          trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i###check
        }
        
        ##If the descendant node is a duplication node, where jump does not take place; then we need to reestimate the trait value as now the ancestral trait value has changed 
        else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="duplication") & (i_node%in%dups_no_jump))
        {
          trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i
          trait_jumped<-recursive_mixed(tree,new_tree,i_node,tip_label_new,trait_jumped,new_trait_i,sig2,events,ntips,dups_no_jump,spe_no_jump)
          
        }
        
        ##If the descendant node is a duplication node, where trait jump takes place; we go to next loop
        else if((i_node > ntips) & (events[which(names(events)%in%i_node)]=="duplication") & (i_node%in%dups_jump))
        {
          trait_jumped[ which(names(trait_jumped) == i_node)]<- new_trait_i###check
        }
        
      }
      
      #changed_values<-which(trait_jumped != trait)	##To identify the nodes where trait changes happened
      
      ##Tip trait values after jump
      trait_oc <-trait_jumped[tip_names]
      
      ##Phenogram to compare
      #fancyTree(tree, type="phenogram95",x=trait_oc)
      
      ##Estimating the Phylogenetic Independent Contrasts (PICs) for time tree, substitution tree and pseudo time tree respectively
      pic_Ttree = pic( trait_oc, tree_painted, var.contrasts=TRUE ) 
      pic_Stree = pic( trait_oc, tree_sim_painted, var.contrasts=TRUE )
      pic_pseudoTtree = pic( trait_oc,calibrated_painted, var.contrasts = TRUE)
      
      ##Logging result
      temp_results <- matrix(ncol=12, nrow=tree_sim_painted$Nnode) ##First 4 columns for Time tree data 
      temp_results[,1] <- c((ntips+1) : (ntips+tree_painted$Nnode))
      temp_results[,2] <- events
      temp_results[,3] <- pic_Ttree[,1]
      temp_results[,4] <- pic_Ttree[,2]
      temp_results[,5] <- c((ntips+1) : (ntips+tree_sim_painted$Nnode)) #Second 4 columns for Substitution tree data
      temp_results[,6] <- events
      temp_results[,7] <- pic_Stree[,1]
      temp_results[,8] <- branching.times(tree_sim_painted)
      temp_results[,9] <- c((ntips+1) : (ntips+calibrated_painted$Nnode)) ## Third 4 columns for Pseudo time tree data
      temp_results[,10] <- events
      temp_results[,11] <- pic_pseudoTtree[,1]
      temp_results[,12] <- pic_pseudoTtree[,2]
      
      colnames(temp_results) <- c("node_Ttree", "event_Ttree", "PIC_Ttree", "Var_Ttree", "node_Stree", "event_Stree", "PIC_Stree", "b_time_Stree", "node_pseudotree", "event_pseudotree", "PIC_pseudotree", "var_pseudotree")
      results[,,i] <- temp_results
    }
  }  
  return(results)
}


print("calling function")

#Calling the function
test_asym_jump_PIC_mixed_null_0.5 <- asymmetric_jump_mixed (0.02, 100, Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.3,0.3) 
test_asym_jump_PIC_mixed_oc_0.5 <- asymmetric_jump_mixed (0.02, 100, Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.5,0.2) 

## Saving data for further analyses and plotting
save.image("Asym_jump_10000_0.5dup_0.02sig2_PIC_all.rda")

## Cleaning memory
rm(dt)
rm(tree1)
rm(Timetree_0.5_painted1)
rm(Substree_0.5_0.02painted_null)
rm(Pseudotree_0.5_0.02painted_null_new)

## Saving data for further analyses and plotting
save.image("Asym_jump_10000_0.5dup_0.02sig2_PIC_plot_data1.rda")

############ Analyses further ###############
name<-c("null","oc")
dff<-list(test_asym_jump_PIC_mixed_null_0.5,test_asym_jump_PIC_mixed_oc_0.5)

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

  ##Collecting the PIC values separately for speciation and duplication events for real time tree, substitution tree and for pseudo time tree
  ##To test, we also are interested in magnitude rather than direction. So, we took the absolute value of the contrast 

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
  cutoff<-as.numeric(0.05/84)
  star<-stars.pval(0.05)
  
  ##Two sided Wilcoxon tests  
  label_p_ot1<-wilcox.test(pic_duplication_Ttree, pic_speciation_Ttree)$p.value
  if ((label_p_ot1 != 0) & (label_p_ot1 > 2.2e-16))
  { 
    label_p_ot1 <- format(label_p_ot1, digits= 3, scientific = TRUE)
    if(label_p_ot1 <= cutoff){label_p_ot<-paste0("P = ",label_p_ot1,star)}
    if(label_p_ot1 > cutoff){label_p_ot<-paste0("P = ",label_p_ot1)}
  }
  else{label_p_ot <- paste0("P < 2.2e-16",star)}
  assign(paste("label_p_ot",p_dup,"_",name[x],sep=""), label_p_ot)

  label_p_st1<-wilcox.test(pic_duplication_Stree, pic_speciation_Stree)$p.value
  if ((label_p_st1 != 0) & (label_p_st1 > 2.2e-16))
  { 
    label_p_st1 <- format(label_p_st1, digits= 3, scientific = TRUE)
    if(label_p_st1 < cutoff){label_p_st<-paste0("P = ",label_p_st1,star)}
    if(label_p_st1 >= cutoff){label_p_st<-paste0("P = ",label_p_st1)}
  }
  else{label_p_st <- paste0("P < 2.2e-16",star)}
  assign(paste("label_p_st",p_dup,"_",name[x],sep=""), label_p_st)

  label_p_pseut1<-wilcox.test(pic_duplication_pseudotree,pic_speciation_pseudotree)$p.value
  if ((label_p_pseut1 != 0) & (label_p_pseut1 > 2.2e-16))
  {
    label_p_pseut1 <- format(label_p_pseut1, digits= 3, scientific = TRUE)
    if(label_p_pseut1 < cutoff){label_p_pseut<-paste0("P = ",label_p_pseut1,star)}
    if(label_p_pseut1 >= cutoff){label_p_pseut<-paste0("P = ",label_p_pseut1)}
  }
  else{label_p_pseut <- paste0("P < 2.2e-16",star)}
  assign(paste("label_p_pseut",p_dup,"_",name[x],sep=""), label_p_pseut)

  ## Generating dataframe for plotting
  simulation.df<-data.frame(pic=NA, Event=NA,label=NA)
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_speciation_Ttree, Event="Speciation", label="Original time trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_duplication_Ttree, Event="Duplication", label="Original time trees"))
  #simulation.df<-rbind(simulation.df,data.frame(pic=pic_speciation_Stree, Event="Speciation", label="Substitution trees"))
  #simulation.df<-rbind(simulation.df,data.frame(pic=pic_duplication_Stree, Event="Duplication", label="Substitution trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_speciation_pseudotree, Event="Speciation", label="Calibrated time trees"))
  simulation.df<-rbind(simulation.df,data.frame(pic=pic_duplication_pseudotree, Event="Duplication", label="Calibrated time trees"))
  simulation.df<-simulation.df[-1,]
  simulation.df$label <- factor(simulation.df$label, levels=c("Original time trees","Substitution trees","Calibrated time trees"))
  simulation.df$Event <- factor( simulation.df$Event, levels=c( "Speciation", "Duplication" ) )

  ## Plotting figure with the simulated data
  dodge <- position_dodge(width = 0.51)

  print("now_plotting")
  plot<-NULL
  plot<-ggplot(simulation.df,aes( x=label, y=pic, fill=Event)) + 
  guides(colour = guide_legend( override.aes = list( shape = 16 ))) +
  geom_boxplot( width=0.5,outlier.colour=NA, position = dodge, notch = T) +
  xlab( NULL ) +
  ylab(expression(bold(paste("PICs of ", tau)))) +
  ylim(0, 1) +
  theme_classic()+
  theme(legend.title=element_blank(),legend.position=c(0.9,0.9)) +
  theme(axis.text=element_text(size=10,face="bold")) +
  theme(legend.text=element_text(size=10,face="bold")) +
  annotate("text", x = 1.0, y = 0.5, label= label_p_ot, fontface = 4,size=4) +
  annotate("text", x = 2.0, y = 0.5, label= label_p_pseut, fontface = 4,size=4)
  assign(paste("jump_model_PIC_plot_",p_dup,"_",name[x],sep=""), plot)

}

## Saving data 
save.image("Asym_jump_10000_0.5dup_0.02sig2_PIC_plot_data2.rda")

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
rm(test_asym_jump_PIC_mixed_null_0.5)
rm(test_asym_jump_PIC_mixed_oc_0.5)

## Saving data for further analyses and plotting
save.image("Asym_jump_10000_0.5dup_0.02sig2_PIC_plot_data_latest.rda")
