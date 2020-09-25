
############## Simulating jump in trait asymetrically by Pairwise comparisons ####################
### This model considers asymmetric trait jumps after 30% of the randomly chosen duplication, and 30% of the speciation nodes under a uniform null scenario
### This model considers asymmetric trait jumps after 50% of the randomly chosen duplication, and 20% of the speciation nodes under an ortholog conjecture scenario
### Hence, we find progressive changes in the trait values of descendant speciation and duplication nodes, where trait jump takes place

##Loading required libraries
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(ggplot2)
library(gnm)
library(gtools)

sig2<-0.02 ##Initial trait evolutionary rate
ntips<-100 ##Number of tips
noise <- 0.00 ## Random noise
pdup<-0.5 ###proportions of duplication


## Loading previously stored data
load("Calibration_tree_painted_pdup0.5_sig0.02_NULL_new.rda")

##functions required
source("functions_TMM.R")

## setting random seed number for reproducibility
set.seed(1234)

##Clearing memory
rm(dt)
rm(Timetree_painted1)

## Parameters setting 
noise <- 0.00 ## Random noise
pdup<-0.5 ### Proportion of duplications
#sig2<-Trait evolutionary rates
#ntips<- Number of tips
#Timetree<-Original simulated time tree
#Substree<-Substitution rates tree
#Pseudotree<- Time calibrated substitution rates tree

## Function for pairwise comparisons using asymmetric trait jump model 
asymmetric_jump_mixed_pw <- function (sig2,ntips,Timetree,Substree,Pseudotree,pdup_jump,pspe_jump) { 
  
  PR_results <- NULL
  tree_calibrated<-0
  
  for(i in 1:length(Timetree[[1]]))
  {
    #calibrated_painted <-Psuedotree_0.5_0.02painted_null_new[[1]][[i]] 
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
      tree_sim_painted <- Substree[[1]][[i]] ## The rates remain same for speciation as well as for duplication
      bt<- branching.times(tree_painted)
      events <- ifelse(names(bt)%in%dup_nodes,"duplication", "speciation") ##Identifying two events
      names(events) <- names(bt)
      
      #Simulating trait values between 0,1 and rate = sig2 across each tree
      trait <- fastBM(tree, sig2=sig2, bounds=c(0,1), internal=T)
      trait_names_tips<- grep("^t.", names(trait), value = T)##Just to get the trait values at the tips 
      tip_names<-trait_names_tips ## Initial tipnames
      trait_tips<-trait[tip_names] ##Trait valus at tips
      
      ##Creating a hash table to match integer tiplabels (5) to alphanumeric tiplabels (like t25)
      tip_label_new <- list()
      for(x in 1:length(tip_names))
      {
        tip_label_new[[x]] <-tip_names[x]
      }
      
      ##Sampling nodes for trait jump
      size_dup<-round(length(dup_nodes)*pdup_jump)
      size_spe<-round(length(spe_nodes)*pspe_jump)
      
      ##Selecting duplication and speciation nodes, where trait jump takes place
      dups_jump<-sample(dup_nodes,size_dup,replace = F)
      spe_jump<-sample(spe_nodes,size_spe,replace = F)
      
      nodes_jump<-NULL
      nodes_jump<-append(dups_jump, spe_jump)
      
      ##Identifying duplication and speciation nodes, where jump do not take place
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
        
        # Selecting the nodes that descend on the i duplication node, and sampling randomly only one of those to be affected by the jump, for asymmetric trait jump
        i_node <-sample(tree$edge[ which(tree$edge[,1] == names(i_bt)) ,2], 1)
        #print(i_node)
        new_tree$edge.length <- tree$edge.length[which(tree$edge[,1] == names(i_bt) & tree$edge[,2]==i_node)]
        
        # Starting from the simulated trait value at the beginning, but assuming the jump model, it changes instantaneously.
        #  so, the trait value is drawn randomly between 0 to 1, following a trait jump. 
        i_trait <- runif(1,min=0,max=1)
        
        # Estimate the new trait based on a Brownian (BM) process, the branch length of the randomly selected branch, and the jumped ancestral state i_trait
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
        
        ##If the descendant node is a duplication node where jump does not take place; then we need to reestimate the trait value as now the ancestral trait value has changed 
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
      
      ##Tip tau values after jump
      trait_oc <-trait_jumped[tip_names]
      
      ##getting age of all nodes
      age_all<-branching.times(tree_painted)
      #plotBranchbyTrait(tree_painted, trait_oc)
      
      ##Identifying tips and edges connecting tips
      tip_edge_nodes<- unique(tree_painted$edge[which(tree_painted$edge[,2] <= ntips), 1])
      tip_edges_all<- unique(tree_painted$edge[which(tree_painted$edge[,2] <= ntips), c(1,2)])
      tip_data <- NULL
      trait.label <- NULL
      tip_edge_event<-NULL
      tip_node_age<-NULL
      trait_value <-NULL
      
      ##Assigning trait data at the tips of the tree
      names(tree_painted$tip.label) <- tree_painted$tip.label
      tree_painted$tip.label <- round(trait_oc,2)
      
      ##Creating a dataframe of tip data for pairwise comparison
      tip_data <-data.frame(tip_edges_all)
      colnames(tip_data)<-c("Node","Tip")
      dup_count=0
      spe_count=0
      
      for(x in 1:nrow(tip_edges_all))
      {
        tip_edge_event <-append(tip_edge_event,events[labels(events)%in%tip_edges_all[x,1]])
        if(tip_edge_event[x]=="duplication") dup_count=dup_count+1
        tip_node_age <-append(tip_node_age,age_all[labels(age_all)%in%tip_edges_all[x,1]])
        trait_value <- append(trait_value, tips(tree_painted,tip_edges_all[x,2]))
      }
      spe_count<-nrow(tip_data)-dup_count
      tip_data$Event <- tip_edge_event
      tip_data$Node.Age <-tip_node_age
      tip_data$Tau <- trait_value 
      tip_data$tip.label <- names(trait_value)
      tip_data <- tip_data[c(1:4,6,5)]
      
      #' Summarize pairwise comparisons between tips of a tree
      ## Creating a data frame with one row for each pairwise combination of tips
      ##making pairwise combination of  tip events
      all_pairwise<-data.frame(t(combn(tip_data$Tip,2)))
      colnames(all_pairwise) <-c("Tip","tip_2")
      all_pairwise<-merge(all_pairwise,tip_data,by="Tip")
      colnames(all_pairwise)[1]<-"tip_1"
      colnames(all_pairwise)[7]<-"Tau_1"
      colnames(all_pairwise)[2]<-"Tip"
      all_pairwise<-merge(all_pairwise,tip_data,by="Tip")
      all_pairwise<-all_pairwise[c(2,1,3,8,4,5,10,6,11,7,12)]
      colnames(all_pairwise)<-c("tip1","tip2","Node1","Node2","Event","Node.Age1","Node.Age2","tip.label1","tip.label2","Tau1","Tau2")
      all_pairwise$mrca<-apply(all_pairwise, 1, function(y) getMRCA(tree_painted, c(as.integer(y['tip1']), as.integer(y['tip2']))))
      MRCA.age<-data.frame(age_all)
      mrcaagedf<-data.frame('mrca.age'=  age_all)
      mrcaagedf$mrca=as.numeric(rownames(mrcaagedf))
      all_pairwise <-merge(all_pairwise,mrcaagedf,by='mrca')
      all_pairwise <-all_pairwise[c(2:12,1,13)]
      MRCA.event <- data.frame(events)
      mrcaeventdf <- data.frame('mrca.event'=  events)
      mrcaeventdf$mrca<-as.numeric(rownames(mrcaeventdf))
      all_pairwise<- merge(all_pairwise,mrcaeventdf,by='mrca')
      pairwise_all <- all_pairwise[c(2:12,1,13:14)]
      
      #Logging results in a data frame for pairwise data
      temp_results<-NULL
      temp_results<-data.frame(Event=pairwise_all$mrca.event,Tau_tip1=pairwise_all$Tau1,Tau_tip2=pairwise_all$Tau2,Time.MRCA=pairwise_all$mrca.age)
      PR_results <- rbind(PR_results,temp_results)
      
    } 
  }
  return (PR_results)
  
} 
print("calling function")

#Calling the function
test_asym_jump_pw_mixed_null_0.5 <- asymmetric_jump_mixed_pw(0.02, 100, Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new, 0.3, 0.3) 
test_asym_jump_pw_mixed_oc_0.5 <- asymmetric_jump_mixed_pw(0.02, 100, Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new, 0.5, 0.5) 

## Saving data for further analyses and plotting
save.image("Asym_jump_10000_0.5dup_0.02sig2_pw_all.rda")

############ Further analyses & plotting ###############
name<-c("null","OC")
dff<-list(test_asym_jump_pw_mixed_null_0.5,test_asym_jump_pw_mixed_oc_0.5)

for(x in 1:length(dff))
{
  test_OC_pairwise<-NULL
  test_OC_pairwise<-dff[[x]] ## Every time change the dataframe for analyses until plotting

  ## Correlation anaysis function
  ## Function to compute correlation in a time interval of 0.5 My for pairwise data
  processed_corr<-NULL
  processed_corr<-corr_analysis_pw(test_OC_pairwise)
  assign(paste("processed_corr_",p_dup,"_",name[x],sep=""), processed_corr)

  ### Best fit model selection by curve fitting for pairwise analyses
  best_fit_model<-NULL
  best_fit_model<-model_selection(processed_corr)
  assign(paste("best_fit_model_",p_dup,"_",name[x],sep=""), best_fit_model)

  print("now_comp_diff")
  ### Writing the output of generalized nonlinear model output  
  Pval<-NULL
  filename<-NULL
  filename<-paste("BM_",p_dup,"_pw_10000_asym_jump_",name[x],".txt",sep="")
  Pval<-Nonlinear_pairwise(processed_corr,filename)
  assign(paste("Pval_",p_dup,"_",name[x],sep=""), Pval)

  ## Before plotting
  label_p1<-NULL
  label_p<-NULL
  
  label_p1<-as.numeric(strsplit(Pval,"p = ")[[1]][2]) ## For pdup<-0.5
  cutoff<-as.numeric(0.05/84)
  star<-stars.pval(0.05)
  
  lp<-NULL

  if((label_p1 != 0) & (label_p1 > 2.2e-16))
  {
    lp <- format(label_p1, digits= 3, scientific = TRUE)
    if(label_p1 < cutoff) {label_p <-paste0("P = ",lp,star)}
    if(label_p1 >= cutoff) {label_p <-paste0("P = ",lp)}
  }
  else{label_p <- paste0("P < 2.2e-16", star)}

  processed_corr$Event[which(processed_corr$Event=="speciation")]<-"Speciation"
  processed_corr$Event[which(processed_corr$Event=="duplication")]<-"Duplication"
  processed_corr$event = factor(processed_corr$Event, levels=c("Speciation", "Duplication" ))
  
  print("now_plotting")
  plot<-NULL
  plot<- ggplot(processed_corr, aes(x=Mya, y=R, color=event)) + 
  geom_point(alpha = 0.2, shape = 16,size=3) +
  geom_smooth(method = "loess",  size = 1,se=F) + 
  xlab(expression(bold("Time in million years "))) +
  ylab(expression(bold("Pearson correlation coeffcient, R"))) +
  xlim(0,15) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.title=element_blank(),legend.position=c(0.9,0.9)) +
  theme(axis.text=element_text(size=10,face="bold")) +
  ggtitle(expression(bold("Asymmetric jump"))) +
  annotate("text", x = 7, y = 0.6, label=label_p, fontface = 4) 
  assign(paste("PW_plot_jump_model",p_dup,"_",name[x],sep=""), plot)

}

##Cleaning memory
rm(test_OC_pairwise)
rm(label_p)
rm(label_p1)
rm(dff)
rm(name)
rm(processed_corr)
rm(plot)
rm(Pval)
rm(filename)
rm(best_fit_model)
rm(test_asym_jump_pw_mixed_null_0.5)
rm(test_asym_jump_pw_mixed_oc_0.5)

## Saving data for further analyses and plotting
save.image("Asym_jump_10000_0.5dup_0.02sig2_pw_plot_data_latest.rda")
