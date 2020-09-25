
## This script is for pairwise comparisons of traits for speciation and duplication events for 3 types of simulated tree models 
##Loading required libraries
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(dplyr)
library(ggplot2)
library(gtools) ## for stars.pval
library(ggrepel)
library(dplyr)
library(gnm)
library(nlme)


## Loading previously stored data
load("Calibration_tree_painted_0.2_0.02two_ways_new.rda")
source("functions_TMM.R")

## setting random seed number to reproduce
set.seed(1234)


## Setting of parameters 
pdup<-0.2 ### Proportion of duplications
# sig2<-Trait evolutionary rates
#ntips <-Number of tips
#dup_adjust<-Times faster following a duplication event
#Timetree<- Original simulated time tree
#Substree<- Substitution rates tree
#Pseudotree<- Time calibrated substitution rates tree

#Function for asymmetric acceleration of sequence and trait evolutionary rates following duplication
asymmetric_function_pw <- function (sig2,ntips,dup_adjust,Timetree,Substree,Pseudotree) {
  
  sigma_vector_oc <- c(sig2, sig2*dup_adjust) 
  names(sigma_vector_oc) <- c("S", "D")
  
  ## Initialization
  PR_results <- NULL
  tree_calibrated<-0
  results <- array(dim=c(ntips-1, 12, length(tree1)))
  
  ## Running loop for all 10000 simulated trees
  for(i in 1:length(Timetree[[1]]))
  {
    #calibrated_painted <-Psuedotree_0.5_0.02painted_null_new[[1]][[i]] 
    calibrated_painted <-Pseudotree[[2]][[i]]
    
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
      bt<- branching.times(tree_painted)
      events <- ifelse(names(bt)%in%dup_nodes,"duplication", "speciation") ##Identifying the two events
      names(events) <- names(bt)
      
      # Simulating trait with values between 0,1 and rate = sig2 for all over the tree
      trait_oc<-sim.rates_bound(tree_painted, sig2=sigma_vector_oc)
      
      ## Obtaining age of all nodes
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
      
      ##Assigning trait data to the tips of the tree
      names(tree_painted$tip.label) <- tree_painted$tip.label
      tree_painted$tip.label <- round(trait_oc,2)
      
      ##Creating a dataframe of tip data for pairwise comparisons
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
      spe_count=nrow(tip_data)-dup_count
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
      #MRCA.age<-age_all[names(age_all) %in% duplicate_pairwise$mrca]
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
test_OC_pairwise_null_0.2<- asymmetric_function_pw(0.02, 100, 1, Timetree_0.2_painted1, Substree_0.2_0.02painted_null, Pseudotree_0.2_0.02painted_null_new)
test_OC_pairwise_frate_0.2<- asymmetric_function_pw(0.02, 100, 5, Timetree_0.2_painted1, Substree_0.2_0.02painted_frate_OC, Pseudotree_0.2_0.02painted_frate_OC_new)
test_OC_pairwise_ftrait_0.2<- asymmetric_function_pw(0.02, 100, 1, Timetree_0.2_painted1, Substree_0.2_0.02painted_ftrait_OC, Pseudotree_0.2_0.02painted_ftrait_OC_new)
test_OC_pairwise_OC1_0.2<- asymmetric_function_pw(0.02, 100, 8, Timetree_0.2_painted1, Substree_0.2_0.02painted_OC1, Pseudotree_0.2_0.02painted_OC1_new)
test_OC_pairwise_OC2_0.2<- asymmetric_function_pw(0.02, 100, 8, Timetree_0.2_painted1, Substree_0.2_0.02painted_OC2, Pseudotree_0.2_0.02painted_OC2_new)
test_OC_pairwise_OC3_0.2<- asymmetric_function_pw(0.02, 100, 2, Timetree_0.2_painted1, Substree_0.2_0.02painted_OC3, Pseudotree_0.2_0.02painted_OC3_new)

## Saving data for further analyses and plotting
save.image("Asym_PW_pdup0.2_sig0.02_all.rda")

## Cleaning memory
rm(dt)
rm(tree1)
rm(Timetree_painted1)
rm(Timetree_0.2_painted1)
rm(Substree_0.2_0.02painted_null)
rm(Substree_0.2_0.02painted_ftrait_OC)
rm(Substree_0.2_0.02painted_frate_OC)
rm(Substree_0.2_0.02painted_OC1)
rm(Substree_0.2_0.02painted_OC2)
rm(Substree_0.2_0.02painted_OC3)
rm(Pseudotree_0.2_0.02painted_null_new)
rm(Pseudotree_0.2_0.02painted_ftrait_OC_new)
rm(Pseudotree_0.2_0.02painted_frate_OC_new)
rm(Pseudotree_0.2_0.02painted_OC1_new)
rm(Pseudotree_0.2_0.02painted_OC2_new)
rm(Pseudotree_0.2_0.02painted_OC3_new)

## Saving data for further analyses and plotting
save.image("Asym_PW_pdup0.2_sig0.02_plot_data1.rda")

############ Analyses further ###############
name<-c("null","ftrait","frate","OC1","OC2","OC3")
dff<-list(test_OC_pairwise_null_0.2,test_OC_pairwise_ftrait_0.2,test_OC_pairwise_frate_0.2,test_OC_pairwise_OC1_0.2,test_OC_pairwise_OC2_0.2,test_OC_pairwise_OC3_0.2)

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
  filename<-paste("model2_",p_dup,"_pw_10000_asym_",name[x],".txt",sep="")
  Pval<-Nonlinear_pairwise(processed_corr,filename)
  assign(paste("Pval_",p_dup,"_",name[x],sep=""), Pval)

  ## Before plotting
  label_p1<-NULL
  label_p<-NULL
  label_p1<-as.numeric(strsplit(Pval,"p = ")[[1]][2]) ## For pdup<-0.2
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
  processed_corr$event = factor(processed_corr$Event, levels=c( "Speciation", "Duplication" ) )
  
  print("now_plotting")
  plot<-NULL
  plot<- ggplot(processed_corr, aes(x=Mya, y=R, color=event)) + 
  geom_point(alpha = 0.2, shape = 16,size=3) +
  #geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,se=F) + 
  geom_smooth(method = "loess",  size = 1,se=F) + 
  #stat_smooth(se=F) + 
  #scale_color_manual(values=c("#00BFC4","#F8766D")) +
  xlab(expression(bold("Time in million years "))) +
  ylab(expression(bold("Pearson correlation coeffcient, R"))) +
  xlim(0,15) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.title=element_blank(),legend.position=c(0.85,0.85)) +
  theme(axis.text=element_text(size=10,face="bold")) +
  ggtitle(expression(bold("Asymmetric acceleration"))) +
  annotate("text", x = 7, y = 0.6, label=label_p, fontface = 4) 
  assign(paste("PW_plot_",p_dup,"_",name[x],sep=""), plot)

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
rm(test_OC_pairwise_null_0.2)
rm(test_OC_pairwise_frate_0.2)
rm(test_OC_pairwise_ftrait_0.2)
rm(test_OC_pairwise_OC1_0.2)
rm(test_OC_pairwise_OC2_0.2)
rm(test_OC_pairwise_OC3_0.2)

## Saving data for further analyses and plotting
save.image("Asym_PW_pdup0.2_sig0.02_plot_data_latest.rda")

