
## This script takes into account trait evolutionary rate heterogeneity after oldest duplication for a specific period for pairwise comparisons
###  This model assumes trait evolutionary rates are same following a speciation ('S'), and a duplication ('D') for a uniform null model
###  After duplication, trait evoluionary rates accelerates for certain period of time like 50% or 60% or 70% or 80% or 90%
###  After that period, trait evolves with a rate of preduplication, i.e., speciation state


##Loading required libraries
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(ggplot2)
library(dplyr)
library(gnm)
library(gtools)

## Loading previously stored data after running "Simulation_tree.R"
load("Calibration_tree_painted_pdup0.5_sig0.02_NULL_new.rda")

## Clearing memory
rm(Timetree_painted1)

print("memory cleaned")

## Functions to use
source( "functions_TMM.R" )
set.seed(1234)

## Setting of parameters 
noise <- 0.00 ## Random noise
pdup<-0.5 ###Proportions of duplication

#sig2<-Trait evolutionary rates
#ntips<-Number of tips
#fast<-acceleration times
#Timetree<-Original time tree
#Substree<-Substitution rates tree
#Pseudotree<-Time calibrated substitution rates tree
# acceleration_time<-the period of trait acceleration

## For this model we only need to consider calibrated null trees
timeacc_pw_model_sym <- function (sig2,ntips,fast,Timetree,Substree,Pseudotree,acceleration_time) { 
  
  ##Acceleration in trait evolutionary rate
  sigma_vector_accelaration <- c(sig2, sig2*fast)
  names(sigma_vector_accelaration) <- c("S", "D")
  PR_results <- NULL
  tree_calibrated<-0
  
  for(i in 1:length(Timetree[[1]]))
  {
    calibrated_painted <-Pseudotree[[1]][[i]] 
    
    if(!is.null(calibrated_painted)) 
    {
      print (i)
      tree<-tree1[[i]]
      Timetree_sym<- Timetree[[1]][[i]]
      tree_painted1<-NULL
      tree_painted1<- Timetree_sym
      dup_nodes<- as.numeric(Timetree[[3]][[i]])
      spe_nodes<- as.numeric(Timetree[[4]][[i]])
      dup_edges_sym<- as.numeric(Timetree[[5]][[i]])
      
      ## Considering edgelength of the duplication branches
      dup_edgelength<-tree_painted1$edge.length[which(tree_painted1$edge[,1] %in% dup_nodes)]
      
      
      tree_sim_painted <- Substree[[1]][[i]] ## as the rates remain same for speciation as well as for duplication
      bt<- branching.times(tree)
      events <- ifelse(names(bt)%in%dup_nodes,"duplication", "speciation") ##Identifying the two events
      names(events) <- names(bt)
      
      dup_desc<-tree_painted1$edge[which(tree_painted1$edge[,1] %in% dup_nodes),2]
      dup_node2edge<-tree_painted1$edge[which(tree_painted1$edge[,1] %in% dup_nodes),c(1,2)]
      test<-NULL
      test<-data.frame(node=dup_node2edge[,1],desc=dup_node2edge[,2])
      test$edge.length<-tree_painted1$edge.length[which(tree_painted1$edge[,1] %in% test$node & tree_painted1$edge[,2] %in% test$desc)]
      eventdf<-data.frame('node.event'=events)
      eventdf$node<-as.integer(row.names(eventdf))
      test<-left_join(test,eventdf,by='node')
      colnames(eventdf)[2]<-"desc" 
      test<-left_join(test,eventdf,by='desc')
      colnames(test)[4]<-'node.event'
      colnames(test)[5]<-'desc.event'
      
      ## Since, the scale of each simulated tree varies, we identified the oldest duplication time of each tree
      ## Then based on the age of oldest duplication age, we fix time so that till that time after duplication, the trait accelerates
      time<-NULL
      oldest_dup_age<-NULL
      oldest_dup_age<-as.numeric(max(bt[which(names(bt) %in%  dup_nodes)]))
      time<-round(oldest_dup_age*acceleration_time,0)
      print(paste0("trait accelerates until : ", time, " My after duplication",sep=""))
 
      ## Identifying long and short edges (based on time scale) to paint the branches accordingly
      dupnodes_long_edge<-test$desc[which(test$edge.length > time)]
      dupnodes_small_edge<-test$desc[!(test$desc %in% dupnodes_long_edge)]
      
      ##For smaller daughter branches, the whole branch will be painted with "D" state following a duplication event
      tree_painted<-NULL
      tree_painted<-paintBranches(tree, dupnodes_small_edge, state = "D",anc.state = "S")
      
      ## The previous painted tree is used to repaint longer duplication branches with two states
      ## For each daughter edgelength, segment length == time should be painted with "D" state 
      ## Rest of the segment should be pained with "S" state
      if(length(dupnodes_long_edge) > 0)
      {
        tree_painted<-paintBranches.tina(tree_painted, dupnodes_long_edge, state = "D", anc.state = "S", time)
      }
      
      ### Now we need to identify the nodes of smaller duplication edgelengths
      ### If the descendant event is a speciation event, we need to identify them to repaint it with two states
      ### the first segment (time - length of ancestral duplication branch) should be painted with "D" state 
      ## Rest of the segment should be pained with "S" state
      ## If the descendant is a tip, since the value is smaller than time, we paint it with "D" state
      test_new.df<-NULL
      test_new.df<-test[(test$desc %in%dupnodes_small_edge) & (test$desc.event=="speciation" | (is.na(test$desc.event))),]
      tree_new<-tree_painted
      for(x in 1:nrow(test_new.df))
      {
        
        i_parent_node<-test_new.df$node[x]
        i_desc_node<-test_new.df$desc[x]
        parent.edgelength<-tree_new$edge.length[which(tree_new$edge[,1]==i_parent_node & tree_new$edge[,2]==i_desc_node)]
        
        ## To check whether this duplication has happend at the terminal node (near tips)
        ## since the value is smaller than time, we paint it with a "D" state
        if(i_desc_node <= ntips) 
        {
          tree_new<-paintBranches(tree_new,i_desc_node,state = "D",anc.state = "S")
        }
        
        ##To check whether this duplication has happend at the internal node whose descendant is a speciation node
        else if((i_desc_node > ntips) & (events[which(names(events)%in%i_desc_node)]=="speciation")) ## This condition is already refined in the dataframe
        {
          # print("I am here in main\n")
          tree_new<-paintBranches(tree_new,i_desc_node,state = "D",anc.state = "S")
          new_parent.node<-i_desc_node
          tree_new<-recursive.painting(tree_new,ntips,new_parent.node,parent.edgelength,time,events)
        }
      }
      tree_painted<-tree_new
      
      
      ## Simulating trait with bounds between 0,1 on real time tree
      trait_oc<-sim.rates_bound(tree_painted, sig2=sigma_vector_accelaration)
      
      ##Getting age of all nodes
      age_all<-branching.times(tree_painted)
      
      ##Identifying tips and edges connecting tips
      tip_edge_nodes<- unique(tree_painted$edge[which(tree_painted$edge[,2] <= ntips), 1])
      tip_edges_all<- unique(tree_painted$edge[which(tree_painted$edge[,2] <= ntips), c(1,2)])
      tip_data <- NULL
      trait.label <- NULL
      tip_edge_event<-NULL
      tip_node_age<-NULL
      trait_value <-NULL
      
      ##Assigning trait data at the tip of the tree
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
      spe_count=nrow(tip_data)-dup_count
      tip_data$Event <- tip_edge_event
      tip_data$Node.Age <-tip_node_age
      tip_data$Tau <- trait_value 
      tip_data$tip.label <- names(trait_value)
      tip_data <- tip_data[c(1:4,6,5)]
      
      #' Summarize pairwise comparisons between tips of a tree
      ## Creating a data frame with one row for each pairwise combination of tips
      ##making pairwise combination of the tip events
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
      
      #Collecting results in a data frame for pairwise comparisons
      temp_results<-NULL
      temp_results<-data.frame(Event=pairwise_all$mrca.event,Tau_tip1=pairwise_all$Tau1,Tau_tip2=pairwise_all$Tau2,Time.MRCA=pairwise_all$mrca.age)
      PR_results <- rbind(PR_results,temp_results)
    } 
  }
  return (PR_results)
} 

#Calling the function
PW_null<- timeacc_pw_model_sym(0.02,100,1,Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.5)
PW_OC_0.5times<- timeacc_pw_model_sym(0.02,100,5,Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.5)
PW_OC_0.6times<- timeacc_pw_model_sym(0.02,100,5,Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.6)
PW_OC_0.7times<- timeacc_pw_model_sym(0.02,100,5,Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.7)
PW_OC_0.8times<- timeacc_pw_model_sym(0.02,100,5,Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.8)
PW_OC_0.9times<- timeacc_pw_model_sym(0.02,100,5,Timetree_0.5_painted1, Substree_0.5_0.02painted_null, Pseudotree_0.5_0.02painted_null_new,0.9) 

## Cleaning memory
rm(dt)
rm(tree1)
rm(Timetree_painted1)
rm(Timetree_0.5_painted1)
rm(Substree_0.5_0.02painted_null)
rm(Pseudotree_0.5_0.02painted_null_new)

## Saving data for further analyses and plotting
save.image("timemodel_PW_pdup0.5_sig0.02_plot_data1.rda")


############ Further analyses and plotting ###############
name<-c("null","oc0.5times","oc0.6times","oc0.7times","0c0.8times","oc0.9times")
dff<-list(PW_null,PW_OC_0.5times,PW_OC_0.6times,PW_OC_0.7times,PW_OC_0.8times,PW_OC_0.9times)

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
  filename<-paste("timemodel_",p_dup,"_pw_10000_asym_",name[x],".txt",sep="")
  Pval<- Nonlinear_pairwise2(processed_corr,filename)
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
  processed_corr$event = factor(processed_corr$Event, levels=c( "Speciation", "Duplication" ) )

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
  theme(legend.title=element_blank(),legend.position=c(0.85,0.85)) +
  theme(axis.text=element_text(size=10,face="bold")) +
  theme(legend.text=element_text(size=10,face="bold")) +
  ggtitle(expression(bold("Time model"))) +
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
rm(PW_null)
rm(PW_OC_0.5times)
rm(PW_OC_0.6times)
rm(PW_OC_0.7times)
rm(PW_OC_0.8times)
rm(PW_OC_0.9times)
rm(lp)


## Saving data for further analyses and plotting
save.image("timemodel_PW_pdup0.5_sig0.02_plot_data_latest.rda")

