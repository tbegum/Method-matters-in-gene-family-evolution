
## This script takes into account trait evolutionary rate heterogeneity after oldest duplication for a specific period
###  This model assumes trait evolutionary rates are same after a speciation ('S'), and a duplication ('D') for a uniform null model
###  After duplication, trait evoluionary rates accelerates for certain period of time like 50% or 60% or 70% or 80% or 90%
###  After that period, trait evolves with a rate of preduplication, i.e., speciation state
## This script is for phylogenetic indeoendent contrasts analyses

## Loading libraries
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(ggplot2)
library(dplyr)
library(gtools)

## Loading previously stored data after running "Simulation_tree.R"
load("Calibration_tree_painted_pdup0.8_sig0.02_NULL_new.rda")

## Clearing memory
rm(Timetree_painted1)

print("memory cleaned")

## Functions to use
source( "functions_TMM.R" )
set.seed(1234)

## Setting of parameters 
noise <- 0.00 ## Random noise
pdup<-0.8 ###Proportions of duplication

#sig2<-Trait evolutionary rates
#ntips<-Number of tips
#fast<-acceleration times
#Timetree<-Original time tree
#Substree<-Substitution rates tree
#Pseudotree<-Time calibrated substitution rates tree
# acceleration_time<-the period of trait acceleration


## For this model we only need to consider calibrated null tress
timeacc_model_sym <- function (sig2,ntips,fast,Timetree,Substree,Pseudotree,acceleration_time) { 
  
   ##Acceleration in trait evolutionary rate
  sigma_vector_accelaration <- c(sig2, sig2*fast)
  names(sigma_vector_accelaration) <- c("S", "D")
  
  results <- array(dim=c(ntips-1, 12, length(tree1)))
  
  for(i in 1:length(Timetree[[1]])) 
  {
    
    ## For this model we need to consider the symmetrically painted tree
    calibrated_painted <-Pseudotree[[1]][[i]]
    
    if(!is.null(calibrated_painted)) 
    {
      print (i)
      tree<-tree1[[i]]
      Timetree_sym<- Timetree[[1]][[i]]
      tree_painted1<-NULL
      tree_painted1 <- Timetree_sym
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

      ## Since, the scale of each simulated tree varies, we identified oldest duplication time
      ## Then based on the age of oldest duplication age, we fix time so that till that time after duplication, the trait accelerates
      time<-NULL
      oldest_dup_age<-as.numeric(max(bt[which(names(bt) %in%  dup_nodes)]))
      time<-round(oldest_dup_age*acceleration_time,0)
      print(paste0("trait accelerates until : ", time, " My after duplication",sep=""))

      
      ## Identifying long and short edges (based on time scale) to paint branches accordingly
      dupnodes_long_edge<-test$desc[which(test$edge.length > time)]
      dupnodes_small_edge<-test$desc[!(test$desc %in% dupnodes_long_edge)]
      
      ##For small (branch length < time) daughter branches, the whole branch will be painted with "D" state following duplication
      tree_painted<-NULL
      tree_painted<-paintBranches(tree, dupnodes_small_edge, state = "D",anc.state = "S")
      
      ## The previous painted tree will be used to repaint longer duplication branches with two states
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
      ## If the descendant is tip, since the value is smaller than time, we paint it with "D" state
      test_new.df<-NULL
      test_new.df<-test[(test$desc %in%dupnodes_small_edge) & (test$desc.event=="speciation" | (is.na(test$desc.event))),]
      tree_new<-tree_painted
      for(x in 1:nrow(test_new.df))
      {
        
        i_parent_node<-test_new.df$node[x]
        
        # Considering the nodes that descend on the duplication node
        i_desc_node<-test_new.df$desc[x]
        parent.edgelength<-tree_new$edge.length[which(tree_new$edge[,1]==i_parent_node & tree_new$edge[,2]==i_desc_node)]
        
        ## To check whether this duplication has happend at the terminal node (near tips)
        ## since the value is smaller than time, we paint it with "D" state
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
      ##Assigning this newly painted tree to tree_painted
      tree_painted<-tree_new
      
      # simulating trait with bounds between 0,1 on real time tree
      trait_oc<-sim.rates_bound(tree_painted, sig2=sigma_vector_accelaration)
      
      ##Estimating the Phylogenetic Independent Contrasts (pics) for time tree, substitution tree and pseudo time tree respectively
      pic_Ttree = pic( trait_oc, tree_painted1, var.contrasts=TRUE ) 
      pic_Stree = pic( trait_oc, tree_sim_painted, var.contrasts=TRUE )
      pic_pseudoTtree = pic(trait_oc,calibrated_painted, var.contrasts = TRUE)
      
      ##Collecting result
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


#Calling the function
#For a uniform null
PIC_null<- timeacc_model_sym(0.02,100,1,Timetree_0.8_painted1, Substree_0.8_0.02painted_null, Pseudotree_0.8_0.02painted_null_new,0.5)

## Trait accelerates until 50% time since oldest duplication 
PIC_OC_0.5times<- timeacc_model_sym (0.02,100,5,Timetree_0.8_painted1, Substree_0.8_0.02painted_null, Pseudotree_0.8_0.02painted_null_new,0.5) 

## Trait accelerates until 60% time since oldest duplication
PIC_OC_0.6times<- timeacc_model_sym (0.02,100,5,Timetree_0.8_painted1, Substree_0.8_0.02painted_null, Pseudotree_0.8_0.02painted_null_new,0.6) 

## Trait accelerates until 70% time since oldest duplication
PIC_OC_0.7times<- timeacc_model_sym (0.02,100,5,Timetree_0.8_painted1, Substree_0.8_0.02painted_null, Pseudotree_0.8_0.02painted_null_new,0.7)

## Trait accelerates until 80% time since oldest duplication
PIC_OC_0.8times<- timeacc_model_sym (0.02,100,5,Timetree_0.8_painted1, Substree_0.8_0.02painted_null, Pseudotree_0.8_0.02painted_null_new,0.8) 

## Trait accelerates until 90% time since oldest duplication
PIC_OC_0.9times<- timeacc_model_sym (0.02,100,5,Timetree_0.8_painted1, Substree_0.8_0.02painted_null, Pseudotree_0.8_0.02painted_null_new,0.9)        

## Cleaning memory
rm(dt)
rm(tree1)
rm(Timetree_0.8_painted1)
rm(Substree_0.8_0.02painted_null)
rm(Pseudotree_0.8_0.02painted_null_new)

## Saving data for further analyses and plotting
save.image("Sym_timemodel_10000_0.8dup_0.02sig2_PIC_all.rda")

############ Further analyses & plotting  ###############
name<-c("null","oc0.5","oc0.6","oc0.7","oc0.8","oc0.9")
dff<-list(PIC_null,PIC_OC_0.5times,PIC_OC_0.6times,PIC_OC_0.7times,PIC_OC_0.8times,PIC_OC_0.9times)

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
   
  ##Two sided Wilcoxon tests   
  label_p_ot1<-wilcox.test(pic_duplication_Ttree, pic_speciation_Ttree)$p.value 
  if ((label_p_ot1 != 0) & (label_p_ot1 > 2.2e-16)) 
  { 
    label_p_ot1 <- format(label_p_ot1, digits= 3, scientific = TRUE) 
    if(label_p_ot1 < cutoff){label_p_ot<-paste0("P = ",label_p_ot1,star)} 
    if(label_p_ot1 >= cutoff){label_p_ot<-paste0("P = ",label_p_ot1)} 
  } 
  else{label_p_ot <- paste0("P < 2.2e-16",star)} 
  #label_p_ot<-two_tailed_wilcox(pic_duplication_Ttree, pic_speciation_Ttree) 
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
    ylim(0, 1.2) +
    theme_classic()+
    theme(legend.title=element_blank(),legend.position=c(0.85,0.85)) +
    theme(axis.text=element_text(size=10,face="bold")) +
    theme(legend.text=element_text(size=10,face="bold")) +
    annotate("text", x = 1.0, y = 0.5, label= label_p_ot, fontface = 4,size=4) +
    annotate("text", x = 2.0, y = 0.5, label= label_p_pseut, fontface = 4,size=4)
  assign(paste("time_model_PIC_plot_",p_dup,"_",name[x],sep=""), plot)
  
}

## Saving data for backup
save.image("time_model_10000_0.8dup_0.02sig2_PIC_plot_data2.rda")

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
rm(PIC_null)
rm(PIC_OC_0.5times)
rm(PIC_OC_0.6times)
rm(PIC_OC_0.7times)
rm(PIC_OC_0.8times)
rm(PIC_OC_0.9times)


## Saving data for backup
save.image("time_model_10000_0.8dup_0.02sig2_PIC_plot_data_latest.rda")
