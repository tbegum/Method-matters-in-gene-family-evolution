
### This program is written to generate painted trees for time dependent trait acceleretion model, which consideres trait evolutionary rate heterogeneity following duplication for a certain period
### The code generates a sample tree model for trees with different proportions of duplications with two different trait acceleration periods
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(ggplot2)
library(dplyr)

## Functions written for this study
source("functions_TMM.R")
set.seed(1234)

## sig2<-Trait evolutionary rates
## ntips<- Number of tips
## fast<-Times trait evolutionary rates faster after duplication than speciation
## Timetree<- Previously stored time tree
## Substree<- Previously stored substitution rates tree
## Pseudotree <- Previously stored time calibrated substitution rates tree
## acceleration_time<-Trait acceleration for proportion of times of oldest duplication age of the tree

## Required function for painting trees in this model
## We considered calibrated null trees for this purpose
timeacc_model_sym_ex <- function (sig2,ntips,fast,Timetree,Substree,Pseudotree,acceleration_time) { 
  
  ##Setting trait evolutionary rates for "speciation" (S) and "duplication" (D)
  sigma_vector_accelaration <- c(sig2, sig2*fast)
  names(sigma_vector_accelaration) <- c("S", "D")
  
   ## For this model we need to consider the symmetrically painted tree
   calibrated_painted <-Pseudotree[[1]][[4]]
    
   if(!is.null(calibrated_painted)) 
   {
      tree<-tree1[[4]]
      Timetree_sym<- Timetree[[1]][[4]]
      tree_painted1<-NULL
      tree_painted1 <- Timetree_sym
      dup_nodes<- as.numeric(Timetree[[3]][[4]])
      spe_nodes<- as.numeric(Timetree[[4]][[4]])
      dup_edges_sym<- as.numeric(Timetree[[5]][[4]])
      
      ## Considering edgelength of the duplication branches
      dup_edgelength<-tree_painted1$edge.length[which(tree_painted1$edge[,1] %in% dup_nodes)]
      
      tree_sim_painted <- Substree[[1]][[4]] ## as the rates remain same for speciation as well as for duplication
      bt<- branching.times(tree)
      events <- ifelse(names(bt)%in%dup_nodes,"duplication", "speciation") 
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
      max_age<-NULL
      max_age<-as.numeric(max(bt))
      
      ## if 50% time after duplication, trait accelerates
      time<-round(oldest_dup_age*acceleration_time,0)
      print(paste0("trait accelerates until : ", time, " My after duplication",sep=""))
      
      ## Identifying long and short edgelength (based on time) to paint branches accordingly
      dupnodes_long_edge<-test$desc[which(test$edge.length > time)]
      dupnodes_small_edge<-test$desc[!(test$desc %in% dupnodes_long_edge)]
      
      
      ##For duplication with small daughter branches the whole branch will be painted with "D" state
      tree_painted<-NULL
      tree_painted<-paintBranches(tree, dupnodes_small_edge, state = "D",anc.state = "S")
      
      ## Now the previous painted tree will be used to repaint longer duplication branches with two states
      ## For each daughter edgelength, segment length = time (fixed) should be painted with "D" state 
      ## Rest of the segment should be pained with "S" state
      if(length(dupnodes_long_edge) > 0)
      {
        tree_painted<-paintBranches.tina(tree_painted, dupnodes_long_edge, state = "D", anc.state = "S", time)
      }
      
      ### Now we need to identify the nodes of smaller duplication edgelengths
      ### If the descendant event is duplication we do not need to consider
      ### If the descendant event is speciation we need to identify them to repaint it with two states
      ### the first segment (time - length of ancestral duplication branch) should be painted with "D" state 
      ## Rest of the segment should be pained with "S" state
      ## If the descendant is tip, since the value is smaller than time, we paint it with "D" state
      test_new.df<-NULL
      test_new.df<-test[(test$desc %in%dupnodes_small_edge) & (test$desc.event=="speciation" | (is.na(test$desc.event))),]
      tree_new<-tree_painted
      for(x in 1:nrow(test_new.df))
      {
        
        i_parent_node<-test_new.df$node[x]
        
        # Considering the nodes that descends from the duplication node
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
      ##Assign this newly painted tree to tree_painted
      tree_painted<-tree_new
    }
  return(tree_painted)
}


## Loading previously stored tree for proportion of duplication 0.2
## We randomly condered i=4 for this example
load("~/Desktop/nwork/Empirical/Manuscript_3/data/data1/PIC0.2_new/NULL/Calibration_tree_painted_pdup0.2_sig0.02_NULL_new.rda") 

## Clearing memory
rm(Substree_0.2_0.02painted_frate_OC)
rm(Substree_0.2_0.02painted_ftrait_OC)
rm(Substree_0.2_0.02painted_OC1)
rm(Substree_0.2_0.02painted_OC2)
rm(Substree_0.2_0.02painted_OC3)
rm(Timetree_painted1)

## Calling function
tree_painted_0.2_a<-timeacc_model_sym_ex(0.02,100,5,Timetree_0.2_painted1,Substree_0.2_0.02painted_null,
                                         Pseudotree_0.2_0.02painted_null_new,0.5)
tree_painted_0.2_b<-timeacc_model_sym_ex(0.02,100,5,Timetree_0.2_painted1,Substree_0.2_0.02painted_null,
                                         Pseudotree_0.2_0.02painted_null_new,0.9)

## Loading previously stored tree for proportion of duplication 0.5
## We randomly condered i=4 for this example
load("~/Desktop/nwork/Empirical/Manuscript_3/data/data1/PIC0.5_new/NULL/Calibration_tree_painted_pdup0.5_sig0.02_NULL_new.rda") 

## Clearing memory
rm(Timetree_painted1)

## Calling function
tree_painted_0.5_c<-timeacc_model_sym_ex(0.02,100,5,Timetree_0.5_painted1,Substree_0.5_0.02painted_null,
                                         Pseudotree_0.5_0.02painted_null_new,0.5)
tree_painted_0.5_d<-timeacc_model_sym_ex(0.02,100,5,Timetree_0.5_painted1,Substree_0.5_0.02painted_null,
                                         Pseudotree_0.5_0.02painted_null_new,0.9)


## Loading previously stored tree for proportion of duplication 0.8
## We randomly condered i=4 for this example
load("~/Desktop/nwork/Empirical/Manuscript_3/data/data1/PIC0.8_new/NULL/Calibration_tree_painted_pdup0.8_sig0.02_NULL_new.rda") 

## Clearing memory
rm(Timetree_painted1)

## Calling function
tree_painted_0.8_e<-timeacc_model_sym_ex(0.02,100,5,Timetree_0.8_painted1,Substree_0.8_0.02painted_null,
                                         Pseudotree_0.8_0.02painted_null_new,0.5)
tree_painted_0.8_f<-timeacc_model_sym_ex(0.02,100,5,Timetree_0.8_painted1,Substree_0.8_0.02painted_null,
                                         Pseudotree_0.8_0.02painted_null_new,0.9)

######### Plotting Calibrated time_tree ##############
## Setting colors
colSimmap <- c("#F8766D", "#00BFC4")
names(colSimmap) <- c("S", "D")

## Plotting

tiff("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Plots/Fig2.tiff",width=700,height=500)
par(mfrow=c(3,2))
## Tree with prop of duplication 0.2
dup_nodes<-NULL
dup_nodes<- as.numeric(Timetree_0.2_painted1[[3]][[4]])
plotSimmap(tree_painted_0.2_a,col=colSimmap,fsize=0.01,mar=c(2,2,2,2))+
  nodelabels(node=dup_nodes,frame="none", pch=16, cex= 1, col="black")+
  title(main = "A. Pdup=0.2, Acceleration time = 50% of old duplication age",cex.main= 1.3)
  #legend("topleft",legend=c("Speciation","Duplication"), ncol=1, pch=19, col=c("#F8766D","#00BFC4"), cex=0.6, text.font=2, bty='n')

## Tree with prop of duplication 0.2
plotSimmap(tree_painted_0.2_b,col=colSimmap,fsize=0.01,mar=c(2,2,2,2))+
nodelabels(node=dup_nodes,frame="none", pch=16, cex= 1, col="black")+
title(main = "B. Pdup=0.2, Acceleration time = 90% of old duplication age",cex.main= 1.3)


## Tree with prop of duplication 0.5
dup_nodes<-NULL
dup_nodes<- as.numeric(Timetree_0.5_painted1[[3]][[4]])
plotSimmap(tree_painted_0.5_c,col=colSimmap,fsize=0.01,mar=c(2,2,2,2))+
  nodelabels(node=dup_nodes,frame="none", pch=16, cex= 1, col="black")+
  title(main = "C. Pdup=0.5, Acceleration time = 50% of old duplication age",cex.main= 1.3)

## for 0.5 proportion of duplication
plotSimmap(tree_painted_0.5_d,col=colSimmap,fsize=0.01,mar=c(2,2,2,2))+
  nodelabels(node=dup_nodes,frame="none", pch=16, cex= 1, col="black")+
  title(main = "D. Pdup=0.5, Acceleration time = 90% of old duplication age",cex.main= 1.3)
 
## Tree with prop of duplication 0.8
dup_nodes<-NULL
dup_nodes<- as.numeric(Timetree_0.8_painted1[[3]][[4]])
plotSimmap(tree_painted_0.8_e,col=colSimmap,fsize=0.01,mar=c(2,2,2,2))+
  nodelabels(node=dup_nodes,frame="none", pch=16, cex= 1, col="black")+
  title(main = "E. Pdup=0.8, Acceleration time = 50% of old duplication age",cex.main= 1.3)

## for 0.5 proportion of duplication
plotSimmap(tree_painted_0.8_f,col=colSimmap,fsize=0.01,mar=c(2,2,2,2))+
  nodelabels(node=dup_nodes,frame="none", pch=16, cex= 1, col="black")+
  title(main = "F. Pdup=0.8, Acceleration time = 90% of old duplication age",cex.main= 1.3)

dev.off()
