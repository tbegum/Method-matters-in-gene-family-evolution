

load("Asym_PIC_pdup0.2_sig0.02_plot_data_latest.rda")
load("Asym_PW_pdup0.2_sig0.02_plot_data_latest.rda")

library(ggplot2)
library(gridExtra)

## Setting path to stored the results
#setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Plots/")

png("dup0.2_model2_S2.png",width=800,height=800)

## Plotting results altogether
p1_pic_model2<-PIC_plot_0.2_ftrait + 
  labs(title=expression(bold(paste("A: Fixed trait ( r"[T]," = 1, r"[K]," = 5)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold"))+
  theme(axis.text.x = element_text(colour = "black",face="bold"))

p2_pw_model2<-PW_plot_0.2_ftrait+ 
  labs(title=expression(bold(paste("B: Fixed trait ( r"[T]," = 1, r"[K]," = 5)"))), subtitle =  "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5))+
  theme(legend.text = element_text(size=10, face="bold"))


p3_pic_model2<-PIC_plot_0.2_OC1 + 
  labs(title=expression(bold(paste("C: Variable trait and sequence rates ( r"[T]," = 8, r"[K]," = 2)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p4_pw_model2<- PW_plot_0.2_OC1 + 
  labs(title=expression(bold(paste("D: Variable trait and sequence rates ( r"[T]," = 8, r"[K]," = 2)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

p5_pic_model2<-PIC_plot_0.2_OC2 + 
  labs(title=expression(bold(paste("E: Variable trait and sequence rates ( r"[T]," = 8, r"[K]," = 8)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p6_pw_model2<- PW_plot_0.2_OC2 + 
  labs(title=expression(bold(paste("F: Variable trait and sequence rates ( r"[T]," = 8, r"[K]," = 8)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

gridExtra::grid.arrange(p1_pic_model2,p2_pw_model2,p3_pic_model2,p4_pw_model2,p5_pic_model2,p6_pw_model2,ncol=2,nrow=3)

dev.off()  
  
