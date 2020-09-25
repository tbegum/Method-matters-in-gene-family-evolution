

load("time_model_10000_0.2dup_0.02sig2_PIC_plot_data_latest.rda")
load("timemodel_PW_pdup0.2_sig0.02_plot_data_latest.rda")

library(ggplot2)
library(gridExtra)

## Setting path to stored the results
setwd("/Users/admin/Desktop/")

png("time_model_0.2.png",width=800,height=600)
par(mfrow=c(3,3))
## Plotting results altogether
p1_pic_model1<-time_model_PIC_plot_0.2_null + 
  labs(title=expression(bold(paste("A: Null model"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold"))+
  theme(axis.text.x = element_text(colour = "black",face="bold"))

p2_pw_model1<-PW_plot_0.2_null+ 
  labs(title=expression(bold(paste("B: Null model"))), subtitle =  "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5)) 
# theme(axis.text.x = element_blank())

p3_pic_model1<-time_model_PIC_plot_0.2_oc0.5 + 
  labs(title=expression(bold(paste("C: Trait acceleration for 50% time of old duplication age"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p4_pw_model1<- PW_plot_0.2_oc0.5times + 
  labs(title=expression(bold(paste("D: Trait acceleration for 50% time of old duplication age"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

p5_pic_model1<-time_model_PIC_plot_0.2_oc0.9 + 
  labs(title=expression(bold(paste("E: Trait acceleration for 90% time of old duplication age"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

p6_pw_model1<-PW_plot_0.2_oc0.9times + 
  labs(title=expression(bold(paste("F: Trait acceleration for 90% time of old duplication age"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

gridExtra::grid.arrange(p1_pic_model1,p2_pw_model1,p3_pic_model1,p4_pw_model1,p5_pic_model1,p6_pw_model1,ncol=2,nrow=3)

dev.off()  
  