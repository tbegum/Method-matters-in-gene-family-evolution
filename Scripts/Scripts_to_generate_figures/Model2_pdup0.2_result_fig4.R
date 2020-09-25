library(ggplot2)
library(gridExtra)

load("Asym_PIC_pdup0.2_sig0.02_plot_data_latest.rda")
load("Asym_PW_pdup0.2_sig0.02_plot_data_latest.rda")


png("/Users/admin/Desktop/dup0.2_model2_fig4.png", width=800, height=800)
par(mfrow=c(3,2))

## Plotting results altogether
p1_pic_model1<-PIC_plot_0.2_null + 
  labs(title=expression(bold(paste("A: Null simulation (r"[T]," = 1, r"[K]," = 1)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold"))+
  theme(axis.text.x = element_text(colour = "black",face="bold"))

p2_pw_model1<-PW_plot_0.2_null+ 
  labs(title=expression(bold(paste("B: Null simulation (r"[T]," = 1, r"[K]," = 1)"))), subtitle =  "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5))+
  theme(legend.text = element_text(size=10, face="bold"))

p3_pic_model1<-PIC_plot_0.2_frate + 
  labs(title=expression(bold(paste("C: Fixed rate (r"[T]," = 5, r"[K]," = 1)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p4_pw_model1<- PW_plot_0.2_frate + 
  labs(title=expression(bold(paste("D: Fixed rate (r"[T]," = 5, r"[K]," = 1)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

p5_pic_model1<-PIC_plot_0.2_OC3 + 
  labs(title=expression(bold(paste("E: Variable trait and sequence rates (r"[T]," = 2, r"[K]," = 8)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p6_pw_model1<- PW_plot_0.2_OC3 + 
  labs(title=expression(bold(paste("F: Variable trait and sequence rates (r"[T]," = 2, r"[K]," = 8)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

gridExtra::grid.arrange(p1_pic_model1,p2_pw_model1,p3_pic_model1,p4_pw_model1,p5_pic_model1,p6_pw_model1,ncol=2,nrow=3)

dev.off()
