
library(ggplot2)
library(gridExtra)

load("Asym_jump_10000_0.8dup_0.02sig2_PIC_plot_data_latest.rda")
load("Asym_jump_10000_0.8dup_0.02sig2_pw_plot_data_latest.rda")


png("all_0.8dup_model3.png", width=800, height=600)
par(mfrow=c(2,2))
## Plotting results altogether
p1_pic_model3<-jump_model_PIC_plot_0.8_null + 
  labs(title=expression(bold(paste("A: Null simulation"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold"))+
  theme(axis.text.x = element_text(colour = "black",face="bold"))

p2_pw_model3<-PW_plot_jump_model0.8_null+ 
  labs(title=expression(bold(paste("B: Null simulation"))), subtitle =  "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5))+
  theme(legend.text = element_text(size=10, face="bold"))

p3_pic_model3<-jump_model_PIC_plot_0.8_oc + 
  labs(title=expression(bold(paste("C: OC simulation"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p4_pw_model3<- PW_plot_jump_model0.8_OC + 
  labs(title=expression(bold(paste("D: OC simulation"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")



gridExtra::grid.arrange(p1_pic_model3,p2_pw_model3,
                        p3_pic_model3,p4_pw_model3,
                        ncol=2,nrow=2)


dev.off()
