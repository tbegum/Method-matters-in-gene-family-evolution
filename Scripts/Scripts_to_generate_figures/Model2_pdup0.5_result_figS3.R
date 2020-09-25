

load("Asym_PIC_pdup0.5_sig0.02_plot_data_latest.rda")
load("Asym_PW_pdup0.5_sig0.02_plot_data_latest.rda")

library(ggplot2)
library(gridExtra)

## Plotting

png("Pdup0.5_all_in_one.png", width=1500, height=1000)
par(mfrow=c(3,4))
## Plotting results altogether
p1_pic_model2<-PIC_plot_0.5_null + 
  labs(title=expression(bold(paste("A: Null simulation (r"[T]," = 1, r"[K]," = 1)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold"))+
  theme(axis.text.x = element_text(colour = "black",face="bold"))

p2_pw_model2<-PW_plot_0.5_null+ 
  labs(title=expression(bold(paste("B: Null simulation (r"[T]," = 1, r"[K]," = 1)"))), subtitle =  "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5))+
  theme(legend.text = element_text(size=10, face="bold"))


p3_pic_model2<-PIC_plot_0.5_ftrait + 
  labs(title=expression(bold(paste("C: Fixed trait (r"[T]," = 1, r"[K]," = 5)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold",hjust = 0.5)) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p4_pw_model2<- PW_plot_0.5_ftrait + 
  labs(title=expression(bold(paste("D: Fixed trait (r"[T]," = 1, r"[K]," = 5)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(legend.position="none")

p5_pic_model2<-PIC_plot_0.5_frate + 
  labs(title=expression(bold(paste("E: Fixed rate (r"[T]," = 5, r"[K]," = 1)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p6_pw_model2<- PW_plot_0.5_frate + 
  labs(title=expression(bold(paste("F: Fixed rate (r"[T]," = 5, r"[K]," = 1)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold"))

p7_pic_model2<-PIC_plot_0.5_OC1 + 
  labs(title=expression(bold(paste("G: Variable trait and sequence rates (r"[T]," = 8, r"[K]," = 2)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p8_pw_model2<- PW_plot_0.5_OC1 + 
  labs(title=expression(bold(paste("H: Variable trait and sequence rates (r"[T]," = 8, r"[K]," = 2)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold"))

p9_pic_model2<-PIC_plot_0.5_OC2 + 
  labs(title=expression(bold(paste("I: Variable trait and sequence rates (r"[T]," = 8, r"[K]," = 8)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p10_pw_model2<- PW_plot_0.5_OC2 + 
  labs(title=expression(bold(paste("J: Variable trait and sequence rates(r"[T]," = 8, r"[K]," = 8)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold")) 

p11_pic_model2<-PIC_plot_0.5_OC3 + 
  labs(title=expression(bold(paste("K: Variable trait and sequence rates (r"[T]," = 2, r"[K]," = 8)"))), subtitle = "PIC method") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(colour = "black",face="bold"))+
  theme(legend.position="none")

p12_pw_model2<- PW_plot_0.5_OC3 + 
  labs(title=expression(bold(paste("L: Variable trait and sequence rates (r"[T]," = 2, r"[K]," = 8)"))), subtitle = "Pairwise comparisons") +
  theme(plot.title = element_text(color="blue", size=12, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 12, face="bold")) 


gridExtra::grid.arrange(p1_pic_model2,p2_pw_model2,p3_pic_model2,
                        p4_pw_model2,p5_pic_model2,p6_pw_model2,
                        p7_pic_model2,p8_pw_model2,p9_pic_model2,
                        p10_pw_model2,p11_pic_model2,p12_pw_model2,ncol=4,nrow=3)


dev.off()

