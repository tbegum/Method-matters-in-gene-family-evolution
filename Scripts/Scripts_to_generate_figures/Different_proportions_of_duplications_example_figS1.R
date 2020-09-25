
##This code is written to illustrate trees with different proportions of duplications

## Loading previously stored data to save time
load("Timetree_painted_pdup0.2_sig0.02_two_ways_latest.rda")
load("Timetree_painted_pdup0.5_sig0.02_two_ways_latest.rda")
load("Timetree_painted_pdup0.8_sig0.02_two_ways_latest.rda")
#fixed_focal_spe_nodes<-c(102, 141, 153, 111, 129, 175, 194, 164, 150, 143)

#  Setting colors
colSimmap <- c("#F8766D", "#00BFC4")
names(colSimmap) <- c("S", "D")

## Plotting
par(mfrow=c(1,3))

## Tree with prop od duplication 0.2

plotSimmap(Timetree_0.2_painted1[[1]][[2]],col=colSimmap,fsize=0.01,mar=c(5,2,2,2))
#nodelabels(node=fixed_focal_spe_nodes,frame="none", pch=16, cex= 1.3, col="blue")
axisPhylo()
mtext("Times in My", side=1, line=3,font = 2, adj=0.5,cex=0.8, col=NA)
title(main = "A. Proportion of duplications: 0.2",cex.main= 1.3)
legend("topleft",legend=c("Speciation","Duplication"), ncol=1, pch=19, col=c("#F8766D","#00BFC4"), cex=1, text.font=2, bty='n')


## Tree with prop od duplication 0.5

plotSimmap(Timetree_0.5_painted1[[1]][[2]],col=colSimmap,fsize=0.01,mar=c(5,2,2,2))
#nodelabels(node=fixed_focal_spe_nodes,frame="none", pch=16, cex= 1.3, col="blue")
axisPhylo()
mtext("Times in My", side=1, line=3,font = 2, adj=0.5, cex=0.8,col=NA)
title(main = "B. Proportion of duplications: 0.5",cex.main= 1.3)


## Tree with prop od duplication 0.8

plotSimmap(Timetree_0.8_painted1[[1]][[2]],col=colSimmap,fsize=0.01, mar=c(5,2,2,2))
#nodelabels(node=fixed_focal_spe_nodes,frame="none", pch=16, cex= 1.3, col="blue")
axisPhylo()
mtext("Times in My", side=1, line=3,font = 2, adj=0.5, cex=0.8,col=NA)
title(main = "C. Proportion of duplications: 0.8",cex.main= 1.3)
