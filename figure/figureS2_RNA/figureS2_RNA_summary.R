### tianliu@genomics.cn
###figS2a Types of standard 20 amino acids decoded by tRNAs across the MAGs. ###
library(ggpubr)
MAG.dset <- read.table("TabS2-MAG.txt", header = T, sep = "\t")
MAG.HQ <- subset(MAG.dset, mimag_quality_level == "high_quality")
MAG.MQ <- subset(MAG.dset, mimag_quality_level == "medium_quality")
length(which(MAG.HQ[,"tRNA_Completeness"] >= 18))
length(which(MAG.MQ[,"tRNA_Completeness"] == 18))

### Distribution of tRNAs' complete set 
MAG.dset$mimag_quality_level <- factor(MAG.dset$mimag_quality_level, levels = c("medium_quality", "high_quality"))
p.MAG <- print(gghistogram(MAG.dset, x = "tRNA_Completeness", y =  "..density..", bins = 21, color = "mimag_quality_level", fill = "mimag_quality_level", palette = c("#00AFBB", "#E7B800")))

ggexport(p.MAG, filename = "tRNA_completeness.pdf", height = 5, width = 8)

### figS2b Three types of rRNA were predicted across the MAGs. ###
n.HQ = length(MAG.HQ$X16S_rRNA)
n.MQ = length(MAG.MQ$X16S_rRNA)

for (rRNA in c("X16S_rRNA", "X5S_rRNA", "X23S_rRNA")){
  rRNA.HQ = length(which(MAG.HQ[,rRNA]!=0))
  rRNA.MQ = length(which(MAG.MQ[,rRNA]!=0))
  print(paste("HQ", rRNA, rRNA.HQ, rRNA.HQ/n.HQ, sep = " "))
  print(paste("MQ", rRNA, rRNA.MQ, rRNA.MQ/n.MQ, sep = " "))
}

### rRNA venn diagram
### http://blog.sciencenet.cn/home.php?mod=space&uid=2985160&do=blog&id=957210
library("VennDiagram")

#### High Quality
HQ.16s <- as.character(subset(MAG.HQ, X16S_rRNA!=0)$bin_id)
HQ.5s <- as.character(subset(MAG.HQ, X5S_rRNA!=0)$bin_id)
HQ.23s <- as.character(subset(MAG.HQ, X23S_rRNA!=0)$bin_id)

venn.diagram(x=list(rRNA_16s = HQ.16s, rRNA_5s = HQ.5s, rRNA_23s = HQ.23s),
             "HQ_rRNA_venn.tiff", imagetype = "tiff",
             height = 500, width = 500, 
             resolution =300, 
             col="white", fill=c(colors()[616], colors()[38], colors()[468]), 
             alpha=c(0.6, 0.6, 0.6), lwd=c(1, 1, 1), 
             #cex=0.4, cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180), cat.cex=0)
             cex=0.4, cat.cex=0)

#### Medium Quality
MQ.16s <- as.character(subset(MAG.MQ, X16S_rRNA!=0)$bin_id)
MQ.5s <- as.character(subset(MAG.MQ, X5S_rRNA!=0)$bin_id)
MQ.23s <- as.character(subset(MAG.MQ, X23S_rRNA!=0)$bin_id)

venn.diagram(x=list(rRNA_16s = MQ.16s, rRNA_5s = MQ.5s, rRNA_23s = MQ.23s),
             "MQ_rRNA_venn.tiff", imagetype = "tiff",
             height = 500, width = 500, 
             resolution =300, 
             col="white", fill=c(colors()[616], colors()[38], colors()[468]), 
             alpha=c(0.6, 0.6, 0.6), lwd=c(1, 1, 1), 
             #cex=0.4, cat.dist=c(-0.07, -0.07, -0.05), cat.pos=c(300, 60, 180), cat.cex=0)
             cex=0.4, cat.cex=0.45)