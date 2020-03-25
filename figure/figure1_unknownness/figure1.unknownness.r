library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

### figure1c Genome numbers distribution of uSGBs and kSGBs
SGB_size <- read.table("01.fig1_summary/mgs_size.summary.csv", header = F)
names(SGB_size) <- c("Genomes", "Type", "SGBs")
SGB_size$Genomes <- factor(SGB_size$Genomes, levels = c("1","2","3","4","5","6-10","11-20","21-50",">51"))
SGB_size$Type <- factor(SGB_size$Type, levels = c( "kMGS_k", "kMGS","uMGS"))
size_p <- ggplot(SGB_size, aes(x = Genomes, y = SGBs, fill = Type)) + geom_bar(stat="identity") + theme_bw(base_size = 16)+ theme(panel.grid =element_blank(), panel.border = element_blank())+ guides(fill=FALSE) 
ggsave(paste0("01.fig1_summary/SGB_size.pdf"), size_p, width = 6, height = 5)

### figure1d Distribution of the fraction of uMAGs in each sample by oral sites and country
#setwd("~/tmp/01.oral/result/figture/")
sample_sgb <- read.delim("00.data/sample_unknown_MAGs.tsv",header = T)
sample_sgb$unknown_ratio <- sample_sgb$uMAGs / sample_sgb$MAGs.HQ.MQ.

p_sample_country <- ggplot(sample_sgb, aes(x = country, y = unknown_ratio, fill = country)) + geom_violin(draw_quantiles = c( 0.5))  + xlab("") + guides(fill=FALSE) + theme_bw(base_size = 16) + coord_flip() + theme(panel.grid =element_blank(), panel.border = element_blank())

p_sample_study <- ggplot(sample_sgb, aes(x = studys, y = unknown_ratio, fill = studys)) + geom_violin(draw_quantiles = c( 0.5))  + theme_bw(base_size = 16) + xlab("") + guides(fill=FALSE) + coord_flip() + theme(panel.grid =element_blank(), panel.border = element_blank())

p_sample_site <- ggplot(sample_sgb, aes(x = oral_site, y = unknown_ratio, fill = oral_site)) + geom_violin(draw_quantiles = c( 0.5)) + theme_bw(base_size = 16) + xlab("") + guides(fill=FALSE) + coord_flip() + theme(panel.grid =element_blank(), panel.border = element_blank())

ggsave("unknown_ratio_country.pdf", p_sample_country, width = 6, height = 6)
ggsave("unknown_ratio_study.pdf", p_sample_study, width = 6, height =7)
ggsave("unknown_ratio_site.pdf", p_sample_site,width = 5.5, height = 3 )
