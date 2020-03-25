### tianliu@genomics.cn
### FigS3a.Human rate in noOral reads ###
library(ggplot2)
library(RColorBrewer)
setwd("~/Desktop/01project/01.oral/result/figure_scripts/github/figure/figureS3_kraken2/")
kk_human <- read.table("00.data/human_ratio_kraken2.txt",  header = T, sep = "\t")
p.S3 <- print(ggplot(kk_human, aes(x = country, y = human_mapping_rate, fill = country))
             + geom_boxplot()
             + theme_bw(base_size = 11)
             + ylab("Homo sapiens relative abundance")
             + scale_fill_brewer(palette="Set2"))
ggsave("figS1_human_noOral.pdf", p.S3, width = 7, height = 5)

L.dset <- subset(kk_human, country == "Luxembourg")
F.dset <- subset(kk_human, country == "France")
summary(L.dset$human_mapping_rate)
summary(F.dset$human_mapping_rate)
