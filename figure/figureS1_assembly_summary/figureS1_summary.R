### tianliu@genomics.cn
### Supplementary Figure 1. Summary of assembly quality, SGBs distribution across 7 studies

####a. Numbers of medium quality and high quality MAGs in samples from 7 studies####
library(ggplot2)
setwd("~/Desktop/01project/01.oral/result/figure_scripts/github/figure/figureS1_assembly_summary/")

### The samples without HQ/MQ MAGs were not show in the sample_unknown_MAGs.tsv
MAG <- read.table("../figure1_unknownness/00.data/sample_unknown_MAGs.tsv", header = T, sep = "\t")
rmhost <- read.table("00.data/rmhsot_reads.txt", header = T, sep = "\t")

a.dset <- merge(MAG, rmhost, by = "sample_id")
a.dset$oral_site <- factor(a.dset$oral_site, levels = c("saliva", "tongue", "dental"))

p_MAG_study <- print(
  ggplot(a.dset, aes( x = MAGs.HQ.MQ., y = rmhost_Gbase, color = studys, shape = oral_site)) 
  #+ geom_violin(draw_quantiles = c( 0.5)) 
  + geom_point(alpha = 0.8, size = 1.2)
  + xlab("Number of high quality & medium quality MAGs reconstructed from a single sample") 
  + ylab("Bases of sample after remove host (Gb)")
  + guides(fill=FALSE) 
  + theme_bw(base_size = 12))
ggsave("S1a_MAGs_in_studys.pdf", p_MAG_study, width = 9, height = 6)

#### c. The number of all SGBs (>0.001 abundance) for each samples. ####
profile_file = "../../../../abundance_profile_strain.jgi.tsv"
meta <- read.table("00.data/rmhsot_reads.txt", header = T, sep = "\t")

pf <- read.table(profile_file, header = T, sep = "\t")
pf_b <- pf
pf_b[pf_b<=0.001] = 0
pf_b[pf_b>0.001]  = 1 
SGB_sum = colSums(pf_b[,-1])
SGB_sum.dset <- data.frame("Sample_ID" = names(SGB_sum), "SGB_num" = SGB_sum)
figc.dset <- merge(SGB_sum.dset, meta[,c(1,2)], by.x = "Sample_ID", by.y = "sample_id")
p.c <- print(ggplot(figc.dset, aes(x = Study, y = SGB_num, fill = Study)) 
             + geom_violin()
             + xlab("number of SGBs per sample")
             + guides(fill=FALSE)
             + theme_bw(base_size = 16) 
             + coord_flip())

ggsave("SGB_number_per_sample.pdf", p.c, width = 6, height = 6)

#### d. uSGB richness (number of uSGB/number of all SGB) for each sample ####
uSGB <- read.delim("00.data/uSGB_id.txt", header = F)
names(uSGB) <- "lineages_strain_new"
uSGB_pf <- merge(uSGB, pf, by = "lineages_strain_new")
uSGB_pf_b <-uSGB_pf
uSGB_pf_b[uSGB_pf_b<=0.001] = 0
uSGB_pf_b[uSGB_pf_b>0.001]  = 1 
uSGB_sum = colSums(uSGB_pf_b[,-1])
uSGB_sum.dset <- data.frame("Sample_ID" = names(uSGB_sum), "uSGB_num" = uSGB_sum)
figd.dset <- merge(uSGB_sum.dset, figc.dset, by = "Sample_ID")
figd.dset$uSGB_richnness = figd.dset$uSGB_num / figd.dset$SGB_num

p.d <- print(ggplot(figd.dset, aes(x = Study, y = uSGB_richnness, fill = Study)) 
             + geom_violin()
             + xlab("uSGBs richness per sample")
             + guides(fill=FALSE)
             + theme_bw(base_size = 16) 
             + coord_flip())

ggsave("Sd_uSGB_richnness_per_sample.pdf", p.d, width = 6, height = 6)

#### e. Sum of all uSGB abundance for each samples. ####
### uSGB_num_ab
usum = colSums(uSGB_pf[,-1])
usum.dset = data.frame("Sample_ID" = names(usum), "uSGB_abundance" = usum)

uSGB_sum.dset <- merge(usum.dset, meta[,c(1,2)], by.x = "Sample_ID", by.y = "sample_id")

p.e <- print(ggplot(uSGB_sum.dset, aes(x = Study, y = uSGB_abundance, fill = Study)) 
                    + geom_violin()
                    + guides(fill=FALSE)
                    + theme_bw(base_size = 16) 
                    + coord_flip())
ggsave("Se.uSGB_abundance.pdf", p.e, width = 6, height = 6)

