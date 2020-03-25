### figure5 Geographical distribution of oral SGBs and strains.
### tianliu@genomics.cn

#### figure 5a. PcoA of oral SGB's profile ####
library(RColorBrewer)
library(ggplot2)
library(vegan)
library(ade4)
library(tibble)

setwd("~/Desktop/01project/01.oral/result/figure_scripts/figure5_geography/")

meta <- read.table("00.data/metadata.txt", header = T, sep = "\t")
profile <- read.table("00.data/abundance_profile_strain.jgi.tsv",header = T, sep = "\t", row.names = 1)

### filter out genomes that are not detected in all samples.
p.filter <- profile[rowSums(profile) != 0, ]
p.filter <- as.data.frame(t(p.filter))

## do pcoa
pcoa=dudi.pco(vegdist(p.filter), scann = FALSE, nf = 2) 
pcoa.dset = as.data.frame(pcoa$li)
pcoa.gg.dset = data.frame(Sample=rownames(pcoa.dset), pco1=pcoa.dset$A1, pco2=pcoa.dset$A2)
pcoa.gg.dset <- merge(pcoa.gg.dset, meta, by.x = "Sample", by.y = "Run_id")
pcoa.gg.dset$Country <- factor(pcoa.gg.dset$Country, 
                               levels = c("China_SZ", "China_YN", "China_BJ",
                                          "Fiji","France","Germany", "Luxembourg", "USA"))
## ggplot
# show_col(hue_pal()(9))
pco_p <- print(ggplot(pcoa.gg.dset, aes(x=pco1, y=pco2, color = Country))
               + geom_point(size=1, alpha=0.6)
               + guides(colour = guide_legend(override.aes = list(size=5, alpha=0.6)))
               + theme_bw()
               + xlab("PCoA1")
               + ylab("PCoA2")
               + theme(axis.title.y = element_text(size=12))
               + theme(axis.text.y = element_text(size=10))
               + theme(axis.title.x = element_text(size=12))
               + theme(axis.text.x = element_text(size=10)))
ggsave("pcoa_city.pdf", pco_p, width = 7, height = 5)

#### figure 5b. High abundance SGBs from each origin populations. ####
library(pheatmap)
profile <- read.table("00.data/abundance_profile_strain.jgi.tsv",header = T, sep = "\t", row.names = 1)
meta <- read.table("00.data/metadata.txt", header = T, sep = "\t", as.is = c(1,2,3)) 

high_ab_sgb = c()
top_num = 5
for (origin in unique(meta$Country)) {
  #origin = "USA"
  p.origin = profile[, subset(meta, Country == origin)$Run_id]
  origin.sgb = names(sort(rowSums(p.origin), decreasing = T)[1:top_num])
  high_ab_sgb = union(high_ab_sgb, origin.sgb)
}

p.high = profile[high_ab_sgb,]

### log10(relative abundance), to avoid log10(0) error, 0 was converted to 1e-8. 
p.high = log10(p.high*100 + 1e-8)

bk = unique(c(seq(-1,1, length=100)))

SGB_ann <- read.table("00.data/SGB_type.txt", header = T, sep = " ", row.names = 1)

orign_ann <- read.table("00.data/metadata.txt", header = T, sep = "\t", row.names = 1)
ann_colors = list(
  type = c(kSGB = "#FDEF96", uSGB = "#F7B71D"),
  Country = c(China_Beijing = "#B0CE66", China_Shenzhen = "#FBADA7", China_Yunnan = "#E1C066",
              Fiji = "#66D8A4", France = "#66D9DC", Germany = "#66CBFF", Luxembourg = "#DDB0FF", USA = "#FFA0E0")
)

### order by country
p.high.order = p.high[, rownames(orign_ann[order(orign_ann$Country),])]

pheatmap(p.high.order, 
         cluster_rows = T,
         cluster_cols = F,
         treeheight_row=0, 
         treeheight_col=0,
         annotation_col = orign_ann[1],
         annotation_row = SGB_ann,
         breaks = bk,
         annotation_colors = ann_colors,
         show_rownames = T, 
         show_colnames = F,
         filename = "fig5b.sgb_heatmap.png", 
         width = 8, 
         height = 6)

#### figure 5d. zhujie ####

