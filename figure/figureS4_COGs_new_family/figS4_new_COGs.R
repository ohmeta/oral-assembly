### tianliu@genomics.cn

### figure S4b COG functional annotation from all MAGs in the new candidatus family.
library(ggplot2)
library(RColorBrewer)

setwd("~/Desktop/01project/01.oral/result/figure_scripts/github/figure/figureS4_COGs_new_family/")
# Information
info_cog <- c("J","A","K","L","B")
info_color <- brewer.pal(5, "BuGn")

# CELLULAR PROCESSES AND SIGNALING
cell_cog <- c("D","Y","V","T","M","N","Z","W","U","O","X")
cols<-brewer.pal(3, "YlOrRd")
pal<-colorRampPalette(cols)
cell_color <- pal(11)

# METABOLISM
meta_cog = c("C","G","E","F","H","I","P","Q")
meta_color = brewer.pal(8, "Purples")

# POORLY CHARACTERIZED
poor_cog <- c("R","S","nohit")
poor_color <- brewer.pal(3, "Greys")

cog_color <- c(info_color, cell_color, meta_color, poor_color)
names(cog_color) <- c(info_cog, cell_cog, meta_cog, poor_cog)

### barplot
dset = read.table("00.data/new_family_cog.txt", header = T, sep = "\t")
dset$COG_id <- factor(dset$COG_id, levels = names(cog_color))

bgi_names = read.table("00.data/SGB.txt", header = F, sep = "\t",stringsAsFactors = F)
bgi_names.sub = bgi_names[,c(1,4)]
names(bgi_names.sub) <- c("SGB_id","bin_id")

sgb_level = c("sgb_2766", "sgb_2064", "sgb_1640", "sgb_1776", "sgb_1001", "sgb_2141", "sgb_1490", "sgb_2021", "sgb_1788", "sgb_1516", "sgb_2790")
bgi_names.sub$SGB_id <- factor(bgi_names.sub$SGB_id, levels = sgb_level)
bgi_names.sub <- bgi_names.sub[order(bgi_names.sub[,1]),]
dset$sample_id <- factor(dset$sample_id, levels = bgi_names.sub$bin_id)
p = print(ggplot(dset,aes(sample_id,gene_num,fill=COG_id)) 
          + geom_bar(stat="identity",position="fill")
          + scale_fill_manual(values = cog_color)
          + theme_bw(base_size = 10)
          + xlab("") + ylab("gene_rate")
          + coord_flip()
          + theme(panel.grid =element_blank(), panel.border = element_blank()))
ggsave("new_family_cog.pdf", p, width = 12, height = 6)