library(pheatmap)
library(RColorBrewer)
setwd("~/tmp/01.oral/04.new/pyani/")

BGI_ani = read.table("BGI_NewOrder_ani_ave.txt", header = T, row.names = 1, sep = "\t")

o = rownames(BGI_ani)
hc = hclust(as.dist(1 - BGI_ani))
BGI_ani = BGI_ani[hc$order, hc$order]
BGI_ani[upper.tri(BGI_ani)] = NA
BGI_ani = BGI_ani[o, o]

pheatmap(BGI_ani, 
         cluster_rows = hc,
         cluster_cols = hc,
         treeheight_row=0, 
         treeheight_col=50,
         show_rownames = T, 
         show_colnames = F,
         display_numbers = TRUE,
         filename = "NewOrder_ani.pdf", 
         width = 8, 
         height = 6)