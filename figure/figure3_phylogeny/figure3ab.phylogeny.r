## Figure 3. Phylogeny of representative oral SGBs
### Figure 3a Taxonomic composition of the 2,313 uSGB
library(reshape2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gdata)
library(Rmisc)

setwd("~/Desktop/01project/01.oral/result/figure_scripts/figure3_phylogeny/")
dset <- read.table("00.data/uSGB_taxon_summary.txt", sep="\t", header = F)
colnames(dset) = c("Rank", "Taxon", "Counts")
dset$Percentage = dset$Counts/dset[dset$Rank=="Kingdom","Counts"]*100

# order datasets
rank_order = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
dset$Rank = reorder.factor(dset$Rank, new.order=rank_order)
dset = dset[order(dset$Rank,dset$Counts),]

# label other taxa
ranks_selected = c("Phylum", "Class", "Order", "Family", "Genus")
rm(new.dset)
for (r in 1:length(ranks_selected)){
  prefix = tolower(substr(ranks_selected[r], 1, 1))
  dset.rank = dset[dset$Rank == ranks_selected[r],]
  dset.rank.others = dset.rank[1:(nrow(dset.rank)-5),]
  dset.newrow = data.frame(Rank=ranks_selected[r], Taxon=paste(prefix,"__Other", sep=""), Counts=sum(dset.rank.others$Counts), 
                           Percentage=sum(dset.rank.others$Percentage), Colour="grey")
  dset.urow = data.frame(Rank=ranks_selected[r], Taxon=paste(prefix,"__Unclassified", sep=""), Counts= 2313 - sum(dset.rank$Counts), 
                         Percentage= 100 - sum(dset.rank$Percentage), Colour="white")
  dset.rank.top = dset.rank[(nrow(dset.rank)-4):nrow(dset.rank),]
  dset.toprow = data.frame(Rank=ranks_selected[r], Taxon=dset.rank.top$Taxon, Counts=dset.rank.top$Counts, 
                           Percentage=dset.rank.top$Percentage, Colour=rev(brewer.pal(5,"Set2"))) 
  if (!exists("new.dset")) {
    new.dset = rbind(dset.urow, dset.newrow, dset.toprow)
  } else {
    new.dset = rbind(new.dset, dset.urow, dset.newrow, dset.toprow)
  }
}
new.dset$Taxon = factor(new.dset$Taxon, levels=new.dset$Taxon)
new.dset$Rank = factor(new.dset$Rank, levels=c("Phylum", "Class", "Order", "Family", "Genus"))
# Phylum counts
#print(dset[which(dset$Rank == "Phylum"),])

# plot stacked plot taxa counts
uSGB_taxon <- print(ggplot(new.dset, aes(x=Rank, y=Percentage, fill=Taxon)) 
      + geom_bar(stat="identity", colour="black", alpha=0.5, size=0.2)
      + theme_bw()
      + ylab("Proportion (%)")
      + scale_fill_manual(values=as.vector(new.dset$Colour))
      #+ scale_x_discrete(limits=ranks_selected)
      + ylim(0,100)
      #+ coord_flip()
      + theme(axis.title.y = element_text(size=16))
      + theme(axis.text.y = element_text(size=16))
      + theme(axis.title.x = element_blank())
      + theme(panel.grid =element_blank(), panel.border = element_blank())
      + theme(axis.text.x = element_text(size=16, hjust=0.5, vjust=0.5)))
ggsave("01.figure/figure3a.uSGB_taxa.pdf", uSGB_taxon, width = 9, height = 8)

### Figure 3b Proportion of the total phylogenetic diversity provided by the uSGB
kSGB_p = read.table("00.data/kSGB_phylum.tsv", header = F, sep = "\t")
uSGB_p = read.table("00.data/uSGB_phylum.tsv", header = F, sep = "\t")
names(kSGB_p) <- c("Type", "Phylum", "count")
names(uSGB_p) <- c("Type", "Phylum", "count")
taxon_merge = merge(uSGB_p, kSGB_p, by = "Phylum", all=TRUE)
taxon_merge[is.na(taxon_merge)] = 0
taxon_merge$sum <- taxon_merge$count.x + taxon_merge$count.y

# phylum whose SGBs' number > 10.
taxon_merge <- subset(taxon_merge, sum > 10)
taxon_merge$percent <- taxon_merge$count.x / (taxon_merge$count.x + taxon_merge$count.y)
taxon_merge$Phylum = factor(taxon_merge$Phylum, levels = rev(taxon_merge$Phylum))

taxon_merge$Phylum <- factor(taxon_merge$Phylum, levels = taxon_merge[order(taxon_merge$sum),]$Phylum)
uSGB_taxon_e_p = print(ggplot(taxon_merge, aes(x = Phylum, y = percent, fill = Phylum)) 
                       + geom_bar(stat="identity")
                       + theme_bw(base_size = 16) 
                       + coord_flip() 
                       + guides(fill=FALSE)
                       #+ scale_fill_brewer(palette="Set3") 
                       + scale_fill_brewer(palette="Paired") 
                       + xlab("") + ylab(""))
ggsave("01.figure/figure3b.uSGB.proportion.pdf", uSGB_taxon_e_p, width = 7, height = 5)