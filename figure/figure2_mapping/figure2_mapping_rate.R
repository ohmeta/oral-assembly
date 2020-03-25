### Figure2 mapping_rate
### tianliu@genomics.cn

map_rate <- unique(read.table("SupTable4_mapping_rate.txt",  header = T, sep = "\t"))

library(ggplot2)

### 2a country
map_rate$database <- factor(map_rate$database, levels = c("eHOMD", "hSGB_Rep", "oralSGB_Rep", "oralMAG", "oralRef"))
p_country <- print(ggplot(subset(map_rate), aes(x=country, y=mapping_rate, fill=database)) 
                    + geom_boxplot(outlier.size = 0.1, alpha=0.6, position=position_dodge(0.6), width=0.5)
                    + scale_y_continuous(breaks=seq(0, 1, 0.2))
                    + theme_bw() 
                    #+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
                    + xlab("") + scale_fill_manual(values=c("steelblue", "darkgreen","yellow","#E69F00"))
                    + theme(axis.text.y = element_text(size=11)) + theme(axis.text.x = element_text(size=11)))
ggsave("fig2a_country.pdf", p_country, width = 8, height = 5)

### 2b diagnosis
map_rate[map_rate == "type 1 diabetes"] = "diabetes"
map_rate[map_rate == "type 2 diabetes"] = "diabetes"
map_rate$diagnosis <- as.factor(map_rate$diagnosis)
sub.map = subset(map_rate, diagnosis %in% c("colorectal cancer", "diabetes", "pregnant", "rheumatoid arthritis"))

p_diagnosis <- print(ggplot(subset(sub.map), aes(x=diagnosis, y= mapping_rate, fill=database)) 
                      + geom_boxplot(outlier.size = 0.1, alpha=0.6, position=position_dodge(0.6), width=0.5)
                      + scale_y_continuous(breaks=seq(0, 1, 0.2))
                      + theme_bw() 
                      #+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
                      + xlab("") + scale_fill_manual(values=c("steelblue", "darkgreen","yellow","#E69F00"))
                      + theme(axis.text.y = element_text(size=11)) + theme(axis.text.x = element_text(size=11)))
ggsave("fig2b_diagnosis.pdf", p_diagnosis, width = 6, height = 5)