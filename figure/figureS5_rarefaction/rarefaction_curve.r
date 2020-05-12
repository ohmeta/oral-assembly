library(tidyverse)
library(tidymodels)
library(ggthemr)
library(ggthemes)
library(RColorBrewer)
library(stats)
library(phyloseq)
library(vegan)
library(ggsci)
library(here)

# figure4b
change_rare <- function(df) {
   # transform (thanks Jie Zhu Ye)
   dt <- as.data.frame(df)
   dt$id <- paste(dt$species, dt$source)
   # aggregate(dt$sample_num,
   #          list(label = dt$id),
   #          function(x)
   #          {
   #             names(table(x))
   #          }
   # )
   vars <- unique(dt$id)
   fdt <- c()
   for (var in vars)
   {
      dt1 <- dt[dt$id == var, ]
      if (length(unique(dt1$sample_num)) >= 3) {
         fm3DNase1 <- nls(gene_count ~ SSasymp(sample_num, Asym, R0, lrc),
                          data = dt1)
         dt1$pred <- predict(fm3DNase1)
         dt1$sample_num <- dt1$sample_num / max(dt1$sample_num)
         head <- dt1[1:2, ]
         head$sample_num <- 0
         head$gene_count <- 0
         # head$sample_num = 1
         # head$gene_count = ?
         head$pred <- 0
         dt1 <- rbind(head, dt1)
         fdt <- rbind(fdt, dt1)
      }
   }
   df <- fdt %>%
      arrange(gene_count) %>%
      filter(gene_count > 0)
   return(df)
}

min_max <- function(x) {
   print(min(x$gene_count))
   print(max(x$gene_count))
}

# http://www.sthda.com/english/wiki/ggplot2-axis-scales-and-transformations
logp3_trans <- trans_new(
   name = "logp",
   trans = function(x) log(x),
   inverse = function(x) exp(x),
   breaks = log_breaks()
)

plot_rare_curve <- function(df, genus_name, outdir) {
   curve <-
      ggplot(df, aes(sample_num / sum(sample_num), gene_count, group = id)) +
      geom_line(
         size = 2,
         aes(
            x = sample_num,
            y = pred,
            linetype = source,
            color = species,
            #fill = species,
            group = id
         )
      ) +
      scale_colour_brewer(palette = "Set3") +
      scale_fill_brewer(palette = "Spectral") +
      scale_y_continuous(trans = logp3_trans) + # ,
      # breaks = trans_breaks("log2", function(x) 2^x),
      # labels = trans_format("log2", math_format(2^.x)))
      # theme(
      #   legend.position = c(.05, .95),
      #   legend.justification = c("left", "top"),
      #   legend.box.just = "left",
      #   legend.margin = margin(3, 3, 3, 3)) +
      xlab("Genomes") +
      ylab("Genes in pangenome") +
      theme(axis.text.x = element_text(colour = "black", size = 6)) +
      # facet_wrap(genus ~ ., scales = 'free') +
      ggtitle(str_c(genus_name, " pangenome rarefaction curve plot")) +
      #scale_color_npg() +
      #scale_fill_npg()
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            #legend.position = 'top'.
            axis.text.x = element_text(size = 12),
            legend.text = element_text(size = 12))

   print(curve)
   ggsave(str_c(outdir, "/ficS6_", genus_name,
                "_pangenome_rarefaction_curve_plot.pdf"),
          width = 12, height = 8)
   return(curve)
}

plot_rare_curve2 <- function(df, genus_name) {
   curve <-
      ggplot(df,
             aes(sample_num, gene_count,
                 group = interaction(species, source))) +
      geom_smooth(
         method = "nls",
         formula = y ~ SSasymp(x, Asym, R0, lrc), se = F, size = 0.5,
         aes(
            linetype = source, color = species, fill = species,
            group = interaction(species, source)
         )
      ) +
      scale_colour_brewer(palette = "Set3") +
      # theme(
      #   legend.position = c(.05, .95),
      #   legend.justification = c("left", "top"),
      #   legend.box.just = "left",
      #   legend.margin = margin(3, 3, 3, 3)) +
      xlab("Genomes") +
      ylab("Genes in pangenome") +
      theme(axis.text.x = element_text(colour = "black", size = 10)) +
      ggtitle(str_c("fic4b_", genus_name, " pangenome rarefaction curve plot")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"))
      #scale_fill_npg() +
      #scale_color_npg()
   print(curve)
   return(curve)
}

# 'g__Prevotella'
# 'g__Neisseria'
# 'g__Streptococcus'
# 'g__Veillonella'
# 'g__Porphyromonas'
# 'g__Fusobacterium'
# 'g__Acinetobacter'
# 'g__Actinomyces'
# 'g__Pauljensenia'
# 'g__Haemophilus_D'
# 'g__Rothia'
# 'g__F0040'

# rare_all.tsv is tool big
rare_df <- read_tsv("rare_all.tsv") %>%
   mutate(
      source = as.factor(source),
      size = MAG + reference_genome_num
   )

rare_prevotella <- change_rare(rare_df %>% filter(genus == "Prevotella"))
rare_neisseria <- change_rare(rare_df %>% filter(genus == "Neisseria"))
rare_streptococcus <- change_rare(rare_df %>% filter(genus == "Streptococcus"))
rare_veillonella <- change_rare(rare_df %>% filter(genus == "Veillonella"))
rare_porphyromonas <- change_rare(rare_df %>% filter(genus == "Porphyromonas"))
rare_fusobacterium <- change_rare(rare_df %>% filter(genus == "Fusobacterium"))
rare_pauljensenia <- change_rare(rare_df %>% filter(genus == "Pauljensenia"))
rare_haemophilus_D <- change_rare(rare_df %>% filter(genus == "Haemophilus_D"))
# rare_acinetobacter <- change_rare(rare_df %>% filter(genus == "Acinetobacter"))
# rare_actinomyces <- change_rare(rare_df %>% filter(genus == "Actinomyces"))

plot_rare_curve(rare_prevotella,
                "Prevotella", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_neisseria,
                "Neisseria", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_streptococcus,
                "Streptococcus", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_veillonella,
                "Veillonella", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_porphyromonas,
                "Porphyromonas", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_fusobacterium,
                "Fusobacterium", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_pauljensenia,
                "Pauljensenia", "figure4/rarefaction_curve_3")
plot_rare_curve(rare_haemophilus_D,
                "Haemophilus_D", "figure4/rarefaction_curve_3")

###
plot_rare_curve2(rare_df %>% filter(genus == "Prevotella"), "Prevotella")
plot_rare_curve2(rare_df %>% filter(genus == "Neisseria"), "Neisseria")
plot_rare_curve2(rare_df %>% filter(genus == "Streptococcus"), "Streptococcus")
plot_rare_curve2(rare_df %>% filter(genus == "Veillonella"), "Veillonella")
plot_rare_curve2(rare_df %>% filter(genus == "Porphyromonas"), "Porphyromonas")
plot_rare_curve2(rare_df %>% filter(genus == "Fusobacterium"), "Fusobacterium")
plot_rare_curve2(rare_df %>% filter(genus == "Pauljensenia"), "Pauljensenia")
plot_rare_curve2(rare_df %>% filter(genus == "Haemophilus_D"), "Haemophilus_D")
