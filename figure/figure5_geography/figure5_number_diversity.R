library(tidyverse)
library(vegan)
library(RColorBrewer)
library(itertools)
library(Hmisc)
library(scales)

fna_summary <- read_tsv("metadata/all_fna_v2_and_metainfo.tsv")

# high quality mgs
hq_genus_info <- read_tsv("metadata/profile_plot.tsv") %>%
    mutate(size = MAG + reference_genome_num) %>%
    filter(size >= 10) %>% # choose size > 10
    mutate(anim_path = str_c("top_genus/g__", capitalize(genus),
                             "/anim/", mgs_id, "/ANIm_percentage_identity.tab"),
           coverage_path = str_c("top_genus/g__", capitalize(genus),
                                 "/anim/", mgs_id, "/ANIm_alignment_coverage.tab"))
                                 
# genomes number
genus_info_num <- hq_genus_info %>%
   filter(variable == "Reference genomes" |
          variable == "Reconstructed genomes") %>%
   mutate(variable = factor(variable,
          levels = c("Reference genomes", "Reconstructed genomes"))) %>%
   mutate(strain = factor(strain))

genus_num_plot <-
   ggplot(genus_info_num,
          aes(strain, value, group = variable, shape = variable)) +
   geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
   coord_flip() +
   facet_wrap(~genus, scales = "free", nrow = 4) +
   geom_text(
      aes(label = value, y = value),
      position = position_dodge(0.9),
      vjust = 0.5,
      hjust = 1,
      size = 2.5
   ) +
   theme_bw() +
   scale_fill_npg() +
   scale_color_npg()
ggsave("fic5b_genus_genome_num_barplot.pdf", width = 26, height = 16)


# color setting
group_color <- c(
  "China_Beijing" = "#93AA00",
  "China_Shenzhen" = "#F8766D",
  "China_Yunnan" = "#D39200",
  "Fiji" = "#00BA38",
  "France" = "#00C19F",
  "Germany" = "#00B9E3",
  "Luxembourg" = "#DB72FB",
  "USA" = "#FF61C3",
  "eHOMD" = "#E31A1C",
  "Ref" = "#000000"
)

mgs_id_list <- c("mgs_3273", "mgs_3467", "mgs_3225")

for (i in 1:nrow(hq_genus_info)) {
  anim_f <- as.character(hq_genus_info[i, "anim_path"])
  cov_f <- as.character(hq_genus_info[i, "coverage_path"])
  if (file.exists(anim_f) && file.exists(cov_f)) {
     mgs_id <- as.character(hq_genus_info[i, "mgs_id"])
     if (!(mgs_id %in% mgs_id_list)) {
       next;
     }
     species <- as.character(hq_genus_info[i, "species"])
     genus <- as.character(hq_genus_info[i, "genus"])
     ani <- read.table(anim_f, row.names = 1, header = T)
     cov <- read.table(cov_f, row.names = 1, header = T)
     ani[cov < 0.3] <- 0
     ani_dist <- 1 - ani
     set.seed(0)
     outdir <- str_c("ani_plot/", genus, "/")
     if (!dir.exists(outdir)) {
       dir.create(outdir)
     }
     pca <- metaMDS(as.dist(ani_dist), k = 2)
     pcs <- as.data.frame(pca$points)

     pcs <- pca$points %>%
              as_tibble(rownames = "id_2") %>%
              left_join(fna_summary)
     pcs$group[is.na(pcs$group)] <- "Ref"
     pcs <- pcs %>%
       mutate(group2 = case_when(
                city == "Shenzhen" ~ "China_Shenzhen",
                city == "Beijing" ~ "China_Beijing",
                city == "Yunnan" ~ "China_Yunnan",
                group == "refseq" ~ "Ref",
                #group == "ehomd" ~ "eHOMD",
                group == "ehomd" ~ "Ref",
                TRUE ~ group
              ))
     pcs_group <- intersect(names(group_color), unique(pcs$group2))
     pcs <- pcs %>%
       mutate(group2 = factor(group2, levels = pcs_group))
     # View(pcs)
     mds_plot <-
         ggplot(pcs, aes(x = MDS1, y = MDS2)) +
         geom_point(aes(colour = group2), size = 2, alpha = 0.6) +
         scale_color_manual(
             name = "Source",
             labels = pcs_group,
             values = group_color[pcs_group]
         ) +
       ylab("MDS2") +
       xlab("MDS1") +
       guides(colour = guide_legend(override.aes = list(size = 5, alpha = 0.6))) +
       ggtitle(str_c(genus, " ", species)) +
       theme(axis.title.y = element_text(size = 14),
             axis.text.y  = element_text(size = 12),
             axis.title.x = element_text(size = 14),
             axis.text.x  = element_text(size = 12)) +
       theme_set(theme_bw())

     save_file <- str_c(outdir, mgs_id, ".pdf")
     ggsave(save_file, mds_plot, width = 8, height = 5)
  }
}
