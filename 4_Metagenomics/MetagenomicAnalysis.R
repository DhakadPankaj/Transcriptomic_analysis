library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(vegan)
library(ggrepel)
library(tidyr)
library(data.table)
library(Polychrome)


setwd("Z:/Fly_infections/metagenomics/")


# Load the file using fread
blast_dt <- fread("All_diamond_results.tsv", sep = "\t", header = FALSE)

colnames(blast_dt) <- c("sample","qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                    "stitle", "sskingdoms", "skingdoms", "species")

name2taxid <- fread("name2taxid.txt", sep = "\t", header = FALSE)
colnames(name2taxid) <- c("species", "taxid")

#merge the two data tables on the "species" column
blast_dt <- merge(blast_dt, name2taxid, by = "species", all.x = TRUE)

# taxonomy table
tax_dt <- fread("taxa_classification.txt", sep = "\t", header = FALSE)
colnames(tax_dt) <- c("taxid", "taxonomy")

#merge the two data tables on the "taxid" column
blast_dt <- merge(blast_dt, tax_dt, by = "taxid", all.x = FALSE)

# remove species column
blast_dt[, species := NULL]

# split taxonomy into separate columns
blast_dt <- blast_dt[, c("kingdom", "phylum", "class", "order", "family", "genus", "species") := tstrsplit(taxonomy, ";", fixed = TRUE)] 

# remove taxonomy column
blast_dt[, taxonomy := NULL]

# get label for each sample (DEFL1P5FN_NODE_9571_length_610_cov_1.702048_g8797_i0_1) get "DEFL1P5FN" and a column "label"
blast_dt[, label := gsub("(_NODE_.*)", "", qseqid)]


# read read counts for each contigs in all species
read_counts <- fread("All_contigs.idxstats.tsv", sep = "\t", header = FALSE)
colnames(read_counts) <- c("ID", "contig", "length", "mapped", "unmapped")

read_counts[, label := gsub("(_NODE_.*)", "", contig)]

# Calculate total mapped reads for each label
read_counts[, total_mapped := sum(mapped), by = label]

# Calculate proportion of mapped reads for each contig within each label
read_counts[, proportion_mapped := mapped / total_mapped]

# drop unnecessary columns "length" "Unmapped" "total_mapped"
read_counts[, c("ID", "length", "unmapped", "total_mapped", "label") := NULL]


# remove orf number from conting name
blast_dt <- blast_dt %>%
  mutate(
    contig = sub("_\\d+$", "", qseqid), 
    orf_id = sub(".*_", "", qseqid)     
  )

#merge blast_dt with read_counts on contig and qseqid
blast_dt <- left_join(blast_dt, read_counts, by = c("contig"))

# Remove Insecta hits (class == Insecta)
#blast_dt <- blast_dt[!(class == "Insecta")]
#blast_dt <- blast_dt[!(kingdom == "unclassified")]

#blast_virus <- blast_dt[(kingdom == "Viruses")]

#blast_bacteria <- blast_dt[(kingdom == "Bacteria")]

# Pie chart of diversity




org <- "Eukaryota"


for (org in c ("Viruses", "Bacteria")){ 
  
blast_org <- blast_dt[(kingdom == org)]
  
# keep unique contig-genus combinations
unique_contig_genus <- blast_org %>%
  group_by(contig, genus) %>%
  summarise(mapped = first(mapped)) %>%
  ungroup()

unique_contig_genus <- as.data.table(unique_contig_genus)

unique_contig_genus[, label := gsub("(_NODE_.*)", "", contig)]

# summarise genus mapped reads per label
genus_summary <- unique_contig_genus %>%
    group_by(label, genus) %>%   # group by sample and genus
    summarise(mapped = sum(mapped)) %>%  # sum mapped reads for that genus in that sample
    ungroup()

fly_mapped_reads <- read.table("fly_mapped_read_counts.tsv", header=FALSE) 
colnames(fly_mapped_reads) <- c("label", "reads_counts_fly")

genus_summary <- genus_summary %>%
  left_join(fly_mapped_reads, by="label")

genus_summary <- genus_summary %>%
  group_by(label) %>%   # within each sample
  mutate(rel_abundance_fly = (mapped / reads_counts_fly) * 1e6) %>%   # compute relative abundance
  ungroup()

genus_summary <- genus_summary %>%
  group_by(label) %>%   # within each sample
  mutate(rel_abundance = (mapped / sum(mapped)) * 100) %>%   # compute relative abundance
  ungroup()

# Get top 10 genera in each sample and then take the union of those genera to get the top genera across all samples. 
top_genera <- genus_summary %>%
  group_by(label) %>%
  arrange(desc(rel_abundance_fly)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(genus) %>%
  distinct()

#top_genera <- genus_summary %>%
 #   group_by(label) %>%
 ##   summarise(total_mapped = sum(mapped), .groups = "drop") %>%
 # ungroup()

#fly_count <- genus_summary %>%
#    group_by(label) %>%
#    summarise(fly_mapped = sum(reads_counts_fly), .groups = "drop") %>%
#    ungroup()

#top_genera <- genus_summary %>%
#  group_by(genus) %>%
#  summarise(total_mapped = sum(rel_abundance_fly), .groups = "drop") %>%
#  arrange(desc(total_mapped)) %>%
#  slice_head(n = 15)


genus_summary <- genus_summary %>%
  mutate(genus_label = ifelse(genus %in% top_genera$genus, genus, "Other"))

#genus_summary <- genus_summary %>%
#  filter(genus_label != "Other")


# Assigne grey colour to "Other" and assign colours to the top 15 genera
set.seed(723451)
palette_colors <- c(
  "Other" = "grey",       # Assign grey to Other
  colorRampPalette(c("#845e9b", "#a3cee3", "#f7f191", "#FF496C", "orange",  "blue", "brown",
                     "#2b80b9", "#54af48", "#6c4141", "#008080", "#EF98AA",
                     "#FC8EAC", "purple", "#CD4A4A", "#1F75FE", "#FB7EFD", "#1CA9C9" ))(length(top_genera$genus))
)

# Assign names to the remaining colors
names(palette_colors) <- c("Other", top_genera$genus)

# Calculate the number of samples each genus is present in
genus_presence <- genus_summary %>%
  group_by(genus_label) %>%
  reframe(
    sample_count = ifelse(genus_label == "Other", 0, n_distinct(label))  # Set sample_count for "Other" to 0
  ) %>%
  arrange(desc(sample_count))  # Arrange by the number of samples (descending)
# Remove duplicates (if any)
genus_presence <- genus_presence %>%
  distinct(genus_label, .keep_all = TRUE)



# Reorder genus_label based on sample_count
genus_summary <- genus_summary %>%
  mutate(genus_label = factor(genus_label, levels = rev(genus_presence$genus_label)))

#genus_summary <- genus_summary %>%
#  mutate(log_rel_abundance_fly = log10(rel_abundance_fly))

genus_summary <- genus_summary %>%
  filter(!is.na(rel_abundance_fly) & rel_abundance_fly >= 1)

# Custom labeling function
label_k <- function(x) {
  ifelse(x >= 1000, paste0(x / 1000000, "M"), x)
}

genus_summary <- genus_summary %>%
  mutate(species = case_when(
    grepl("^CAME", label) ~ "Hcam",
    grepl("^CONF", label) ~ "Hconf",
    grepl("^DEFL", label) ~ "Sdef",
    TRUE ~ NA_character_
  ))

genus_summary <- genus_summary %>%
  mutate(treatment = ifelse(grepl("I$", label), "Infected", "Uninfected"))

# Stacked barplot showing the relative abundance of genus in each sample
bar_plot <- ggplot(genus_summary, aes(x = label, y = rel_abundance_fly, fill = genus_label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_y_continuous(labels = label_k) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Samples", y = "Reads per million host-mapped reads", fill = paste0("Genus (", org, ")")) +
  scale_fill_manual(values = palette_colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "right"
  ) 

bar_plot

box_plot <- ggplot(genus_summary, aes(x = reorder(genus_label, -rel_abundance_fly, FUN = median), y = rel_abundance_fly, fill = genus_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, shape = 1) +
  theme_classic() +
  scale_y_continuous(labels = label_k) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  labs(x = "Genus", y = "Reads per million host-mapped reads", fill = paste0("Genus (", org,")")) +
  scale_fill_manual(values = palette_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none")

box_plot

box_plot + theme(legend.position = "none")

#ggsave(paste0(org,"_genus_abundances_A.emf"), plot = bar_plot + theme(legend.position = "none"), width = 8.27, height = 11.69/2.2, units = "in")
#ggsave(paste0(org,"_genus_abundances_B.svg"), plot = box_plot, width = 8.27, height = 11.69/2.2, units = "in")
#ggsave(paste0(org,"_genus_abundances_legends.svg"), plot = get_legend(bar_plot), width = 8.27, height = 11.69/2, units = "in")

}

####
# Wilcoxon test
####

sp <- "Hcam"

for (sp in c("Hcam", "Hconf", "Sdef")) {
  # Filter data for the specific species
  genus_summary_sp <- genus_summary %>% filter(species == sp)

  # Create a new column for treatment (Infected or Uninfected)
  genus_summary_sp <- genus_summary_sp %>%
    mutate(treatment = ifelse(grepl("I$", label), "Infected", "Uninfected"))

  # Get unique genera for testing
  genera <- unique(genus_summary_sp$genus_label)

  # Run Wilcoxon test for each genus
  wilcox_results <- lapply(genera, function(gen) {
    test_data <- genus_summary_sp %>% filter(genus_label == gen)
    if(length(unique(test_data$treatment)) == 2) {
      test <- wilcox.test(rel_abundance_fly ~ treatment, data = test_data)
      data.frame(
        genus = gen,
        p_value = test$p.value
      )
    } else {
      data.frame(
        genus = gen,
        p_value = NA
      )
    }
  }) %>% bind_rows()

  # Adjust for multiple testing (FDR)
  wilcox_results <- wilcox_results %>%
    mutate(adj_p = p.adjust(p_value, method = "fdr"))

  # View results
  wilcox_results %>% arrange(adj_p)

  write.table(wilcox_results %>% arrange(adj_p), file = paste0("wilcox_results_", sp, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

# Get unique genera for testing
genera <- unique(genus_summary$genus_label)

# Run Wilcoxon test for each genus
wilcox_results <- lapply(genera, function(gen) {
  test_data <- genus_summary %>% filter(genus_label == gen)
  if(length(unique(test_data$treatment)) == 2) {
    test <- wilcox.test(rel_abundance ~ treatment, data = test_data)
    data.frame(
      genus = gen,
      p_value = test$p.value
    )
  } else {
    data.frame(
      genus = gen,
      p_value = NA
    )
  }
}) %>% bind_rows()


# Adjust for multiple testing (FDR)
wilcox_results <- wilcox_results %>%
  mutate(adj_p = p.adjust(p_value, method = "fdr"))

# View results
wilcox_results %>% arrange(adj_p)

write.table(wilcox_results %>% arrange(adj_p), file = "wilcox_results_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)

####
# Alpha-Diversity
####
#genus_summary <- genus_summary %>%
#  mutate(treatment = ifelse(grepl("I$", label), "Infected", "Uninfected"))

for (sp in c("Hcam", "Hconf", "Sdef")) {
  alpha_data_sp <- genus_summary %>% filter(species == sp) %>%
    select(label, genus, rel_abundance_fly, species, treatment) %>%
    pivot_wider(names_from = genus, values_from = rel_abundance_fly, values_fill = 0)

  sample_labels_sp <- alpha_data_sp$label
  alpha_counts_sp <- alpha_data_sp %>% select(-label, -species, -treatment)

  alpha_metrics_sp <- data.frame(
    sample = sample_labels_sp,
    shannon = diversity(alpha_counts_sp, index = "shannon"),
    simpson = diversity(alpha_counts_sp, index = "simpson")
  )

  alpha_metrics_sp <- alpha_metrics_sp %>%
    left_join(unique(alpha_data_sp[, c("label", "treatment")]), by = c("sample" = "label"))
  
  test <- wilcox.test(shannon ~ treatment, data = alpha_metrics_sp)
  wilcox_results <- data.frame(species = sp,
                               p_value = test$p.value)
  #print(wilcox.test(shannon ~ treatment, data = alpha_metrics_sp))
  
  fileConn <- file(paste0("alpha_metrics_", sp, ".txt"), open = "w")
  writeLines(capture.output(wilcox_results), fileConn)
  writeLines(capture.output(alpha_metrics_sp), fileConn)
  close(fileConn)

  #write.table(alpha_metrics_sp, file = paste0("alpha_metrics_", sp, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# x axis 
  alpha_plot_sp <- ggplot(alpha_metrics_sp, aes(x = treatment, y = shannon, 
                                                fill = treatment, group = treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, size = 0.8, width = 0.3) +
    geom_jitter(width = 0.05, size = 2, alpha = 0.3) +
    #geom_text_repel(size = 3, max.overlaps = 1, show.legend = FALSE) +
    scale_color_manual(values=c("#029fc9", "#c65120")) +
    theme_classic() +
    labs(title = "", y = "Alpha Diversity (Shannon)", x = "") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 12))

  #gsave(paste0("alpha_", sp, ".svg"), plot = alpha_plot_sp, width = 8.27/2, height = 11.69/3, units = "in")
alpha_plot_sp
# PCoA analysis for each species
  
  sample_metadata_sp <- alpha_data_sp %>%
    select(label, treatment)

  genus_matrix_sp <- alpha_data_sp %>%
    select(-label, -treatment, -species) %>%
    as.data.frame()

  # Bray-Curtis dissimilarity
  bray_dist_sp <- vegdist(genus_matrix_sp, method = "bray")

  pcoa_res_sp <- cmdscale(bray_dist_sp, k = 2, eig = TRUE)

  pcoa_df_sp <- data.frame(
    Sample = sample_metadata_sp$label,
    Treatment = sample_metadata_sp$treatment,
    PCoA1 = pcoa_res_sp$points[,1],
    PCoA2 = pcoa_res_sp$points[,2]
  )

  PCoA_plot_sp <- ggplot(pcoa_df_sp, aes(x = PCoA1, y = PCoA2, color = Treatment, label=Sample)) +
    geom_point(size=2) +
    scale_color_manual(values=c("#c65120", "#029fc9")) +
    geom_text_repel(size = 3, max.overlaps = 6, show.legend = FALSE) +
    #stat_ellipse(level = 0.95, linetype = "dotted") +
    theme_classic() +
    labs(
      title = paste("PCoA: ", sp),
      x = paste0("PCoA1 (", round(pcoa_res_sp$eig[1] / sum(pcoa_res_sp$eig) * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(pcoa_res_sp$eig[2] / sum(pcoa_res_sp$eig) * 100, 1), "%)")
    ) +
    theme(
      text = element_text(size = 12),
      legend.position = "none"
    )

  #ggsave(paste0("PCoA_", sp, ".svg"), plot = PCoA_plot_sp, width = 8.27/2, height = 11.69/3, units = "in")

  #comb_plot <- plot_grid(alpha_plot_sp, PCoA_plot_sp, ncol = 2))
  ggsave(paste0("Alpha_", sp, ".svg"), plot = alpha_plot_sp, width = 8.27/2, height = 11.69/3, units = "in")

  ggsave(paste0("PCoA_", sp, ".svg"), plot = PCoA_plot_sp, width = 8.27/2, height = 11.69/3, units = "in")

}


alpha_data <- genus_summary %>%
  select(label, genus, mapped) %>%
  pivot_wider(names_from = genus, values_from = mapped, values_fill = 0)

sample_labels <- alpha_data$label
alpha_counts <- alpha_data %>% select(-label)


# Calculate Shannon and Simpson diversity
alpha_metrics <- data.frame(
  sample = sample_labels,
  shannon = diversity(alpha_counts, index = "shannon"),
  simpson = diversity(alpha_counts, index = "simpson")
)

alpha_metrics <- alpha_metrics %>%
  left_join(unique(genus_summary[, c("label", "treatment")]), by = c("sample" = "label"))

head(alpha_metrics)
write.table(alpha_metrics, file = "alpha_metrics_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)


alpha_plot <- ggplot(alpha_metrics, aes(x = treatment, y = shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
  scale_color_manual(values=c("#029fc9", "#c65120")) +
  theme_classic() +
  labs(title = "Shannon diversity across treatments", y = "Shannon diversity", x = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

ggsave("alpha_all_samples.svg", plot = alpha_plot, width = 8.27/2, height = 11.69/3, units = "in")


#########
## PCoA analysis
#########

genus_abundance_wide <- genus_summary %>% 
  select(label, genus, rel_abundance, treatment) %>%
  pivot_wider(names_from = genus, values_from = rel_abundance, values_fill = 0)

sample_metadata <- genus_abundance_wide %>%
  select(label, treatment)

genus_matrix <- genus_abundance_wide %>%
  select(-label, -treatment) %>%
  as.data.frame()

# Bray-Curtis dissimilarity
bray_dist <- vegdist(genus_matrix, method = "bray")

pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = sample_metadata$label,
  Treatment = sample_metadata$treatment,
  PCoA1 = pcoa_res$points[,1],
  PCoA2 = pcoa_res$points[,2]
)

PCoA_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment, label=Sample)) +
  geom_point(size=2) +
  scale_color_manual(values=c("#c65120", "#029fc9")) +
  geom_text_repel(size = 3, max.overlaps = 6, show.legend = FALSE) +
  #stat_ellipse(level = 0.95, linetype = "dotted") +
  theme_classic() +
  labs(
    title = "Beta diversity (PCoA) based on Bray-Curtis dissimilarity",
    x = paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)")
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "none"
  )

ggsave("PCoA_all_samples.svg", plot = PCoA_plot, width = 8.27/2, height = 11.69/3, units = "in")

