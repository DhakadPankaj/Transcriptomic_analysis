library(readxl)
library(ggplot2)


setwd("Z:/Fly_infections")

# Read excel workbook
sp <- "Sdef"

# empty data frame to store the top genes
top_genes_dt <- data.frame()

for (sp in c("Hcam", "Hconf", "Sdef")) {
  # Read the excel file
  DE_dt <- read_excel("DE_Dmel_association.xlsx", sheet = sp)
  
  DE_dt <- as.data.frame(DE_dt)
  DE_dt <- DE_dt[DE_dt$significance == "Upregulated", ]
  
  # If Name is NA then replace it with Gene
  DE_dt[DE_dt$Name=="NA",]$Name <- DE_dt[DE_dt$Name=="NA",]$Gene
  
  top1 <- 30
  top2 <- 5
  if(sp == "Sdef") {
     top1 <- 10
     top2 <- 5
  }
  # Select top 10 genes based on log2FoldChange
  top_genes <- head(DE_dt[order(-DE_dt$log2FoldChange), ], top1)
  
  # select top 10 genes based on adjusted p-value
  top_genes_adj <- head(DE_dt[order(DE_dt$padj), ], top2)
  
  # get the union of the two top genes
  top_genes <- unique(rbind(top_genes, top_genes_adj))
  
  top_genes$Name <- factor(top_genes$Name, levels = top_genes$Name[order(top_genes$log2FoldChange, decreasing = FALSE)])
  
  # Add species column
     top_genes$species <- sp
 # Add to the top_genes_dt data frame
  top_genes_dt <- rbind(top_genes_dt, top_genes[,c("Gene", "Name", "species", "log2FoldChange", "padj")])
}

DE_dt <- read_excel("DE_Dmel_association.xlsx", sheet = sp)


DE_dt <- as.data.frame(DE_dt)
DE_dt <- DE_dt[DE_dt$significance == "Upregulated", ]

# If Name is NA then replace it with Gene
DE_dt[DE_dt$Name=="NA",]$Name <- DE_dt[DE_dt$Name=="NA",]$Gene

# Select top 10 genes based on log2FoldChange
top_genes <- head(DE_dt[order(-DE_dt$log2FoldChange), ], 10)

# select top 10 genes based on adjusted p-value
top_genes_adj <- head(DE_dt[order(DE_dt$padj), ], 5)

# get the union of the two top genes
top_genes <- unique(rbind(top_genes, top_genes_adj))

top_genes$Name <- factor(top_genes$Name, levels = top_genes$Name[order(top_genes$log2FoldChange, decreasing = FALSE)])


top_genes <- top_genes %>%
  mutate(xmin = 0.8,  # center on x = 1, make tile narrower
         xmax = 0.85,
         ymin = as.numeric(factor(Name)) - 0.2,
         ymax = as.numeric(factor(Name)) + 0.2)

# save plot as plot_sp 
assign(
  paste0("plot_", sp), ggplot(top_genes) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = log2FoldChange), color = "grey", linewidth = 0.001) +
  scale_fill_gradient2(low = "white", high = "red") +
  theme_minimal() +
  labs(title = paste("Top DE genes for", sp),
       y = "", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "left"
  ) +
  scale_y_continuous(
    breaks = as.numeric(factor(top_genes$Name)),
    labels = top_genes$Name,
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(0.8, 1.2)) 
)

plot_Hcam
plot_Hconf
plot_Sdef



###

top_genes <- top_genes_dt

gene_counts <- top_genes %>%
  distinct(Name, species) %>%
  count(Name, name = "n_species")

# Create a full set of combinations
full_grid <- expand.grid(
  Name = unique(top_genes$Name),
  species = unique(top_genes$species)
)

# Join with original data and gene_counts
top_genes_complete <- full_grid %>%
  left_join(top_genes, by = c("Name", "species")) %>%
  left_join(gene_counts, by = "Name")

top_genes_complete <- top_genes_complete %>%
  arrange(desc(n_species), Name) %>%
  mutate(Name = factor(Name, levels = unique(Name)))

# Use ymin/ymax for geom_rect; limit y-scale to max_genes
plot <- ggplot(top_genes_complete, aes(x = species, y = Name, fill = log2FoldChange)) +
  geom_tile(color = "grey", width = 0.15, height = 1.2) +
  scale_fill_gradient2(
    low = "white", high = "red",
    na.value = "white",  # Blank tiles for missing values
    midpoint = 0
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_flip()+
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size=8, angle = 90, vjust = 0.3),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )
plot
ggsave("Top_DE_genes.svg", plot = plot, width = 8.27, height = 11.7/3, units = "in")
