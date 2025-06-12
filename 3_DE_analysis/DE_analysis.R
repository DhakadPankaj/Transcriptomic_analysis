## RNA seq analysis with DESeq2

library(DESeq2)
library(ggplot2)
library(draw)
library(cowplot)
library(pheatmap)
library(ggrepel)
setwd("Z:/Fly_infections")

# 
species_list <- c("Hcam", "Hconf")

for (species in species_list) {
  cat("Processing species:", species, "\n")
  filename <- paste(species, "_combined_counts.txt", sep="")
  dat <- read.table(filename, header=TRUE, row.names=1)
  dat <- dat[, order(grepl("FN$", colnames(dat)), decreasing = TRUE)]
  # Convert to matrix
  dat <- as.matrix(dat)
  head(dat)
  
  # Read and filter the condition data for the current species
  condition <- read.table("Fly_infection_Prettgeri_WithSubmissionNames.txt", header=TRUE, row.names="label_submitted")
  condition <- condition[condition$Label == species,]
  
  # Ensure the row names of condition match the column names of dat
  condition <- condition[match(colnames(dat), rownames(condition)), ]
  
  coldata <- data.frame(row.names=row.names(condition), treatment=as.factor(condition$Treatment))
  coldata$treatment <- relevel(coldata$treatment, ref="naive")
  # DESEq: normalizing to model fitting
  dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design=~ treatment)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)
  
  plotDispEsts(dds, main=paste("Dispersion plot for", species))
  
  # Regularized log transformation for clustering/heatmaps
  rld <- rlogTransformation(dds)
  
  custom_theme <- theme(plot.title = element_text(size=20, hjust = 0.5, vjust = -2),
                        plot.margin = margin(30, 30, 30, 30),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=12),
                        panel.background = element_rect(fill='transparent'),
                        plot.background = element_rect(fill='transparent', color=NA),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.position = "none",
                        axis.line = element_line(colour = "black" )) 
  
  # PCA plot with customized theme and labels
  pcaData <- plotPCA(rld, intgroup="treatment", returnData=TRUE, ntop=500)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, label=name)) +
    geom_point(size=2) +
    scale_color_manual(values=c("#029fc9", "#c65120")) +
    geom_text_repel(size = 3, max.overlaps = 6, show.legend = FALSE) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    custom_theme 
  print(pca)
  
  #ggsave("pca_legends.svg", plot = pca, width = 8.27, height = 8.27, units = "in")
  
  # Get differential expression results
  #res <- results(dds)
  res <- results(dds, name="treatment_infected_vs_naive")
  table(res$padj<0.05)
  
  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  # Merge with normalized count data
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  
  # Write results
  write.table(resdata, file=paste0(species, "_diffexpr_results.tsv"), quote = FALSE, row.names = F, sep = "\t")
  
  # Filter out genes with very low baseMean
  resdata <- resdata[resdata$baseMean >= 10, ]
  resdata <- na.omit(resdata)
  # Add a column to indicate significance
  resdata$significance <- "Not Significant"
  resdata$significance[resdata$padj < 0.05 & resdata$log2FoldChange >= 1] <- "Upregulated"
  resdata$significance[resdata$padj < 0.05 & resdata$log2FoldChange <= -1] <- "Downregulated"
  
  # Volcano plot
  volcano <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
    geom_point(size = 2, alpha=0.5) +
    scale_color_manual(values=c("blue", "grey", "red")) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
    custom_theme +
    labs(x="Log2 Fold Change", y="-Log10 (Adj. P-value)") 
  
  print(volcano)
  #ggsave("Volcano_legends.svg", plot = volcano, width = 8.27, height = 8.27, units = "in")
  
  # Clustered heatmap
  sig <- resdata[(abs(resdata$log2FoldChange) > 1 & resdata$padj < 0.05), ]
  write.table(sig, file=paste0(species, "_SDE_genes.tsv"), quote = FALSE, row.names = F, sep = "\t")
  
  
  sig <- sig[order(sig$padj), ]
  
  mat <- assay(rld)[sig$Gene, ]
  mat <- mat - rowMeans(mat)
  heatmap <- pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize=10, border_color=NA, color = colorRampPalette(c("blue", "white", "red"))(100))

  library(draw)
  heatmap_grob <- grid::grid.grabExpr(grid::grid.draw(heatmap$gtable))

 #combined_plot <- plot_grid(pca, volcano, heatmap_grob, labels = c("A", "B", "C"), ncol = 2, nrow = 2)

 #savefilename <- paste0(species,"_PCA_volcano_heatmap.svg")
  ggsave(paste0(species, "_PCA.svg"), plot = pca, width = 8.27/2, height = 11.7/3, units = "in")
  ggsave(paste0(species, "_Volcano.svg"), plot = volcano, width = 8.27/2, height = 11.7/3, units = "in")
  ggsave(paste0(species, "_Heatmap.svg"), plot = heatmap_grob, width = 8.27/2, height = 11.7/3, units = "in")
 
}

dat <- read.table("Sdef_combined_counts.txt", header=TRUE, row.names=1)
dat <- dat[, order(grepl("FN$", colnames(dat)), decreasing = TRUE)]
# Convert to matrix
dat <- as.matrix(dat)
head(dat)

# Create a data frame with the sample information
condition <- read.table("Fly_infection_Prettgeri_WithSubmissionNames.txt", header=TRUE, row.names="label_submitted")

# Filter for sp
condition <- condition[condition$Label == "Sdef",]

# Ensure the row names of condition match the column names of dat
condition <- condition[match(colnames(dat), rownames(condition)), ]

# Filter out male samples
condition <- condition[condition$Sex == "F", ]
dat <- dat[, colnames(dat) %in% rownames(condition)]
condition

coldata <- data.frame(row.names=row.names(condition), treatment=as.factor(condition$Treatment), sex = as.factor(condition$Sex))

coldata$treatment <- relevel(coldata$treatment, ref="naive")

coldata
# DESEq: normalizing to model fitting
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design =~ treatment)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

plotDispEsts(dds, main="Dispersion plot")

# Regularized log transformation for clustering/heatmaps

rld <- rlogTransformation(dds)
head(assay(rld))

custom_theme <- theme(plot.title = element_text(size=20, hjust = 0.5, vjust = -2),
                      plot.margin = margin(30, 30, 30, 30),
                      axis.title = element_text(size=16),
                      axis.text = element_text(size=12),
                      panel.background = element_rect(fill='transparent'),
                      plot.background = element_rect(fill='transparent', color=NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position = "none",
                      axis.line = element_line(colour = "black" )) 

# PCA plot with customized theme and labels
pcaData <- plotPCA(rld, intgroup=c("treatment", "sex"), returnData=TRUE, ntop=500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=sex, label=name)) +
    geom_point(size=2) +
    scale_color_manual(values=c("#029fc9", "#c65120")) +
    geom_text_repel(size = 3, max.overlaps = 6, show.legend = FALSE) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    custom_theme
pca


# Get differential expression results
#res <- results(dds)
res <- results(dds, name="treatment_infected_vs_naive")
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.table(resdata, file="Sdef_Fonly_diffexpr_results.tsv",quote = FALSE,row.names = F, sep = "\t")

# Filter out genes with very low baseMean
resdata <- resdata[resdata$baseMean >= 10, ]
resdata <- na.omit(resdata)
# Add a column to indicate significance
resdata$significance <- "Not Significant"
resdata$significance[resdata$padj < 0.05 & resdata$log2FoldChange >= 1] <- "Upregulated"
resdata$significance[resdata$padj < 0.05 & resdata$log2FoldChange <= -1] <- "Downregulated"

# Volcano plot
volcano <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
  geom_point(size = 2, alpha=0.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  custom_theme +
  labs( x="Log2 Fold Change", y="-Log10 (Adj. P-value)") 

print(volcano)

# Clustered heatmap

# extract significant genes (FC > 2 or FC < 0.5 and p-adj < 0.05)
sig <- resdata[(abs(resdata$log2FoldChange) > 1 & resdata$padj < 0.05), ]

write.table(sig, file="Sdef_Fonly_SDE_genes.tsv", quote = FALSE, row.names = F, sep = "\t")

sig <- sig[order(sig$padj), ]

mat <- assay(rld)[sig$Gene, ]
mat <- mat - rowMeans(mat)
heatmap <- pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize=10, border_color=NA, color = colorRampPalette(c("blue", "white", "red"))(100))

heatmap_grob <- grid::grid.grabExpr(grid::grid.draw(heatmap$gtable))

#combined_plot <- plot_grid(pca, volcano, heatmap_grob, labels = c("A", "B", "C"), ncol = 2, nrow = 2)

#savefilename <- paste0("Sdef_Fonly","_PCA_volcano_heatmap.svg")
#ggsave(savefilename, plot = combined_plot, width = 8.27, height = 8.27, units = "in")
ggsave(paste0("Sdef_PCA.svg"), plot = pca, width = 8.27/2, height = 11.7/3, units = "in")
ggsave(paste0("Sdef_Volcano.svg"), plot = volcano, width = 8.27/2, height = 11.7/3, units = "in")
ggsave(paste0("Sdef_Heatmap.svg"), plot = heatmap_grob, width = 8.27/2, height = 11.7/3, units = "in")
