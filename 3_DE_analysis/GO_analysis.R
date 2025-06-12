library(clusterProfiler)
library(enrichplot)
library(data.table)
library(topGO)
library("Rgraphviz")
library(readxl)
library(dplyr)
library(ggplot2)

setwd("Z:/Fly_infections")

for(sp in c("Hconf", "Hcam", "Sdef")){
#sp<-"Sdef"
eggnog_data <- read_excel(paste0("eggNog/", sp, ".emapper.annotations.xlsx"))

gene2go <- eggnog_data[, c("query", "GOs")]
gene2go <- gene2go %>% filter(GOs != "-")

write.table(gene2go, file=paste0("eggNog/", sp, "_GO_universe.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = FALSE)

gene2go <- readMappings(file = paste0("eggNog/", sp, "_GO_universe.tsv"), sep = "\t")

geneUniverse <- names(gene2go)

DE_genes <- read_excel("DE_genes.xlsx", sheet=sp)
DE_up <- DE_genes %>% filter(significance == "Upregulated")
DE_down <- DE_genes %>% filter(significance == "Downregulated")

genesOfInterest_up <- as.character(DE_up$Transcript)
genesOfInterest_down <- as.character(DE_down$Transcript)

geneList_up <- factor(as.integer(geneUniverse %in% genesOfInterest_up))
names(geneList_up) <- geneUniverse

geneList_down <- factor(as.integer(geneUniverse %in% genesOfInterest_down))
names(geneList_down) <- geneUniverse

for (category in c("BP", "MF", "CC")){
#category<-"BP"
myGOdata_up <- new("topGOdata", description="GO enrichment UP", ontology= category, allGenes=geneList_up,  annot = annFUN.gene2GO, gene2GO = gene2go)
myGOdata_down <- new("topGOdata", description="GO enrichment DOWN", ontology= category, allGenes=geneList_down,  annot = annFUN.gene2GO, gene2GO = gene2go)

#myGOdata

resultFisher_up <- runTest(myGOdata_up, algorithm="weight01", statistic="fisher")
resultFisher_down <- runTest(myGOdata_down, algorithm="weight01", statistic="fisher")


# see how many results we get where weight01 gives a P-value <= 0.05:
mysummary_up <- summary(attributes(resultFisher_up)$score <= 0.05)
mysummary_up

numsignif_up <- as.integer(mysummary_up[[3]]) # how many terms is it true that P <= 0.001
#numsignif_up
allRes_up <- GenTable(myGOdata_up, topgoFisher = resultFisher_up, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_up)
#allRes_up
write.table(allRes_up, file=paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = TRUE)

mysummary_down <- summary(attributes(resultFisher_down)$score <= 0.05)
numsignif_down <- as.integer(mysummary_down[[3]]) # how many terms is it true that P <= 0.001
allRes_down <- GenTable(myGOdata_down, topgoFisher = resultFisher_down, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_down)
write.table(allRes_down, file=paste0("GO_enrichment/", sp,"_", category , "_DOWN_GOenrich.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = TRUE)


myterms_up <- allRes_up$GO.ID
mygenes_up <- genesInTerm(myGOdata_up, myterms_up)

myterms_down <- allRes_down$GO.ID
mygenes_down <- genesInTerm(myGOdata_down, myterms_down)

fileConn1 <- file(paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich_genes.tsv"), open = "w")
for (i in 1:length(myterms_up)){
  myterm <- myterms_up[i]
  mygenesforterm <- mygenes_up[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest_up # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  # write each line to the file
  writeLines(paste(myterm, mygenesforterm2, sep = "\t"), fileConn1)
}
close(fileConn1)

# Combine the allRes_up and file with genes
mygenes_up <- read.table(paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich_genes.tsv"), header = FALSE, sep = "\t")
colnames(mygenes_up) <- c("GO.ID", "genes")
allRes_up$GeneRatio <- allRes_up$Significant / allRes_up$Annotated
allRes_up <- merge(allRes_up, mygenes_up, by = "GO.ID")
#sort the table by the topgoFisher column
allRes_up$topgoFisher <- as.numeric(as.character(allRes_up$topgoFisher))
allRes_up <- allRes_up[order(allRes_up$topgoFisher),]
write.table(allRes_up, file=paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = TRUE)

# delete the file with genes
file.remove(paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich_genes.tsv"))


fileConn2 <- file(paste0("GO_enrichment/", sp, "_", category , "_DOWN_GOenrich_genes.tsv"), open = "w")

for (i in 1:length(myterms_down)){
  myterm <- myterms_down[i]
  mygenesforterm <- mygenes_down[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest_down # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  # write each line to the file
  writeLines(paste(myterm, mygenesforterm2, sep = "\t"), fileConn2)
}
close(fileConn2)

# Combine the allRes_down and file with genes
mygenes_down <- read.table(paste0("GO_enrichment/", sp, "_", category , "_DOWN_GOenrich_genes.tsv"), header = FALSE, sep = "\t")
colnames(mygenes_down) <- c("GO.ID", "genes")
allRes_down$GeneRatio <- allRes_down$Significant / allRes_down$Annotated
allRes_down <- merge(allRes_down, mygenes_down, by = "GO.ID")
#sort the table by the topgoFisher column (NUMERIC)   some values are e-10, e-20, etc
allRes_down$topgoFisher <- as.numeric(as.character(allRes_down$topgoFisher))
allRes_down <- allRes_down[order(allRes_down$topgoFisher),]
write.table(allRes_down, file=paste0("GO_enrichment/", sp, "_", category , "_DOWN_GOenrich.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = TRUE)

# delete the file with genes
file.remove(paste0("GO_enrichment/", sp, "_", category , "_DOWN_GOenrich_genes.tsv"))

 }
 }

library(clusterProfiler)
library(enrichplot)
library(data.table)
library(topGO)
library(dplyr)
library(readxl)
library(ggplot2)

for(sp in c("Hconf", "Hcam", "Sdef")) {
  # Read the Excel file
  eggnog_data <- read_excel(paste0("eggNog/", sp, ".emapper.annotations.xlsx"))
  
  # Extract the gene-to-GO mappings
  gene2go <- eggnog_data[, c("query", "GOs")]
  gene2go <- gene2go %>% filter(GOs != "-")
  
  # Write the gene2go data to a file
  write.table(gene2go, file=paste0("eggNog/", sp, "_GO_universe.tsv"), quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)
  
  # Read the gene2go mappings
  gene2go <- readMappings(file = paste0("eggNog/", sp, "_GO_universe.tsv"), sep = "\t")
  
  # Create the gene universe
  geneUniverse <- names(gene2go)
  
  # Load the DE genes
  DE_genes <- read_excel("DE_genes.xlsx", sheet=sp)
  DE_up <- DE_genes %>% filter(significance == "Upregulated")
  DE_down <- DE_genes %>% filter(significance == "Downregulated")
  
  genesOfInterest_up <- as.character(DE_up$Transcript)
  genesOfInterest_down <- as.character(DE_down$Transcript)
  
  geneList_up <- factor(as.integer(geneUniverse %in% genesOfInterest_up))
  names(geneList_up) <- geneUniverse
  
  geneList_down <- factor(as.integer(geneUniverse %in% genesOfInterest_down))
  names(geneList_down) <- geneUniverse
  
  for (category in c("BP", "MF","CC")) {
    myGOdata_up <- new("topGOdata", description="GO enrichment UP", ontology= category, allGenes=geneList_up,  annot = annFUN.gene2GO, gene2GO = gene2go)
    myGOdata_down <- new("topGOdata", description="GO enrichment DOWN", ontology= category, allGenes=geneList_down,  annot = annFUN.gene2GO, gene2GO = gene2go)
    
    resultFisher_up <- runTest(myGOdata_up, algorithm="weight01", statistic="fisher")
    resultFisher_down <- runTest(myGOdata_down, algorithm="weight01", statistic="fisher")
    
    mysummary_up <- summary(attributes(resultFisher_up)$score <= 0.05)
    numsignif_up <- as.integer(mysummary_up[[3]])
    allRes_up <- GenTable(myGOdata_up, topgoFisher = resultFisher_up, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_up)
    write.table(allRes_up, file=paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = TRUE)
    
    mysummary_down <- summary(attributes(resultFisher_down)$score <= 0.05)
    numsignif_down <- as.integer(mysummary_down[[3]])
    allRes_down <- GenTable(myGOdata_down, topgoFisher = resultFisher_down, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif_down)
    write.table(allRes_down, file=paste0("GO_enrichment/", sp,"_", category , "_DOWN_GOenrich.tsv"), quote = FALSE, row.names = F, sep = "\t", col.names = TRUE)
    
    # Open a file to write the results to
    fileConn1 <- file(paste0("GO_enrichment/", sp, "_", category , "_UP_GOenrich_genes.tsv"), open = "w")
    myterms_up <- allRes_up$GO.ID
    mygenes_up <- genesInTerm(myGOdata_up, myterms_up)
    for (i in 1:length(myterms_up)) {
      myterm <- myterms_up[i]
      mygenesforterm <- mygenes_up[[myterm]]
      myfactor <- mygenesforterm %in% genesOfInterest_up
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
      mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
      writeLines(paste(myterm, mygenesforterm2, sep = "\t"), fileConn1)
    }
    close(fileConn1)
    
    fileConn2 <- file(paste0("GO_enrichment/", sp, "_", category , "_DOWN_GOenrich_genes.tsv"), open = "w")
    myterms_down <- allRes_down$GO.ID
    mygenes_down <- genesInTerm(myGOdata_down, myterms_down)
    for (i in 1:length(myterms_down)) {
      myterm <- myterms_down[i]
      mygenesforterm <- mygenes_down[[myterm]]
      myfactor <- mygenesforterm %in% genesOfInterest_down
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
      mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
      writeLines(paste(myterm, mygenesforterm2, sep = "\t"), fileConn2)
    }
    close(fileConn2)
  }
}

# Plot a dotplot for BP terms
category <- "BP"
bp_files <- list.files("GO_enrichment/", pattern = paste0(category, "_UP_GOenrich.tsv"), full.names = TRUE)
bp_files

go_data <- lapply(bp_files, function(file) {
  cat("Reading file:", file, "\n") 
  data <- read.table(file, header = TRUE, sep = "\t", fill = TRUE, quote = "", comment.char = "")
  data$Species <- gsub("_.*", "", basename(file)) # Extract species name
  return(data)
}) %>% bind_rows()


top_terms <- go_data %>%
  group_by(Species) %>%
  arrange(topgoFisher) %>%
  slice_head(n = 25) %>%
  ungroup()

top_terms <- top_terms %>%
  mutate(GO_ID_Term = paste0(GO.ID, ": ", Term))

top_terms$GO_ID_Term <- reorder(top_terms$GO_ID_Term, -log10(top_terms$topgoFisher))


dot_plot <- ggplot(top_terms, aes(x = GeneRatio, y = GO_ID_Term)) +
  geom_point(aes(size = Significant, color = -log10(topgoFisher))) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)) + 
  labs(x = "Gene Ratio", y = "GO Term (BP)", 
       color = "-log10(p-value)", size = "Count") +
  theme(axis.text.y = element_text(size = 12) ,
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1) ,
        axis.title.y = element_blank(),
        legend.position = "none") +
  facet_wrap(~ Species) # if you have multiple species

dot_plot

ggsave(paste0("GO_enrichment/", category, "_dotplot.svg"), plot = dot_plot, width = 8.27, height = 11.7/1.2, units = "in")
ggsave(paste0("GO_enrichment/", category, "_dotplot_legends.svg"), plot = get_legend(dot_plot), width = 8.27/2, height = 11.7/2, units = "in")


####################