# Script: GO-ORA-analysis.R
# Description:  Gene Ontology overrepresentation analysis, visualization of enrichment results
#               Example dataset: Chron's disease - rectum samples from https://ibdmdb.org/
# Version: 1.0
# Last updated: 2024-05-28
# Author: mkutmon

# #############################################
# SET UP WORKSPACE
# #############################################
# Let's install the required packages if needed
if (!("BiocManager" %in% installed.packages())) { install.packages("BiocManager") }
if (!("rstudioapi" %in% installed.packages())) { BiocManager::install("rstudioapi") }
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db") }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr") }
if (!("EnhancedVolcano" %in% installed.packages())) { BiocManager::install("EnhancedVolcano") }
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler") }
if (!("enrichplot" %in% installed.packages())) { BiocManager::install("enrichplot") }
if (!("Rgraphviz" %in% installed.packages())) { BiocManager::install("Rgraphviz") }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr") }

# Loading packages
library(rstudioapi)
library(org.Hs.eg.db)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(Rgraphviz)
library(readr)

# #############################################
# SETUP
# #############################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
out.folder <- "output/"
dir.create(out.folder)

# #############################################
# READ DATASET
# #############################################
options(scipen=999)
# importing CD gene expression dataset
data <- read.table("CD-dataset.tsv", header=TRUE, sep = "\t", dec = ".")

# #############################################
# DIFFERENTIALLY EXPRESSED GENES
# #############################################

# Let's get the differentially expressed genes
# TODO: discuss cutoffs!
logfc.cutoff <- 0.58
pvalue.cutoff <- 0.05
degs <- data[abs(data$log2FC) > logfc.cutoff & data$pvalue < pvalue.cutoff,]
write.table(degs, file=paste0(out.folder,"degs.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

# Let's create the volcano plots (ignore warning if folders already exist)
# you might need to run the code to generate the volcano plots alone if the output remains white
filename <- paste0(out.folder,"CD-volcano.png")
png(filename , width = 2000, height = 1500, res = 150)
EnhancedVolcano::EnhancedVolcano(data, title = paste0(nrow(degs), " DEGs"), lab = data$GeneSymbol, x = "log2FC", y = "pvalue", pCutoff = pvalue.cutoff, FCcutoff = logfc.cutoff, labSize = 3, xlim = c(-3.5,3.5), ylim=c(0,10))
dev.off()


# ==========================================
# PATHWAY ENRICHMENT ANALYSIS
# ==========================================
# overrepresentation analysis
# identification of over-represented GO terms in the differentially
# expressed genes 
# focus on biological processes

# this might take a little bit
res.go <- clusterProfiler::enrichGO(degs$ENTREZ.ID, OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, minGSSize = 20, maxGSSize = 400)
res.go.df <- as.data.frame(res.go)

# Q: Check the enrichment result - what are the most significant pathways?
# What information does the GeneRatio and BgRatio columns contain?

num <- nrow(as.data.frame(res.go))
print(paste0(num, " significant biological processes"))
  
write.table(as.data.frame(res.go), file=paste0(out.folder,"GO-ora.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ==========================================
# VISUALIZING THE ENRICHMENT RESULT
# ==========================================
res.go.sim <- enrichplot::pairwise_termsim(res.go)

# for the visualization there are too many enriched processes
# only the top 200 are included in the visualizations
ep <- emapplot(res.go.sim, showCategory = 200)
tp <- treeplot(res.go.sim, label_format = 0.5, cluster.params = list(n = 15), showCategory = 200)

filename <- paste0(out.folder,"GO_ORA_Emappplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(ep)
dev.off()
filename <- paste0(out.folder,"GO_ORA_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(tp)
dev.off()

# Look at the two figures in your output folder
# Q: How would you interpret them?
