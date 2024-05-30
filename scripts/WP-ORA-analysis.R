# Script: WP-ORA-analysis.R
# Description:  WikiPathways overrepresentation analysis, visualization of enrichment results
#               Example dataset: Chron's disease - rectum samples from https://ibdmdb.org/
# Version: 1.0
# Last updated: 2024-05-28
# Author: mkutmon

# #############################################
# SET UP WORKSPACE
# #############################################
# Let's install the required packages if needed
# it will only install the packages if they are not yet installed
if (!("BiocManager" %in% installed.packages())) { install.packages("BiocManager") }
if (!("rstudioapi" %in% installed.packages())) { BiocManager::install("rstudioapi") }
# if you get a non-zero exit status for org.Hs.eg.db - run the following line twice
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db") }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr") }
if (!("EnhancedVolcano" %in% installed.packages())) { BiocManager::install("EnhancedVolcano") }
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler") }
if (!("enrichplot" %in% installed.packages())) { BiocManager::install("enrichplot") }
if (!("Rgraphviz" %in% installed.packages())) { BiocManager::install("Rgraphviz") }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr") }
if (!("RColorBrewer" %in% installed.packages())) { BiocManager::install("RColorBrewer") }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3") }

Sys.setenv(LANG = "en")

# Loading packages
library(rstudioapi)
library(org.Hs.eg.db)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(Rgraphviz)
library(readr)
library(RColorBrewer)
library(RCy3)


# #############################################
# SETUP
# #############################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
out.folder <- "output/" #it created a folder named "output" in the directory we are working in
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
logfc.cutoff <- 0.58 #2^0.58 = 1.50 and 2^-0.58 = 0.67
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
# identification of sign. altered pathways from WikiPathways 

res.wp <- clusterProfiler::enrichWP(degs$ENTREZ.ID, organism = "Homo sapiens", pAdjustMethod = "fdr", pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 400)
res.wp.df <- as.data.frame(res.wp)

# Q: Check the enrichment result - what are the most significant pathways?
# What information does the GeneRatio and BgRatio columns contain?

num <- nrow(as.data.frame(res.wp))
print(paste0(num, " significant pathways"))
  
write.table(as.data.frame(res.wp), file=paste0(out.folder,"wikipathways-ora.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ==========================================
# VISUALIZING THE ENRICHMENT RESULT
# ==========================================
res.wp.sim <- enrichplot::pairwise_termsim(res.wp)
  
ep <- emapplot(res.wp.sim, showCategory = num)
tp <- treeplot(res.wp.sim, label_format = 0.5, cluster.params = list(n = 15), showCategory = num)

filename <- paste0(out.folder,"WikiPathways_ORA_Emappplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(ep)
dev.off()
filename <- paste0(out.folder,"WikiPathways_ORA_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(tp)
dev.off()

# Look at the two figures in your output folder
# Q: How would you interpret them?

# ==========================================
# VISUALIZING PATWHAY DIAGRAMS
# ==========================================

# Make sure you have Cytoscape started
RCy3::cytoscapePing()

# Install the WikiPathways app the first time - afterwards add a # in front
RCy3::installApp("wikipathways")

# You can select any of the significant pathways to visualize the data on the diagrams
# Some are quite big but you can explore them in Cytoscape
# You simply need to change the WP identifier in the next line from WP545 to any other significant pathway identifier
RCy3::commandsRun('wikipathways import-as-pathway id=WP545')

toggleGraphicsDetails()
loadTableData(data, data.key.column = "Ensembl.ID", table.key.column = "Ensembl")
control.points <- c (-1.0, 0.0, 1.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FC", control.points, colors, style.name = "WikiPathways")

RCy3::commandsRun('wikipathways import-as-pathway id=WP5115')

# Q: How can you read this pathway diagram? 
# Is the pathway more up- or down-regulated in Chron's disease patients?