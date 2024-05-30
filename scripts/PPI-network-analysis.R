# Script: PPI-network-analysis.R
# Description:  Create a PPI network with STRING and analyse the network topology
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
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr") }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr") }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3") }
if (!("randomcoloR" %in% installed.packages())) { BiocManager::install("randomcoloR") }

# Loading packages
library(rstudioapi)
library(dplyr)
library(readr)
library(RCy3)
library(randomcoloR)


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

# ==================================================================
# PPI network creation with the stringApp for Cytoscape
# ==================================================================
# make sure Cytoscape is running
# whenever the call connects to Cytoscape you see things changing live
RCy3::cytoscapePing()

# the first time install the stringApp
installApp('stringApp') 

# now let's query the database and retrieve all known interactions between the DEGs
# I am using the high confidence cutoff of 0.7 
query <- format_csv(as.data.frame(degs$ENTREZ.ID), col_names=F, escape = "double", eol =",")
commandsRun(paste0('string protein query cutoff=0.7 newNetName="PPI network" query="',query,'" limit=0'))

# Look at the network - what do you see? 
# Do you see any clusters/patterns?


# ==================================================================
# Let's analyse the network to find the most important nodes
# ==================================================================

# let's extract the largest connected compontent
RCy3::commandsRun('network select subnetwork createSubnetwork=TRUE')
RCy3::analyzeNetwork()
RCy3::createVisualStyle("log2FC vis")
RCy3::setNodeLabelMapping("display name", style.name = "log2FC vis")
RCy3::setNodeSizeMapping("Degree", table.column.values = c(1,90), sizes = c(30,150), mapping.type = "c", default.size = 30, style.name = "log2FC vis")
RCy3::setNodeColorDefault("#E5F5E0", style.name = "log2FC vis")
RCy3::setVisualStyle("log2FC vis")
RCy3::clearSelection()
RCy3::exportPNG(paste0(out.folder,"PPI-degree.png"), allGraphicsDetails = TRUE, zoom = 500)
# Look at the network. Do you see any hub nodes? Where are they located? 

# What happens if you change the node size mapping based on the betweenness? how would you interpret this result?
RCy3::setNodeSizeMapping("BetweennessCentrality", table.column.values = c(0,0.5), sizes = c(30,150), mapping.type = "c", default.size = 30, style.name = "log2FC vis")
RCy3::exportPNG(paste0(out.folder,"PPI-betweenness.png"),  allGraphicsDetails = TRUE, zoom = 500)


# What happens if you change the node size mapping based on the clustering coefficient? how would you interpret this result?
RCy3::setNodeSizeMapping("ClusteringCoefficient", table.column.values = c(0,1), sizes = c(30,150), mapping.type = "c", default.size = 30, style.name = "log2FC vis")
RCy3::exportPNG(paste0(out.folder,"PPI-clusteringcoefficient.png"),  allGraphicsDetails = TRUE, zoom = 500)


# ==========================================
# VISUALIZING THE EXPRESSION DATA
# ==========================================

loadTableData(data, data.key.column = "ENTREZ.ID", table.key.column = "query term")
control.points <- c (-2.0, 0.0, 2.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FC", control.points, colors, style.name = "log2FC vis")

RCy3::exportPNG(paste0(out.folder,"PPI-log2FC.png"), zoom = 500)
# Q: How would you interpret this figure? do they up-regulated / down-regulated genes cluster together?


# ==========================================
# PERFORM COMMUNITY CLUSTERING
# ==========================================

# let's perform some network clustering to find subnetworks to see what genes are in the network

RCy3::installApp("clustermaker2")
suid <- RCy3::getNetworkSuid()
# this algorithm has a random component so it might look slightly different if you rerun it again
RCy3::commandsRun(paste0("cluster glay network=",suid))

# now let's visualize the different clusters - did you expect these clusters?
clusters <- RCy3::getTableColumns(columns=c("SUID", "__glayCluster"))
RCy3::copyVisualStyle("log2FC vis", "clusters")
num.c <- length(unique(clusters$`__glayCluster`))
palette <- distinctColorPalette(num.c)
setNodeColorMapping("__glayCluster",table.column.values = c(1:num.c), colors = palette, mapping.type = "d", style.name = "clusters")
setVisualStyle("clusters")
RCy3::exportPNG(paste0(out.folder,"PPI-clustered.png"), zoom = 500)


# let's extract the cluster with number 1 (you can pick any that you like)
nodes.c1 <- clusters[clusters$`__glayCluster`==1,]
RCy3::selectNodes(nodes.c1$SUID)
RCy3::createSubnetwork(nodes.c1$SUID)
RCy3::setVisualStyle("log2FC vis")

RCy3::exportPNG(paste0(out.folder,"PPI-cluster1.png"), zoom = 500)

RCy3::commandsRun(paste0('string retrieve enrichment allNetSpecies="Homo sapiens"'))
RCy3::commandsRun('string show enrichment')

# check the enrichment result for this subnetwork 
# are there any processes you also found in the pathway enrichment analysis?

# repeat with other clusters if you have time


# ==========================================
# SAVE CYTOSCAPE SESSION
# ==========================================

RCy3::saveSession(filename=paste0(out.folder,"PPI.cys"))
