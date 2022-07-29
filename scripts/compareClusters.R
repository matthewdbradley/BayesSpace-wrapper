#!/usr/bin/env Rscript
library(dplyr)
library(BayesSpace)
library(ggplot2)
library(Seurat)
library(patchwork)
args <- commandArgs(trailingOnly=TRUE)
s <- readRDS(paste0("outs/dataObjects/",args[1], "_seurat.Rds"))
# Parse arg string
clusters <- args[2] %>% strsplit(split=";") %>% unlist() %>% strsplit(split=":")
# Create Labeling function
label_spot <- Vectorize(function(cluster) {
  if (cluster == clusters[[1]][1]) {
    clusters[[1]][2]
  } else if (cluster == clusters[[2]][1]) {
    clusters[[2]][2]
  } else {
    "Other"
  }
})
DE.labels <- label_spot(s$spatial.cluster) %>% factor()
print("Created DE Labels")
# Add the region labels to each cell, and then subset out the other cells
s <- Seurat::AddMetaData(s,DE.labels, col.name="DE.label") %>% Seurat::SetIdent(value="DE.label") %>% subset(DE.label!="Other")
# Scale up data
s@assays$Spatial@scale.data <- s@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
# Find markers for scaled cluster data
print("Finding markers...")
top_markers <- Seurat::FindAllMarkers(s, assay='Spatial', slot='data',group.by='DE.label', only.pos=TRUE)
# Calculate percent difference between the clusters
top_markers$pct.diff <- top_markers$pct.1 - top_markers$pct.2

print("Writing out top markers...")
write.csv(top_markers, paste0("outs/markers/", args[1], "_",clusters[[1]][2], "_vs_",clusters[[2]][2], ".csv"))

# Find top 5 genes from each cluster based on avg_log2FC and pct.diff respectively
top5_markers_log2FC<-top_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top5_markers_pct.diff <- top_markers %>% group_by(cluster) %>% top_n(5,pct.diff)

# Print out heatmaps
print("Creating heatmaps...")
hm.log2FC <- Seurat::DoHeatmap(s, features = top5_markers_log2FC$gene, slot='scale.data', group.by = "DE.label", angle=0, size=4, raster=FALSE, label=FALSE)
hm.pct.diff <- Seurat::DoHeatmap(s, features = top5_markers_pct.diff$gene, slot='scale.data', group.by = "DE.label", angle=0, size=4, raster=FALSE, label=FALSE)
ggsave(hm.log2FC, file=paste0("outs/markers/", args[1], "_",clusters[[1]][2], "_vs_",clusters[[2]][2], "_log2FCHeatMap.pdf"), dpi=300, width=7, height=7,units="in")
ggsave(hm.pct.diff, file=paste0("outs/markers/", args[1], "_",clusters[[1]][2], "_vs_",clusters[[2]][2], "_pctdiffHeatMap.pdf"), dpi=300, width=7, height=7,units="in")
