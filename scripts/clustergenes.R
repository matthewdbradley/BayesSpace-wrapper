library(dplyr)
library(BayesSpace)
library(ggplot2)
library(Seurat)
library(patchwork)

## Read in data
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
s <- readRDS(paste0("outs/dataObjects/",id, ".Rds"))
s.en <- readRDS(paste0("outs/dataObjects/",id,"_enhanced.Rds"))

## ID highly variable genes, remove ribosomal genes
#dec <- scran::modelGeneVar(s)
#top <- scran::getTopHVGs(dec, n = )
top <- rownames(s)
hvgs <- top[grep("^RP[LS]", top, invert=TRUE)]
hvgs <- hvgs[grep("MT-", hvgs, invert=TRUE)]

## Create subplot identifiers with novel gene expression rows
s.en <- enhanceFeatures(s.en, s, model="xgboost", feature_names=hvgs,nrounds=0)
saveRDS(s.en, file=paste0("outs/dataObjects/", id, "_enhancedFeatures.Rds"))

## Convert SCE to Seurat object and use BayesSpace cluster as identifier
sobj <- Seurat::CreateSeuratObject(counts=logcounts(s.en),
                                   assay='Spatial',
                                   meta.data=as.data.frame(colData(s.en)))
sobj <- Seurat::SetIdent(sobj, value = "spatial.cluster")

## Scale data
sobj@assays$Spatial@scale.data <-
  sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
sobj@assays$Spatial@var.features <- hvgs

saveRDS(sobj, paste0("outs/dataObjects/", id, "_seurat.Rds"))
## Select top n markers from each cluster (by log fold change)
top_markers <- Seurat::FindAllMarkers(sobj, assay='Spatial', slot='data',
                                      group.by='spatial.cluster',
                                      only.pos=TRUE) %>% 
  group_by(cluster) %>% mutate(pct.diff=pct.1-pct.2) %>% arrange(cluster)

saveRDS(top_markers, file=paste0("outs/markers/", id, "_top_markers.Rds"))
write.csv(top_markers, file=paste0("outs/markers/", id, "_top_markers.csv"))

top10_markers <- top_n(top_markers, 10, avg_log2FC) %>% arrange(cluster)

p1 <- Seurat::DoHeatmap(sobj, features = top10_markers$gene, slot='scale.data',
                  group.by = "spatial.cluster", 
                  angle=0) + guides(col = FALSE)
ggsave(p1, file=paste0("outs/markers/", id, "_heatmap.pdf"), dpi=300, width=7, height=7,units="in")
