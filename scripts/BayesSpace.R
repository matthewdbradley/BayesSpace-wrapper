library(BayesSpace)
library(Seurat)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
id <- args[1]
n <- as.numeric(args[2])
jp <- as.numeric(args[3])
js <- as.numeric(args[4])
s <- readRDS(paste0("outs/dataObjects/", id, ".Rds"))
s.en <- spatialEnhance(s, q=n, platform="Visium", d=15,
                                   jitter_scale=js,jitter_prior=jp)
saveRDS(file=paste0("outs/dataObjects/", id, "_enhanced.Rds"), s.en)
p1 <- clusterPlot(s.en)
ggsave(paste0("outs/clusterPlots/", id,"_subspots.pdf"), plot=p1, width=7, height=7, dpi=300, units="in")
