library(BayesSpace)
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 2){
	s <- readVisium(args[1])
	b <- read.csv(args[3])
	s <- s[,b$Barcode] %>% spatialPreprocess(platform="Visium", n.PCs=40)
}else{
	s <- readVisium(args[1]) %>% spatialPreprocess(platform="Visium", n.PCs=40)
}

s <- qTune(s, qs=seq(3,20))
plt <- qPlot(s)
ggsave(plt, file=paste0("outs/qPlot/", args[2],".pdf"), dpi=300, width=7, height=7,units="in")
saveRDS(s, file=paste0("outs/dataObjects/", args[2], ".Rds"))
