library(BayesSpace)
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
nclus <- as.numeric(args[2])
jp <- as.numeric(args[3])
s <- readRDS(paste0("outs/dataObjects/", id, ".Rds"))

# Cluster the sample

s <- spatialCluster(s, q=nclus, platform="Visium", d=15, init.method="mclust", model="t") 
saveRDS(s, file = paste0("outs/dataObjects/", id, ".Rds"))

# Save the basic res clusterPlot

p1 <- clusterPlot(s, color="black") +theme_bw() + xlab("Column") + ylab("Row") + labs(fill="BayesSpace\ncluster", title=paste0("Spatial clustering", id))
ggsave(file = paste0("outs/clusterPlots/",id, "_normal_spot.pdf"), plot=p1, width=7, height=7, dpi=300, units="in") 

# Enhance the res iteratively to find optimal jitter scale

for(i in seq(1,5,by=0.25)){
  s.en <- spatialEnhance(s, q=nclus, platform="Visium", d=15, nrep=7000, burn.in=250,jitter_scale=i,save.chain=T,verbose=TRUE)
  chain <- mcmcChain(s.en, "Ychange")
  saveRDS(chain, file = paste0("outs/chains/", id, "/mcmc-", as.character(i) ,".rds"))
}

# Melt MCMC chain data together to plot on one graph

x <- data.frame(lapply(seq(1,5,by=0.25), function(x) readRDS(paste0("outs/chains/", id, "/mcmc-", as.character(x), ".rds"))))[2:71,]
colnames(x) <- seq(1,5,by=0.25)
x$iter <- seq(1,70)
x.melt <- melt(x, id.vars="iter")
p2 <- ggplot(x.melt, mapping=aes(x=iter, y=value, col=variable)) + geom_line() + xlab("Iterations (100x)") + ylab("Ychange") + labs(color="Jitter scale")
ggsave(file = paste0("outs/chains/", id, "/chains-plot.pdf"), plot = p2, width = 7, height =7, dpi = 300, units = "in")
