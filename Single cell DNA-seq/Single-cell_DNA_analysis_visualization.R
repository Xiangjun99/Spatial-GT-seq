setwd("/path/to/single-cell_dna/")

library(copykit)
library(SummarizedExperiment)
library(magrittr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(cowplot)
library(GenomicRanges)
library(BiocParallel)

register(MulticoreParam(workers = 32, progressbar = T),default = T)

copykit <- runVarbin('bamFile',genome='mm10',resolution = "1Mb",remove_Y=T,remove_X = TRUE,is_paired_end=T,min_bincount=1)
save(copykit, file = "copykit_V1.rdata")


col_fun=circlize::colorRamp2(breaks = c(-0.5,0,0.5), 
                             c("#4575b4","white","#d73027"))

copykit <- runMetrics(copykit)
copykit <- findAneuploidCells(copykit, remove_XY = T, simul = F)

########%%%%%%%% step3: cluster analysis %%%%%%%%########
copykit <- runUmap(copykit)

subclone_pal = c(
  'c1' = "#EB545C",
  'c2' = "#5085C4",
  'c3' = "#F6EC1B",
  'c4' = '#87C55F',
  'c5' = '#CC66FF',
  'c6' = '#59A14F')

copykit  <- findClusters(copykit, k_subclones = 20, method ="leiden") #"hdbscan", "leiden", "louvain"

p2 <- plotUmap(copykit, label = 'subclones') + 
  scale_fill_manual(values = subclone_pal, name = "subclones")
p2
ggsave("single-cell_UMAP.pdf", plot = p2, width = 8, height = 8, dpi = 300)


# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
copykit <- calcConsensus(copykit)
copykit <- runConsensusPhylo(copykit)

pdf("single-cell-heatmap.pdf", width = 16, height = 8)
plotHeatmap(
  copykit,
  genome = "mm10",
  #consensus = TRUE,
  label = c('subclones'),
  #genes = c("TP53", "BRAF", "MYC"),
  order_cells = 'consensus_tree',
  col = col_fun,
  label_colors = list(subclones = subclone_pal),
  n_threads = 100
)
dev.off()

save(copykit, file = "copykit_V2.rdata")
