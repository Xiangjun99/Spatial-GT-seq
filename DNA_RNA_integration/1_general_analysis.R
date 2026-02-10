setwd("/gpfs/ycga/work/liu_yang/xd97/@@2026_coprofile_revised/8_human_1708_72/1_GENO1708_resequence_merged3/4-Integration_RNA_CNVfilter_ADT")

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
library(aCGH)
library(tidyverse)


#------add RNA cluster information to CNV copykit----#####
#load RNA object
load("RNA1708_celltype_V2.rdata")
obj_rna

#load CNV object
load("1708_copykit.rdata")
copykit

#load ADT object
load("ADT_celltype_V2.rdata")
obj_adt

#extract CNV metadata and match to RNA based on barcode AxB
name_dna <- copykit@colData %>% 
  as.data.frame() %>% 
  dplyr::select(c("subclones","sample"))  
head(name_dna)

#add RNA_ to RNA object metadata to avoid mix with DNA metadata
colnames(obj_rna@meta.data) <- ifelse(
  stringr::str_detect(colnames(obj_rna@meta.data), "RNA"), 
  colnames(obj_rna@meta.data), 
  paste0("RNA_", colnames(obj_rna@meta.data))
)

DNA_rna_meta <- obj_rna@meta.data %>%
  rownames_to_column(var = "sample") %>%  
  left_join(name_dna, by = "sample") %>% 
  dplyr::select(
    sample,
    nCount_RNA, nFeature_RNA,  RNA_nCount_SCT, RNA_nFeature_SCT,
    subclones, RNA_celltype
  )
head(DNA_rna_meta)

#add ADT_ to protin object metadata to avoid mix with DNA metadata
colnames(obj_adt@meta.data) <- ifelse(
  stringr::str_detect(colnames(obj_adt@meta.data), "ADT"), 
  colnames(obj_adt@meta.data), 
  paste0("ADT_", colnames(obj_adt@meta.data))
)
name_adt = data.frame(obj_adt@active.ident)
head(name_adt)
colnames(name_adt)[1] <- "adt_cluster"

#add adt_seurat clusters to dna_rna_meta
DNA_rna_adt_meta <- DNA_rna_meta %>%
  left_join(
    name_adt %>% rownames_to_column("sample"),
    by = "sample"
  )

#add RNA metadata to CNV copykit
#remove row with NA 
DNA_rna_adt_meta1 <- DNA_rna_adt_meta %>% 
  drop_na()
rownames(DNA_rna_adt_meta1) <- DNA_rna_adt_meta1$sample
DNA_rna_adt_meta1_aligned <- DNA_rna_adt_meta1[
  match(colnames(copykit), rownames(DNA_rna_adt_meta1)),
]

#check
identical(
  colnames(copykit),
  rownames(DNA_rna_adt_meta1_aligned)
)

copykit@colData <- cbind(
  copykit@colData,
  DNA_rna_adt_meta1_aligned
)
save(copykit, file = "./results_0122/1708_copykit_addRNAcelltype_adt.rdata")

#CNV heatmap plot
load("/gpfs/gibbs/project/liu_yang/xd97/CNV_ref_genome/copykit_environment_human.RData")

col_fun=circlize::colorRamp2(breaks = c(-1,0,1), 
                             c("#74add1","white","#fd7f64"))

#set color for RNA cluster
rna_pal <- c(
  "Tumor1"="#F9D367",
  "CAF1"="#63B5B7",
  "Tumor2"="#C2412D",
  "CAF2"="#4C6A9A",
  "Immune cells"="#8C6BB1",
  "Tumor3"="#b15928"  ##FB9A99
)

subclone_pal = c(
  'Diploid' = "#5085C4",
  'C1' = "#EB545C",
  'C2' = "#F6EC1B",
  'C3' = '#87C55F',
  'C4' = '#CC66FF')

adt_pal =c( "Tumor" = "#fb9a99",
            "CAF" = "#2ec4b6",
            "Immune cells" = "#D4A2D9"
            )

pdf("./results_0122/Integration_CNV_heatmap_RNA_ADT_annotation.pdf", width = 16, height = 8)
plotHeatmap(
  copykit,
  genome = "hg38",
  #consensus = TRUE,
  label = c('subclones','adt_cluster','RNA_celltype'),
  order_cells = 'consensus_tree',
  #order_cells = "hclust",
  col = col_fun,
  label_colors = c(list(subclones = subclone_pal,adt_cluster=adt_pal,RNA_celltype=rna_pal)),
  n_threads = 100
)
dev.off()

#spatial——plot
#add RNA UMAP info
# 对齐行名再合并
umap_mat <- obj_rna@reductions$umap@cell.embeddings
umap_mat <- umap_mat[DNA_rna_adt_meta$sample, ]  
DNA_rna_adt_meta2 <- cbind(DNA_rna_adt_meta, umap_mat)

head(DNA_rna_adt_meta2)
colnames(DNA_rna_adt_meta2)[
  colnames(DNA_rna_adt_meta2) %in% c("umap_1", "umap_2")
] <- c("RNA_umap_1", "RNA_umap_2")

#add adt UMAP info
umap_mat <- obj_adt@reductions$umap@cell.embeddings
umap_mat <- umap_mat[DNA_rna_adt_meta$sample, ]  
DNA_rna_adt_meta3 <- cbind(DNA_rna_adt_meta2, umap_mat)

head(DNA_rna_adt_meta3)
colnames(DNA_rna_adt_meta3)[
  colnames(DNA_rna_adt_meta3) %in% c("umap_1", "umap_2")
] <- c("ADT_umap_1", "ADT_umap_2")

###spatial plot
test <- DNA_rna_adt_meta3 %>% separate(sample, c("A", "B"),  sep = "x")

#CNV plot
pdf(file = paste("./results_0122/Integrated_CNV_Spatial_clustering.pdf",sep =""), width=11.6, height=11.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=subclones)) + 
  scale_color_manual(values = subclone_pal)  + 
  ggtitle("UMAP") +
  geom_point(shape = 16, size = 4.0)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  coord_equal(xlim=c(0,73),ylim=c(73,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        #axis.text=element_text(size=20),
        #axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +   
  theme(plot.background = element_rect(fill = "black"))
dev.off()

pdf(file = paste("./results_0122/Integrated_RNA_Spatial_clustering.pdf",sep =""), width=11.6, height=11.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=RNA_celltype)) + 
  scale_color_manual(values = rna_pal)  + 
  ggtitle("UMAP") +
  geom_point(shape = 16, size = 4.0)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  coord_equal(xlim=c(0,73),ylim=c(73,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        #axis.text=element_text(size=20),
        #axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +   
  theme(plot.background = element_rect(fill = "black"))
dev.off()

pdf(file = paste("./results_0122/Integrated_ADT_Spatial_clustering.pdf",sep =""), width=11.6, height=11.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=adt_cluster)) + 
  scale_color_manual(values = adt_pal)  + 
  ggtitle("UMAP") +
  geom_point(shape = 16, size = 4.0)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  coord_equal(xlim=c(0,73),ylim=c(73,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        #axis.text=element_text(size=20),
        #axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +   
  theme(plot.background = element_rect(fill = "black"))
dev.off()

#UMAP plot
p1 <- ggplot(DNA_rna_adt_meta3,aes(x = RNA_umap_1, y = RNA_umap_2, fill = RNA_celltype)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = rna_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
ggsave("./results_0122/1708_RNAUMAP_all_cells_color_by_rna_clusters.pdf", plot = p1, width = 8, height = 6, dpi = 300)

p2 <- ggplot(DNA_rna_adt_meta3,aes(x = RNA_umap_1, y = RNA_umap_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = subclone_pal, na.translate =F) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
ggsave("./results_0122/1708_RNAUMAP_all_cells_color_by_dna_clusters.pdf", plot = p2, width = 8, height = 6, dpi = 300)

p3 <- ggplot(DNA_rna_adt_meta3,aes(x = RNA_umap_1, y = RNA_umap_2, fill = adt_cluster)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = adt_pal, na.translate =F) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p3
ggsave("./results_0122/1708_RNAUMAP_all_cells_color_by_adt_clusters.pdf", plot = p3, width = 8, height = 6, dpi = 300)

write.csv(DNA_rna_adt_meta3, file="./results_0122/1708_DNA_rna_adt_meta3_containUMAP.csv")

#add adt_cluster and subclone info to RNA_cluster
#check order
all(colnames(obj_rna) %in% rownames(DNA_rna_adt_meta3))

obj_rna$subclones <- DNA_rna_adt_meta3$subclones
obj_rna$adt_cluster <- DNA_rna_adt_meta3$adt_cluster
save(obj_rna,file="./results_0122/1708_RNAobj_addCNVclone_adt_V3.rdata")

#add RNA_celltype and subclone info to adt
#check order
all(colnames(obj_adt) %in% rownames(DNA_rna_adt_meta3))

obj_adt$subclones <- DNA_rna_adt_meta3$subclones
obj_adt$RNA_celltype <- DNA_rna_adt_meta3$RNA_celltype
save(obj_adt,file="./results_0122/1708_ADTobj_addCNVclone_RNAcelltype_V3.rdata")

#plots
#CNV UMAP as based
p4 <- plotUmap(copykit, label = 'subclones') + 
  scale_fill_manual(values = subclone_pal, name = "subclones")
p4
ggsave("./results_0122/CNVUMAP_singlePlot/CNV_UMAP.pdf", plot = p4, width = 8, height = 6, dpi = 300)

p5 <- plotUmap(copykit, label = 'RNA_celltype') + 
  scale_fill_manual(values = rna_pal, name = "RNA_celltype")
p5
ggsave("UMAP-3clones-renamed_260122.pdf", plot = p2, width = 8, height = 6, dpi = 300)

#plot RNA celltype

outdir <- "./results_0122/CNVUMAP_singlePlot"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

celltypes <- names(rna_pal)

umap_df <- as.data.frame(copykit@int_colData@listData[["reducedDims"]])
umap_df$cell <- rownames(umap_df)
umap_df$RNA_celltype <- copykit$RNA_celltype

for (ct in celltypes) {
  
  df_bg <- umap_df[umap_df$RNA_celltype != ct, ]
  df_fg <- umap_df[umap_df$RNA_celltype == ct, ]  # 高亮点
  
  p <- ggplot() +
    geom_point(
      data = df_bg,
      aes(x = umap.V1, y = umap.V2),
      shape = 21,
      size = 2.5,
      stroke = 0.03,
      fill = "grey85",
      color = "grey85"
    ) +
    geom_point(
      data = df_fg,
      aes(x = umap.V1, y = umap.V2),
      shape = 21,
      size = 2.5,
      stroke = 0.03,
      fill = rna_pal[ct],
      color = "black"
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle(ct)
  
  ggsave(
    filename = file.path(
      outdir,
      paste0("CNVUMAP_RNAhighlight_", gsub(" ", "_", ct), ".pdf")
    ),
    plot = p,
    width = 4,
    height = 3
  )
}


p6 <- plotUmap(copykit, label = 'adt_cluster') + 
  scale_fill_manual(values = adt_pal, name = "adt_cluster")
p6
ggsave("UMAP-3clones-renamed_260122.pdf", plot = p2, width = 8, height = 6, dpi = 300)

#plot adt cluster
celltypes <- names(adt_pal)

umap_df <- as.data.frame(copykit@int_colData@listData[["reducedDims"]])
umap_df$cell <- rownames(umap_df)
umap_df$adt_cluster <- copykit$adt_cluster

for (ct in celltypes) {

  df_bg <- umap_df[umap_df$adt_cluster != ct, ]
  df_fg <- umap_df[umap_df$adt_cluster == ct, ] 
  
  p <- ggplot() +
  
    geom_point(
      data = df_bg,
      aes(x = umap.V1, y = umap.V2),
      shape = 21,
      size = 2.5,
      stroke = 0.03,
      fill = "grey85",
      color = "grey85"
    ) +
   
    geom_point(
      data = df_fg,
      aes(x = umap.V1, y = umap.V2),
      shape = 21,
      size = 2.5,
      stroke = 0.03,
      fill = adt_pal[ct],
      color = "black"
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle(ct)
  
  ggsave(
    filename = file.path(
      outdir,
      paste0("CNVUMAP_ADThighlight_", gsub(" ", "_", ct), ".pdf")
    ),
    plot = p,
    width = 4,
    height = 3
  )
}

###CNV polytree plot
#Visualiza
pdf("./results_0122/Phylotree—3clones_color.pdf", width = 8, height = 6)
plotPhylo(copykit, label = 'subclones', consensus = TRUE, label_colors=subclone_pal)
dev.off()

###each modality UMAP plot
copykit

obj_rna

p1 <- ggplot(DNA_rna_adt_meta3,aes(x = RNA_umap_1, y = RNA_umap_2, fill = RNA_celltype)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = rna_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
ggsave("./results_0122/RNA_UMAP.pdf", plot = p1, width = 8, height = 6, dpi = 300)


obj_adt
p2 <- ggplot(DNA_rna_adt_meta3,aes(x = ADT_umap_1, y = ADT_umap_2, fill = adt_cluster)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = adt_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
ggsave("./results_0122/ADT_UMAP.pdf", plot = p2, width = 8, height = 6, dpi = 300)
