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

# Add RNA cluster information to CNV copykit
# Load RNA object
load("RNA1708_celltype_V2.rdata")
obj_rna

# Load CNV object
load("1708_copykit.rdata")
copykit

# Load ADT object
load("ADT_celltype_V2.rdata")
obj_adt

# Extract CNV metadata and match to RNA based on barcode AxB
name_dna <- copykit@colData %>%
  as.data.frame() %>%
  dplyr::select(c("subclones", "sample"))
head(name_dna)

# Add RNA_ prefix to RNA metadata to avoid mixing with DNA metadata
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
    nCount_RNA, nFeature_RNA, RNA_nCount_SCT, RNA_nFeature_SCT,
    subclones, RNA_celltype
  )
head(DNA_rna_meta)

# Add ADT_ prefix to protein metadata
colnames(obj_adt@meta.data) <- ifelse(
  stringr::str_detect(colnames(obj_adt@meta.data), "ADT"),
  colnames(obj_adt@meta.data),
  paste0("ADT_", colnames(obj_adt@meta.data))
)

name_adt <- data.frame(obj_adt@active.ident)
colnames(name_adt)[1] <- "adt_cluster"

# Add ADT clusters to DNA/RNA metadata
DNA_rna_adt_meta <- DNA_rna_meta %>%
  left_join(
    name_adt %>% rownames_to_column("sample"),
    by = "sample"
  )

# Remove rows with NA
DNA_rna_adt_meta1 <- DNA_rna_adt_meta %>%
  drop_na()
rownames(DNA_rna_adt_meta1) <- DNA_rna_adt_meta1$sample

DNA_rna_adt_meta1_aligned <- DNA_rna_adt_meta1[
  match(colnames(copykit), rownames(DNA_rna_adt_meta1)),
]

# Check alignment
identical(
  colnames(copykit),
  rownames(DNA_rna_adt_meta1_aligned)
)

# Add metadata to copykit object
copykit@colData <- cbind(
  copykit@colData,
  DNA_rna_adt_meta1_aligned
)
save(copykit, file = "./results_0122/1708_copykit_addRNAcelltype_adt.rdata")

# CNV heatmap
load("/gpfs/gibbs/project/liu_yang/xd97/CNV_ref_genome/copykit_environment_human.RData")

col_fun <- circlize::colorRamp2(
  breaks = c(-1, 0, 1),
  c("#74add1", "white", "#fd7f64")
)

# Color palettes
rna_pal <- c(
  "Tumor1" = "#F9D367",
  "CAF1" = "#63B5B7",
  "Tumor2" = "#C2412D",
  "CAF2" = "#4C6A9A",
  "Immune cells" = "#8C6BB1",
  "Tumor3" = "#b15928"
)

subclone_pal <- c(
  "Diploid" = "#5085C4",
  "C1" = "#EB545C",
  "C2" = "#F6EC1B",
  "C3" = "#87C55F",
  "C4" = "#CC66FF"
)

adt_pal <- c(
  "Tumor" = "#fb9a99",
  "CAF" = "#2ec4b6",
  "Immune cells" = "#D4A2D9"
)

pdf("./results_0122/Integration_CNV_heatmap_RNA_ADT_annotation.pdf", width = 16, height = 8)
plotHeatmap(
  copykit,
  genome = "hg38",
  label = c("subclones", "adt_cluster", "RNA_celltype"),
  order_cells = "consensus_tree",
  col = col_fun,
  label_colors = list(
    subclones = subclone_pal,
    adt_cluster = adt_pal,
    RNA_celltype = rna_pal
  ),
  n_threads = 100
)
dev.off()

# Spatial plots
# Add RNA UMAP coordinates
umap_mat <- obj_rna@reductions$umap@cell.embeddings
umap_mat <- umap_mat[DNA_rna_adt_meta$sample, ]
DNA_rna_adt_meta2 <- cbind(DNA_rna_adt_meta, umap_mat)

colnames(DNA_rna_adt_meta2)[
  colnames(DNA_rna_adt_meta2) %in% c("umap_1", "umap_2")
] <- c("RNA_umap_1", "RNA_umap_2")

# Add ADT UMAP coordinates
umap_mat <- obj_adt@reductions$umap@cell.embeddings
umap_mat <- umap_mat[DNA_rna_adt_meta$sample, ]
DNA_rna_adt_meta3 <- cbind(DNA_rna_adt_meta2, umap_mat)

colnames(DNA_rna_adt_meta3)[
  colnames(DNA_rna_adt_meta3) %in% c("umap_1", "umap_2")
] <- c("ADT_umap_1", "ADT_umap_2")

# Spatial coordinate plot
test <- DNA_rna_adt_meta3 %>%
  separate(sample, c("A", "B"), sep = "x")

