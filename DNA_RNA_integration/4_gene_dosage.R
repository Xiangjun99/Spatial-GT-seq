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
library(aCGH)

load("/gpfs/gibbs/project/liu_yang/xd97/CNV_ref_genome/copykit_environment_mouse.RData")

load("rna_obj.rdata")
pbmc

load("copykit.rdata")
copykit

#----combine with DNA meta data 
dna_meta_temp <- copykit@colData %>% as.data.frame() %>% dplyr::select("subclones","sample")
library(tibble)
rna_meta_temp <- pbmc@meta.data %>% as.data.frame() %>% rownames_to_column() %>% 
  left_join(dna_meta_temp, by = c("rowname" = "sample")) %>% column_to_rownames() 
rna_meta_temp2 <- cbind(rna_meta_temp, pbmc@reductions$umap@cell.embeddings)

#add CNV subclones information to RNA object
pbmc$subclones <- rna_meta_temp2$subclones
table(pbmc$subclones)
head(pbmc@meta.data)

rna_col <- c( 'normal cells' = "#5085C4",
              'c1' = "#EB545C",
                'Epithelial cells' = "#73A1BE",
                'Tumor1' = '#D87070',
                'Tumor2' = '#FCFF0C',
                'Tumor3' = "#9896f1",
                'Muscle cells' = "#80FF08"
              )
              
p <- DimPlot(pbmc, group.by = "subclones", label = TRUE) + 
  NoLegend() +
  scale_color_manual(
    values = rna_col[1:length(unique(pbmc$subclones))]
  )
p

###gene dosage-calculate mean of RNA expression
#log normalize RNA data for compare
pbmc <- NormalizeData(pbmc, assay = "RNA")

subclones <- as.vector(pbmc@meta.data$subclones)
object_bulk <- as.data.frame(t(GetAssayData(pbmc, assay = "RNA", slot = "data")))
object_bulk$subclone <- subclones
all_genes <- colnames(object_bulk)

object_bulk_mean <- object_bulk %>%
  dplyr:: group_by(subclone) %>%
  dplyr:: summarise(across(everything(),   ~mean(., na.rm = TRUE)))
object_bulk_mean <- as.data.frame(t(object_bulk_mean))

colnames(object_bulk_mean) <- object_bulk_mean[1, ]
object_bulk_mean <- object_bulk_mean[-1, ]
object_bulk_mean$gene <- rownames(object_bulk_mean)
object_bulk_mean


copykit
varbin_mtx_all_log2 <- copykit

varbin_mtx_all_log2  <- calcConsensus(varbin_mtx_all_log2 )
varbin_mtx_all_log2  <- runConsensusPhylo(varbin_mtx_all_log2 )

raw_mtx_cnv <- as.data.frame(varbin_mtx_all_log2@consensus)
#add bin name
chr_ranges <-as.data.frame(SummarizedExperiment::rowRanges(varbin_mtx_all_log2))
raw_mtx_cnv$bin <- chr_ranges$seqnames
raw_mtx_cnv$abspos <- chr_ranges$abspos
head(raw_mtx_cnv)

raw_mtx_cnv <- raw_mtx_cnv %>%
  arrange(bin, abspos) %>%                      
  mutate(bin = paste0(bin, ".", row_number()))  

raw_mtx_cnv <- as.data.frame(raw_mtx_cnv)
rownames(raw_mtx_cnv) <- raw_mtx_cnv$bin
raw_mtx_cnv$bin <- NULL
head(raw_mtx_cnv)

colnames(raw_mtx_cnv)[colnames(raw_mtx_cnv) == "normal cells"] <- "Diploid"

cnv_changes_diploid <- raw_mtx_cnv

clone_cols <- setdiff(colnames(raw_mtx_cnv), c("Diploid", "abspos"))

for (clone in clone_cols) {
  new_col_name <- paste0("cnv_", clone, "_vsDiploid")
  cnv_changes_diploid[[new_col_name]] <- raw_mtx_cnv[[clone]] - raw_mtx_cnv$Diploid
}

cnv_changes_diploid

gene_bin <- read.table("mm10_binNumber_protein_coding_map.tsv", header = 1)

#add bins information
cnv_changes_diploid$bins <- as.numeric(gsub(".*\\.", "", rownames(cnv_changes_diploid)))

gene_bin_sel <- gene_bin %>% 
  dplyr::filter(gene %in% all_genes ) #all_genes from RNA expression data:object_bulk

names(gene_bin_sel) <- c("gene", "bins")

cnv_genes_mtx_diploid <- left_join(cnv_changes_diploid, gene_bin_sel, by = "bins")

cnv_genes_mtx_diploid

cnv_genes_mean_exp_mtx_diploid <- left_join(cnv_genes_mtx_diploid,object_bulk_mean, by = "gene" )
names(cnv_genes_mean_exp_mtx_diploid) <- gsub(".y$", "_exp", names(cnv_genes_mean_exp_mtx_diploid))
names(cnv_genes_mean_exp_mtx_diploid) <- gsub(".x$", "_cnv", names(cnv_genes_mean_exp_mtx_diploid))
head(cnv_genes_mean_exp_mtx_diploid)

colnames(cnv_genes_mean_exp_mtx_diploid)[colnames(cnv_genes_mean_exp_mtx_diploid) == "normal cells"] <- "Diploid_exp"

exp_columns <- grepl("exp$",colnames(cnv_genes_mean_exp_mtx_diploid))
cnv_genes_mean_exp_mtx_diploid[exp_columns] <- lapply(cnv_genes_mean_exp_mtx_diploid[exp_columns], as.numeric)

cnv_genes_mean_exp_mtx_diploid <- cnv_genes_mean_exp_mtx_diploid %>%
  dplyr::mutate(across(ends_with("_exp"), ~ ifelse(. > 0.05, ., NA))) 

gene_per_bin <- cnv_genes_mean_exp_mtx_diploid %>%
  dplyr::group_by(bins) %>%
  dplyr::summarise(count=n()) 

median(gene_per_bin$count)

clone_exp_cols <- grep("_exp$", colnames(cnv_genes_mean_exp_mtx_diploid), value = TRUE)
clone_exp_cols <- setdiff(clone_exp_cols, "Diploid_exp") 

for (col in clone_exp_cols) {
  new_col_name <- paste0("sub_", col)  
  cnv_genes_mean_exp_mtx_diploid[[new_col_name]] <- cnv_genes_mean_exp_mtx_diploid[[col]] - cnv_genes_mean_exp_mtx_diploid$Diploid_exp
}

head(cnv_genes_mean_exp_mtx_diploid)

clone_name <- c("c1")

for (i in 1:length(clone_name)) {
  
  cnv_col <- paste0(clone_name[i], "_cnv")         
  group_col <- paste0("group_", clone_name[i])     
  
  group <- 1  
  cnv_genes_mean_exp_mtx_diploid[[group_col]] <- NA  

  for (j in 1:(nrow(cnv_genes_mean_exp_mtx_diploid)-1)) {

    cnv_genes_mean_exp_mtx_diploid[[group_col]][j] <- group
    
    if (cnv_genes_mean_exp_mtx_diploid[[cnv_col]][j+1] == cnv_genes_mean_exp_mtx_diploid[[cnv_col]][j]) {
      cnv_genes_mean_exp_mtx_diploid[[group_col]][j+1] <- group
    } else {
      group <- group + 1 
      cnv_genes_mean_exp_mtx_diploid[[group_col]][j+1] <- group
    }
  }
}

sample_mean_mtx_name <- vector("character", length = 0)
clone_name <- c("c1") 
for (i in 1:length(clone_name)) {
  group_col <- paste0("group_", clone_name[i])   
  exp_col <- paste0("sub_", clone_name[i], "_exp") 
  cnv_col <- paste0(clone_name[i], "_cnv")     
  diff_cnv_col <- paste0("cnv_",clone_name[i],"_vsDiploid" ) 
  name <- paste0("cnv_genes_exp_mean_mtx_bin_", clone_name[i])
  
  mtx <- data.frame(
    group = cnv_genes_mean_exp_mtx_diploid[[group_col]],
    cnv = cnv_genes_mean_exp_mtx_diploid[[cnv_col]], 
    sub_bin = cnv_genes_mean_exp_mtx_diploid[[exp_col]],
    bins = cnv_genes_mean_exp_mtx_diploid$bins,
    #bin = cnv_genes_mean_exp_mtx_diploid$bin,
    #gene_start = cnv_genes_mean_exp_mtx_diploid$gene_start,
    diff_cnv = cnv_genes_mean_exp_mtx_diploid[[diff_cnv_col]]
  )
  
  mtx_group <- mtx %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(sub_group = median(sub_bin, na.rm = TRUE))
  
  mtx <- left_join(mtx, mtx_group, by = "group")
  
  write.csv(mtx, file = paste0(name, "_2210bins_251126.csv"), row.names = FALSE)
  
  assign(name, mtx)
  sample_mean_mtx_name <- append(sample_mean_mtx_name, name)
}

head(cnv_changes_diploid)

chr_sizes <- read.csv("mm10_chr_bin_coords.csv", header = FALSE)
head(chr_sizes)

head(cnv_genes_exp_mean_mtx_bin_c1)

library(tidyverse)
library(cowplot)

colnames(chr_sizes) <- c("chr", "chr_end")
chr_labels <- chr_sizes$chr
chr_breaks <- chr_sizes$chr_end
chr_mid <- c(0, chr_breaks[-length(chr_breaks)]) + diff(c(0, chr_breaks))/2  # 用于在图上标 chr 名
names(chr_mid) <- chr_labels

p_cnv <- ggplot(cnv_changes_diploid, aes(x = bins, y = cnv_c1_vsDiploid)) +
  geom_line(color = "steelblue", size = 0.7, na.rm = TRUE) +
  #geom_smooth(aes(y = cnv_C4_vsDiploid), method = "loess", se = FALSE, color = "darkblue") +
  theme_cowplot() +
  labs(y = "CNV (C4 vs Diploid)", x = NULL) +
  # 染色体分界线
  geom_vline(xintercept = chr_breaks, linetype = "dashed", color = "grey70") +
  # 染色体标签
  scale_x_continuous(breaks = chr_mid, labels = chr_labels,limits = c(1, 2210),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  theme(axis.text.x = element_blank(),  # 上 panel 不显示 x
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=12)) +
  ggtitle("CNV differences")

p_cnv

p2_rna <- ggplot(cnv_genes_exp_mean_mtx_bin_c1, aes(x = bins, y = sub_group)) + 
  geom_point(aes(y = sub_bin), alpha = 0.9, size = 0.1, color = "black") +
  geom_line(color = "red", size = 0.7, na.rm = TRUE)+
  theme_cowplot() +
  labs(y = "RNA diff", x = "Genomic coordinates") +
  geom_vline(xintercept = chr_breaks, linetype = "dashed", color = "grey70") +
  scale_x_continuous(breaks = chr_mid, labels = chr_labels,expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1)) +
  ggtitle("RNA expression differences")

p2_rna


