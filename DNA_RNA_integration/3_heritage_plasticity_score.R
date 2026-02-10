library(Seurat)
library(tidyverse)
library(forcats)
library(stringr)
library(fs)
library(copykit)
library(VennDiagram)
library(SummarizedExperiment)
library(ggpubr)
library(ruok)
library(scales)
library(ComplexHeatmap)
library(UCell)
library(ggtree)
library(phylogram)
library(dendextend)
library(phytools)
library(ape)
library(dendsort)
library(rstatix)
library(patchwork)
library(msigdbr)

obj_rna
tumor <- obj_rna[,  obj_rna$subclones %in% c("C1","C2","C3")]
tumor

DefaultAssay(tumor) <- "SCT"
msig_h <- msigdbr(species="Homo sapiens", category="H") %>%
  dplyr::select(gs_name, gene_symbol)
db <- split(msig_h$gene_symbol, msig_h$gs_name)
str(db)

tumor <- AddModuleScore(
  tumor,
  features = db,
  name = "MSIGDB"
)
for (i in seq_along(db)) {
  tumor@meta.data[[names(db)[i]]] <- tumor@meta.data[[paste0("MSIGDB", i)]]
  tumor@meta.data[[paste0("MSIGDB", i)]] <- NULL
}
df_signature_res <- tumor@meta.data[, names(db)]
df_meta <- tumor@meta.data
mat_pheno_viz <- df_meta[, c(names(db), "subclones")]

mat_pheno_viz <- aggregate(
  mat_pheno_viz[, names(db)],
  by = list(subclones = mat_pheno_viz$subclones),
  FUN = mean,
  na.rm = TRUE
)

rownames(mat_pheno_viz) <- mat_pheno_viz$subclones
mat_pheno_viz$subclones <- NULL
mat_pheno_viz <- as.matrix(mat_pheno_viz)
mat_pheno_viz

#viz module score by Vln Plot or Feature plot
VlnPlot(
  tumor, 
  feature = names(db)[22], 
  group.by = "subclones", #seurat_clusters
  pt.size = 0.5
)

mat_pheno_viz <- df_meta %>%
  dplyr::select(all_of(c(phenotypes_opts, z))) %>%
  dplyr::group_by_at(z) %>%
  dplyr::summarise_all(mean, na.rm = TRUE) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames(var = z) %>% #z是subclones
  as.matrix()
mat_pheno_viz_z <- scale(mat_pheno_viz)

#############################################################
# CNA distance
cna_dist_mat <- tmp_cna
stopifnot(identical(rownames(cna_dist_mat), colnames(cna_dist_mat)))

cna_dist <- as.dist(cna_dist_mat)
head(cna_dist)
#############################################################

mat_pheno_viz_z <- mat_pheno_viz_z[rownames(cna_dist_mat), phenotypes_opts]

calc_heritability_score <- function(d1, d2) {
lb <- labels(d1)
m2 <- as.matrix(d2)
m2 <- m2[lb, lb]
d2 <- as.dist(m2)

o <- cor.test(d1, d2, method = "pearson")

c(
  r = as.numeric(o$estimate),
  pval = as.numeric(o$p.value)
)
}

res <- sapply(
  phenotypes_opts,
  function(xx) {
    calc_heritability_score(
      cna_dist,
      dist(mat_pheno_viz_z[, xx, drop = FALSE], method = "euclidean")
    )
  }
)

res <- t(res)
res <- as.data.frame(res)
res$module <- rownames(res)
res
res$plasticity <- 1 - res$r 


###lollipop plot

r_vec <- res$r
names(r_vec) <- res$module

r_vec <- r_vec[colnames(module_score_mat)]

r_mat <- matrix(r_vec, nrow = 1)
colnames(r_mat) <- colnames(module_score_mat)
rownames(r_mat) <- "r"

df <- data.frame(
  pathway = colnames(r_mat),
  r = as.numeric(r_mat["r", ])
)

df$pathway <- sub("^HALLMARK_", "", df$pathway)

df <- df[order(df$r), ]

p1 <- ggplot(df, aes(y = reorder(pathway, r), x = r)) +
  geom_segment(
    aes(x = 0, xend = r, yend = pathway),
    linewidth = 0.8,
    color = "grey40"
  ) +
  geom_point(
    aes(color = r > 0),
    size = 3
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(
    values = c(`TRUE` = "#C2412D", `FALSE` = "#2F6DAE"),
    guide = "none"
  ) +
  labs(
    x = "Pearson correlation (r)",
    y = NULL
  ) +
  theme_classic(base_size = 12)
p1
ggsave("pathway_heritage_score.pdf", plot = p1,width = 6,height = 10)




