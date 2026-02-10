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

load("/gpfs/gibbs/project/liu_yang/xd97/CNV_ref_genome/copykit_environment_human.RData")

load("1708_copykit_addRNAcelltype_adt_V3_260122.rdata")

copykit<- runPhylo(copykit, metric = 'manhattan')

#Visualiza
plotPhylo(copykit, label = 'subclones', consensus = TRUE, label_colors=subclone_pal)

###===============load function=================###
get_tree_tips_order <- function(phylo) {
  # Ref: https://github.com/Puriney/copykit/blob/faf47072f34130dda263037fda637520899ca013/R/plotHeatmap.R#L359-L362
  is_tip <- phylo$edge[, 2] <= length(phylo$tip.label)
  ordered_tips_index <- phylo$edge[is_tip, 2]
  phylo_tips_order <- phylo$tip.label[ordered_tips_index] %>% rev()
  return(phylo_tips_order)
}

calc_cophenetic_cor <- function(phylo, hc) { 
  a <- cophenetic(phylo) ## a matrix
  b <- cophenetic(hc) # a dist
  a <- a[labels(b), labels(b)]
  a <- as.dist(a)
  cor(a, b)
}

get_hclust_tip_order <- function(hc){ hc$labels[hc$order] }
get_phylo_tip_order <- function(tree) {
  # ref: https://stackoverflow.com/a/34364914/1608734
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  return(tree$tip.label[ordered_tips])
}

###==============step1: calculate DNA distance=================
css_mat <- copykit@consensus 
dist_mat <- as.matrix(dist(t(css_mat), method = "manhattan"))

library(ape)
dist_obj <- as.dist(dist_mat)

me_tree <- ape::fastme.bal(dist_obj)

me_tree <- ape::root(me_tree, outgroup = "Diploid", resolve.root = TRUE)

library(ggtree)

subclones <- factor(colnames(css_mat))
list_samples <- split(subclones, subclones)  

treeplt <- ggtree(me_tree) %>%
  ggtree::groupOTU(list_samples) %<+% 
  data.frame(group=subclones) +  
  geom_tiplab(aes(color=group), size=5, hjust=0.9, vjust=1.8) +
  geom_tippoint(aes(color=group), size=5) +
  scale_color_manual(values = subclone_pal) +
  theme_tree2() +
  theme(legend.position="right")
treeplt

library(cowplot)
cowplot::ggsave2(
  filename = "1708_ME_tree.pdf",
  plot = treeplt,
  width = 6,
  height = 5
)

write_rds(me_tree, file.path(dir_res, 'Genome146_copykit_fastme_tree.rds'))
write_rds(me_tree, file='Genome1708_copykit_fastme_tree.rds')
write_rds(
  dist_mat, 
  file='Genome1708_copykit_dist_manhattan.mat.rds')

write_rds(treeplt, file='Genome1708_copykit_ggtree_me_tree.rds')

all_clones <- colnames(css_mat)
other_clones <- setdiff(all_clones, "Diploid")
most_diploid_nb <- other_clones[which.min(dist_mat["Diploid", other_clones])]

me_tree <- ladderize(me_tree)
me_tree_reroot <- phytools::reroot(me_tree, node.number = which(me_tree$tip.label == most_diploid_nb))
me_tree_nodipoid <- ape::drop.tip(me_tree_reroot, tip="Diploid")

subclones <- factor(other_clones)
levels(subclones)

me_tree_nodipoid$tip.label

ggtree_misc_meta <- data.frame(
  Newick_label = me_tree_nodipoid$tip.label,
  group = me_tree_nodipoid$tip.label 
)

list_samples <- as.list(me_tree_nodipoid$tip.label)
names(list_samples) <- me_tree_nodipoid$tip.label

treeplt_reroot <- ggtree(me_tree_nodipoid, size = 0.2) +
  geom_tiplab(aes(color = label), size = 5, hjust = -0.2, vjust = 0.4) +
  geom_tippoint(aes(color = label), size = 5) +
  scale_color_manual(values = subclone_pal) +
  theme_tree2() +
  theme(legend.position = "right")

treeplt_reroot

cowplot::ggsave2(
  filename = "1708_ME_tree_rerooted_noDiploid.pdf",
  plot = treeplt_reroot,
  width = 6, height = 5
)

write_rds(me_tree_reroot, file='Genome1708_copykit_fastme_tree_rerooted.rds')
write_rds(treeplt_reroot, file='Genome1708_copykit_ggtree_me_tree_rerooted.rds')

cna_dist_val_range <- range(dist_mat)


###==============step2: calculate RNA distance=================
#load RNA object with CNV sunclone info
load("1708_RNAobj_addCNVclone_adt_V3_260122.rdata")
pbmc=obj_rna

DefaultAssay(pbmc) <- "SCT"
n_pcs <- ncol(Embeddings(pbmc, reduction = "pca"))

pc_idx <- 1:30
rna_events <- Embeddings(pbmc, reduction = "pca")[, pc_idx]

dict_cell2clone <- pbmc@meta.data$subclones
names(dict_cell2clone) <- Cells(pbmc)
dict_cell2clone <- forcats::fct_drop(dict_cell2clone)
dict_cell2clone

rna_events_clone <- apply(rna_events, 2, function(vv) {
  tapply(vv, dict_cell2clone, median)
})

rna_dist <- dist(rna_events_clone)
rna_dist <- as.matrix(rna_dist)

ident_vec=subclones
ident_vec
rna_dist <- rna_dist[ident_vec, ident_vec]
dist_mat <- dist_mat[ident_vec, ident_vec]

#save RNA_distance
print(rna_dist)
print(range(rna_dist))
print(dist_mat)

dir=getwd()

write_csv(
  as.data.frame(rna_dist),
  file.path(dir, "rna_pairwise_distance.csv")
)

write_csv(
  as.data.frame(dist_mat),
  file.path(dir, "cna_pairwise_distance.csv")
)

#rescale to CNA distance level
rna_dist <- scales::rescale(
  x  = as.matrix(as.dist(rna_dist)),
  to = cna_dist_val_range
)

#------ concordance ------
cna_dist_mat = dist_mat
stopifnot(all.equal(rownames(cna_dist_mat), rownames(rna_dist)))
stopifnot(all.equal(colnames(cna_dist_mat), colnames(rna_dist)))

#Method1: Pearson correlation
score_GP_concordance_obj <- cor.test(
  as.dist(cna_dist_mat), 
  as.dist(rna_dist), method = 'pearson')
score_GP_concordance <- as.numeric(score_GP_concordance_obj$estimate)
cor.test(as.dist(cna_dist_mat), 
         as.dist(rna_dist), 
         method = 'pearson')
sample_name = "GENO1708_copykit"
cat('score_GP_concordance = ', signif(score_GP_concordance, digits = 2), 'in sample', sample_name, '.\n')
write_lines(score_GP_concordance,
            file = 'pearsonscore_GP_concordance.txt')


library(ade4)
run_mantel_rtest <- function(m1, m2){
  lbs <- rownames(m1)
  stopifnot(identical(intersect(lbs, rownames(m1)), lbs))
  stopifnot(identical(intersect(lbs, rownames(m2)), lbs))
  m1 <- m1[lbs, lbs]
  m2 <- m2[lbs, lbs]
  set.seed(42)
  mantel_res <- ade4::mantel.rtest(as.dist(m1), as.dist(m2))
  # print(mantel_res)
  return(c('pvalue'=mantel_res$pvalue, 'r'=mantel_res$obs))
}
score_GP_concordance_pval <- as.numeric(
  score_GP_concordance_obj$p.value)

#------ viz: concordance ------------------------
genotree_tip_orders <- get_phylo_tip_order(me_tree_nodipoid)

print(genotree_tip_orders)
write_lines(genotree_tip_orders, file='tree_tip_order.txt')

#check
stopifnot(all(genotree_tip_orders %in% rownames(cna_dist_mat)))
stopifnot(all(genotree_tip_orders %in% rownames(rna_dist)))

#Step 1
clone_order <- rev(genotree_tip_orders)  
clone_order

tmp_cna  <- cna_dist_mat[clone_order, clone_order]
tmp_rna <- rna_dist[clone_order, clone_order]

tmp_cna
tmp_rna

#Step 2
library(circlize)
library(viridis)

col_fun <- circlize::colorRamp2(
  quantile(c(c(tmp_cna), c(tmp_rna)), c(0.01, 0.99)),
  hcl.colors(2, "Peach")
)

#Step 3
library(ComplexHeatmap)
library(grid)

pt_cna <- Heatmap(
  tmp_cna,
  name = "CNA dist",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_names_side = "left",
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (i >= j) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
      if (tmp_cna[i, j] > 0) {
        grid.text(tmp_cna[i, j], x, y, gp = gpar(fontsize = 9))
      }
    }
  }
)

pt_cna

#Step 4
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

pt_rna <- Heatmap(
  tmp_rna,
  name = "rna dist",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_names_side = "left",
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (i >= j) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
      if (tmp_rna[i, j] > 0) {
        lab <- if (is.wholenumber(tmp_rna[i, j])) {
          as.integer(tmp_rna[i, j])
        } else {
          sprintf("%.2f", tmp_rna[i, j])
        }
        grid.text(lab, x, y, gp = gpar(fontsize = 9))
      }
    }
  }
)

pt_rna

#Step 5
tmp_diff <- abs(tmp_cna - tmp_rna)

#Step 6
pt_diff <- Heatmap(
  tmp_diff,
  name = "abs diff",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = circlize::colorRamp2(
    c(cna_dist_val_range[1],
      mean(cna_dist_val_range),
      cna_dist_val_range[2]),
    hcl.colors(3, "Cyan-Magenta")
    #hcl.colors(3, "YlOrBr")
  ),
  row_names_side = "left",
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (i >= j) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
      if (tmp_diff[i, j] > 0) {
        grid.text(sprintf("%.2f", tmp_diff[i, j]), x, y, gp = gpar(fontsize = 9))
      }
    }
  }
)

pt_diff

#Step 7
pdf(file="heatmap.GP_concordance_Pearson_color.pdf",
    width = 9,
    height = 3)

draw(
  pt_cna + pt_rna + pt_diff,
  column_title = sprintf(
    "concordance = %.2f (p = %.3g)",
    score_GP_concordance,
    score_GP_concordance_pval
  )
)
dev.off()
