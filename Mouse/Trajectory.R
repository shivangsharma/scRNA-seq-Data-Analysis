library(Seurat)
library(DropletUtils)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(pheatmap)
library(readxl)
library(RColorBrewer)
library(BiocManager)
library(minet)
library(kableExtra)
library(Matrix)
library(writexl)
library(scales)
library(RColorBrewer)
library(htmlwidgets)
library(ggrepel)
library(plotly)
library(Hmisc)
library(slingshot)
library(tradeSeq)
library(monocle3)
library(SeuratWrappers)



obj <- readRDS("~/Desktop/Mouse/Step 2 - Analysis/D10A vs H/Comparison/Urothelial/Combined (U+P)/Combined_Combined_SeuratObject.rds")
obj@assays$integrated@counts <- obj@assays$integrated@data
# ----------------------------------Slingshot-----------------------------------
sdsP <- slingshot(Embeddings(obj,"umap"), clusterLabels = obj$seurat_clusters, start.clus = 8, stretch = 0)
sds <- as.SlingshotDataSet(sdsP)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
cell_colors <- cell_pal(obj$seurat_clusters, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(obj$seurat_clusters, hue_pal())
plot(Embeddings(obj, "umap"), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
plot(Embeddings(obj,"umap"), col = brewer.pal(9,"Set1")[obj$seurat_clusters], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
lin1 <- getLineages(Embeddings(obj,"umap"), obj$seurat_clusters, start.clus = '8')
crv1 <- getCurves(lin1)
df <- as.data.frame(sdsP@assays@data@listData[["pseudotime"]])
df[is.na(df)] <- 0
df$`Avg P_Time` <- pmax(df$Lineage1,df$Lineage2)
df <- df[order(df$`Avg P_Time`),]
# ------------------------------------------------------------------------------

# -----------------------------------Monocle3-----------------------------------
obj_monocle3 <- as.cell_data_set(obj)
obj_monocle3 <- cluster_cells(cds = obj_monocle3, reduction_method = "UMAP", cluster_method = "louvain", num_iter = 20)
obj_monocle3 <- learn_graph(obj_monocle3, use_partition = TRUE)
plot_cells(obj_monocle3, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE, cell_size = 0.75)
start_cells <- WhichCells(obj, ident = "0")
obj_monocle3 <- order_cells(obj_monocle3, reduction_method = "UMAP", root_cells = start_cells)
# Urothelial_UMAP_Combined_Monocle3_Leiden
plot_cells(cds = obj_monocle3, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 0.75, label_roots = FALSE, label_branch_points = TRUE, label_leaves = TRUE, trajectory_graph_color = "green")
# ------------------------------------------------------------------------------


# ------------Create Monocle3 object from Seurat integrated dataset-------------
gene_annotations <- as.data.frame(rownames(obj), row.names = rownames(obj))
colnames(gene_annotations) <- "gene_short_name"
cell_metadata <- as.data.frame(colnames(obj), row.names = colnames(obj))
colnames(cell_metadata) <- "barcode"
expression_matrix <- obj@assays$RNA@counts
cds <- new_cell_data_set(expression_data = expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotations)
cds <- preprocess_cds(cds, num_dim = 30, pseudo_count = 10)
cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP", umap.metric = "cosine", umap.min_dist = 0.3, umap.n_neighbors = 30L, umap.nn_method = "annoy")
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 20, cluster_method = "leiden", num_iter = 20)
recreate.partitions <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData$UMAP$partitions <- recreate.partitions
list_cluster <- obj@meta.data$CellAnn
names(list_cluster) <- colnames(obj)
cds@clusters@listData$UMAP$clusters <- list_cluster
cds@clusters@listData$UMAP$louvain_res <- "NA"
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- obj@reductions$umap@cell.embeddings
cds@int_colData@listData[["reducedDims"]]@listData[["PCA"]] <- obj@reductions$pca@cell.embeddings
cds@preprocess_aux@listData[["gene_loadings"]] <- obj@reductions$pca@feature.loadings
cds <- learn_graph(cds)
plot_cells(cds, label_principal_points = TRUE)
Idents(obj) <- "CellAnn"
root_cells <- WhichCells(obj, ident = "Basal")
cds <- order_cells(cds, root_pr_nodes = "Y_42")
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.75, trajectory_graph_color = "green")
test <- graph_test(cds, neighbor_graph = "principal_graph")
test <- subset(test, q_value < 0.05)
gene_module_df <- find_gene_modules(cds[row.names(test), ])
cell_group_df <- tibble::tibble(cell = row.names(colData(cds)), cell_group = cds@clusters@listData[["UMAP"]][["clusters"]])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale = "column", clustering_method = "ward.D2")
# ------------------------------------------------------------------------------

obj_sce <- as.SingleCellExperiment(obj)


counts = obj@assays[["H_Urothelial"]]@scale.data
pseudotime = sdsP@assays@data@listData[["pseudotime"]]
cellWeights = sdsP@assays@data@listData[["weights"]]
pseudotime[is.na(pseudotime)] <- 0
obj.markers <- read_excel("~/Desktop/Mouse/Step 2 - Analysis/D10A vs H/Comparison/Urothelial/Combined/Combined_Combined_Markers.xlsx")
genes <- unique(obj.markers$gene)
obj_sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, genes = genes)
ATres <- associationTest(obj_sce)
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
heatdata <- obj@assays[["H_Urothelial"]]@scale.data[topgenes, rownames(df)]
df2 <- as.data.frame(obj@meta.data)
heatclust <- df2[rownames(df), ]$seurat_clusters