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
library(ggrepel)

#------------------Clustering and Annotation of H_Urothelial--------------------
obj_H_Urothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/H_Urothelial_SeuratObject.rds") %>% 
  CreateSeuratObject(assay = "RNA") %>% 
  SCTransform(assay = "RNA", variable.features.n = 5000, return.only.var.genes = FALSE) %>% 
  RunPCA(npcs = 50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(algorithm = 4) %>% 
  RunUMAP(dims = 1:50)
plotUMAP_H_Urothelial <- DimPlot(obj_H_Urothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_H_Urothelial
DefaultAssay(obj_H_Urothelial) <- "RNA"
obj_H_Urothelial <- NormalizeData(obj_H_Urothelial, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj_H_Urothelial)
obj_H_Urothelial <- ScaleData(obj_H_Urothelial, features = all.genes)
obj_H_Urothelial.markers <- FindAllMarkers(obj_H_Urothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_H_Urothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/H_Urothelial_markers.xlsx")
new.cluster.ids <- c("Intermediate", "Superficial", "Intermediate", "Basal", "Intermediate")
names(new.cluster.ids) <- levels(obj_H_Urothelial)
obj_H_Urothelial <- RenameIdents(obj_H_Urothelial, new.cluster.ids)
obj_H_Urothelial$CellAnn <- obj_H_Urothelial@active.ident
Idents(obj_H_Urothelial) <- "CellAnn"
obj_H_Urothelial.markers <- FindAllMarkers(obj_H_Urothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_H_Urothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/H_Urothelial_Annotated_markers.xlsx")
plotUMAP_H_Urothelial <- DimPlot(obj_H_Urothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_H_Urothelial
saveRDS(obj_H_Urothelial, "~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/H_Urothelial_Annotated_SeuratObject.rds")
gc()

#------------------Clustering and Annotation of D10A_Urothelial-----------------
obj_D10A_Urothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/D10A_Urothelial_SeuratObject.rds") %>% 
  CreateSeuratObject(assay = "RNA") %>% 
  SCTransform(assay = "RNA", variable.features.n = 5000, return.only.var.genes = FALSE) %>% 
  RunPCA(npcs = 50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(algorithm = 4) %>% 
  RunUMAP(dims = 1:50)
plotUMAP_D10A_Urothelial <- DimPlot(obj_D10A_Urothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_D10A_Urothelial
DefaultAssay(obj_D10A_Urothelial) <- "RNA"
obj_D10A_Urothelial <- NormalizeData(obj_D10A_Urothelial, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj_D10A_Urothelial)
obj_D10A_Urothelial <- ScaleData(obj_D10A_Urothelial, features = all.genes)
obj_D10A_Urothelial.markers <- FindAllMarkers(obj_D10A_Urothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_D10A_Urothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/D10A_Urothelial_markers.xlsx")
new.cluster.ids <- c("Intermediate", "Superficial", "Basal")
names(new.cluster.ids) <- levels(obj_D10A_Urothelial)
obj_D10A_Urothelial <- RenameIdents(obj_D10A_Urothelial, new.cluster.ids)
obj_D10A_Urothelial$CellAnn <- obj_D10A_Urothelial@active.ident
Idents(obj_D10A_Urothelial) <- "CellAnn"
obj_D10A_Urothelial.markers <- FindAllMarkers(obj_D10A_Urothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_D10A_Urothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/D10A_Urothelial_Annotated_markers.xlsx")
plotUMAP_D10A_Urothelial <- DimPlot(obj_D10A_Urothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_D10A_Urothelial
saveRDS(obj_D10A_Urothelial, "~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/D10A_Urothelial_Annotated_SeuratObject.rds")
gc()

#----------------Integration of H_Urothelial and D10A_Urothelial----------------
obj_H_Urothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/H_Urothelial_Annotated_SeuratObject.rds")
obj_D10A_Urothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/D10A_Urothelial_Annotated_SeuratObject.rds")
all.genes <- intersect(rownames(obj_H_Urothelial), rownames(obj_D10A_Urothelial))
list <- c(obj_H_Urothelial, obj_D10A_Urothelial)d
features <- SelectIntegrationFeatures(object.list = list)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors, features.to.integrate = all.genes)
orig.ident <- c(rep("H", dim(obj_H_Urothelial)[2]), rep("D10A", dim(obj_D10A_Urothelial)[2]))
orig.ident <- as.factor(orig.ident)
combined@meta.data[["orig.ident"]] <- orig.ident
combined@meta.data[["old.ident"]] <- orig.ident
obj <- combined
DefaultAssay(obj) <- "integrated"
remove(combined, list, obj_D10A_Urothelial, obj_H_Urothelial, anchors)
gc()

obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca")
plotUMAP_IntegrationPlot <- DimPlot(obj, reduction = "umap", group.by = "orig.ident",pt.size=0.5)
plotUMAP <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size=0.5)
plotUMAP
plotUMAP_IntegrationPlot
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj)
obj <- ScaleData(obj, features = all.genes)
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/Integrated_Urothelial_markers.xlsx")
Idents(obj) <- "orig.ident"
obj.markers <- FindAllMarkers(obj, on346ly.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,paste("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/02 - Integration","/Integrated_Urothelial_OrigIdent_markers.xlsx",sep=""))

new.cluster.ids <- c("Intermediate (H)", "Urothelial", "Intermediate (D10A)", "Basal")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj$CellAnnCoarse <- obj@active.ident
Idents(obj) -> "CellAnn"
plotUMAP_Annotated <- DimPlot(obj, reduction = "umap", pt.size=0.5)
plotUMAP_Annotated
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,paste("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/02 - Integration","/Combined_annotated_markers.xlsx",sep=""))

saveRDS(obj, "~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/Integrated_Urothelial_Annotated_SeuratObject.rds")

obj_Int_Urothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/Integrated_Urothelial_Annotated_SeuratObject.rds")
Idents(obj_Int_Urothelial) <- "orig.ident"
obj_Int_Urothelial.markers <- FindMarkers(obj_Int_Urothelial, ident.1 = "H", ident.2 = "D10A",assay = "RNA", slot = "data", logfc.threshold = 0, only.pos = FALSE, test.use = "wilcox", min.pct = 0, min.cells.feature = 0, min.cells.group = 0, random.seed = 5)
max.limit = ceiling(max(temp$avg_log2FC))
min.limit = floor(min(temp$avg_log2FC))
limit = max(c(abs(max.limit),abs(min.limit)))
temp$gene_symbol <- rownames(temp)
temp$diffExp <- "NO"
temp$diffExp[temp$avg_log2FC > 0.6 & temp$p_val_adj < 0.05] <- "UP"
temp$diffExp[temp$avg_log2FC < -0.6 & temp$p_val_adj < 0.05] <- "DOWN"
temp$delabel[temp$diffExp != "NO"] <- temp$gene_symbol[temp$diffExp != "NO"]
p <- ggplot(data=temp, aes(x=avg_log2FC, y=log10(1-log10(p_val_adj)), col=diffExp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=log10(1-log10(0.05)), col="red")
#--------------------------Intermediate cells comparison------------------------
Idents(obj) <- "seurat_clusters"
temp <- subset(obj, idents = c("0", "2"))
Idents(temp) <- "orig.ident"
temp <- FindMarkers(temp, ident.1 = "H", ident.2 = "D10A",assay = "RNA", slot = "data", logfc.threshold = 0, only.pos = FALSE, test.use = "wilcox", min.pct = 0, min.cells.feature = 0, min.cells.group = 0, random.seed = 5)
max.limit = ceiling(max(temp$avg_log2FC))
min.limit = floor(min(temp$avg_log2FC))
limit = max(c(abs(max.limit),abs(min.limit)))
temp$gene_symbol <- rownames(temp)
temp$diffExp <- "NO"
temp$diffExp[temp$avg_log2FC > 0.6 & temp$p_val_adj < 0.05] <- "UP"
temp$diffExp[temp$avg_log2FC < -0.6 & temp$p_val_adj < 0.05] <- "DOWN"
temp$delabel[temp$diffExp != "NO"] <- temp$gene_symbol[temp$diffExp != "NO"]
p <- ggplot(data=temp, aes(x=avg_log2FC, y=log10(1-log10(p_val_adj)), col=diffExp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=log10(1-log10(0.05)), col="red")

#---------------------------Superficial cells comparison------------------------
Idents(obj) <- "seurat_clusters"
temp <- subset(obj, idents = c("1"))
Idents(temp) <- "orig.ident"
temp <- FindMarkers(temp, ident.1 = "H", ident.2 = "D10A",assay = "RNA", slot = "data", logfc.threshold = 0, only.pos = FALSE, test.use = "wilcox", min.pct = 0, min.cells.feature = 0, min.cells.group = 0, random.seed = 5)
max.limit = ceiling(max(temp$avg_log2FC))
min.limit = floor(min(temp$avg_log2FC))
limit = max(c(abs(max.limit),abs(min.limit)))
temp$gene_symbol <- rownames(temp)
temp$diffExp <- "NO"
temp$diffExp[temp$avg_log2FC > 0.6 & temp$p_val_adj < 0.05] <- "UP"
temp$diffExp[temp$avg_log2FC < -0.6 & temp$p_val_adj < 0.05] <- "DOWN"
temp$delabel[temp$diffExp != "NO"] <- temp$gene_symbol[temp$diffExp != "NO"]
p <- ggplot(data=temp, aes(x=avg_log2FC, y=log10(1-log10(p_val_adj)), col=diffExp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=log10(1-log10(0.05)), col="red")

#------------------------------Basal cells comparison---------------------------
Idents(obj) <- "seurat_clusters"
temp <- subset(obj, idents = c("3"))
Idents(temp) <- "orig.ident"
temp <- FindMarkers(temp, ident.1 = "H", ident.2 = "D10A",assay = "RNA", slot = "data", logfc.threshold = 0, only.pos = FALSE, test.use = "wilcox", min.pct = 0, min.cells.feature = 0, min.cells.group = 0, random.seed = 5)
max.limit = ceiling(max(temp$avg_log2FC))
min.limit = floor(min(temp$avg_log2FC))
limit = max(c(abs(max.limit),abs(min.limit)))
temp$gene_symbol <- rownames(temp)
temp$diffExp <- "NO"
temp$diffExp[temp$avg_log2FC > 0.6 & temp$p_val_adj < 0.05] <- "UP"
temp$diffExp[temp$avg_log2FC < -0.6 & temp$p_val_adj < 0.05] <- "DOWN"
temp$delabel[temp$diffExp != "NO"] <- temp$gene_symbol[temp$diffExp != "NO"]
p <- ggplot(data=temp, aes(x=avg_log2FC, y=log10(1-log10(p_val_adj)), col=diffExp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=log10(1-log10(0.05)), col="red")

#-----------------------Common Urothelial subtype markers-----------------------
obj_H_Urothelial.markers_shortlisted <- obj_H_Urothelial.markers[obj_H_Urothelial.markers$p_val_adj < 0.001, ]
DotPlot(obj_D10A_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[1:30]))
DotPlot(obj_H_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[1:30]))
DotPlot(obj_H_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[31:50]))
DotPlot(obj_D10A_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[31:50]))
DotPlot(obj_H_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[31:50]))
DotPlot(obj_H_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[51:67]))
DotPlot(obj_D10A_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[51:67]))
DotPlot(obj_H_Urothelial, features = unique(obj_H_Urothelial.markers_shortlisted$gene[51:67]))

obj_D10A_Urothelial.markers_shortlisted <- obj_D10A_Urothelial.markers[obj_D10A_Urothelial.markers$p_val_adj < 0.053, ]
DotPlot(obj_D10A_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[1:30]))
DotPlot(obj_H_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[1:30]))
DotPlot(obj_H_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[31:50]))
DotPlot(obj_D10A_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[31:50]))
DotPlot(obj_H_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[31:50]))
DotPlot(obj_H_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[51:67]))
DotPlot(obj_D10A_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[51:67]))
DotPlot(obj_H_Urothelial, features = unique(obj_D10A_Urothelial.markers_shortlisted$gene[51:67]))

DotPlot(obj_H_Urothelial, features = c("Upk1b", "Jun", "Gata3", "Stox2", "Fosb", "Gm48099", "Dock8", "Fmo5", "Jund", "S100a6", "Chka", "Flnb", "Trp63", "Abhd2", "Itga2", "Lmo7", "Dsp", "Rbm39", "Hivep2", "Dhx40", "Parm1", "Nedd4l", "Krt23", "Ube2k", "Rock2", "Mir6236", "Plec", "Scarna2", "Atg4a", "Gpc6", "Akap12"))
DotPlot(obj_D10A_Urothelial, features = c("Upk1b", "Jun", "Gata3", "Stox2", "Fosb", "Gm48099", "Dock8", "Fmo5", "Jund", "S100a6", "Chka", "Flnb", "Trp63", "Abhd2", "Itga2", "Lmo7", "Dsp", "Rbm39", "Hivep2", "Dhx40", "Parm1", "Nedd4l", "Krt23", "Ube2k", "Rock2", "Mir6236", "Plec", "Scarna2", "Atg4a", "Gpc6", "Akap12"))

#----------------Create Monocle3 object from H Urothelial cells-----------------
obj <- obj_H_Urothelial
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

#----------------Compute Pseudotrajectory for H Urothelial cells----------------
cds <- learn_graph(cds)
plot_cells(cds, label_principal_points = TRUE)
Idents(obj) <- "CellAnn"
root_cells <- WhichCells(obj, ident = "Basal")
cds <- order_cells(cds, root_pr_nodes = "Y_42")
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.75, trajectory_graph_color = "green")

#--------------Create Monocle3 object from D10A Urothelial cells----------------
obj <- obj_D10A_Urothelial
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

#---------------Compute Pseudotrajectory for D10A Urothelial cells--------------
cds <- learn_graph(cds)
plot_cells(cds, label_principal_points = TRUE)
Idents(obj) <- "CellAnn"
root_cells <- WhichCells(obj, ident = "Basal")
cds <- order_cells(cds, root_pr_nodes = "Y_42")
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.75, trajectory_graph_color = "green")