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
library(monocle3)

#------------------Clustering and Annotation of H_Endothelial-------------------
obj_H_Endothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/H_Endothelial_SeuratObject.rds") %>% 
  CreateSeuratObject(assay = "RNA") %>% 
  SCTransform(assay = "RNA", variable.features.n = 5000, return.only.var.genes = FALSE) %>% 
  RunPCA(npcs = 50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(algorithm = 4) %>% 
  RunUMAP(dims = 1:50)
plotUMAP_H_Endothelial <- DimPlot(obj_H_Endothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_H_Endothelial
DefaultAssay(obj_H_Endothelial) <- "RNA"
obj_H_Endothelial <- NormalizeData(obj_H_Endothelial, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj_H_Endothelial)
obj_H_Endothelial <- ScaleData(obj_H_Endothelial, features = all.genes)
obj_H_Endothelial.markers <- FindAllMarkers(obj_H_Endothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_H_Endothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/H_Endothelial_markers.xlsx")
obj_H_Endothelial.markers_shortlisted <- obj_H_Endothelial.markers[obj_H_Endothelial.markers$pct.1 > 0.085 & obj_H_Endothelial.markers$pct.1/obj_H_Endothelial.markers$pct.2 > 2, ]
new.cluster.ids <- c("1", "2", "3", "4")
names(new.cluster.ids) <- levels(obj_H_Endothelial)
obj_H_Endothelial <- RenameIdents(obj_H_Endothelial, new.cluster.ids)
obj_H_Endothelial$CellAnn <- obj_H_Endothelial@active.ident
Idents(obj_H_Endothelial) <- "CellAnn"
obj_H_Endothelial.markers <- FindAllMarkers(obj_H_Endothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_H_Endothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/H_Endothelial_Annotated_markers.xlsx")
plotUMAP_H_Endothelial <- DimPlot(obj_H_Endothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_H_Endothelial
saveRDS(obj_H_Endothelial, "~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/H_Endothelial_Annotated_SeuratObject.rds")
gc()

#------------------Clustering and Annotation of D10A_Endothelial-----------------
obj_D10A_Endothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/D10A_Endothelial_SeuratObject.rds") %>% 
  CreateSeuratObject(assay = "RNA") %>% 
  SCTransform(assay = "RNA", variable.features.n = 5000, return.only.var.genes = FALSE) %>% 
  RunPCA(npcs = 50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(algorithm = 4) %>% 
  RunUMAP(dims = 1:50)
plotUMAP_D10A_Endothelial <- DimPlot(obj_D10A_Endothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_D10A_Endothelial
DefaultAssay(obj_D10A_Endothelial) <- "RNA"
obj_D10A_Endothelial <- NormalizeData(obj_D10A_Endothelial, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj_D10A_Endothelial)
obj_D10A_Endothelial <- ScaleData(obj_D10A_Endothelial, features = all.genes)
obj_D10A_Endothelial.markers <- FindAllMarkers(obj_D10A_Endothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_D10A_Endothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/D10A_Endothelial_markers.xlsx")
new.cluster.ids <- c("1")
names(new.cluster.ids) <- levels(obj_D10A_Endothelial)
obj_D10A_Endothelial <- RenameIdents(obj_D10A_Endothelial, new.cluster.ids)
obj_D10A_Endothelial$CellAnn <- obj_D10A_Endothelial@active.ident
Idents(obj_D10A_Endothelial) <- "CellAnn"
obj_D10A_Endothelial.markers <- FindAllMarkers(obj_D10A_Endothelial, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj_D10A_Endothelial.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/D10A_Endothelial_Annotated_markers.xlsx")
plotUMAP_D10A_Endothelial <- DimPlot(obj_D10A_Endothelial, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP_D10A_Endothelial
saveRDS(obj_D10A_Endothelial, "~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/D10A_Endothelial_Annotated_SeuratObject.rds")
gc()

#----------------Integration of H_Endothelial and D10A_Endothelial----------------
obj_H_Endothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/H_Endothelial_Annotated_SeuratObject.rds")
obj_D10A_Endothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/D10A_Endothelial_Annotated_SeuratObject.rds")
all.genes <- intersect(rownames(obj_H_Endothelial), rownames(obj_D10A_Endothelial))
list <- c(obj_H_Endothelial, obj_D10A_Endothelial)
features <- SelectIntegrationFeatures(object.list = list)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors, features.to.integrate = all.genes)
orig.ident <- c(rep("H", dim(obj_H_Endothelial)[2]), rep("D10A", dim(obj_D10A_Endothelial)[2]))
orig.ident <- as.factor(orig.ident)
combined@meta.data[["orig.ident"]] <- orig.ident
combined@meta.data[["old.ident"]] <- orig.ident
obj <- combined
DefaultAssay(obj) <- "integrated"
remove(combined, list, obj_D10A_Endothelial, obj_H_Endothelial, anchors)
gc()

obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca")
plotUMAP_IntegrationPlot <- DimPlot(obj, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size=0.5)
plotUMAP <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size=0.5)
plotUMAP
plotUMAP_IntegrationPlot
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj)
obj <- ScaleData(obj, features = all.genes)
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,"~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/Integrated_Endothelial_markers.xlsx")
Idents(obj) <- "orig.ident"
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,paste("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/02 - Integration","/Integrated_Endothelial_OrigIdent_markers.xlsx",sep=""))

group_temp_H <- WhichCells(obj_H_Endothelial, idents = "MF")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("S-MF")
group_temp_H <- WhichCells(obj_H_Endothelial, idents = "DMF")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("D-MF")
group_temp_H <- WhichCells(obj_H_Endothelial, idents = "MIF")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("MIF")
group_temp_H <- WhichCells(obj_H_Endothelial, idents = "MC")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("MC")

group_temp_H <- WhichCells(obj_D10A_Endothelial, idents = "MF1")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("S-MF")
group_temp_H <- WhichCells(obj_D10A_Endothelial, idents = "MF2")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("P-MF")
group_temp_H <- WhichCells(obj_D10A_Endothelial, idents = "MyF")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("MyF")
group_temp_H <- WhichCells(obj_D10A_Endothelial, idents = "InF")
DimPlot(obj, cells.highlight = c(group_temp_H), group.by = "orig.ident", cols.highlight = c("red"), cols = c("grey")) + NoLegend() + ggtitle("InF")

# Idents(obj) <- "seurat_clusters"
# new.cluster.ids <- c("MF", "Endothelial", "MIF", "Basal", "MC")
# names(new.cluster.ids) <- levels(obj)
# obj <- RenameIdents(obj, new.cluster.ids)
# obj$CellAnnCoarse <- obj@active.ident
# Idents(obj) -> "CellAnn"
# plotUMAP_Annotated <- DimPlot(obj, reduction = "umap", pt.size=0.5)
# plotUMAP_Annotated
# obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
# write_xlsx(obj.markers,paste("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/02 - Integration","/Combined_annotated_markers.xlsx",sep=""))

saveRDS(obj, "~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/Integrated_Endothelial_SeuratObject.rds")

obj_Int_Endothelial <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/Integrated_Endothelial_SeuratObject.rds")
Idents(obj_Int_Endothelial) <- "orig.ident"
obj_Int_Endothelial.markers <- FindMarkers(obj_Int_Endothelial, ident.1 = "H", ident.2 = "D10A",assay = "RNA", slot = "data", logfc.threshold = 0, only.pos = FALSE, test.use = "wilcox", min.pct = 0, min.cells.feature = 0, min.cells.group = 0, random.seed = 5)
temp <- obj_Int_Endothelial.markers
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

#-----------------------Common Endothelial subtype markers----------------------
Idents(obj_H_Endothelial) <- "CellAnn"

H_MC_geneset <- c("Msln", "Upk3b", "Gpm6a", "Igfbp5", "C3", "Plxna4", "Abi1", "Sox6", "Syne2", "Sntg1", "Spock2")
H_MF_geneset <- c("Col14a1", "Pi16", "Cygb", "Dcn", "Gsn", "Fbln1")
H_MIF_geneset <- c("Tnc", "Spon1", "Runx1", "Rad51b", "Kdm6b", "Nfkb1", "Nup210l", "Thsd4", "Rnf213", "Col6a3", "Col12a1")
H_DMF_geneset <- c("Jarid2", "Filip1l", "Il31ra", "Mir6236", "Hexb", "Scarna2", "Lars2", "Col1a1", "Col1a2", "Col3a1")
DotPlot(obj_H_Endothelial, features = unique(c(H_MF_geneset, H_MIF_geneset, H_DMF_geneset, H_MC_geneset)))

temp <- obj_H_Endothelial.markers_shortlisted[obj_H_Endothelial.markers_shortlisted$cluster == "MC",]
temp <- temp[temp$avg_log2FC > 1.5, ]
temp <- temp[temp$pct.1 > 0.2, ]
temp$gene
H_MC_geneset <- c("Col24a1", "Ncam1", "Tango6", "Il18r1", "Piezo2", "Taco1", "Cald1", "Plxna4", "Mical2", "C3", "Schip1", "Kalrn", "Gm20186", "Fbln2", "Col5a3", "Hmcn1", "Tmem108", "Actn1", "Ppp2r2b", "Thbs1", "Map1b", "Smad3", "Smad7", "Tnc", "Wt1", "Tpm1", "Flnb", "Fn1", "Myh9", "Aebp1")

temp <- obj_H_Endothelial.markers_shortlisted[obj_H_Endothelial.markers_shortlisted$cluster == "MF",]
temp <- temp[temp$avg_log2FC > 1, ]
temp$gene
H_MF_geneset <- c("Serpina3n", "Fbln1", "Pi16", "Ly6a", "Gsn", "Mgp", "Clec3b", "C1ra", "Sparcl1", "Aldh2", "Id3", "Rps12", "Rplp1", "Bgn", "C4b", "Dcn", "Cd34", "Plac9a", "Junb", "Egr1", "Timp2", "Ubc", "Scara5", "C1s1", "Mmp14", "Zfp36", "Eef2", "Fos", "Mt2", "Ccn1", "Actb", "Ier2")

temp <- obj_H_Endothelial.markers_shortlisted[obj_H_Endothelial.markers_shortlisted$cluster == "MIF",]
temp <- temp[temp$avg_log2FC > 1.5, ]
temp$gene
H_MIF_geneset <- c("Pid1", "Sgcz", "Adamtsl1", "Cntn4", "Frmpd4", "Arhgap6", "Frmd5", "Gm45837", "Fmnl2", "Gm3336", "Kcnq5", "9530026P05Rik", "Rad51b", "Myo1f", "Lsamp", "Lrrk1", "Lama2")

temp <- obj_H_Endothelial.markers_shortlisted[obj_H_Endothelial.markers_shortlisted$cluster == "DMF",]
temp <- temp[temp$avg_log2FC > 1, ]
temp$gene
H_DMF_geneset <- c("Ebf1", "Tnfaip6", "Man1a", "Marchf3", "Nfkb1", "Slc4a4", "Ebf2", "Creb5", "Rhoj", "Nav3", "Nova1", "Kcnq1ot1", "Sulf1", "Fbn2", "Mapk9"  )

DotPlot(obj_H_Endothelial, features = unique(c(H_MC_geneset, H_MF_geneset, H_MIF_geneset, H_DMF_geneset)))

Idents(obj_D10A_Endothelial) <- "CellAnn"

temp <- obj_D10A_Endothelial.markers_shortlisted[obj_D10A_Endothelial.markers_shortlisted$cluster == "MyF",]
temp <- temp[temp$avg_log2FC >1.47, ]
temp$gene
D10A_MyF_geneset <- c("Col24a1", "Ncam1", "Tango6", "Il18r1", "Piezo2", "Taco1", "Cald1", "Plxna4", "Mical2", "C3", "Schip1", "Kalrn", "Gm20186", "Fbln2", "Col5a3", "Hmcn1", "Tmem108", "Actn1", "Ppp2r2b", "Thbs1", "Map1b", "Smad3", "Smad7", "Tnc", "Wt1", "Tpm1", "Flnb", "Fn1", "Myh9", "Aebp1")

temp <- obj_D10A_Endothelial.markers_shortlisted[obj_D10A_Endothelial.markers_shortlisted$cluster == "MF1",]
temp <- temp[temp$avg_log2FC >1.5, ]
temp <- temp[temp$pct.1 > 0.2, ]
temp$gene
D10A_MF1_geneset <- c("Serpina3n", "Fbln1", "Pi16", "Ly6a", "Gsn", "Mgp", "Clec3b", "C1ra", "Sparcl1", "Aldh2", "Id3", "Rps12", "Rplp1", "Bgn", "C4b", "Dcn", "Cd34", "Plac9a", "Junb", "Egr1", "Timp2", "Ubc", "Scara5", "C1s1", "Mmp14", "Zfp36", "Eef2", "Fos", "Mt2", "Ccn1", "Actb", "Ier2")

temp <- obj_D10A_Endothelial.markers_shortlisted[obj_D10A_Endothelial.markers_shortlisted$cluster == "MF2",]
temp <- temp[temp$avg_log2FC >1.5, ]
temp$gene
D10A_MF2_geneset <- c("Pid1", "Sgcz", "Adamtsl1", "Cntn4", "Frmpd4", "Arhgap6", "Frmd5", "Gm45837", "Fmnl2", "Gm3336", "Kcnq5", "9530026P05Rik", "Rad51b", "Myo1f", "Lsamp", "Lrrk1", "Lama2")

temp <- obj_D10A_Endothelial.markers_shortlisted[obj_D10A_Endothelial.markers_shortlisted$cluster == "InF",]
temp <- temp[temp$avg_log2FC >1.5, ]
temp$gene
D10A_InF_geneset <- c("Ebf1", "Tnfaip6", "Man1a", "Marchf3", "Nfkb1", "Slc4a4", "Ebf2", "Creb5", "Rhoj", "Nav3", "Nova1", "Kcnq1ot1", "Sulf1", "Fbn2", "Mapk9"  )

DotPlot(obj_D10A_Endothelial, features = unique(c(D10A_MyF_geneset, D10A_MF1_geneset, D10A_MF2_geneset, D10A_InF_geneset)))

#----------------Create Monocle3 object from H Endothelial cells-----------------
obj <- obj_H_Endothelial
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

#----------------Compute Pseudotrajectory for H Endothelial cells----------------
cds <- learn_graph(cds)
plot_cells(cds, label_principal_points = TRUE)
Idents(obj) <- "CellAnn"
root_cells <- WhichCells(obj, ident = "Basal")
cds <- order_cells(cds, root_pr_nodes = "Y_65")
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.75, trajectory_graph_color = "green")

#--------------Create Monocle3 object from D10A Endothelial cells----------------
obj <- obj_D10A_Endothelial
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

#---------------Compute Pseudotrajectory for D10A Endothelial cells--------------
cds <- learn_graph(cds)
plot_cells(cds, label_principal_points = TRUE)
Idents(obj) <- "CellAnn"
root_cells <- WhichCells(obj, ident = "Basal")
cds <- order_cells(cds, root_pr_nodes = "Y_26")
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.75, trajectory_graph_color = "green")