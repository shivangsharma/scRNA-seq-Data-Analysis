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

base <- "~/Desktop/Mouse"
old <- "/Step 1 - Convert ParseBio format to 10X format"
new <- "/Step 2 - Analysis (Manuscript)/01 - Clustering"
sample <- "/D10A"
if(!dir.exists(paste(base,new,sample,sep="")))
{
  dir.create(paste(base,new,sample,sep=""))
}
obj.data <- Read10X(data.dir = paste(base,old,sample,sep=""))
obj <- CreateSeuratObject(counts = obj.data,assay = "RNA")

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-") + PercentageFeatureSet(obj, pattern = "^Mt-") + PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")+PercentageFeatureSet(obj, pattern = "^RN18S")+PercentageFeatureSet(obj, pattern = "^RNA[1234567890]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]")
obj[["log10GenesPerUMI"]] <- log10(obj[["nFeature_RNA"]])/log10(obj[["nCount_RNA"]])

selected_c <- WhichCells(obj, expression = nFeature_RNA > 100)
selected_f <- rownames(obj)[Matrix::rowSums(obj) > 3]
obj <- subset(obj, features = selected_f, cells = selected_c, subset = percent.mt < 10)
obj <- SCTransform(obj, assay = "RNA", variable.features.n = 5000, return.only.var.genes = FALSE)
obj <- RunPCA(obj, npcs = 50)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, algorithm = 4)
obj <- RunUMAP(obj, dims = 1:50)
plotUMAP <- DimPlot(obj, reduction = "umap",label = TRUE,pt.size=0.5)
plotUMAP
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj)
obj <- ScaleData(obj, features = all.genes)
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,paste(base,new,sample,"/D10A_markers.xlsx",sep=""))
gc()

new.cluster.ids <- c("Urothelial", "Mesenchymal", "Mesenchymal", "Immune", "EMT", "Urothelial", "Endothelial", "Mesenchymal")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj$CellAnnCoarse <- obj@active.ident
Idents(obj) <- "CellAnnCoarse"
plotUMAP <- DimPlot(obj, reduction = "umap", pt.size=0.5)
plotUMAP
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(x = obj)
obj <- ScaleData(obj, features = all.genes)
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
write_xlsx(obj.markers,paste(base,new,sample,"/D10A_annotated_markers.xlsx",sep=""))
obj_D10A@meta.data$SampleID <- c(rep("D10A", times = dim(obj_D10A)[2]))

saveRDS(obj_D10A,paste(base,new,sample,"/D10A_SeuratObject.rds",sep=""))


#--------------------------Revisiting the data----------------------------------
# base <- "~/Desktop/Mouse"
# new <- "/Step 2 - Analysis (Manuscript)/01 - Clustering"
# sample <- "/D10A"
# obj <- readRDS(paste(base,new,sample,"/D10A_SeuratObject.rds",sep=""))
# plotUMAP <- DimPlot(obj, reduction = "umap",label = TRUE,pt.size=0.5)
# plotUMAP
# obj.markers <- read_excel(paste(base,new,sample,"/D10A_markers.xlsx",sep=""))
# View(obj.markers)