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
library(scater)

base <- "/home/shiva/Desktop/Mouse"
old <- "/Step 1 - Convert ParseBio format to 10X format"
new <- "/Step 2 - Analysis"
sample <- "/D10A"
obj.data <- Read10X(data.dir = paste(base,old,sample,sep=""))
obj <- CreateSeuratObject(counts = obj.data,assay = "RNA")
remove(obj.data)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-") + PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")+PercentageFeatureSet(obj, pattern = "^Rn18s")+PercentageFeatureSet(obj, pattern = "^Rna[1234567890]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^(p)]")
obj[["log10GenesPerUMI"]] <- log10(obj[["nFeature_RNA"]])/log10(obj[["nCount_RNA"]])
selected_c <- WhichCells(obj, expression = nFeature_RNA > 100)
selected_f <- rownames(obj)[Matrix::rowSums(obj) > 3]
obj <- subset(obj, features = selected_f, cells = selected_c, subset = percent.mt < 10)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes <- str_to_title((cc.genes.updated.2019$g2m.genes))
obj <- CellCycleScoring(object = obj, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(obj)
top50 <- head(VariableFeatures(obj), 50)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
plot2
obj_sce <- as.SingleCellExperiment(obj)
pltHighExpr <- plotHighestExprs(obj_sce, exprs_values = "counts", colour_cells_by = "ident")
gc()
obj_D10A <- obj
remove(obj)

base <- "/home/shiva/Desktop/Mouse"
old <- "/Step 1 - Convert ParseBio format to 10X format"
new <- "/Step 2 - Analysis"
sample <- "/H"
obj.data <- Read10X(data.dir = paste(base,old,sample,sep=""))
obj <- CreateSeuratObject(counts = obj.data,assay = "RNA")
remove(obj.data)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-") + PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")+PercentageFeatureSet(obj, pattern = "^Rn18s")+PercentageFeatureSet(obj, pattern = "^Rna[1234567890]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^(p)]")
obj[["log10GenesPerUMI"]] <- log10(obj[["nFeature_RNA"]])/log10(obj[["nCount_RNA"]])
selected_c <- WhichCells(obj, expression = nFeature_RNA > 100)
selected_f <- rownames(obj)[Matrix::rowSums(obj) > 3]
obj <- subset(obj, features = selected_f, cells = selected_c, subset = percent.mt < 10)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes <- str_to_title((cc.genes.updated.2019$g2m.genes))
obj <- CellCycleScoring(object = obj, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(obj)
top50 <- head(VariableFeatures(obj), 50)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
plot2
obj_sce <- as.SingleCellExperiment(obj)
pltHighExpr <- plotHighestExprs(obj_sce, exprs_values = "counts", colour_cells_by = "ident")

gc()
obj_H <- obj
remove(obj)

original_ident <- rep(c("H", "D10A"), times = c(dim(obj_H)[2], dim(obj_D10A)[2]))
obj <- merge(obj_H, y = obj_D10A)
obj@meta.data[["orig.ident"]] <- original_ident
Idents(obj) <- "orig.ident"
#remove(obj_D10A, obj_H, original_ident, base, g2m.genes, new, old, s.genes, sample, selected_c, selected_f)
saveRDS(obj, "~/Desktop/Mouse/Step 2 - Analysis/D10A vs H/QC/Combined_SeuratObject.rds")
gc()

obj@meta.data[["nCount_RNA"]] <- log10(obj@meta.data[["nCount_RNA"]])
obj@meta.data[["nFeature_RNA"]] <- log10(obj@meta.data[["nFeature_RNA"]])
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "nCount_RNA")
VlnPlot(obj, features = "percent.mt")
VlnPlot(obj, features = "percent.ribo")
VlnPlot(obj, features = "percent.hb")
VlnPlot(obj, features = c("S.Score", "G2M.Score"))
obj@meta.data[["nCount_RNA"]] <- 10^(obj@meta.data[["nCount_RNA"]])
obj@meta.data[["nFeature_RNA"]] <- 10^(obj@meta.data[["nFeature_RNA"]])
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.ribo")

obj_sce <- as.SingleCellExperiment(obj)
pltHighExpr <- plotHighestExprs(obj_sce, exprs_values = "counts", colour_cells_by = "ident")

cells_D10A <- colnames(subset(x = obj,idents = "D10A"))
cells_H <- colnames(subset(x = obj,idents = "H"))
volcano_markers <- FindMarkers(object = obj,ident.1 = cells_D10A,ident.2 = cells_H,assay = "RNA", slot = "data", logfc.threshold = 0, only.pos = FALSE, test.use = "wilcox", min.pct = 0, min.cells.feature = 0, min.cells.group = 0, random.seed = 5)
max.limit = ceiling(max(volcano_markers$avg_log2FC))
min.limit = floor(min(volcano_markers$avg_log2FC))
limit = max(c(abs(max.limit),abs(min.limit)))
volcano_markers$gene_symbol <- rownames(volcano_markers)
volcano_markers$diffExp <- "NO"
volcano_markers$diffExp[volcano_markers$avg_log2FC > 0.6 & volcano_markers$p_val_adj < 0.05] <- "UP"
volcano_markers$diffExp[volcano_markers$avg_log2FC < -0.6 & volcano_markers$p_val_adj < 0.05] <- "DOWN"
volcano_markers$delabel[volcano_markers$diffExp != "NO"] <- volcano_markers$gene_symbol[volcano_markers$diffExp != "NO"]
p <- ggplot(data=volcano_markers, aes(x=avg_log2FC, y=log10(1 - log10(p_val_adj)), col=diffExp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=log10(1 - log10(0.05)), col="red")

volc <- volc[(volc$diffExp == "UP" | volc$diffExp == "DOWN") & (volc$pct.1 + volc$pct.2 > 0.08), ]
varF <- c(VariableFeatures(obj_H), VariableFeatures(obj_D10A))
finvarF <- intersect(varF, volc$gene_symbol)
volc <- volc[volc$gene_symbol %in% finvarF, ]


percentage_D10A <- data.frame(Cell_Types = c("Urothelial (S)", "Urothelial (I)", "Urothelial (B)", "Urothelial (P)", "Fibrocytes", "Fibroblasts", "Myofibroblasts", "Endothelial", "Immune (C)"), Count = c(172, 206, 504, 37, 439, 74, 307, 151, 290))
percentage_D10A$Percentage <- percentage_D10A$Count/sum(percentage_D10A$Count)*100
percentage_D10A$ymax <- cumsum(percentage_D10A$Percentage)
percentage_D10A$ymin <- c(0, head(percentage_D10A$ymax, n = -1))
ggplot(percentage_D10A, aes(x = "", y = Percentage, fill = Cell_Types)) + geom_bar(stat = "identity", width = 1, color = "white") + coord_polar(theta = "y", start = 0) + theme_void() + scale_fill_brewer(palette = "Set1")
percentage_H <- data.frame(Cell_Types = c("Urothelial (S)", "Urothelial (I)", "Urothelial (B)", "Fibroblasts", "Endothelial", "SMC", "VSMC", "Neurons", "Immune (C)"), Percentage = c(1007, 705, 357, 2949, 593, 218, 132, 108, 41))
percentage_H$Percentage <- percentage_H$Count/sum(percentage_H$Count)*100
percentage_H$ymax <- cumsum(percentage_H$Percentage)
percentage_H$ymin <- c(0, head(percentage_H$ymax, n = -1))
ggplot(percentage_H, aes(x = "", y = Percentage, fill = Cell_Types)) + geom_bar(stat = "identity", width = 1, color = "white") + coord_polar(theta = "y", start = 0) + theme_void() + scale_fill_brewer(palette = "Set1")
