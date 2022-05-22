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

#------------------------------Reading data for H and D10A----------------------
obj_H <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/01 - Clustering/H/H_SeuratObject.rds")
plotUMAP_H <- DimPlot(obj_H, reduction = "umap", pt.size = 0.5, label = TRUE)
plotUMAP_H
H_markers <- read_excel("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/01 - Clustering/H/H_markers.xlsx")
obj_D10A <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/01 - Clustering/D10A/D10A_SeuratObject.rds")
plotUMAP_D10A <- DimPlot(obj_D10A, reduction = "umap", pt.size = 0.5, label = TRUE)
plotUMAP_D10A
D10A_markers <- read_excel("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/01 - Clustering/D10A/D10A_markers.xlsx")

#-----------------Finding common markers for H and D10A-------------------------
H_markers_temp <- H_markers[H_markers$p_val_adj < 0.05 & (H_markers$pct.1 > 0.20), ]
H_markers_temp <- H_markers_temp[((H_markers_temp$pct.1)/(H_markers_temp$pct.2) > 5), ]
H_diff_genes <- unique(H_markers_temp$gene)
D10A_markers_temp <- D10A_markers[(D10A_markers$p_val_adj < 0.05) & (D10A_markers$pct.1 > 0.20), ]
D10A_markers_temp <- D10A_markers_temp[((D10A_markers_temp$pct.1)/(D10A_markers_temp$pct.2) > 5), ]
D10A_diff_genes <- unique(D10A_markers_temp$gene)
diff_genes <- intersect(H_diff_genes,D10A_diff_genes)
H_markers_common <- H_markers_temp[H_markers_temp$gene %in% diff_genes, ]
D10A_markers_common <- D10A_markers_temp[D10A_markers_temp$gene %in% diff_genes, ]

#-----------Eyeball the featureplot to decide the final markers in H------------
FeaturePlot(obj_H, features = H_markers_common$gene[1:10])
# Mesenchymal markers
# "Fbln1", "Pi16", "Cygb"
# Endothelial markers
# "Flt1", "Cyyr1", "Adamts9", "Ptprb", "Cdh5", "Adgrf5"
FeaturePlot(obj_H, features = H_markers_common$gene[10:20])
# Endothelial markers
# "Erg", "Pecam1", "Iigp1", "Ptprm", "Ldb2", "Ebf1", "Prkch", "Cdh13", "Fli1"
# "Iigp1", "Adamts9" is every blood vessel cells
# "Erg". "Adgrf5", "Flt1", "Ptprb", "Cdh5" is strictly endothelial cells
FeaturePlot(obj_H, features = H_markers_common$gene[20:30])
# Endothelial markers
# "Schip1"
# Urothelial markers
# "Dsp", "Pkhd1", "Rab27b", "Las1l"
FeaturePlot(obj_H, features = H_markers_common$gene[30:40])
# Mesenchymal markers
# "C3"
# VSMC, Myofib markers
# "Cald1", "Sorbs1"
# Urothelial markers
# "Parm1", "Cdh1", "9530036O11Rik", "Atp10b"
FeaturePlot(obj_H, features = H_markers_common$gene[40:50])
# Urothelial markers
# "Lmo7", "Tmprss2", "Ctnnd2", "Wwc1"
FeaturePlot(obj_H, features = H_markers_common$gene[50:56])
# Mesenchymal/Endthelial markers
# "Pkrg1", "Gm20532", "Arhgap31"
# Immune cells
# "Dock2"
# Neural cells
# Ncam1
FeaturePlot(obj_H, features = c("Baiap2la", "Rbm47", "Dnmpbp", "Dhx40", "Acacb", "Shank3", "Meg3", "Dcn"))
# Urothelial markers
# "Rbm47", "Dhx40"
# Immune markers
# "Acacb"
# Endothelial markers
# "Shank3"
# Mesenchymal markers
# "Meg3"
FeaturePlot(obj_H, features = c("Syne", "Fbn1", "Tjp2", "Aqp3", "Gata3", "Myh11", "Tagln", "Acta2", "Mylk", "Smtn"))
# Urothelial markers
# "Gata3"
# SMC markers
# "Myh11", "Tagln", "Mylk"
FeaturePlot(obj_H, features = c("Des", "Csmd1", "Dok6", "Cdh19", "Kcna6", "Cadm2", "Ncam1", "Cd74", "H2-Ab1", "H2-Eb1", "Il1b"))
# Neuronal
# "Csmd1"
FeaturePlot(obj_H, features = c("Col1a1", "Col1a2", "Col3a1", "Flt1", "Ptprb", "Ptprc", "Cd84", "Il1b", "Ccr5"))
# Mesenchymal markers
# "Col1a1", "Col1a2", "Col3a1"
# Endothelial markers
# "Flt1", "Ptprb"
# Immune markers
# "Ptprc", "Ccr5", "Il1b"

new.cluster.ids <- c("Mesenchymal", "Urothelial", "Mesenchymal", "Endothelial", "Urothelial", "Urothelial", "Mesenchymal", "SMC", "Urothelial", "Urothelial", "VSMC", "Immune (TR)", "Neurons")
names(new.cluster.ids) <- levels(obj_H)
obj_H <- RenameIdents(obj_H, new.cluster.ids)
obj_H$CellAnnCoarse <- obj_H@active.ident
DimPlot(obj_H, reduction = "umap", pt.size = 0.5)
Idents(obj_H) <- "CellAnnCoarse"

#---------Eyeball the featureplot to decide the final markers in D10A-----------
FeaturePlot(obj_D10A, features = D10A_markers_common$gene[1:10])
# Urothelial markers
# "Tmprss2", "Ctnnd2", "Dsp", "Cdh1", "Parm1", "Rab27b", "Lmo7", "Atp10b", "9530036O11Rik", "Wwc1"
FeaturePlot(obj_D10A, features = H_markers_common$gene[10:20])
# Endothelial markers
# "Pecam1", "Ptprm", "Ldb2", "Ebf1", "Prkch", "Cdh13", "Fli1"
FeaturePlot(obj_D10A, features = H_markers_common$gene[20:30])
# Urothelial markers
# "Schip1", "Ablim1", "Dsp", "Pkhd1", "Rab27b", "Las1l"
FeaturePlot(obj_D10A, features = H_markers_common$gene[30:40])
# Mesenchymal markers
# "C3", "Plxna4"
# VSMC, Myofib markers
# "Cald1"
# Urothelial markers
# "Pkhd1", "9530036O11Rik", "Parm1", "Cdh1", "Ikzf2", "Atp10b"
FeaturePlot(obj_D10A, features = H_markers_common$gene[40:49])
# Urothelial markers
# "Rab27b","AY036118", "Lmo7", "Tmprss2", "Ctnnd2", "Wwc1"
FeaturePlot(obj_D10A, features = c("Baiap2la", "Rbm47", "Dnmpbp", "Dhx40", "Acacb", "Shank3", "Meg3", "Dcn"))
# Urothelial markers
# "Rbm47", "Dhx40"
# Endothelial markers
# "Shank3"
# Mesenchymal markers
# "Meg3"
FeaturePlot(obj_D10A, features = c("Syne", "Fbn1", "Tjp2", "Aqp3", "Gata3", "Myh11", "Tagln", "Acta2", "Mylk", "Smtn"))
# Urothelial markers
# "Gata3", "Tjp2"
# SMC markers
# "Myh11", "Tagln", "Mylk"
# Mesenchymal markers
# "Fbn1"
FeaturePlot(obj_D10A, features = c("Des", "Csmd1", "Dok6", "Cdh19", "Kcna6", "Cadm2", "Ncam1", "Cd74", "H2-Ab1", "H2-Eb1", "Il1b"))
# Neuronal
# "Ncam1"
FeaturePlot(obj_D10A, features = c("Col1a1", "Col1a2", "Col3a1", "Flt1", "Ptprb", "Ptprc", "Cd84", "Il1b", "Ccr5"))
# Mesenchymal markers
# "Col1a1", "Col1a2", "Col3a1"
# Endothelial markers
# "Flt1", "Ptprb"
# Immune markers
# "Ptprc", "Cd84", "Ccr5", "Il1b"

new.cluster.ids <- c("Urothelial", "Mesenchymal", "Mesenchymal", "Immune", "Urothelial", "Urothelial", "Endothelial", "Mesenchymal")
names(new.cluster.ids) <- levels(obj_D10A)
obj_D10A <- RenameIdents(obj_D10A, new.cluster.ids)
obj_D10A$CellAnnCoarse <- obj_D10A@active.ident
DimPlot(obj_D10A, reduction = "umap", pt.size = 0.5)
Idents(obj_D10A) <- "CellAnnCoarse"

H_scale_data <- t(obj_H@assays$H@scale.data)
H_scale_data <- as.data.frame((H_scale_data))
colnames(H_scale_data) <- obj_H@assays[["H"]]@data@Dimnames[[1]]
rownames(H_scale_data) <- obj_H@assays[["H"]]@data@Dimnames[[2]]
H_scale_data <- H_scale_data[, H_pearson_genes_short]
H_cor_pearson <- cor(H_scale_data)
heatmap(H_cor_pearson, scale = "none", col = col)
pheatmap(H_cor_pearson, cluster_rows = FALSE, cluster_cols = FALSE)

temp_scale_data <- t(obj_D10A@assays$RNA@counts)
temp_scale_data <- as.data.frame(temp_scale_data)
colnames(temp_scale_data) <- temp@assays$RNA@data@Dimnames[[1]]
rownames(temp_scale_data) <- temp@assays$RNA@data@Dimnames[[2]]
temp_scale_data <- temp_scale_data[, cc_features$G2M.Score]
temp_cor_pearson <- cor(temp_scale_data)
temp_cor_pearson[is.na(temp_cor_pearson)] <- 0
cor_rowSum <- temp_cor_pearson^2
cor_rowSum <- Matrix::rowSums(temp_cor_pearson)
hist(cor_rowSum, breaks = 100)
pheatmap(temp_cor_pearson, cluster_rows = FALSE, cluster_cols = FALSE)


D10A_scale_data <- t(obj_D10A@assays$S1@scale.data)
D10A_scale_data <- as.data.frame((D10A_scale_data))
colnames(D10A_scale_data) <- obj_D10A@assays[["S1"]]@data@Dimnames[[1]]
rownames(D10A_scale_data) <- obj_D10A@assays[["S1"]]@data@Dimnames[[2]]
D10A_scale_data <- D10A_scale_data[, D10A_pearson_genes_short]
D10A_cor_pearson <- cor(D10A_scale_data)
heatmap(D10A_cor_pearson, scale = "none", col = col)
pheatmap(D10A_cor_pearson, cluster_rows = FALSE, cluster_cols = FALSE)

# Pan urothelial cell markers
FeaturePlot(obj_D10A, features = c("Pard3", "Cdh1", "Tjp2", "Gata3", "Dsp", "Aqp3"))

# Superficial urothelial cell markers
FeaturePlot(obj_D10A, features = c("Patj", "Cldn4", "Rab11fip1", "Ocln", "Sh3", "Itch", "Gata3", "Trpv4", "Aqp3", "Sytl2"))
# "Sytl2", "Myo5b", "Arhgef26", "Rhoa", "Snx31", "Dnm2" supports uroplakin on the surface of umbrella cells
# Uroplakin related genes are obtained from https://journals.physiology.org/doi/full/10.1152/physrev.00041.2019
# Nedd4l overexpression in 5 might have increased Na+ ion channel conductance
# Presence of "Snap25", "Stx3", "Cep290", "Rab11fip1", "Actn4" in 5 shows formation of cilia
# Cilia related genes obtained from https://cshperspectives.cshlp.org/content/10/1/a027813.full

# Intermediate urothelial cell markers
FeaturePlot(obj_D10A, features = c("Cgn", "Gata3"))

# Basal cell markers
FeaturePlot(obj_D10A, features = c("Ezr", "Itgb4", "Asic"))

# 8
"Ybx3"

# Fibrocytes = Gpc6, Bgn, Ncam1, Vcan
# Fibroblasts = Podc, Lum, Bgn, Col12a1, Col15a1
# Piezo2 upregulation in 2 might be a reason for activation of fibroblast

# Pan mesenchymal cell markers
FeaturePlot(obj_D10A, features = c("Col1a1", "Col3a1", "Syne1", "Fbn1"))

# Pan endothelial cell markers
FeaturePlot(obj_D10A, features = c("Pecam1", "Ptprm", "Flt1", "Adamts9", "Cyyr1", "Fbxl7", "Col4a1"))

# p21 (Cdkn1a), B2m, Cd36, Cd97
# Cd9, Cd93
# Cd80, Cd86

# Pan Immune cell markers
FeaturePlot(obj_D10A, features = c("Ptprj", "Zeb2", "Cd44", "Il1b"))

#Pan Tissue resident immune cell markers
FeaturePlot(obj_H, features = c("H2-Ab1", "H2-Eb1", "Cd74"))

#Pan Circulatory immune cell markers
FeaturePlot(obj_H, features = c("Ptprc", "Cd84", "Ccr5", "Cd86", "Il1b"))

# Pan SMC markers
FeaturePlot(obj_H, features = c("Myh11", "Tagln", "Cald1", "Acta2", "Des", "Mylk", "Smtn"))

#Ptprc, 