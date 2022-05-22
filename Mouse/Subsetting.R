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
obj_H <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/03 - Annotation/H_SeuratObject.rds")
plotUMAP_H <- DimPlot(obj_H, reduction = "umap", pt.size = 0.5, label = TRUE)
plotUMAP_H
H_markers <- read_excel("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/03 - Annotation/H_annotated_markers.xlsx")
obj_D10A <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/03 - Annotation/D10A_SeuratObject.rds")
plotUMAP_D10A <- DimPlot(obj_D10A, reduction = "umap", pt.size = 0.5, label = TRUE)
plotUMAP_D10A
D10A_markers <- read_excel("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/03 - Annotation/D10A_annotated_markers.xlsx")

#------------------------------Subsetting Urothelial cells----------------------
GetAssayData(obj_H, assay = "RNA", slot = "counts")[, WhichCells(obj_H, idents = "Urothelial")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/H_Urothelial_SeuratObject.rds")
GetAssayData(obj_D10A, assay = "RNA", slot = "counts")[, WhichCells(obj_D10A, idents = "Urothelial")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Urothelial/D10A_Urothelial_SeuratObject.rds")

#------------------------------Subsetting Mesenchymal cells---------------------
GetAssayData(obj_H, assay = "RNA", slot = "counts")[, WhichCells(obj_H, idents = "Mesenchymal")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Mesenchymal/H_Mesenchymal_SeuratObject.rds")
GetAssayData(obj_D10A, assay = "RNA", slot = "counts")[, WhichCells(obj_D10A, idents = "Mesenchymal")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Mesenchymal/D10A_Mesenchymal_SeuratObject.rds")

#------------------------------Subsetting Endothelial cells---------------------
GetAssayData(obj_H, assay = "RNA", slot = "counts")[, WhichCells(obj_H, idents = "Endothelial")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/H_Endothelial_SeuratObject.rds")
GetAssayData(obj_D10A, assay = "RNA", slot = "counts")[, WhichCells(obj_D10A, idents = "Endothelial")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Endothelial/D10A_Endothelial_SeuratObject.rds")

#------------------------------Subsetting Immune cells---------------------
GetAssayData(obj_H, assay = "RNA", slot = "counts")[, WhichCells(obj_H, idents = "Immune")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Immune/H_Immune_SeuratObject.rds")
GetAssayData(obj_D10A, assay = "RNA", slot = "counts")[, WhichCells(obj_D10A, idents = "Immune")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/Immune/D10A_Immune_SeuratObject.rds")

#------------------------------Subsetting EMT cells---------------------
GetAssayData(obj_H, assay = "RNA", slot = "counts")[, WhichCells(obj_H, idents = "EMT")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/EMT/H_EMT_SeuratObject.rds")
GetAssayData(obj_D10A, assay = "RNA", slot = "counts")[, WhichCells(obj_D10A, idents = "EMT")] %>% 
  saveRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/05 - Subclustering/EMT/D10A_EMT_SeuratObject.rds")
