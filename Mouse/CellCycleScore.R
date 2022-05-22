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
library(Revelio)

obj_D10A <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/03 - Annotation/D10A_SeuratObject.rds")
obj_D10A <- readRDS("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/03 - Annotation/H_SeuratObject.rds")
cc_genes <- read_excel("~/Desktop/Mouse/Step 2 - Analysis (Manuscript)/04 - CellCycle/Mouse_CellCycle_geneset, Whitfield, 2002.xlsx")
cc_features <- list(G1S = cc_genes$G1.S[cc_genes$G1.S %in% rownames(obj_D10A)], 
                    S = cc_genes$S[cc_genes$S %in% rownames(obj_D10A)], 
                    G2 = cc_genes$G2[cc_genes$G2 %in% rownames(obj_D10A)], 
                    G2M = cc_genes$G2.M[cc_genes$G2.M %in% rownames(obj_D10A)], 
                    MG1 = cc_genes$M.G1[cc_genes$M.G1 %in% rownames(obj_D10A)]
                    )

D10A_counts <- t(obj_D10A@assays$RNA@counts)
D10A_counts <- as.data.frame(D10A_counts)
colnames(D10A_counts) <- obj_D10A@assays$RNA@data@Dimnames[[1]]
rownames(D10A_counts) <- obj_D10A@assays$RNA@data@Dimnames[[2]]
temp_scale_data <- list(G1S = D10A_counts[, cc_features$G1S], 
                        S = D10A_counts[, cc_features$S], 
                        G2 = D10A_counts[, cc_features$G2], 
                        G2M = D10A_counts[, cc_features$G2M], 
                        MG1 = D10A_counts[, cc_features$MG1]
                        )
avg_temp_scale_data <- list(G1S = Matrix::rowSums(temp_scale_data$G1S), 
                            S = Matrix::rowSums(temp_scale_data$S), 
                            G2 = Matrix::rowSums(temp_scale_data$G2), 
                            G2M = Matrix::rowSums(temp_scale_data$G2M),
                            MG1 = Matrix::rowSums(temp_scale_data$MG1)
                            )
cc_features_cor <- list(G1S = data.frame(R = rep(NA, times = length(cc_features$G1S)), row.names = cc_features$G1S), 
                        S = data.frame(R = rep(NA, times = length(cc_features$S)), row.names = cc_features$S), 
                        G2 = data.frame(R = rep(NA, times = length(cc_features$G2)), row.names = cc_features$G2), 
                        G2M = data.frame(R = rep(NA, times = length(cc_features$G2M)), row.names = cc_features$G2M), 
                        MG1 = data.frame(R = rep(NA, times = length(cc_features$MG1)), row.names = cc_features$MG1)
                        )
for (i in c("G1S", "S", "G2", "G2M", "MG1")) {
  genes <- cc_features[[i]]
  for (j in genes) {
    cc_features_cor[[i]][j, ] <- cor(avg_temp_scale_data[[i]],temp_scale_data[[i]][[j]], method = "pearson")
  }
  cc_features_cor[[i]] <- arrange(cc_features_cor[[i]], desc(R))
  cc_features_cor[[i]] <- cc_features_cor[[i]] %>% 
    rownames_to_column("gene_names") %>% 
    filter(R > 0.4) %>% 
    column_to_rownames("gene_names")
}
new_cc_features <- list(G1S = rownames(cc_features_cor$G1S), 
                        S = rownames(cc_features_cor$S), 
                        G2 = rownames(cc_features_cor$G2), 
                        G2M = rownames(cc_features_cor$G2M), 
                        MG1 = rownames(cc_features_cor$MG1)
                        )
new_cc_genes_max_l <- max(sapply(new_cc_features, length))
new_cc_genes <- data.frame(G1.S = c(new_cc_features$G1S, rep(NA, times = new_cc_genes_max_l - length(new_cc_features$G1S))), 
                 S = c(new_cc_features$S, rep(NA, times = new_cc_genes_max_l - length(new_cc_features$S))), 
                 G2 = c(new_cc_features$G2, rep(NA, times = new_cc_genes_max_l - length(new_cc_features$G2))), 
                 G2.M = c(new_cc_features$G2M, rep(NA, times = new_cc_genes_max_l - length(new_cc_features$G2M))), 
                 M.G1 = c(new_cc_features$MG1, rep(NA, times = new_cc_genes_max_l - length(new_cc_features$MG1)))
                 )
new_cc_genes[nrow(new_cc_genes) + 1, ] <- rep(NA, times = ncol(new_cc_genes))
obj <- obj_D10A
obj_revelio <- createRevelioObject(rawData = obj@assays$RNA@counts, cyclicGenes = cc_genes, lowernGeneCutoff = 0) #cell filtering for nUMI > 1200 removed in source code
obj_revelio <- getCellCyclePhaseAssignInformation(dataList = obj_revelio) #cell filtering for outliers removed in source code (Line 67 in getCellCyclePhaseAssignInformation() using 'trace(getCellCyclePhaseAssignInformation, edit = TRUE)' command )
obj@meta.data[["S.Score.Revelio"]] <- obj_revelio@cellInfo[["S_zScore"]]
obj@meta.data[["G2M.Score.Revelio"]] <- obj_revelio@cellInfo[["G2.M_zScore"]]
obj@meta.data[["G1S.Score.Revelio"]] <- obj_revelio@cellInfo[["G1.S_zScore"]]
obj@meta.data[["G2.Score.Revelio"]] <- obj_revelio@cellInfo[["G2_zScore"]]
obj@meta.data[["MG1.Score.Revelio"]] <- obj_revelio@cellInfo[["M.G1_zScore"]]
obj@meta.data[["Phase.Revelio"]] <- as.character(obj_revelio@cellInfo[["ccPhase"]])
obj_D10A <- obj
FeatureScatter(obj_D10A, feature1 = "G1S.Score.Revelio", feature2 = "G2M.Score.Revelio", group.by = "CellAnnCoarse")

obj <- obj_H
obj_revelio <- createRevelioObject(rawData = obj@assays$RNA@counts, cyclicGenes = cc_genes, lowernGeneCutoff = 0) #cell filtering for nUMI > 1200 removed in source code
obj_revelio <- getCellCyclePhaseAssignInformation(dataList = obj_revelio) #cell filtering for outliers removed in source code (Line 67 in getCellCyclePhaseAssignInformation() using 'trace(getCellCyclePhaseAssignInformation, edit = TRUE)' command )
obj@meta.data[["S.Score.Revelio"]] <- obj_revelio@cellInfo[["S_zScore"]]
obj@meta.data[["G2M.Score.Revelio"]] <- obj_revelio@cellInfo[["G2.M_zScore"]]
obj@meta.data[["G1S.Score.Revelio"]] <- obj_revelio@cellInfo[["G1.S_zScore"]]
obj@meta.data[["G2.Score.Revelio"]] <- obj_revelio@cellInfo[["G2_zScore"]]
obj@meta.data[["MG1.Score.Revelio"]] <- obj_revelio@cellInfo[["M.G1_zScore"]]
obj@meta.data[["Phase.Revelio"]] <- as.character(obj_revelio@cellInfo[["ccPhase"]])
obj_H <- obj
FeatureScatter(obj_H, feature1 = "G1S.Score.Revelio", feature2 = "G2M.Score.Revelio", group.by = "CellAnnCoarse")

remove(new_cc_genes)


FeaturePlot(obj, features = rownames(obj[grep("^Ccn[abde][0-9]$", rownames(obj)), ]))
FeaturePlot(obj, features = rownames(obj[grep("^Cdk", rownames(obj)), ]))