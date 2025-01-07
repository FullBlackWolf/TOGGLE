---
title: "Rawdata Preprocessing: RNA Map of Myocardial"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Preprocessing & Quality Control
---

Load required packages
---
```R
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(future)
library(harmony)
library(RColorBrewer)
```


Use highly variable genes
---
```R
# Create a vector to read files
setwd("C:/GEOANALYSIS/GSE253768")
# Load data
load("C:/GEOANALYSIS/GSE253768/MI-FibroblastCell.Rdata")
# Downgrade the matrix to an older version if necessary
testAB.integrated[["RNA"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
# Extract fibroblasts from the myocardial infarction group
testAB.integrated_MI <- subset(testAB.integrated,idents=c("Fibroblasts"),invert = FALSE)
Idents(testAB.integrated) <- "Group"
testAB.integrated_MI <- subset(testAB.integrated,idents=c("MI"),invert = FALSE)
# Recalculate highly variable genes
testAB.integrated_MI <- NormalizeData(testAB.integrated_MI) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
```

Extract the expression matrix of highly variable genes
---
```R
# Extract the expression matrix of highly variable genes from the data slot
hvg_genes <- VariableFeatures(testAB.integrated_MI)  # Get the names of highly variable genes
hvg_matrix <- testAB.integrated_MI@assays$RNA@data[hvg_genes, ]  # Extract the expression data of these highly variable genes
# Transpose the matrix
hvg_matrix_transposed <- t(hvg_matrix)
# Create a new Seurat object and save the transposed matrix into it
testAB.integrated_MI <- CreateSeuratObject(counts = hvg_matrix_transposed)

# Select highly variable genes and perform dimensionality reduction
testAB.integrated_MI <- NormalizeData(testAB.integrated_MI) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
```

Perform UMAP dimensionality reduction
---
```R
# Smaller n.neighbors: tighter clusters; larger n.neighbors: more uniform distribution
# Smaller min.dist: tighter clusters; larger min.dist: more uniform distribution
# Smaller spread: more concentrated plot; larger spread: better separation between clusters
testAB.integrated_MI <- RunUMAP(testAB.integrated_MI, dims = 1:15, 
                                spread = 2, n.neighbors = 100, min.dist = 0.8)# This separates gene groups into 3 distinct clusters
UMAPPlot(testAB.integrated_MI,label=T)
# Save
save(testAB.integrated_MI, file = "FibroblastCellMatrix Transposition of MI.Rdata")
### Ensure the matrix includes all cells
DefaultAssay(testAB.integrated_MI) <- "RNA"
testAB.integrated_MI[["RNA"]] <- as(object = testAB.integrated_MI[["RNA"]], Class = "Assay")# Convert to version 4 matrix
sceasy::convertFormat(
  testAB.integrated_MI,
  from = "seurat",
  to = "anndata",
  outFile = "FibroblastCellMatrix Transposition of MI.h5ad"
)
```

Import the results
---
```R
index_result <- read.csv("result_RNA.csv")
## Ensure the 'Index' in the table matches the cell names in the Seurat object
## Set 'Index' as row names for easier operations
rownames(index_result) <- index_result$Index
## Map the 'Result' column to the metadata of the Seurat object using 'Index'
metadata <- testAB.integrated_MI@meta.data  # Get the metadata of the Seurat object
metadata$Result <- index_result[rownames(metadata), "Result"]
## Update the metadata of the Seurat object
testAB.integrated_MI@meta.data <- metadata
## Check the updated metadata
head(testAB.integrated_MI@meta.data)
# Get the metadata of the Seurat object
metadata <- testAB.integrated_MI@meta.data
# Create a 'Fenqun' column and group based on 'Result' column values
metadata$Fenqun <- NA  # Initialize the 'Fenqun' column
# Use conditional statements to group 'Result' column values by ranges
metadata$Fenqun[metadata$Result >= 1 & metadata$Result <= 16] <- "RNA-G1"
metadata$Fenqun[metadata$Result >= 17 & metadata$Result <= 20] <- "RNA-G2"
metadata$Fenqun[metadata$Result >= 21 & metadata$Result <= 35] <- "RNA-G3"
metadata$Fenqun[metadata$Result >= 36 & metadata$Result <= 41] <- "RNA-G4"
metadata$Fenqun[metadata$Result >= 42 & metadata$Result <= 61] <- "RNA-G5"
metadata$Fenqun[metadata$Result >= 62 & metadata$Result <= 65] <- "RNA-G6"
metadata$Fenqun[metadata$Result >= 66 & metadata$Result <= 72] <- "RNA-G7"
metadata$Fenqun[metadata$Result >= 73 & metadata$Result <= 92] <- "RNA-G8"
# Update the modified metadata into the Seurat object
testAB.integrated_MI@meta.data <- metadata
# View results
head(testAB.integrated_MI@meta.data$Fenqun)
```

Plotting
---
```R
cell_type_cols <- c("#B383B9", "#EE934E","#F5D2A8","#7DBFA7","#fced82","#D2EBC8","#AECDE1","#3c77af")
p1 <- DimPlot(testAB.integrated_MI, reduction = "umap", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "RNA clustering UMAP of MI.pdf", plot = p1, device = 'pdf', width = 21, height = 18, units = 'cm')
```
<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/RNA-1.png" 
     alt="RNA-1.png" 
     title="RNA-1.png">
