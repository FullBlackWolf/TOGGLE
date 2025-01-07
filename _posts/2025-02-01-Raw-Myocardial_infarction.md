---
title: "Rawdata Preprocessing: Myocardial Infarction"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Preprocessing & Quality Control
---


Load required R packages
---
```R
# Load required packages
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
Read file names and set the working directory
---
```R
# Create a vector to read files
setwd("C:/GEOANALYSIS/GSE253768")
## Save the file names in the folder to 'dir_name'
dir_name <- list.files(pattern = "\\.csv$") # Match only CSV files
## View 'dir_name'
dir_name
#[1] "MI1.csv"   "MI2.csv"   "Sham1.csv" "Sham2.csv"
## Assign names to 'dir_name' (use file name without extension)
names(dir_name) <- gsub("\\.csv$", "", dir_name)
## View the renamed 'dir_name'
dir_name
##      MI1         MI2       Sham1       Sham2 
##"MI1.csv"   "MI2.csv" "Sham1.csv" "Sham2.csv" 
```

Batch read data and create Seurat objects
---
```R
## Batch data processing
scRNAlist <- list()
for (i in 1:length(dir_name)) {
  # Read CSV file
  counts <- read.csv(file = dir_name[i], row.names = 1) # Use the first column as row names
  counts <- as.matrix(counts) # Convert to matrix format
  
  # Create Seurat object and add file name as a label
  scRNAlist[[i]] <- CreateSeuratObject(
    counts = counts,
    min.cells = 3,
    min.features = 300,
    project = names(dir_name)[i]
  )
}
```

Calculate mitochondrial and red blood cell proportions
---
```R
# Check the Seurat object list
scRNAlist

# Calculate mitochondrial and red blood cell proportions in batch
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]
  # Calculate mitochondrial proportion
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^Mt-")
  # Calculate red blood cell proportion
  HB_genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz")
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m] 
  HB_genes <- HB_genes[!is.na(HB_genes)] 
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, pattern = "^Hb") 
  # Assign 'sc' back to scRNAlist[[i]]
  scRNAlist[[i]] <- sc
  # Remove 'sc'
  rm(sc)
}
```

Quality control and preliminary merging
---
```R
# Perform a simple merge and then plot quality control (QC)
CI <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]]))
head(colnames(CI))
unique(sapply(X = strsplit(colnames(CI), split = "_"), FUN = "[", 1))

plot1 <- FeatureScatter(CI, feature1 = "nFeature_RNA", feature2 = "mt_percent")
plot2 <- FeatureScatter(CI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))# QC plot 1
VlnPlot(CI, features = c("mt_percent", "nFeature_RNA", "nCount_RNA", "HB_percent"), ncol = 4, pt.size=0)# QC plot 2
VlnPlot(CI, features = c("mt_percent", "nFeature_RNA", "nCount_RNA", "HB_percent"), ncol = 4, pt.size=0.5)# QC plot 3
```

Filter cells
---
```R
# Filter cells in batch
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
                mt_percent < 10 & 
                HB_percent < 5 & 
                nCount_RNA < quantile(nCount_RNA,0.97))})
```

Data normalization, feature selection, and dimensionality reduction
--- 
```R
# Merge Seurat objects
scRNAlist <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]]))
# Select highly variable genes and perform dimensionality reduction
scRNAlist <- NormalizeData(scRNAlist) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
```

Harmony integration analysis
---
```R
# Integrate using Harmony
testAB.integrated <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
# Copy 'orig.ident' to 'Sample'
testAB.integrated@meta.data$Sample <- testAB.integrated@meta.data$orig.ident
# Copy 'orig.ident' to 'Group' and remove numbers
testAB.integrated@meta.data$Group <- gsub("[0-9]", "", testAB.integrated@meta.data$orig.ident)
```


Clustering and dimensionality reduction visualization
---
```R
# Check the updated metadata
head(testAB.integrated@meta.data)

# Add grouping information after integration
metadata <- testAB.integrated@meta.data
write.csv(metadata, file="meta.data.csv")# Export and save
#testAB.integrated@meta.data <- metadata

# Perform clustering
testAB.integrated <- FindNeighbors(testAB.integrated, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.18)#15群
# Perform UMAP/tSNE dimensionality reduction
testAB.integrated <- RunTSNE(testAB.integrated, reduction = "harmony", dims = 1:15)
testAB.integrated <- RunUMAP(testAB.integrated, reduction = "harmony", dims = 1:25)
# Save
save(testAB.integrated, metadata, file = "MI Cell-15 clusters.Rdata")
# Export markers
testAB.integrated <- JoinLayers(testAB.integrated)
CI.markers <- FindAllMarkers(testAB.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CI.markers, file="MI Cell marker.csv")
```

Annotation and cluster labeling
---
```R
# Naming the 15 clusters
new.cluster.ids <- c("Fibroblasts", "Cardiomyocytes", "Endothelial cells",
                     "Cardiomyocytes", "Smooth muscle cells", "Macrophages",
                     "T cells", "Endothelial cells", "Endothelial cells",
                     "Mesothelial cells", "Macrophages", "Endothelial cells",
                     "Schwann cells", "Myofibroblast", "Endothelial cells")
names(new.cluster.ids) <- levels(testAB.integrated)
testAB.integrated <- RenameIdents(testAB.integrated, new.cluster.ids)
testAB.integrated$clusters2 <- testAB.integrated@active.ident

save(testAB.integrated, metadata, file = "MI Cell-15 clusters.Rdata")
```
Export results and visualization
---
```R
# Export the count of each cluster
Table1 <- table(testAB.integrated$Group, testAB.integrated$clusters2)
Table2 <- table(testAB.integrated$Sample, testAB.integrated$clusters2)
write.table(Table1, file = "Cell counts in each group.txt", sep ="\t")
write.table(Table2, file = "Cell counts in each sample.txt", sep ="\t")

# Export UMAP images from preliminary analysis
cell_type_cols <- c("#B383B9", "#F5CFE4","#EE934E","#F5D2A8","#fced82","#D2EBC8","#7DBFA7","#AECDE1","#3c77af")
p1 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "clusters2", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Preliminary grouping of MI - overall.pdf", plot = p1, device = 'pdf', width = 21, height = 18, units = 'cm')

```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-1.png" 
     alt="Myocardial-1.png" 
     title="Myocardial-1.png">

Generate the UMAP Plot
---
```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", split.by = "Group", group.by = "clusters2", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Preliminary grouping of MI - split by group.pdf", plot = p2, device = 'pdf', width = 38, height = 18, units = 'cm')
```

Extract Fibroblast Cells
---

```R
# Extract and save fibroblasts
testAB.integrated <- subset(testAB.integrated,idents=c("Fibroblasts","Myofibroblast"),invert = FALSE)
testAB.integrated$RNA_snn_res.0.18 <- NULL
testAB.integrated$clusters1 <- NULL
testAB.integrated$clusters2 <- NULL
testAB.integrated$seurat_clusters <- NULL
```

Re-normalize and Identify Highly Variable Genes
---

```R
#Re-finding highly variable genes
testAB.integrated <- NormalizeData(testAB.integrated) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
# Save
save(testAB.integrated, file = "MI-FibroblastCell.Rdata")
```

Export as AnnData Format
---

```R
# Export as h5ad version
high_var_genes <- VariableFeatures(testAB.integrated)  
testAB.integrated[["RNA"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
data_high_var <- testAB.integrated@assays$RNA@data[high_var_genes, ]
# Create a new Seurat object containing only highly variable genes
testAB_high_var <- subset(
  x = testAB.integrated,
  features = high_var_genes
)
#Convert to AnnData format and save
sceasy::convertFormat(
  testAB_high_var,
  from = "seurat",
  to = "anndata",
  outFile = "MI-FibroblastCell.h5ad"
)
```

Load and Integrate Additional Data
---

```R

#reload
load("C:/GEOANALYSIS/GSE253768/MI-FibroblastCell.Rdata")
# Import the results
index_result <- read.csv("result_DEG.csv")
## Ensure that the table's Index is consistent with the Cell name of the Seurat object
## Set the Index to the row name to facilitate subsequent operations
rownames(index_result) <- index_result$Index
##Match the Result column to the metadata of the Seurat object according to the Index
metadata <- testAB.integrated@meta.data # Get the metadata of the Seurat object
metadata$Result <- index_result[rownames(metadata), "Result"]
##Update the metadata of the Seurat object
testAB.integrated@meta.data <- metadata
```

Filter Cells Based on Metadata
---

```R
##Check the updated metadata
head(testAB.integrated@meta.data)
##Filter the cells in the metadata with Result column from 1 to 57
cells_to_keep <- rownames(testAB.integrated@meta.data[testAB.integrated@meta.data$Result >= 1 &
                                                        testAB.integrated@meta.data$Result <= 57, ])
##Extract these cells and form a new Seurat object
testAB.integrated <- subset(testAB.integrated, cells = cells_to_keep)
## Save
save(testAB.integrated, file = "MI-FibroblastCell-8000.Rdata")
```

Group Cells into Categories
---

```R
#Add new columns in metadata to group these
testAB.integrated@meta.data <- testAB.integrated@meta.data %>%
  mutate(Fenqun = case_when(
    Result >= 1 & Result <= 2 ~ "FibR1-G1",
    Result >= 3 & Result <= 7 ~ "FibR1-G2",
    Result >= 8 & Result <= 15 ~ "FibR1-G3",
    Result >= 16 & Result <= 26 ~ "FibR1-G4",
    Result >= 27 & Result <= 40 ~ "FibR1-G5",
    Result >= 41 & Result <= 48 ~ "FibR1-G6",
    Result >= 49 & Result <= 57 ~ "FibR1-G7",
    TRUE ~ NA_character_ 
  ))
```

Reprocess and Visualize
---

```R
# Redraw UMAP
testAB.integrated <- SCTransform(testAB.integrated,assay = 'RNA')
testAB.integrated <- RunPCA(testAB.integrated)
ElbowPlot(testAB.integrated)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:10)
UMAPPlot(testAB.integrated,group.by='Fenqun',label=T)
## Save
save(testAB.integrated, file = "MI-FibroblastCell-8000.Rdata")
```

Export Cluster Markers
---

```R
# Export markers for different groups
Idents(testAB.integrated) <- "Fenqun"
CI.markers <- FindAllMarkers(testAB.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CI.markers, file="FibroblastCell markers.csv")
```

Save and Export UMAP Plot
---

```R

## Visualize and export
cell_type_cols <- c("#5a5098","#6693b1","#a3caa9","#deedad","#ffffcc","#efd695","#dd9667","#bd5c56","#842844")
p1 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "FibroblastCell-1.pdf", plot = p1, device = 'pdf', width = 15, height = 12, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-2.png" 
     alt="Myocardial-2.png" 
     title="Myocardial-2.png">

```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", split.by = "Group", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "FibroblastCell-2.pdf", plot = p2, device = 'pdf', width = 24, height = 12, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-3.png" 
     alt="Myocardial-3.png" 
     title="Myocardial-3.png">


<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-4.png" 
     alt="Myocardial-4.png" 
     title="Myocardial-4.png">

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-5.png" 
     alt="Myocardial-5.png" 
     title="Myocardial-5.png">

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-6.png" 
     alt="Myocardial-6.png" 
     title="Myocardial-6.png">

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-7.png" 
     alt="Myocardial-7.png" 
     title="Myocardial-7.png">

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-8.png" 
     alt="Myocardial-8.png" 
     title="Myocardial-8.png">
     

Save
---
```R
## Save 
save(seurat_object, file = "MI-FibroblastCell-mixed with whole transcriptome data.Rdata")
# Export markers
Idents(seurat_object) <- "Fenqun"
CI.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CI.markers, file="8000 Cells mixed with normal transcriptome markers.csv")
```
