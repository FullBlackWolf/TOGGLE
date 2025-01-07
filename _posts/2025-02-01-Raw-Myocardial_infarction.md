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

UMAP
---
```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", split.by = "Group", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "FibroblastCell-2.pdf", plot = p2, device = 'pdf', width = 24, height = 12, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-3.png" 
     alt="Myocardial-3.png" 
     title="Myocardial-3.png">
     
UMAP
---
```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", split.by = "Group", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "FibroblastCell-2.pdf", plot = p2, device = 'pdf', width = 24, height = 12, units = 'cm')
```


<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-4.png" 
     alt="Myocardial-4.png" 
     title="Myocardial-4.png">


𝓣𝓪𝓴𝓮 𝓜𝓘'𝓼 𝓘𝓶𝓶𝓾𝓷𝓮𝓒𝓮𝓵𝓵 𝓪𝓷𝓭 𝓶𝓮𝓻𝓰𝓮 𝓲𝓽 𝔀𝓲𝓽𝓱 𝓕𝓲𝓫𝓻𝓸𝓫𝓵𝓪𝓼𝓽𝓒𝓮𝓵𝓵 𝓽𝓸 𝓶𝓪𝓴𝓮 𝓒𝓮𝓵𝓵 𝓬𝓸𝓶𝓶𝓾𝓷𝓲𝓬𝓪𝓽𝓲𝓸𝓷

Load Data
---
```R
load("C:/GEOANALYSIS/GSE253768/MI-FibroblastCell-8000.Rdata")
```
Extract Fibroblasts
---
```R
#Take out the FibroblastCell from the myocardial infarction group
Idents(testAB.integrated) <- "Fenqun"
Fib_seurat <- subset(testAB.integrated,idents=c("FibR1-G5", "FibR1-G6", "FibR1-G7"),invert = FALSE)
Idents(Fib_seurat) <- "Group"
Fib_seurat <- subset(Fib_seurat,idents=c("MI"),invert = FALSE)
Fib_seurat[["RNA"]] <- as(object = Fib_seurat[["RNA"]], Class = "Assay")
#Take out the ImmuneCell from the previous Cell's MI group
```
Extract Immune Cells
---
```R
load("C:/GEOANALYSIS/GSE253768/MI Cell-15 clusters.Rdata")
Idents(testAB.integrated) <- "Group"
Immu_seurat <- subset(testAB.integrated,idents=c("MI"),invert = FALSE)
Idents(Immu_seurat) <- "clusters2"
Immu_seurat <- subset(Immu_seurat,idents=c("Macrophages", "T cells"),invert = FALSE)
Immu_seurat[["RNA"]] <- as(object = Immu_seurat[["RNA"]], Class = "Assay")
```

Merge the two matrices and keep the Cell type data
---
```R
# Get the expression matrix of Fib_seurat and Immu_seurat
Fib_expr <- Fib_seurat@assays$RNA@counts
Immu_expr <- Immu_seurat@assays$RNA@counts
# Get metadata
Fib_meta <- Fib_seurat@meta.data
Immu_meta <- Immu_seurat@meta.data
# Add a new column Source to metadata to mark the source of data
Fib_meta$Source <- "Fib"
Immu_meta$Source <- "Immu"
# Merge expression matrix
combined_expr <- cbind(Fib_expr, Immu_expr)
# Merge metadata
combined_meta <- bind_rows(
  mutate(Fib_meta, Cell_Barcode = rownames(Fib_meta)),
  mutate(Immu_meta, Cell_Barcode = rownames(Immu_meta))
)
# Ensure that the Cell barcode and expression matrix column names are consistent
combined_meta <- combined_meta %>% filter(Cell_Barcode %in% colnames(combined_expr))
# Create a new Seurat object
combined_seurat <- CreateSeuratObject(
  counts = combined_expr,
  meta.data = combined_meta
)
```

Add Cell Type Information
---
```R
# View the column names in metadata
colnames(combined_seurat@meta.data)
# View the first few rows of metadata to ensure that the correct columns are included
head(combined_seurat@meta.data)
# Ensure that the Source column has been created correctly
print(head(combined_seurat$Source))
# Make sure the Fenqun and clusters2 columns have data
print(head(combined_seurat$Fenqun))
print(head(combined_seurat$clusters2))
# Create a cell_type column and view the results
combined_seurat$cell_type <- dplyr::case_when(
  combined_seurat$Source == "Fib" ~ as.character(combined_seurat$Fenqun),
  combined_seurat$Source == "Immu" ~ as.character(combined_seurat$clusters2),
  TRUE ~ NA_character_
)
```

Clean Metadata
---
```R
# View the newly created cell_type column
print(head(combined_seurat$cell_type))
# View the combined Seurat object
combined_seurat
combined_seurat$Sample <- NULL
combined_seurat$Group <- NULL
combined_seurat$Result <- NULL
combined_seurat$Fenqun <- NULL
combined_seurat$RNA_snn_res.0.18 <- NULL
combined_seurat$seurat_clusters <- NULL
combined_seurat$clusters1 <- NULL
combined_seurat$clusters2 <- NULL
```

Perform UMAP Analysis
---
```R
#do umap
combined_seurat <- SCTransform(combined_seurat,assay = 'RNA')
combined_seurat <- RunPCA(combined_seurat)
ElbowPlot(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)
```

Visualize and Save UMAP
---

```R
UMAPPlot(combined_seurat,group.by='cell_type',label=T)
cell_type_cols <- c("#6693b1","#a3caa9","#efd695","#dd9667","#bd5c56")
p1 <- DimPlot(combined_seurat, reduction = "umap", group.by = "cell_type", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "UMAP of ImmuneCell and grouped FibroblastCell.pdf", plot = p1, device = 'pdf', width = 18, height = 15, units = 'cm')
```


<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-5.png" 
     alt="Myocardial-5.png" 
     title="Myocardial-5.png">


Save Seurat Object
---
```R
# Save
save(combined_seurat, file = "ImmuneCell and grouped FibroblastCell files.Rdata")
```

Load CellChat Library and Preprocess Data
---
```R

#####Start CCC analysis
library(CellChat)
DefaultAssay(combined_seurat) <- "RNA"
combined_seurat <- NormalizeData(combined_seurat, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
```
Prepare Data for CellChat
---

```R
#Propose the required data
data.input  <- combined_seurat@assays$RNA$data
identity = data.frame(group =combined_seurat$cell_type, row.names = names(combined_seurat$cell_type)) 
unique(identity$group) # check the cell labels
```
Create and Configure CellChat Object
---
```R
#Create a cellchat object
cellchat <- createCellChat(object = data.input)

#Add metadata information to the CellChat object
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))
```

Load CellChat Database
---
```R
#Load and set the required CellChatDB database
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB # set the used database in the object
```
Preprocess and Infer Interactions & Compute Communication Probability
---
```R
#Preprocess expression data for cell-to-cell interaction analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
#Infer the interaction network between cells and analyze
cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)#Filter cells with less communication
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
Get all ligand-receptor pairs and their communication probabilities
df.net <- subsetCommunication(cellchat)
#Extract communication information by pathway
df.pathway = subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net, "Fib and ImmuneCell Interaction.csv", quote = F, sep = ',')
write.csv(df.pathway, "Fib and ImmuneCell Interactions - Characterized by Pathways.csv",quote = F,sep = ',')
#Analyze intercellular communication network
cellchat <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat, file = "cellchat-Fib and ImmuneCell.rds")
#Drawing
Fib = c("FibR1-G7", "FibR1-G5", "FibR1-G6")
Immu = c("Macrophages", "T cells")
library(patchwork)
```

Save graph function
---
```R
par(mfrow=c(1,2),xpd=T)
#netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
#                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Fibroblasts as Targets", 
                 sources.use = Immu,
                 targets.use = Fib)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Immune cells as Targets",
                 sources.use = Fib,
                 targets.use = Immu)
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-6.png" 
     alt="Myocardial-6.png" 
     title="Myocardial-6.png">

```R
#Make Cell communication bubble chart
p1 <- netVisual_bubble(cellchat, sources.use = Immu, targets.use = Fib, title.name = "Fibroblasts as Targets", remove.isolate = T) 
p2 <- netVisual_bubble(cellchat, sources.use = Fib, targets.use = Immu,title.name = "Immune cells as Targets", remove.isolate = T)
p1 + p2
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-7.png" 
     alt="Myocardial-7.png" 
     title="Myocardial-7.png">

Perform mixed analysis of the common transcriptome
---
```R
library(SeuratDisk)
# Set the input file path and output file path
input_file <- "tissue_and_sc.h5ad" # Replace with your h5ad file path
output_file <- "tissue_and_sc.h5seurat" # Convert to temporary h5seurat file
# Convert .h5ad to .h5seurat format
Convert(input_file, dest = "h5seurat", overwrite = TRUE)
# Load .h5seurat file to Seurat object
seurat_object <- LoadH5Seurat(output_file)
# View Seurat object information
print(seurat_object)
```


Regenerate meta.data
---
```R
load("C:/GEOANALYSIS/GSE253768/MI-FibroblastCell.Rdata")
# Extract meta.data of testAB.integrated and seurat_object
meta_testAB <- testAB.integrated@meta.data
meta_seurat <- seurat_object@meta.data
# Extract "orig.ident", "Group", "Sample" from testAB.integrated
selected_metadata <- meta_testAB[, c("orig.ident", "Group", "Sample")]
# Ensure row names are aligned
rownames(selected_metadata) <- rownames(meta_testAB)
# Initialize new metadata for seurat_object
meta_seurat$orig.ident <- NA
meta_seurat$Group <- NA
meta_seurat$Sample <- NA
# Align metadata from testAB.integrated to seurat_object
common_cells <- intersect(rownames(meta_seurat), rownames(selected_metadata))
meta_seurat[common_cells, c("orig.ident", "Group", "Sample")] <- selected_metadata[common_cells, ]
# Add new labels for extra cells in seurat_object
additional_cells <- setdiff(rownames(meta_seurat), rownames(selected_metadata))
# Set new Group and Sample names
meta_seurat[additional_cells, "orig.ident"] <- additional_cells # Set orig.ident as the cell label
meta_seurat[additional_cells, "Group"] <- additional_cells # Set Group as the cell label
meta_seurat[additional_cells, "Sample"] <- additional_cells # Set Sample as the cell label
# Update meta.data of seurat_object
seurat_object@meta.data <- meta_seurat
# View the updated meta.data
head(seurat_object@meta.data)
```
To incorporate a new grouping (e.g., "Result") into the metadata of a Seurat object and update it with the grouping results from an external file
---
```R
##Import the new grouping results
index_result <- read.csv("result_mixed mapping.csv")
##Make sure the table's Index is consistent with the Cell name of the Seurat object
##Set the Index as the row name to facilitate subsequent operations
rownames(index_result) <- index_result$Var1
##Match the Result column to the metadata of the Seurat object according to the Index
metadata <- seurat_object@meta.data # Get the metadata of the Seurat object
metadata$Result <- index_result[rownames(metadata), "Result"]
##Update the metadata of the Seurat object
seurat_object@meta.data <- metadata
##Check the updated metadata
head(seurat_object@meta.data)
#Add Fenqun3
```

Specify Cell grouping
---
```R
mi_cells <- c("MI.1", "MI.2", "MI.3", "MI.4")
#Update the meta.data of the Seurat object and create a new Fenqun3 column
seurat_object@meta.data <- seurat_object@meta.data %>%
  mutate(Fenqun4 = case_when(
    Cell_Barcode %in% mi_cells ~ "MI", 
    Result >= 1 & Result <= 9 ~ "FibR2-G1",
    Result >= 10 & Result <= 12 ~ "FibR2-G2",
    Result >= 13 & Result <= 25 ~ "FibR2-G3",
    Result >= 26 & Result <= 33 ~ "FibR2-G4",
    Result >= 34 & Result <= 41 ~ "FibR2-G5",
    Result >= 42 & Result <= 50 ~ "FibR2-G6",
    Result >= 51 & Result <= 59 ~ "FibR2-G7",
    TRUE ~ NA_character_ 
  ))
# View updated meta.data
head(seurat_object@meta.data)
```

Redraw UMAP
---
```R
seurat_object <- SCTransform(seurat_object,assay = 'RNA')
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
```

Make a plot highlighting MI and Sham
---
```R
# Create a new column Fenqun3's UMAP plot
cell_type_cols <- c("#ffffcc","#5a5098","#a3caa9","#deedad","#bd5c56","#efd695","#dd9667","#6693b1","#842844")
# Create a new column to define the point size based on the "Fenqun4" column
seurat_object@meta.data$point_size <- ifelse(seurat_object@meta.data$Fenqun4 %in% c("MI"), 3, 0.1)
# Use DimPlot to make a UMAP plot, keeping the original colors and labels, but resizing the points
p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "Fenqun4", pt.size = seurat_object@meta.data$point_size, label = TRUE, repel = TRUE, raster = FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme( panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "UMAP diagram highlighting MI and Sham points.pdf", plot = p1, device = 'pdf', width = 21, height = 18, units = 'cm') 
```
     

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
