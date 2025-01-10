---
title: "Mature Functions: Function Difference of Myocardial"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---


1.Preprocessing
---

The default data file storage directory is `C:/GEOANALYSIS/GSE253768`

1.1.Load required R packages (R)
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

1.2.Read file names and set the working directory (R)
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

1.3.Batch read data and create Seurat objects (R)
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

1.4.Calculate mitochondrial and red blood cell proportions (R)
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

1.5.Quality control and preliminary merging (R)
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

1.6.Filter cells (R)
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

1.7.Data normalization, feature selection, and dimensionality reduction (R)
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

1.8.Harmony integration analysis (R)
---

```R
# Integrate using Harmony
testAB.integrated <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
# Copy 'orig.ident' to 'Sample'
testAB.integrated@meta.data$Sample <- testAB.integrated@meta.data$orig.ident
# Copy 'orig.ident' to 'Group' and remove numbers
testAB.integrated@meta.data$Group <- gsub("[0-9]", "", testAB.integrated@meta.data$orig.ident)
```


1.9.Clustering and dimensionality reduction visualization (R)
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

1.10.Annotation and cluster labeling (R)
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

1.11.Export results and visualization (R)
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

1.12.Generate the UMAP Plot (R)
---

```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", split.by = "Group", group.by = "clusters2", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Preliminary grouping of MI - split by group.pdf", plot = p2, device = 'pdf', width = 38, height = 18, units = 'cm')
```

2.Calculating similarity between cells
---

2.1.Exporting Data (R)
---

```R
# Extract and save fibroblasts
testAB.integrated <- subset(testAB.integrated,idents=c("Fibroblasts","Myofibroblast"),invert = FALSE)
testAB.integrated$RNA_snn_res.0.18 <- NULL
testAB.integrated$clusters1 <- NULL
testAB.integrated$clusters2 <- NULL
testAB.integrated$seurat_clusters <- NULL
#Re-finding highly variable genes
testAB.integrated <- NormalizeData(testAB.integrated) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
# Save
save(testAB.integrated, file = "MI-FibroblastCell.Rdata")
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
  outFile = "Myocardial_infarction_hypermutable_gene.h5ad" #心梗-高变基因.h5ad
) 
```
Save the R variable table as `variables.Rdata`.   

2.2.Calculating similarity between cells (Python)
---

Move the generated `Myocardial_infarction_hypermutable_gene.h5ad` in the folder to `[LittleSnowFox library file address]\Lib\site-packages\kailin\database\Clustering_sample\fibroblasts\data\`
```html
#link
https://biocomputing.cowtransfer.com/s/b7f5aa9cc9ee4e

#Password  
fl1n09
```


```python
import numpy as np
import os
import kailin as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)
kl.kl_initialize(0)
parent_directory_origin = kl.kl_settings.parent_directory_origin

current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Clustering',1)


h5ad_filename = "Myocardial_infarction_hypermutable_gene.h5ad"
num_samples = 8000


current_folder_filter = os.path.join(current_folder, 'fibroblasts')
current_folder_filter_1 = os.path.join(current_folder_filter, 'data')
filter_sample = os.path.join(current_folder_filter_1, h5ad_filename)
print(filter_sample)

import anndata
adata_filter = anndata.read_h5ad(filter_sample)

total_samples = adata_filter.shape[0]
if num_samples > total_samples:
    raise ValueError(f"Requested {num_samples} samples, but only {total_samples} available.")
random_indices = np.random.choice(total_samples, num_samples, replace=False)

# 根据索引过滤数据
adata_random_subset = adata_filter[random_indices, :]

# 打印选取后的数据维度
print(f"Subset data shape: {adata_random_subset.shape}")
print(adata_random_subset)
print(current_folder_filter_1)

file_name = "fibroblasts_random_8000.h5ad"
file_path = os.path.join(current_folder_filter_1, file_name)
adata_random_subset.write(file_path)


adata_random_subset


# 查看前5个行标签和列标签
print(adata_random_subset)

#选择要使用哪个样本
choosen_sample = "fibroblasts"
#选择.h5ad文件
h5ad_filename = "fibroblasts_random_8000.h5ad"

#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
orig_adata,loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input,1,13000,0.1,0.001,True)
#orig_adata,loading_directory,distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)


#current_folder_input = current_folder
#loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input)
print(loading_directory)
print(choosen_sample)


import os
import scipy.io
current_folder_result = os.path.join(loading_directory, 'result')
print(current_folder_result)
mat_name = "distance_matrix.mat"
mat_path = os.path.join(current_folder_result, mat_name)
scipy.io.savemat(mat_path, {"distance_matrix": distance_matrix})
print(f"Distance matrix saved as MAT file to: {file_path}")


import pandas as pd

obs_names_list = orig_adata.obs_names.tolist()
var_names_list = orig_adata.var_names.tolist()

print(orig_adata.obs['orig.ident'])
index = orig_adata.obs.index
values = orig_adata.obs['orig.ident']

df_orig_adata = pd.DataFrame({
    "Index": index,
    "Value": values
})

print(df_orig_adata)

csv_name = "orig_ident.csv"
csv_path = os.path.join(current_folder_result, csv_name)
print(csv_path)
df_orig_adata.to_csv(csv_path, index=False)
```


2.3.Unsupervised learning process (Matlab)
---

Afterward, execute the following file: `[LittleSnowFox's Anaconda installation directory]\database\Clustering_sample\fibroblasts\main_v3_matlab_run_me_difference_gene.m`


```matlab


%% Load data and Split to compute
%% Load data and Split to compute
MM0 = load('./result/distance_matrix.mat');
MM0 = MM0.distance_matrix;

%% 读取要排序的对象
count_=readtable('./result/orig_ident.csv');

%% 得到边界划分点
%[p,splitlist] = binary_corr_sorting(MM0,20,125,5,5);
[p,splitlist] = binary_corr_sorting(MM0,20,200,5,5);

%% 对划分点去重
[uniqueList, ~, ~] = unique(splitlist, 'stable');

%% 对相似度矩阵排序
MM=MM0(p,p);
split=[];

%% 重排count_result
count_result=count_(p,:);
split_simple=uniqueList;

%% 第一个起始位点置为1
split_simple(1)=1;
split_simple=[split_simple,length(MM0)];

%% 计算均值矩阵
[simple_matrix]=sample_computing(count_result,split_simple,MM,"mean");

%% 合并成小矩阵
ClusterReslut=cluster_map(split_simple,simple_matrix,0,0.0002,0);
count_result.Result = ClusterReslut;

%重排小矩阵
[cluster_map_matrix] = genetic_encoder( ...
    simple_matrix, ...
    60, ...% nPop = 50;  % 种群规模大小为30
    1, ...% nPc = 1; % 子代规模的比例0.8
    200, ...% maxIt = 200; % 最大迭代次数
    5 ...% cycletimes = 200; % 循环计算次数
    );

%重拍小矩阵方案2
% 创建行和列标签（示例）
%row_labels = cluster_map_label;
%column_labels = cluster_map_label;
% 使用 heatmap 函数并传递相应参数
h = heatmap(cluster_map_matrix);
%h.YDisplayLabels = row_labels; % 设置行标签
%h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0,0.0001]%

%if you want to save that
writetable(count_result, './result/result_diff.csv'); #导入R


%% 临近法激活
corr_matrix = relevance_generate(0.00011,1,cluster_map_matrix);
hi = heatmap(corr_matrix);


%% 编码
encode_result = encoder_corr_matrix(0.000111,0.000109,10,1,cluster_map_matrix);
figure(2)
hj = heatmap(encode_result);

%% 解码
figure(3)
[weighting_decode,decode_result] = decoder_corr_matrix(encode_result);
weighting_result = decode_result;
hk = heatmap(weighting_result);
hk.ColorLimits = [30,31]
```

2.4.Use the group information (R)
---

Open the R variable table as `variables.Rdata`.


```R
#reload
load("C:/GEOANALYSIS/GSE253768/MI-FibroblastCell.Rdata")
```

The original location of `result_diff.csv` is `[Package path of LittleSnowFox]\database\Clustering_sample\fibroblasts\result\result_diff.csv`.

```R
# Import the results
index_result <- read.csv("result_diff.csv")
```

```R
## Ensure that the table's Index is consistent with the Cell name of the Seurat object
## Set the Index to the row name to facilitate subsequent operations
rownames(index_result) <- index_result$Index
##Match the Result column to the metadata of the Seurat object according to the Index
metadata <- testAB.integrated@meta.data # Get the metadata of the Seurat object
metadata$Result <- index_result[rownames(metadata), "Result"]
##Update the metadata of the Seurat object
testAB.integrated@meta.data <- metadata
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

2.5.Generate groups (R)
---

```R
#
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
# Redraw UMAP
testAB.integrated <- SCTransform(testAB.integrated,assay = 'RNA')
testAB.integrated <- RunPCA(testAB.integrated)
ElbowPlot(testAB.integrated)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:10)
UMAPPlot(testAB.integrated,group.by='Fenqun',label=T)
## Save
save(testAB.integrated, file = "MI-FibroblastCell-8000.Rdata")
# Export markers for different groups
Idents(testAB.integrated) <- "Fenqun"
CI.markers <- FindAllMarkers(testAB.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CI.markers, file="FibroblastCell markers.csv")
## Visualize and export
cell_type_cols <- c("#5a5098","#6693b1","#a3caa9","#deedad","#ffffcc","#efd695","#dd9667","#bd5c56","#842844")
p1 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "FibroblastCell-1.pdf", plot = p1, device = 'pdf', width = 15, height = 12, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-3.png" 
     alt="Myocardial-3.png" 
     title="Myocardial-3.png">

```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", split.by = "Group", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "FibroblastCell-2.pdf", plot = p2, device = 'pdf', width = 24, height = 12, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-4.png" 
     alt="Myocardial-4.png" 
     title="Myocardial-4.png">


3.Take MI's ImmuneCell and merge it with FibroblastCell to make Cell communication (R)
---

```R
load("C:/GEOANALYSIS/GSE253768/MI-FibroblastCell-8000.Rdata")
#Take out the FibroblastCell from the myocardial infarction group
Idents(testAB.integrated) <- "Fenqun"
Fib_seurat <- subset(testAB.integrated,idents=c("FibR1-G5", "FibR1-G6", "FibR1-G7"),invert = FALSE)
Idents(Fib_seurat) <- "Group"
Fib_seurat <- subset(Fib_seurat,idents=c("MI"),invert = FALSE)
Fib_seurat[["RNA"]] <- as(object = Fib_seurat[["RNA"]], Class = "Assay")
```

3.1.Take out the ImmuneCell from the previous Cell's MI group (R)
---

```R
load("C:/GEOANALYSIS/GSE253768/MI Cell-15 clusters.Rdata")
Idents(testAB.integrated) <- "Group"
Immu_seurat <- subset(testAB.integrated,idents=c("MI"),invert = FALSE)
Idents(Immu_seurat) <- "clusters2"
Immu_seurat <- subset(Immu_seurat,idents=c("Macrophages", "T cells"),invert = FALSE)
Immu_seurat[["RNA"]] <- as(object = Immu_seurat[["RNA"]], Class = "Assay")
```

3.2.Merge the two matrices and keep the Cell type data (R)
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

3.3.View the combined Seurat (R)
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
#do umap
combined_seurat <- SCTransform(combined_seurat,assay = 'RNA')
combined_seurat <- RunPCA(combined_seurat)
ElbowPlot(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)
UMAPPlot(combined_seurat,group.by='cell_type',label=T)
cell_type_cols <- c("#6693b1","#a3caa9","#efd695","#dd9667","#bd5c56")
```

3.4.View the result (R)
---

```R
p1 <- DimPlot(combined_seurat, reduction = "umap", group.by = "cell_type", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "UMAP of ImmuneCell and grouped FibroblastCell.pdf", plot = p1, device = 'pdf', width = 18, height = 15, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Myocardial-5.png" 
     alt="Myocardial-5.png" 
     title="Myocardial-5.png">

4.Cell signaling (R)
---

```R
# Save
save(combined_seurat, file = "ImmuneCell and grouped FibroblastCell files.Rdata")


#####Start CCC analysis
library(CellChat)
DefaultAssay(combined_seurat) <- "RNA"
combined_seurat <- NormalizeData(combined_seurat, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
#Propose the required data
data.input  <- combined_seurat@assays$RNA$data
identity = data.frame(group =combined_seurat$cell_type, row.names = names(combined_seurat$cell_type)) 
unique(identity$group) # check the cell labels

#Create a cellchat object
cellchat <- createCellChat(object = data.input)

#Add metadata information to the CellChat object
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))

#Load and set the required CellChatDB database
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB # set the used database in the object

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
#Get all ligand-receptor pairs and their communication probabilities
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
#Save graph function
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

