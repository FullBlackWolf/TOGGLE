---
title: "Mature Functions: RNA map of Myocardial"
date: 2024-12-04T13:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---


1.Generate an h5ad file. (R)
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

Set the working directory in R to `C:/GEOANALYSIS/GSE253768`.

```R
# Use highly variable genes
# Create a vector to read files
setwd("C:/GEOANALYSIS/GSE253768")
```

```R
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

# Extract the expression matrix of highly variable genes
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
# Perform UMAP dimensionality reduction
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
  outFile = "MI_fiberRNA.h5ad"
)
```

2.Calculate functional mapping.
---

```python
import numpy as np
import os
import LittleSnowFox as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)

kl.kl_initialize(0)

parent_directory_origin = kl.kl_settings.parent_directory_origin

current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Clustering',1)
```

```python
#选择要使用哪个样本
choosen_sample = "fibroblasts"

#选择.h5ad文件
h5ad_filename = "MI-高变基因.h5ad"


#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder


orig_adata,loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix_sample(
    choosen_sample,
    h5ad_filename,
    "draw",
    current_folder_input,
    round_of_smooth=1,
    neighbor_N=13000,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
    )
#orig_adata,loading_directory,distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)
```

```python
import os
import scipy.io
current_folder_result = os.path.join(loading_directory, 'result')
print(current_folder_result)
mat_name = "distance_matrix_RNA.mat"
mat_path = os.path.join(current_folder_result, mat_name)
scipy.io.savemat(mat_path, {"distance_matrix": distance_matrix})
import pandas as pd

obs_names_list = orig_adata.obs_names.tolist()
var_names_list = orig_adata.var_names.tolist()

index = orig_adata.obs.index
df_orig_adata = pd.DataFrame({
    "Index": index,
})

print(df_orig_adata)
csv_name = "orig_ident_RNA.csv"
csv_path = os.path.join(current_folder_result, csv_name)
print(csv_path)
df_orig_adata.to_csv(csv_path, index=False)
```
3.Unsupervised learning
---

```matlab


%% Load data and Split to compute
%% Load data and Split to compute
MM0 = load('./result/distance_matrix_RNA.mat');
MM0 = MM0.distance_matrix;

%% 读取要排序的对象
count_=readtable('./result/orig_ident_RNA.csv');

%% 得到边界划分点
%[p,splitlist] = binary_corr_sorting(MM0,20,125,5,5);
[p,splitlist] = binary_corr_sorting(MM0,20,50,5,5);

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
h.ColorLimits = [0,0.001]%
% writetable(count_result, './result/result_RNA.csv');

%% 临近法激活
corr_matrix = relevance_generate(0.0001,1,cluster_map_matrix);
hi = heatmap(corr_matrix);


%% 编码
% encode_result = encoder_corr_matrix(0.0011,0.001,10,10,cluster_map_matrix);
encode_result = encoder_corr_matrix(0.00011,0.00010,10,1,cluster_map_matrix);
figure(2)
hj = heatmap(encode_result);

%% 解码
figure(3)
[weighting_decode,decode_result] = decoder_corr_matrix(encode_result);
weighting_result = decode_result + 2*weighting_decode;
hk = heatmap(weighting_result);
hk.ColorLimits = [60,65]
```


4.Generate RNA map
---

```R
# Import the results
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

# Plotting
cell_type_cols <- c("#B383B9", "#EE934E","#F5D2A8","#7DBFA7","#fced82","#D2EBC8","#AECDE1","#3c77af")
p1 <- DimPlot(testAB.integrated_MI, reduction = "umap", group.by = "Fenqun", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "RNA clustering UMAP of MI.pdf", plot = p1, device = 'pdf', width = 21, height = 18, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/RNA-1.png" 
     alt="RNA-1.png" 
     title="RNA-1.png">


Save
---

```R
metadata <- testAB.integrated_MI@meta.data
write.csv(metadata, file="Myocardial Fibroblast Cell RNA Clustering Results.csv")
```

