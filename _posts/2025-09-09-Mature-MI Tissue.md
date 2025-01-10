---
title: "Mature Functions: Function Map of MI Tissue"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---

The prerequisite analysis `Mature Functions: Function Difference of Myocardial` needs to be completed first.

5.Perform mixed analysis of the common transcriptome
---

5.1.Prepare the single-cell maps needed for use.
---

```python

import numpy as np
import os
import LittleSnowFox as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)
#初始化函数，将kailin转至工作目录。如果此前初始化过，那么在再次运行def kl_initialize(0)时，
#则拒绝初始化，避免套娃。运行def kl_initialize(1)时，强制重新初始化。
kl.kl_initialize(0)
#获取kailin工作的根目录
parent_directory_origin = kl.kl_settings.parent_directory_origin

#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Clustering',1)
```

```python
h5ad_filename = "心梗-成纤维细胞-所有基因.h5ad"
file_name = "抽样_全基_8000.h5ad"
sample_choose = 'fibroblasts'
num_samples = 8000


current_folder_filter = os.path.join(current_folder, sample_choose)
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

file_path = os.path.join(current_folder_filter_1, file_name)
adata_random_subset.write(file_path)


adata_random_subset
```
```python
#选择要使用哪个样本
choosen_sample = "fibroblasts"
#选择.h5ad文件
h5ad_filename = "抽样_全基_8000.h5ad"

#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
orig_adata,loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input,1,13000,0.1,0.001,True)
#orig_adata,loading_directory,distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)
#运行自带的示例，并获取非稀疏矩阵
#这里需要做非示例的函数进去
#current_folder_input = current_folder
#loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input)
print(loading_directory)
print(choosen_sample)
```

```python
import os
import scipy.io
current_folder_result = os.path.join(loading_directory, 'result')
print(current_folder_result)
mat_name = "distance_matrix_all.mat"
mat_path = os.path.join(current_folder_result, mat_name)
scipy.io.savemat(mat_path, {"distance_matrix": distance_matrix})
print(f"Distance matrix saved as MAT file to: {file_path}")
```

```python
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

csv_name = "orig_ident_all.csv"
csv_path = os.path.join(current_folder_result, csv_name)
print(csv_path)
df_orig_adata.to_csv(csv_path, index=False)
```

5.2.Projecting organizations onto a single-cell functional map


5.2.1.Calculate the cell similarity of the whole genome.
---


5.2.2.Read data

```python
import numpy as np
import os
import LittleSnowFox as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)
#初始化函数，将kailin转至工作目录。如果此前初始化过，那么在再次运行def kl_initialize(0)时，
#则拒绝初始化，避免套娃。运行def kl_initialize(1)时，强制重新初始化。
kl.kl_initialize(0)
#获取kailin工作的根目录
parent_directory_origin = kl.kl_settings.parent_directory_origin

#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Clustering',1)
```

5.2.3.Combine

```python
h5ad_filename = "fibroblasts_all_gene.h5ad" #心梗-成纤维细胞-所有基因.h5ad
csv_filename = "心肌梗死普通转录组.csv" #心梗-成纤维细胞-所有基因.h5ad
sample_choose = 'fibroblasts'
num_samples = 8000

adata_merged = kl.preprocessing.merge_tissue_with_singlecell(
    current_folder,
    h5ad_filename,
    csv_filename,
    sample_choose,
    num_samples
    )


```

5.2.4.Perform the mapping from tissue to single cells.
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

## Regenerate meta.data
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
#Specify Cell grouping
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

# Redraw UMAP
seurat_object <- SCTransform(seurat_object,assay = 'RNA')
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

## Make a plot highlighting MI and Sham
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

Save the result
---

```R
save(seurat_object, file = "MI-FibroblastCell-mixed with whole transcriptome data.Rdata")
# Export markers
Idents(seurat_object) <- "Fenqun"
CI.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CI.markers, file="8000 Cells mixed with normal transcriptome markers.csv")
```
