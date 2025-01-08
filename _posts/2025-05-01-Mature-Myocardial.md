---
title: "Mature Functions: Function Difference of Myocardial"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---

```python
import cospar as cs
import numpy as np
import os
import kailin as kl
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
h5ad_filename = "心梗-成纤维细胞-高变基因.h5ad"
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

file_name = "抽样_高变_8000.h5ad"
file_path = os.path.join(current_folder_filter_1, file_name)
adata_random_subset.write(file_path)


adata_random_subset
```

```python
#选择要使用哪个样本
choosen_sample = "fibroblasts"
#选择.h5ad文件
h5ad_filename = "抽样_高变_8000.h5ad"

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
mat_name = "distance_matrix.mat"
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

csv_name = "orig_ident.csv"
csv_path = os.path.join(current_folder_result, csv_name)
print(csv_path)
df_orig_adata.to_csv(csv_path, index=False)
```
