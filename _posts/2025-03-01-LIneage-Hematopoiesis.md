---
title: "LIneage Tracing: Hematopoiesis"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---

Import external libraries.
---

```python
import numpy as np
import os
import kailin as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)


kl.kl_initialize(0)

#Get the work path
parent_directory_origin = kl.kl_settings.parent_directory_origin
print(parent_directory_origin)
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)
```

Generate similarity matrix
```python
#选择要使用哪个样本
choosen_sample = "Hematopoiesis"

#选择.h5ad文件
h5ad_filename = "Hematopoiesis_progenitor.h5ad"


#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
orig_adata,loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)
#orig_adata,loading_directory,distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)

#运行自带的示例，并获取非稀疏矩阵
#这里需要做非示例的函数进去
#current_folder_input = current_folder
#loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input)
```


Save .csv and .mat

```python

#需要区分dense和sparase
save_list = ["orig_adata.obsm['X_emb']", "orig_adata.obs['label']"]

#将要计算的文件保存到/result
merged_csv,result_directory = kl.workcatalogue.kl_save(loading_directory,choosen_sample,distance_matrix,save_list,orig_adata)
```


Unsupervised learning for only progenitor cells
