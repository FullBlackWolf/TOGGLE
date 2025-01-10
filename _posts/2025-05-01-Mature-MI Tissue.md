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

5.1.Projecting organizations onto a single-cell functional map

5.1.1.Calculate the cell similarity of the whole genome.
---


5.1.2.Read data

```python
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

5.1.3.Combine

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

