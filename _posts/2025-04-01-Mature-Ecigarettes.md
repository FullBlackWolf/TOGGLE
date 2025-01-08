---
title: "Mature Functions: RNA Pathways Changes of Ecigarettes"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---

Adult mice-Electronic cigarette
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
print(parent_directory_origin)
#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Cluster',1)
print(current_folder)
```

```python
#选择要使用哪个样本
choosen_sample = "e-cigarettes_difference expression"
#选择.h5ad文件
h5ad_filename = "成年小鼠-电子烟.xlsx"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
current_folder_input = os.path.join(current_folder_input, 'database')
current_folder_input = os.path.join(current_folder_input, 'Clustering_sample')
current_folder_input = os.path.join(current_folder_input, 'e-cigarettes_difference_expression')
current_save_path = os.path.join(current_folder_input, 'result')
current_folder_input = os.path.join(current_folder_input, 'data')
current_folder_input = os.path.join(current_folder_input, h5ad_filename)

print(current_folder_input)


```

```python
import pandas as pd
import anndata as ad
data = pd.read_excel(current_folder_input)
print(cigarette_data)
# Selecting all rows and all columns after the first in the original data
cigarette_data_X = data.iloc[:, 1:]
# Displaying the first few rows to confirm the selection
print(cigarette_data_X)
cigarette_adata = ad.AnnData(X=cigarette_data_X.values)
print(cigarette_adata)
from kailin.omic_functions import matrix                           
cigarette_matrix = matrix.distance_matrix(
    cigarette_adata,
    round_of_smooth=1,
    neighbor_N=20,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
)
```

```python
import scipy.sparse

# 假设 cigarette_matrix 是一个稀疏矩阵
cigarette_matrix_dense = cigarette_matrix.toarray()  # 转换为非稀疏 numpy 数组

print(cigarette_matrix_dense)
```
```python
import matplotlib
matplotlib.use('TkAgg')  # 或 'Qt5Agg' 等

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

# 设置后端
matplotlib.use('TkAgg')

# 绘制热图
sns.heatmap(cigarette_matrix_dense, cmap='Blues', vmin=0, vmax=0.00001)
plt.title('Distance Matrix')
plt.xlabel('Samples')
plt.ylabel('Samples')
plt.show()

```

```python
import scipy.io as sio
current_save_path_mat = os.path.join(current_save_path, 'Adult_mice-Electronic_cigarette.mat')
sio.savemat(current_save_path_mat, {'cigarette_matrix_dense': cigarette_matrix_dense})
```

```python
label = data.iloc[:, 0]
print(label)
```

```python
import os
import pandas as pd
# 提取 data 的第一列
label = data.iloc[:, 0]  # 如果 data 是 numpy 数组
# 如果 data 是 pandas DataFrame，则用 label = data.iloc[:, 0]
# 确保保存路径存在
save_folder = 'current_save_path'
if not os.path.exists(current_save_path):
    os.makedirs(current_save_path)
# 设置保存文件的完整路径
current_save_path_csv = os.path.join(current_save_path, 'Adult_mice-Electronic_cigarette.csv')
# 保存为 CSV 文件
pd.DataFrame(label).to_csv(current_save_path_csv, index=False)
print(f"文件已保存到 {current_save_path_csv}")
```


Adult_mice-Electronic_freshair
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
print(parent_directory_origin)
#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Cluster',1)
print(current_folder)

#选择要使用哪个样本
choosen_sample = "e-cigarettes_difference expression"
#选择.h5ad文件
h5ad_filename = "成年小鼠-新鲜空气.xlsx"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
current_folder_input = os.path.join(current_folder_input, 'database')
current_folder_input = os.path.join(current_folder_input, 'Clustering_sample')
current_folder_input = os.path.join(current_folder_input, 'e-cigarettes_difference_expression')
current_save_path = os.path.join(current_folder_input, 'result')
current_folder_input = os.path.join(current_folder_input, 'data')
current_folder_input = os.path.join(current_folder_input, h5ad_filename)

print(current_folder_input)

import pandas as pd
import anndata as ad
data = pd.read_excel(current_folder_input)
print(cigarette_data)
# Selecting all rows and all columns after the first in the original data
cigarette_data_X = data.iloc[:, 1:]
# Displaying the first few rows to confirm the selection
print(cigarette_data_X)
cigarette_adata = ad.AnnData(X=cigarette_data_X.values)
print(cigarette_adata)
from kailin.omic_functions import matrix                           
cigarette_matrix = matrix.distance_matrix(
    cigarette_adata,
    round_of_smooth=1,
    neighbor_N=20,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
)


import scipy.sparse

# 假设 cigarette_matrix 是一个稀疏矩阵
cigarette_matrix_dense = cigarette_matrix.toarray()  # 转换为非稀疏 numpy 数组

print(cigarette_matrix_dense)

import scipy.io as sio
current_save_path_mat = os.path.join(current_save_path, 'Adult_mice-Electronic_freshair.mat')
sio.savemat(current_save_path_mat, {'cigarette_matrix_dense': cigarette_matrix_dense})

label = data.iloc[:, 0]
print(label)

import os
import pandas as pd
# 提取 data 的第一列
label = data.iloc[:, 0]  # 如果 data 是 numpy 数组
# 如果 data 是 pandas DataFrame，则用 label = data.iloc[:, 0]
# 确保保存路径存在
save_folder = 'current_save_path'
if not os.path.exists(current_save_path):
    os.makedirs(current_save_path)
# 设置保存文件的完整路径
current_save_path_csv = os.path.join(current_save_path, 'Adult_mice-Electronic_freshair.csv')
# 保存为 CSV 文件
pd.DataFrame(label).to_csv(current_save_path_csv, index=False)
print(f"文件已保存到 {current_save_path_csv}")
```

Adult mice-e-cigarette-only e-cigarette and air-only difference in genes
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
print(parent_directory_origin)
#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Cluster',1)
print(current_folder)

#选择要使用哪个样本
choosen_sample = "e-cigarettes_difference expression"
#选择.h5ad文件
h5ad_filename = "成年小鼠-只电子烟和空气.xlsx"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
current_folder_input = os.path.join(current_folder_input, 'database')
current_folder_input = os.path.join(current_folder_input, 'Clustering_sample')
current_folder_input = os.path.join(current_folder_input, 'e-cigarettes_difference_expression')
current_save_path = os.path.join(current_folder_input, 'result')
current_folder_input = os.path.join(current_folder_input, 'data')
current_folder_input = os.path.join(current_folder_input, h5ad_filename)

print(current_folder_input)

import pandas as pd
import anndata as ad
data = pd.read_excel(current_folder_input)
print(data)
# Selecting all rows and all columns after the first in the original data
cigarette_data_X = data.iloc[:, 1:]
# Displaying the first few rows to confirm the selection
print(cigarette_data_X)
cigarette_adata = ad.AnnData(X=cigarette_data_X.values)
print(cigarette_adata)
from kailin.omic_functions import matrix                           
cigarette_matrix = matrix.distance_matrix(
    cigarette_adata,
    round_of_smooth=1,
    neighbor_N=20,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
)


import scipy.sparse

# 假设 cigarette_matrix 是一个稀疏矩阵
cigarette_matrix_dense = cigarette_matrix.toarray()  # 转换为非稀疏 numpy 数组

print(cigarette_matrix_dense)

import scipy.io as sio
current_save_path_mat = os.path.join(current_save_path, 'Adult mice - e-cigarette and air only.mat')
sio.savemat(current_save_path_mat, {'cigarette_matrix_dense': cigarette_matrix_dense})

label = data.iloc[:, 0]
print(label)

import os
import pandas as pd
# 提取 data 的第一列
label = data.iloc[:, 0]  # 如果 data 是 numpy 数组
# 如果 data 是 pandas DataFrame，则用 label = data.iloc[:, 0]
# 确保保存路径存在
save_folder = 'current_save_path'
if not os.path.exists(current_save_path):
    os.makedirs(current_save_path)
# 设置保存文件的完整路径
current_save_path_csv = os.path.join(current_save_path, 'Adult mice - e-cigarette and air only.csv')
# 保存为 CSV 文件
pd.DataFrame(label).to_csv(current_save_path_csv, index=False)
print(f"文件已保存到 {current_save_path_csv}")
```

Fresh air in children vs adults - Child mice - E-cigarettes
---
Child mice - aVc_c
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
print(parent_directory_origin)
#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Cluster',1)
print(current_folder)

#选择要使用哪个样本
choosen_sample = "e-cigarettes_difference expression"
#选择.h5ad文件
h5ad_filename = "幼年vs成年的新鲜空气-幼年小鼠-电子烟.xlsx"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
current_folder_input = os.path.join(current_folder_input, 'database')
current_folder_input = os.path.join(current_folder_input, 'Clustering_sample')
current_folder_input = os.path.join(current_folder_input, 'e-cigarettes_difference_expression')
current_save_path = os.path.join(current_folder_input, 'result')
current_folder_input = os.path.join(current_folder_input, 'data')
current_folder_input = os.path.join(current_folder_input, h5ad_filename)

print(current_folder_input)

import pandas as pd
import anndata as ad
data = pd.read_excel(current_folder_input)
print(data)
# Selecting all rows and all columns after the first in the original data
cigarette_data_X = data.iloc[:, 1:]
# Displaying the first few rows to confirm the selection
print(cigarette_data_X)
cigarette_adata = ad.AnnData(X=cigarette_data_X.values)
print(cigarette_adata)
from kailin.omic_functions import matrix                           
cigarette_matrix = matrix.distance_matrix(
    cigarette_adata,
    round_of_smooth=1,
    neighbor_N=2000,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
)


import scipy.sparse

# 假设 cigarette_matrix 是一个稀疏矩阵
cigarette_matrix_dense = cigarette_matrix.toarray()  # 转换为非稀疏 numpy 数组

print(cigarette_matrix_dense)

import scipy.io as sio
current_save_path_mat = os.path.join(current_save_path, 'Child mice - aVc_c.mat')
sio.savemat(current_save_path_mat, {'cigarette_matrix_dense': cigarette_matrix_dense})

label = data.iloc[:, 0]
print(label)

import os
import pandas as pd
# 提取 data 的第一列
label = data.iloc[:, 0]  # 如果 data 是 numpy 数组
# 如果 data 是 pandas DataFrame，则用 label = data.iloc[:, 0]
# 确保保存路径存在
save_folder = 'current_save_path'
if not os.path.exists(current_save_path):
    os.makedirs(current_save_path)
# 设置保存文件的完整路径
current_save_path_csv = os.path.join(current_save_path, 'Child mice - aVc_c.csv')
# 保存为 CSV 文件
pd.DataFrame(label).to_csv(current_save_path_csv, index=False)
print(f"文件已保存到 {current_save_path_csv}")
```


Child vs adult electronic cigarettes - Child mice - electronic cigarettes
---
Child mice - c_c
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
print(parent_directory_origin)
#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Cluster',1)
print(current_folder)

#选择要使用哪个样本
choosen_sample = "e-cigarettes_difference expression"
#选择.h5ad文件
h5ad_filename = "幼年vs成年的电子烟-幼年小鼠-电子烟.xlsx"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
current_folder_input = os.path.join(current_folder_input, 'database')
current_folder_input = os.path.join(current_folder_input, 'Clustering_sample')
current_folder_input = os.path.join(current_folder_input, 'e-cigarettes_difference_expression')
current_save_path = os.path.join(current_folder_input, 'result')
current_folder_input = os.path.join(current_folder_input, 'data')
current_folder_input = os.path.join(current_folder_input, h5ad_filename)

print(current_folder_input)

import pandas as pd
import anndata as ad
data = pd.read_excel(current_folder_input)
print(data)
# Selecting all rows and all columns after the first in the original data
cigarette_data_X = data.iloc[:, 1:]
# Displaying the first few rows to confirm the selection
print(cigarette_data_X)
cigarette_adata = ad.AnnData(X=cigarette_data_X.values)
print(cigarette_adata)
from kailin.omic_functions import matrix                           
cigarette_matrix = matrix.distance_matrix(
    cigarette_adata,
    round_of_smooth=1,
    neighbor_N=2000,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
)


import scipy.sparse

# 假设 cigarette_matrix 是一个稀疏矩阵
cigarette_matrix_dense = cigarette_matrix.toarray()  # 转换为非稀疏 numpy 数组

print(cigarette_matrix_dense)

import scipy.io as sio
current_save_path_mat = os.path.join(current_save_path, 'Child mice - c_c.mat')
sio.savemat(current_save_path_mat, {'cigarette_matrix_dense': cigarette_matrix_dense})

label = data.iloc[:, 0]
print(label)

import os
import pandas as pd
# 提取 data 的第一列
label = data.iloc[:, 0]  # 如果 data 是 numpy 数组
# 如果 data 是 pandas DataFrame，则用 label = data.iloc[:, 0]
# 确保保存路径存在
save_folder = 'current_save_path'
if not os.path.exists(current_save_path):
    os.makedirs(current_save_path)
# 设置保存文件的完整路径
current_save_path_csv = os.path.join(current_save_path, 'Child mice - c_c.csv')
# 保存为 CSV 文件
pd.DataFrame(label).to_csv(current_save_path_csv, index=False)
print(f"文件已保存到 {current_save_path_csv}")
```

Electronic cigarettes vs. fresh air - young mice - electronic cigarettes
---
aVc_c
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
print(parent_directory_origin)
#改进：
#添加一个cluster模式
#选择进行Lineage Tracing还是Cluster，并给出可用的列表
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Cluster',1)
print(current_folder)

#选择要使用哪个样本
choosen_sample = "e-cigarettes_difference expression"
#选择.h5ad文件
h5ad_filename = "电子烟vs新鲜空气-幼年小鼠-电子烟.xlsx"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
current_folder_input = os.path.join(current_folder_input, 'database')
current_folder_input = os.path.join(current_folder_input, 'Clustering_sample')
current_folder_input = os.path.join(current_folder_input, 'e-cigarettes_difference_expression')
current_save_path = os.path.join(current_folder_input, 'result')
current_folder_input = os.path.join(current_folder_input, 'data')
current_folder_input = os.path.join(current_folder_input, h5ad_filename)

print(current_folder_input)

import pandas as pd
import anndata as ad
data = pd.read_excel(current_folder_input)
print(data)
# Selecting all rows and all columns after the first in the original data
cigarette_data_X = data.iloc[:, 1:]
# Displaying the first few rows to confirm the selection
print(cigarette_data_X)
cigarette_adata = ad.AnnData(X=cigarette_data_X.values)
print(cigarette_adata)
from kailin.omic_functions import matrix                           
cigarette_matrix = matrix.distance_matrix(
    cigarette_adata,
    round_of_smooth=1,
    neighbor_N=2000,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
)


import scipy.sparse

# 假设 cigarette_matrix 是一个稀疏矩阵
cigarette_matrix_dense = cigarette_matrix.toarray()  # 转换为非稀疏 numpy 数组

print(cigarette_matrix_dense)

import scipy.io as sio
current_save_path_mat = os.path.join(current_save_path, 'aVc_c.mat')
sio.savemat(current_save_path_mat, {'cigarette_matrix_dense': cigarette_matrix_dense})

label = data.iloc[:, 0]
print(label)

import os
import pandas as pd
# 提取 data 的第一列
label = data.iloc[:, 0]  # 如果 data 是 numpy 数组
# 如果 data 是 pandas DataFrame，则用 label = data.iloc[:, 0]
# 确保保存路径存在
save_folder = 'current_save_path'
if not os.path.exists(current_save_path):
    os.makedirs(current_save_path)
# 设置保存文件的完整路径
current_save_path_csv = os.path.join(current_save_path, 'aVc_c.csv')
# 保存为 CSV 文件
pd.DataFrame(label).to_csv(current_save_path_csv, index=False)
print(f"文件已保存到 {current_save_path_csv}")
```
