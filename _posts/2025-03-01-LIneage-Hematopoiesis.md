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


print(kl.__version__)


kl.kl_initialize(0)


parent_directory_origin = kl.kl_settings.parent_directory_origin
print(parent_directory_origin)
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)
```

Generate similarity matrix
```python
choosen_sample = "Hematopoiesis"
h5ad_filename = "Hematopoiesis_progenitor.h5ad"


current_folder_input = current_folder
orig_adata,loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)
#orig_adata,loading_directory,distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)


#current_folder_input = current_folder
#loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input)
```


Save .csv and .mat
---

```python

save_list = ["orig_adata.obsm['X_emb']", "orig_adata.obs['label']"]
merged_csv,result_directory = kl.workcatalogue.kl_save(loading_directory,choosen_sample,distance_matrix,save_list,orig_adata)
```
The files are saved in `[LittleSnowFox's Anaconda installation directory]\database\Tracing_sample\Hematopoiesis\result\` as `merged_data.csv` and `distance_matrix.mat`.

Unsupervised learning for only progenitor cells
---

Run `[LittleSnowFox's Anaconda installation directory]\database\Tracing_sample\Hematopoiesis\main_v3_matlab_run_only_prog.m`

```matlab
clear;clc;
%% Parameters



%% Load data  and  Split to compute
%% Load data and Split to compute
MM0 = load('result/r1n30distance_matrix.mat');
%Input the Umap XY and Label
count_=readtable(['result/merged_data.csv']);
%count_=readtable('result/combined_monocle2.csv');

%computing
MM0 = MM0.distance_matrix;
MM0=MM0(contains(count_.Var4,'_prog'),contains(count_.Var4,'_prog'));
count_=count_(contains(count_.Var4,'_prog'),:);


%% MM0矩阵
%% Iterations_Number求解次数
%% Cell_Resolutio：length(M)<Cell_Resolution ，每个簇细胞数目不能小于Cell_Resolution，否则就不分了
%% Min_Wrong：split<Min_Wrong ，split是负相关性的数目，负相关性不能少于Min_Wrong个，否则就不分了
%% Min_Right：length(M)-split<5 正相关性不能少于5个，否则就不分了
%binary_corr_sorting(M,Iterations_Number,Cell_Resolution,Min_Wrong,Min_Right)
[p,splitlist] = binary_corr_sorting(MM0,20,100,5,5);
[uniqueList, ~, ~] = unique(splitlist, 'stable');
MM=MM0(p,p);
split=[];
count_result=count_(p,:);
split_simple=uniqueList;
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]
[simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM);
[acc]=acc_computing(split_simple,count_result);
[predict_result]=predict_computing(count_,simple_label,split_simple);
count_result_out=[count_result,predict_result];
writetable(count_result_out,'result/map_draw_blood.csv')

[p,splitlist] = binary_corr_sorting(MM0,3,300,5,5);
[uniqueList, ~, ~] = unique(splitlist, 'stable');
MM=MM0(p,p);
split=[];
count_result=count_(p,:);
split_simple=uniqueList;
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]
[simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM);


[predict_result]=predict_computing(count_,simple_label,split_simple);
count_result_out=[count_result,predict_result];
%[simple_label_result,simple_matrix_result]=cluster_map(simple_label,simple_matrix);
%writecell(simple_label_result,'result/Figure_c_label_blood.csv')
%writematrix(simple_matrix_result,'result/Figure_c_matrix_blood.csv')



%重排小矩阵
[cluster_map_label,cluster_map_matrix] = genetic_encoder( ...
    simple_label, ...
    simple_matrix, ...
    60, ...% nPop = 50;  % 种群规模大小为30
    1, ...% nPc = 1; % 子代规模的比例0.8
    200, ...% maxIt = 200; % 最大迭代次数
    5 ...% cycletimes = 200; % 循环计算次数
    );
%genetic_encoder(simple_label,simple_matrix,nPop,nPc,maxIt,cycletimes)
% nVar = 100; % x的长度


%重拍小矩阵方案2
% 创建行和列标签（示例）
row_labels = cluster_map_label;
column_labels = cluster_map_label;
% 使用 heatmap 函数并传递相应参数
h = heatmap(cluster_map_matrix);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0, 0.0007]

%对小矩阵进行排序
%计算pesudotime，两种计算模式，mean和median
[pesudotime_info] = pesudotime_combine(split_simple,count_.Pst,"mean",cluster_map_label)
%使用sigmoid函数处理伪时间
pesudotime_info_sigmoid = sigmoid(pesudotime_info,45,12,1000);
% 使用 heatmap 函数并传递相应参数

column_labels = pesudotime_info_sigmoid;
row_labels = cluster_map_label;
h = heatmap(cluster_map_matrix);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0, 0.0007]


% %% 处理得到时间矩阵
% cluster_map_matrix_debug = zeros(length(cluster_map_matrix),length(cluster_map_matrix));
% for itimes = 1:1:length(pesudotime_info_sigmoid)
%     cluster_map_matrix_debug(itimes,:) = pesudotime_info_sigmoid(itimes).*cluster_map_matrix(itimes,:);
%     cluster_map_matrix_debug(:,itimes) = pesudotime_info_sigmoid(itimes).*cluster_map_matrix(:,itimes);
% end
% column_labels = pesudotime_info_sigmoid;
% row_labels = cluster_map_label;
% h = heatmap(cluster_map_matrix_debug);
% h.YDisplayLabels = row_labels; % 设置行标签
% h.XDisplayLabels = column_labels; % 设置列标签
% h.ColorLimits = [0, 0.0007]
% 
% simple_label_str_result = cluster_map_label;


%% 临近法激活
corr_matrix = relevance_generate(0.00029,2,cluster_map_matrix);
hi = heatmap(corr_matrix);
hi.YDisplayLabels = row_labels; % 设置行标签
hi.XDisplayLabels = column_labels; % 设置列标签


%% 编码
encode_result = encoder_corr_matrix(0.00026,0.00030,10,3,cluster_map_matrix);
figure(2)
hj = heatmap(encode_result);
hj.YDisplayLabels = row_labels; % 设置行标签
hj.XDisplayLabels = column_labels; % 设置列标签

%% 解码
figure(3)
[weighting_decode,decode_result] = decoder_corr_matrix(encode_result);
weighting_result = weighting_decode + decode_result;
hk = heatmap(weighting_result);
hk.ColorLimits = [16,17]
hk.YDisplayLabels = row_labels; % 设置行标签
hk.XDisplayLabels = column_labels; % 设置列标签
```
Use `.\[LittleSnowFox's Anaconda installation directory]\R_processing\Hematopoiesis_prog.R` to generate the picture.



```R
library(ggplot2)

if (!exists("first_run_flag")) {
  setwd("..")
  current_dir <- getwd()
  print("Switched to the parent directory.")
  current_dir
  first_run_flag <- TRUE
} else {
  print("Not the first run, skipping setwd.")
}

print(current_dir)
database_dir <- file.path(current_dir, "database")
Tracing_dir <- file.path(database_dir, "Tracing_sample")
Hematopoiesis_dir <- file.path(Tracing_dir, "Hematopoiesis")
Hematopoiesis_result_dir <- file.path(Hematopoiesis_dir, "result")
Hematopoiesis_map <- file.path(Hematopoiesis_result_dir, "map_draw_blood.csv")

repro <- read.csv(Hematopoiesis_map)

ggplot(repro,aes(x=Var1,y=Var2,color=Var5))+geom_point()
ggplot(repro,aes(x=Var1,y=Var2,color=Var5))+geom_point()
```
<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Hematopoiesis_prog.png" 
     alt="Hematopoiesis_prog.png" 
     title="Hematopoiesis_prog.png">

Unsupervised learning for whole cells
---

```matlab
clear;clc;
%% Parameters



%% Load data  and  Split to compute
%% Load data and Split to compute
MM0 = load('result/r1n30distance_matrix.mat');
%Input the Umap XY and Label
count_=readtable(['result/merged_data.csv']);
%count_=readtable('result/combined_monocle2.csv');

%computing
MM0 = MM0.distance_matrix;
%MM0=MM0(contains(count_.Var4,'_prog'),contains(count_.Var4,'_prog'));
%count_=count_(contains(count_.Var4,'_prog'),:);


%% MM0矩阵
%% Iterations_Number求解次数
%% Cell_Resolutio：length(M)<Cell_Resolution ，每个簇细胞数目不能小于Cell_Resolution，否则就不分了
%% Min_Wrong：split<Min_Wrong ，split是负相关性的数目，负相关性不能少于Min_Wrong个，否则就不分了
%% Min_Right：length(M)-split<5 正相关性不能少于5个，否则就不分了
%binary_corr_sorting(M,Iterations_Number,Cell_Resolution,Min_Wrong,Min_Right)
[p,splitlist] = binary_corr_sorting(MM0,20,100,5,5);
[uniqueList, ~, ~] = unique(splitlist, 'stable');
MM=MM0(p,p);
split=[];
count_result=count_(p,:);
split_simple=uniqueList;
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]
[simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM);
[acc]=acc_computing(split_simple,count_result);
[predict_result]=predict_computing(count_,simple_label,split_simple);
count_result_out=[count_result,predict_result];
writetable(count_result_out,'result/all_map_blood.csv')

[p,splitlist] = binary_corr_sorting(MM0,3,300,5,5);
[uniqueList, ~, ~] = unique(splitlist, 'stable');
MM=MM0(p,p);
split=[];
count_result=count_(p,:);
split_simple=uniqueList;
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]
[simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM);


[predict_result]=predict_computing(count_,simple_label,split_simple);
count_result_out=[count_result,predict_result];
%[simple_label_result,simple_matrix_result]=cluster_map(simple_label,simple_matrix);
%writecell(simple_label_result,'result/Figure_c_label_blood.csv')
%writematrix(simple_matrix_result,'result/Figure_c_matrix_blood.csv')
%writetable(count_result_out,'result/umap_blood.csv')


%重排小矩阵
[cluster_map_label,cluster_map_matrix] = genetic_encoder( ...
    simple_label, ...
    simple_matrix, ...
    60, ...% nPop = 50;  % 种群规模大小为30
    1, ...% nPc = 1; % 子代规模的比例0.8
    200, ...% maxIt = 200; % 最大迭代次数
    5 ...% cycletimes = 200; % 循环计算次数
    );
%genetic_encoder(simple_label,simple_matrix,nPop,nPc,maxIt,cycletimes)
% nVar = 100; % x的长度


%重拍小矩阵方案2
% 创建行和列标签（示例）
row_labels = cluster_map_label;
column_labels = cluster_map_label;
% 使用 heatmap 函数并传递相应参数
h = heatmap(cluster_map_matrix);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0, 0.0007]

%对小矩阵进行排序
%计算pesudotime，两种计算模式，mean和median
[pesudotime_info] = pesudotime_combine(split_simple,count_.Pst,"mean",cluster_map_label)
%使用sigmoid函数处理伪时间
pesudotime_info_sigmoid = sigmoid(pesudotime_info,45,12,1000);
% 使用 heatmap 函数并传递相应参数
row_labels = pesudotime_info_sigmoid;
column_labels = cluster_map_label;
h = heatmap(cluster_map_matrix);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0, 0.0007]


%% 处理得到时间矩阵
cluster_map_matrix_debug = zeros(length(cluster_map_matrix),length(cluster_map_matrix));
for itimes = 1:1:length(pesudotime_info_sigmoid)
    cluster_map_matrix_debug(itimes,:) = pesudotime_info_sigmoid(itimes).*cluster_map_matrix(itimes,:);
    cluster_map_matrix_debug(:,itimes) = pesudotime_info_sigmoid(itimes).*cluster_map_matrix(:,itimes);
end
row_labels = pesudotime_info_sigmoid;
column_labels = cluster_map_label;
h = heatmap(cluster_map_matrix_debug);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0.000005,0.00001]

simple_label_str_result = cluster_map_label;


%% 临近法激活
corr_matrix = relevance_generate(0.0003,3,cluster_map_matrix);
hi = heatmap(corr_matrix);
hi.YDisplayLabels = row_labels; % 设置行标签
hi.XDisplayLabels = column_labels; % 设置列标签


%% 编码
encode_result = encoder_corr_matrix(0.00029,0.00031,10,3,cluster_map_matrix);
figure(2)
hj = heatmap(encode_result);
hj.YDisplayLabels = row_labels; % 设置行标签
hj.XDisplayLabels = column_labels; % 设置列标签

%% 解码
figure(3)
[weighting_decode,decode_result] = decoder_corr_matrix(encode_result);
weighting_result = weighting_decode + decode_result;
hk = heatmap(weighting_result);
hk.ColorLimits = [19,20]
hk.YDisplayLabels = row_labels; % 设置行标签
hk.XDisplayLabels = column_labels; % 设置列标签
```
Use `.\[LittleSnowFox's Anaconda installation directory]\R_processing\Hematopoiesis_all.R` to generate the picture.

```R
library(ggplot2)

if (!exists("first_run_flag")) {
  setwd("..")
  current_dir <- getwd()
  print("Switched to the parent directory.")
  current_dir
  first_run_flag <- TRUE
} else {
  print("Not the first run, skipping setwd.")
}

print(current_dir)
database_dir <- file.path(current_dir, "database")
Tracing_dir <- file.path(database_dir, "Tracing_sample")
Hematopoiesis_dir <- file.path(Tracing_dir, "Hematopoiesis")
Hematopoiesis_result_dir <- file.path(Hematopoiesis_dir, "result")
Hematopoiesis_map <- file.path(Hematopoiesis_result_dir, "all_map_blood.csv")

repro <- read.csv(Hematopoiesis_map)

ggplot(repro,aes(x=Var1,y=Var2,color=Var5))+geom_point()
ggplot(repro,aes(x=Var1,y=Var2,color=Var5))+geom_point()
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Hematopoiesis_all.png" 
     alt="Hematopoiesis_all.png" 
     title="Hematopoiesis_all.png">
