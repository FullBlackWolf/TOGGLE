---
title: "Mature Functions: Programmed Ferroptosis"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---

Identify cells undergoing programmed ferroptosis.
---

```Python
import cospar as cs
import numpy as np
import os
import kailin as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)
# Initialization function to switch Kailin to the working directory. 
# If it has been initialized before, running def kl_initialize(0) again 
# will reject reinitialization to avoid recursion. 
# Running def kl_initialize(1) will forcibly reinitialize.
kl.kl_initialize(0)

# Retrieve the root directory where Kailin works
parent_directory_origin = kl.kl_settings.parent_directory_origin
print(parent_directory_origin)

current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)

# Select which sample to use
choosen_sample = "Nerveferroptosis"
# Select the .h5ad file
h5ad_filename = "GSE232429deatd_粗糙过滤_testAB.integrated"
# Run the built-in example and obtain the sparse matrix
# Here, a non-example function needs to be included
current_folder_input = current_folder
orig_adata, loading_directory, distance_matrix = kl.preprocessing.kl_dense_matrix(
    choosen_sample, 
    h5ad_filename, 
    "draw", 
    current_folder_input, 
    1, 
    13000, 
    0.1, 
    0.001, 
    True
)

print(loading_directory)
print(choosen_sample)
```

Load required packages
---
```R
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(future)
library(RColorBrewer)
library(SummarizedExperiment)
library(scater)
library(cowplot)
library(harmony)
library(monocle3)
library(uwot)
library(ComplexHeatmap)
library(ggrepel)
```

DEG analysis for GSE232429
---

```R
library(SingleCellExperiment)
library(DEsingle)
testAB.integrated[["RNA"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
```

Set active.ident to ranse
---

```R
Idents(testAB.integrated) <- "ranse"
DefaultAssay(testAB.integrated) <- "RNA"
```

Perform pairwise comparisons
---

```R
s0 <- subset(testAB.integrated,idents=c("Group R1-1", "Group R1-2"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-1 vs Group R1-2.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-1", "Group R1-3"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-1 vs Group R1-3.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-1", "Group R1-4"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-1 vs Group R1-4.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-1", "Group R1-5"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-1 vs Group R1-5.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-2", "Group R1-3"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-2 vs Group R1-3.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-2", "Group R1-4"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-2 vs Group R1-4.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-2", "Group R1-5"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-2 vs Group R1-5.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-3", "Group R1-4"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-3 vs Group R1-4.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-3", "Group R1-5"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-3 vs Group R1-5.csv")

s0 <- subset(testAB.integrated,idents=c("Group R1-4", "Group R1-5"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$ranse)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R1-4 vs Group R1-5.csv")
```


Based on all the results, groups 1, 2, and 5 were included in the subsequent analysis
---

```R
load("C:/GEOANALYSIS/GSE232429/GSE232429 Neuron.Rdata")
Idents(testAB.integrated) <- "Biaoqian"
testAB.integrated <- subset(testAB.integrated,idents=c("R1-1","R1-2","R1-5"),invert = FALSE)
cell_type_cols <- c("#6693b1","#a3caa9","#bd5c56")
p2 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "ranse", split.by = "Group", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 3E-1.pdf", plot = p2, device = 'pdf', width = 26, height = 14, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-4.png" 
     alt="Neuron-4.png" 
     title="Neuron-4.png">


Visualization
---

```R
p3 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Biaoqian", split.by = "Group", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 3E-2.pdf", plot = p3, device = 'pdf', width = 26, height = 14, units = 'cm')
save(testAB.integrated,file = 'GSE232429 after removing 3 and 4.Rdata')
```
<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-5.png" 
     alt="Neuron-5.png" 
     title="Neuron-5.png">


Distinguish the stages of ferroptosis in cells. Due to shared RNA pathways, many cells undergoing apoptosis are mixed in.
---

```Python
import cospar as cs
import numpy as np
import os
import kailin as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)
# Initialization function to switch Kailin to the working directory. 
# If it has been initialized before, running def kl_initialize(0) again 
# will reject reinitialization to avoid recursion. 
# Running def kl_initialize(1) will forcibly reinitialize.
kl.kl_initialize(0)

# Retrieve the root directory where Kailin works
parent_directory_origin = kl.kl_settings.parent_directory_origin
print(parent_directory_origin)

current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)

# Select which sample to use
choosen_sample = "Nerveferroptosis_15_21"
# Select the .h5ad file
h5ad_filename = "2024.10.28的15-21数据.h5ad"
# Run the built-in example and obtain the sparse matrix
# Here, a non-example function needs to be included
current_folder_input = current_folder
orig_adata, loading_directory, distance_matrix = kl.preprocessing.kl_dense_matrix(
    choosen_sample, 
    h5ad_filename, 
    "draw", 
    current_folder_input, 
    1, 
    13000, 
    0.1, 
    0.001, 
    True
)
# orig_adata, loading_directory, distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(
#     choosen_sample, 
#     h5ad_filename, 
#     "draw", 
#     current_folder_input
# )
# Run the built-in example and obtain the non-sparse matrix
# Here, a non-example function needs to be included
# current_folder_input = current_folder
# loading_directory, distance_matrix = kl.preprocessing.kl_dense_matrix(
#     choosen_sample, 
#     h5ad_filename, 
#     "draw", 
#     current_folder_input
# )
print(loading_directory)
print(choosen_sample)
```


Perform pseudo-time inference on groups 1, 2, and 5
---

```R
load("C:/GEOANALYSIS/GSE232429/GSE232429 after removing 3 and 4.Rdata")
testAB.integrated[["RNA4"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
#re-do umap
UMAPPlot(testAB.integrated,group.by='ranse',label=T)
testAB.integrated = SCTransform(testAB.integrated,assay = 'RNA')
testAB.integrated = RunPCA(testAB.integrated)
ElbowPlot(testAB.integrated)
testAB.integrated = RunUMAP(testAB.integrated,dims = 1:10)
UMAPPlot(testAB.integrated,group.by='ranse',label=T)
```

Extract the matrix to make pseudo time
---

```R
DefaultAssay(testAB.integrated) <- "RNA4"
data <- as(as.matrix(testAB.integrated@assays$RNA4@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame',
          data = testAB.integrated@meta.data)
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
pd <- as(pd, "data.frame")
fd <- as(fd, "data.frame")
#Start make pseudo time
cds <- new_cell_data_set(data,
                         cell_metadata = pd
                         ,gene_metadata = fd)
cds <- preprocess_cds(cds, num_dim = 10)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, resolution = 0.000001)
plot_cells(cds, color_cells_by = "ranse",show_trajectory_graph = F,group_label_size=12)
#Replace the original umap data
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,1] = testAB.integrated@reductions$umap@cell.embeddings[,1]
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,2] = testAB.integrated@reductions$umap@cell.embeddings[,2]
#Continue to the next step
cds <- learn_graph(cds)
cds <- order_cells(cds)
p3 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = T,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 1.5) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 4A.pdf", plot = p3, device = 'pdf', width = 18, height = 16, units = 'cm')
save(testAB.integrated,file = 'GSE232429 after removing 3 and 4.Rdata')
```



<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-6.png" 
     alt="Neuron-6.png" 
     title="Neuron-6.png">

     
Regroup 1, 2, and 5
---

```R
testAB.integrated=get(load(file = 'GSE232429 after removing 3 and 4.Rdata'))
#Import the new grouping results
# Read CSV file
result_data <- read.csv("pseudotime_map.csv", stringsAsFactors = FALSE)
# Ensure the file contains columns 'Var1' and 'Result'; check file content
head(result_data)
# Check if all 'Var1' values exist in Seurat object's cell names
common_cells <- intersect(result_data$Var1, rownames(testAB.integrated@meta.data))
if (length(common_cells) < nrow(result_data)) {
  warning("Some cells in 'pseudotime_map.csv' are not found in testAB.integrated metadata!")
}
```

Map 'Result' values to Seurat object's metadata based on 'Var1'
---

```R
# First, create a new column 'Result' and set it to NA
testAB.integrated@meta.data$Result <- NA
# Use match() to merge corresponding values
matching_indices <- match(rownames(testAB.integrated@meta.data), result_data$Var1)
testAB.integrated@meta.data$Result <- result_data$Result[matching_indices]
## Group based on 'Result'
# get metadata
metadata <- testAB.integrated@meta.data
# Create a new column shijian based on the value of the Result column
metadata$shijian <- with(metadata, 
                     ifelse(Result >= 1 & Result <= 14, "Group R2-1",
                     ifelse(Result >= 15 & Result <= 18, "Group R2-2",
                      ifelse(Result >= 19 & Result <= 21, "Group R2-3",
                      ifelse(Result == 22, "Group R2-4",
                      ifelse(Result == 23, "Group R2-5",
                      ifelse(Result == 24, "Group R2-6",
                      ifelse(Result == 25, "Group R2-7",
                      ifelse(Result >= 26 & Result <= 44, "Group R2-8",
                      ifelse(Result == 45, "Group R2-9", NA)))))))))
)
```

Build the Biaoqian column, remove the "Group" in shijian
---

```R
metadata$Biaoqian <- gsub("^Group ", "", metadata$shijian)
```

Assign updated metadata back to the Seurat object
---

```R
testAB.integrated@meta.data <- metadata
```

Check results
---

```R
head(testAB.integrated@meta.data)
#Save
save(testAB.integrated,file = 'GSE232429 after removing 3 and 4.Rdata')
```


Plotting
---

```R
#Figure 4C
cell_type_cols <- c("#5a5098","#6693b1","#a3caa9","#deedad","#ffffcc","#efd695","#dd9667","#bd5c56","#842844")
p4 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "shijian", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 4C-1.pdf", plot = p4, device = 'pdf', width = 21, height = 18, units = 'cm')
```
     
<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-7.png" 
     alt="Neuron-7.png" 
     title="Neuron-7.png">

Visualization
---

```R
p5 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Biaoqian", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 4C-2.pdf", plot = p5, device = 'pdf', width = 21, height = 18, units = 'cm')
```


<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-8.png" 
     alt="Neuron-8.png" 
     title="Neuron-8.png">

Export cell proportions
---

```R

Table1 <- table(testAB.integrated$newresults)
write.table(Table1, file = "The number of cells per cluster.txt", sep ="\t")
cell_type_cols <- c("#5a5098","#6693b1","#a3caa9","#deedad","#ffffcc","#efd695","#dd9667","#bd5c56","#842844")
p5 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "shijian", split.by = "Group", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 4E-1.pdf", plot = p5, device = 'pdf', width = 26, height = 14, units = 'cm')
```



<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-9.png" 
     alt="Neuron-9.png" 
     title="Neuron-9.png">

Visualization
---

```R
p6 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Biaoqian", split.by = "Group", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 4E-2.pdf", plot = p6, device = 'pdf', width = 26, height = 14, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-10.png" 
     alt="Neuron-10.png" 
     title="Neuron-10.png">

Group Difference
---

```R
Table1 <- table(testAB.integrated$Group, testAB.integrated$shijian2)
write.table(Table1, file = "The cell proportion of the new grouping result-group.txt", sep ="\t")
# Plot cells elevated compared to MCAO group
tb <- data.frame(table(testAB.integrated$shijian2,testAB.integrated$Sample, testAB.integrated$Group))
tb=tb[,c(1,3,4)]

tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var3 == x[2],3]))
tb<- tb %>% mutate(Percentage = round(Freq/Total,3) * 100)
tb=tb[,c(1,2,5)]
tb$Var1=as.factor(tb$Var1)
tb$Var3=as.factor(tb$Var3)
head(tb)
df= do.call(rbind,
            lapply(split(tb,tb$Var1), function(x){
              # x= split(tb,tb$Var1)[[1]]
              tmp = t.test(x$Percentage ~ x$Var3)
              return(c(tmp$p.value, tmp$estimate[1]-tmp$estimate[2]))
            }))

colnames(df) = c("pval","Difference")
df = as.data.frame(df)
df$threshold = factor(ifelse(df$Difference > 0 ,'Down','Up'))

ggplot(df,aes(x=Difference,y=-log10(pval),color=threshold))+
  geom_point()+
  geom_text_repel(
    aes(label = rownames(df)),
    size = 4,
    segment.color = "black", show.legend = FALSE ) +
  theme_bw()+# Modify plot background
  theme(legend.title = element_blank()) +  # Hide legend title
  ylab('-log10(pval)')+  # Update Y-axis label
  xlab('Difference')+  # Update X-axis label
  geom_vline(xintercept=c(0),lty=3,col="black",lwd=0.5)
ggsave("Figure 4G.pdf",width = 5,height = 3.8)
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-11.png" 
     alt="Neuron-11.png" 
     title="Neuron-11.png">
     

Do DEG analysis for regrouping results
---

```R
testAB.integrated[["RNA"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
#set active.ident to shijian
Idents(testAB.integrated) <- "shijian"
DefaultAssay(testAB.integrated) <- "RNA"
#Do DEG analysis
s0 <- subset(testAB.integrated,idents=c("Group R2-1", "Group R2-5"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-1 vs Group R2-5.csv")

s0 <- subset(testAB.integrated,idents=c("Group R2-1", "Group R2-6"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-1 vs Group R2-6.csv")

s0 <- subset(testAB.integrated,idents=c("Group R2-1", "Group R2-7"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-1 vs Group R2-7.csv")

s0 <- subset(testAB.integrated,idents=c("Group R2-3", "Group R2-5"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-3 vs Group R2-5.csv")

s0 <- subset(testAB.integrated,idents=c("Group R2-3", "Group R2-6"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-3 vs Group R2-6.csv")

s0 <- subset(testAB.integrated,idents=c("Group R2-3", "Group R2-7"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-3 vs Group R2-7.csv")

s0 <- subset(testAB.integrated,idents=c("Group R2-2", "Group R2-4"),invert = FALSE)
s0 <- as.SingleCellExperiment(s0)
group0 <- factor(s0$shijian2)
results0 <- DEsingle(counts = s0, group = group0, parallel = TRUE)
write.csv(results0, file="Group R2-2 vs Group R2-4.csv")
```

After determining the head and tail, we will find the differential genes of the head and tail cell groups
---

```R
testAB.integrated=get(load(file = 'GSE232429 after removing 3 and 4.Rdata'))
DefaultAssay(testAB.integrated) <- "RNA" 
testAB.integrated <- JoinLayers(testAB.integrated)
Idents(testAB.integrated) <- "shijian"
chayi1 <- FindMarkers(testAB.integrated, ident.1 = "Group R2-3", ident.2 = "Group R2-9",
                      assay = "RNA", slot="counts",
                      only.pos = F, min.pct = 0, logfc.threshold = 0)
write.csv(chayi1, file="Differential genes between Group R2-3 and Group R2-9.csv")

```





Verify the ferroptosis ratio
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
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)

#选择要使用哪个样本
choosen_sample = "Nerveferroptosis_19_21_Group8"
#选择.h5ad文件
h5ad_filename = "重画矩阵GSE232429_testAB.integrated.h5add"
#运行自带的示例，并获取稀疏矩阵
#这里需要做非示例的函数进去
current_folder_input = current_folder
updated_folder = os.path.join(current_folder, "Nerveferroptosis_19_21_Group8/data")
h5ad_path = os.path.join(updated_folder, "2024.10.30-有Health和GROUP8细分群的数据.h5ad")

print(h5ad_path)

```
```python
import anndata as ad
adata = ad.read_h5ad(h5ad_path)
adata.obs['nFeature_RNA']
adata.obs['fenqun1']
adata_Health_RNA = adata[adata.obs['fenqun1'] == 'Health'].var.index.tolist()
# 获取特定层的形状
length_of_adata_Healt_RNA = len(adata_Health_RNA)
print(adata.X.shape)
print(length_of_adata_Healt_RNA)

print(tarc_column)
adata_Death_RNA = adata[adata.obs['fenqun1'] != 'Health'].var.index.tolist()
# 获取特定层的形状
length_of_adata_Death_RNA= len(adata_Death_RNA)
print(adata.X.shape)
print(length_of_adata_Death_RNA)

#
import pandas as pd

# 初始化 check_list 和一个空的 DataFrame
check_list = ['Smad7', 
              'Pex2', 
              'Far1', 
              'Mtch1', 
              'Lpin1', 
              'Nras', 
              'Agps', 
              'Wipi1', 
              'Hmgb1', 
              'Mapk3', 
              'Cd82', 
              'Elovl5', 
              'Scp2', 
              'Lgmn', 
              'Adam23', 
              'Emc2', 
              'Ulk2', 
              'Hddc3', 
              'Gstz1', 
              'Map3k11', 
              'Cirbp']

df = pd.DataFrame(columns=['Value', 'Average', 'Rate'])  # 定义一个空的 DataFrame

# 读取数据有多少个
nras_expression = adata[adata.obs['fenqun1'] != 'Health'][:, check_list[1]].X
total_num = adata[adata.obs['fenqun1'] != 'Health'].X.shape[0]
print('表达数据的行数为:', total_num)



num_rows = adata[adata.obs['fenqun1'] != 'Health'].n_obs
# 打印行数
print('adata_Death_RNA 的行数为:', num_rows)
#读取adata中

num_rows_1 = adata[adata.obs['fenqun1'] == 'Health'].n_obs
# 打印行数
print('adata_Health_RNA 的行数为:', num_rows_1)
#读取adata中


# 循环遍历 check_list，将每个值保存到 DataFrame 中
for i, value in enumerate(check_list):
    df.loc[i, 'Value'] = value  # 只给 'Value' 列赋值
    
    #取check_list中的Gene
    Chosen_computing = check_list[i]
    #print(Chosen_computing)
    
    #选取健康细胞作为参照值
    #adata[adata.obs['fenqun1'] != 'Health']
    tarc_column = adata[:, adata.var.index == Chosen_computing].X
    
    #取出平均值
    tarc_mean = np.median(tarc_column.toarray())
    df.loc[i, 'Average'] = tarc_mean
    
    #把要计算的基因样本抽出来
    compare_expression = adata[adata.obs['fenqun1'] != 'Health'][:, check_list[1]].X.toarray()
    #print('compare_expression:',compare_expression )
    average_standa = df.loc[i, 'Average']
    
    #计算比较总数
    count_number = sum(value >= average_standa for value in compare_expression)
    #count_number = sum(compare_expression >= df.loc[i, 'Average'])
    #print('Sum:',count_number)
    
    #计算比率
    rate = count_number/total_num
    
    #保存比率
    df.loc[i, 'Rate'] = rate
    


# 打印最终的 DataFrame
print(df)

```

















Distinguish between ferroptosis and apoptosis
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
current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)


#选择要使用哪个样本
choosen_sample = "Nerveferroptosis_remove_R1_3_4"
#选择.h5ad文件
h5ad_filename = "重画矩阵GSE232429_testAB.integrated.h5ad"
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

orig_adata.obs['shijian2']


#需要区分dense和sparase
save_list = ["orig_adata.obsm['X_umap']", "orig_adata.obs['shijian2']"]

#将要计算的文件保存到/result
merged_csv,result_directory = kl.workcatalogue.kl_save(loading_directory,choosen_sample,distance_matrix,save_list,orig_adata)

```



Cells from Group R2-2 and Group R2-3 of the MCAO group were taken for further analysis
---

```R
testAB.integrated=get(load(file = 'GSE232429 after removing 3 and 4.Rdata'))
Idents(testAB.integrated) <- "Biaoqian"
testAB.integrated <- subset(testAB.integrated,idents=c("R2-2","R2-3"),invert = FALSE)
Idents(testAB.integrated) <- "Group"
testAB.integrated <- subset(testAB.integrated,idents=c("MCAO"),invert = FALSE)
DefaultAssay(testAB.integrated) <- "RNA4"
###First generate h5ad for further analysis###
# Make sure you select a matrix that contains all genes
sceasy::convertFormat(
  testAB.integrated,
  from = "seurat",
  to = "anndata",
  outFile = "Cells from Group R2-2 and Group R2-3.h5ad"
)
####Then I got 15-21-result.csv, so I imported it
# Read CSV file
result_data <- read.csv("15-21-result.csv", stringsAsFactors = FALSE)
# Make sure the columns in the file are named "Var3" and "Result", check the file contents
head(result_data)
# Check if all 'Var1' values exist in Seurat object's cell names
common_cells <- intersect(result_data$Var3, rownames(testAB.integrated@meta.data))
if (length(common_cells) < nrow(result_data)) {
  warning("Some cells in '15-21-result.csv' are not found in testAB.integrated metadata!")
}
# Map 'Result' values to Seurat object's metadata based on 'Var1'
# First, create a new column 'Result' and set it to NA
testAB.integrated@meta.data$Result <- NA
# Use match() to merge corresponding values
matching_indices <- match(rownames(testAB.integrated@meta.data), result_data$Var3)
testAB.integrated@meta.data$Result <- result_data$Result[matching_indices]
## Group based on 'Result'
# get metadata
metadata <- testAB.integrated@meta.data
# Create a new column fenqun1 based on the value of the Result column
metadata$fenqun1 <- with(metadata, 
                         ifelse(Result >= 1 & Result <= 6, "Group R3-1",
                                ifelse(Result >= 7 & Result <= 12, "Group R3-2",
                                       ifelse(Result >= 13 & Result <= 19, "Group R3-3",
                                              ifelse(Result >= 20 & Result <= 27, "Group R3-4",
                                                     ifelse(Result >= 28 & Result <= 34, "Group R3-5", NA))))))
# Build the Biaoqian column and remove "Group" from fenqun1
metadata$Biaoqian <- gsub("^Group ", "", metadata$fenqun1)

# Assign updated metadata back to the Seurat object
testAB.integrated@meta.data <- metadata
# Check results
head(testAB.integrated@meta.data)
#Save
save(testAB.integrated,file = 'Cells from Group R2-2 and Group R2-3.Rdata')
```

Drawing
---

```R
cell_type_cols <- c("#6693b1","#a3caa9","#efd695","#dd9667","#bd5c56")
testAB.integrated = RunPCA(testAB.integrated)
ElbowPlot(testAB.integrated)
testAB.integrated = RunUMAP(testAB.integrated,dims = 1:20)
p5 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "fenqun1", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 6A-1.pdf", plot = p5, device = 'pdf', width = 15, height = 12, units = 'cm')

```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-12.png" 
     alt="Neuron-12.png" 
     title="Neuron-12.png">

Visualization
---

```R
p6 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Biaoqian", pt.size=0.5, label = T,repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 6A-2.pdf", plot = p6, device = 'pdf', width = 15, height = 12, units = 'cm')

```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-13.png" 
     alt="Neuron-13.png" 
     title="Neuron-13.png">
