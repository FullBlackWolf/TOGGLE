---
title: "Rawdata Preprocessing: Neuron"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Preprocessing & Quality Control
---

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

Create vector to read files
---
```R
setwd("C:/GEOANALYSIS/GSE232429")
```

Read data
---
```R
Sham1<- Read10X(data.dir = "Sham1")
MCAO1<- Read10X(data.dir = "MCAO1")
MCAO2<- Read10X(data.dir = "MCAO2")
```
Create Seurat object and filter. Add code to filter out cells with fewer than 200 genes (min.features = 200) and genes covered by fewer than 3 cells (min.cells = 3)
---
```R
Sham1<- CreateSeuratObject(counts =Sham1, project = "Sham1", min.features = 200, min.cells = 3)
MCAO1<- CreateSeuratObject(counts =MCAO1, project = "MCAO1", min.features = 200, min.cells = 3)
MCAO2<- CreateSeuratObject(counts =MCAO2, project = "MCAO2", min.features = 200, min.cells = 3)

# Calculate mitochondrial DNA
Sham1[["percent.mt"]] <- PercentageFeatureSet(Sham1, pattern = "^mt-")
MCAO1[["percent.mt"]] <- PercentageFeatureSet(MCAO1, pattern = "^mt-")
MCAO2[["percent.mt"]] <- PercentageFeatureSet(MCAO2, pattern = "^mt-")

# Merge data and plot quality control (QC)
CI <- merge(Sham1, y=c(MCAO1, MCAO2), add.cell.ids = c("Sham1", "MCAO1", "MCAO2"), project = "all")
head(colnames(CI))
unique(sapply(X = strsplit(colnames(CI), split = "_"), FUN = "[", 1))

plot1 <- FeatureScatter(CI, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))#Quality control plot1
VlnPlot(CI, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size=0)#Quality control plot2
VlnPlot(CI, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size=0.5)#Quality control plot3
```

Remove cells with high mitochondrial gene expression or extreme values
---
```R
# Parameters referenced from https://cloud.tencent.com/developer/article/2195816 and https://zhuanlan.zhihu.com/p/484804392
# nFeature_RNA: Genes expressed in each cell > 300 and < 7000;
# mt_percent: Mitochondrial gene expression < 25% of total gene expression;
Sham1<- subset(Sham1, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & 
                 percent.mt < 25)
MCAO1<- subset(MCAO1,subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & 
                 percent.mt < 25)
MCAO2<- subset(MCAO2, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & 
                 percent.mt < 25)
```

Perform CCA integration
---
```R
myfunction1 <- function(testA.seu){
  testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)
  return(testA.seu)
}
Sham1<- myfunction1(Sham1)
MCAO1<- myfunction1(MCAO1)
MCAO2<- myfunction1(MCAO2)
```

Integration
---
```R
list <- list (Sham1, MCAO1, MCAO2)
testAB.anchors <- FindIntegrationAnchors(object.list = list, dims = 1:20)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
```

Add sample and group information
---
```R
# Retrieve metadata
metadata <- testAB.integrated@meta.data
# Copy 'orig.ident' to new column 'Sample'
metadata$Sample <- metadata$orig.ident
# Create new column 'Group' based on 'orig.ident'
metadata$Group <- ifelse(grepl("Sham", metadata$orig.ident), "Sham",
                         ifelse(grepl("MCAO", metadata$orig.ident), "MCAO", NA))
# Ensure updated metadata is reassigned to Seurat object
testAB.integrated@meta.data <- metadata
# Check results
head(testAB.integrated@meta.data)
```
As per documentation, use 'integrated' for finding cluster markers and 'RNA' (normalized data) for differential analysis
---
```R
# Set default matrix to 'integrated' for subsequent steps
DefaultAssay(testAB.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.1)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:10)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)
save(testAB.integrated, file = "GSE232429 Neuron.Rdata")
```


Save the file as h5ad for further analysis in Python
---
```R
library(sceasy)
# Ensure default assay is RNA
DefaultAssay(testAB.integrated) <- "RNA"
# Convert RNA assay to version 4 matrix format
testAB.integrated[["RNA"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
# Use FindVariableFeatures to select highly variable genes
testAB.integrated <- FindVariableFeatures(
  object = testAB.integrated,
  selection.method = "vst", 
  nfeatures = 2000          
)
# Check if highly variable genes were correctly selected
variable_genes <- VariableFeatures(testAB.integrated)
cat("Number of variable genes selected:", length(variable_genes), "\n")
head(variable_genes)
# Export as h5ad file, ensuring inclusion of highly variable gene information
sceasy::convertFormat(
  testAB.integrated,
  from = "seurat",
  to = "anndata",
  outFile = "GSE232429 Neuron.h5ad"
)
```

Save the file as h5ad for further analysis in Python
---
```R
library(sceasy)
# Ensure default assay is RNA
DefaultAssay(testAB.integrated) <- "RNA"
# Convert RNA assay to version 4 matrix format
testAB.integrated[["RNA"]] <- as(object = testAB.integrated[["RNA"]], Class = "Assay")
# Use FindVariableFeatures to select highly variable genes
testAB.integrated <- FindVariableFeatures(
  object = testAB.integrated,
  selection.method = "vst", 
  nfeatures = 2000          
)
# Check if highly variable genes were correctly selected
variable_genes <- VariableFeatures(testAB.integrated)
cat("Number of variable genes selected:", length(variable_genes), "\n")
head(variable_genes)
# Export as h5ad file, ensuring inclusion of highly variable gene information
sceasy::convertFormat(
  testAB.integrated,
  from = "seurat",
  to = "anndata",
  outFile = "GSE232429 Neuron.h5ad"
)
```


Visualization
---
```R
cell_type_cols <- c("#6693b1","#a3caa9","#deedad","#dd9667","#bd5c56")
p1 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "ranse", split.by = "Group", pt.size=0.5, label = T, repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 3A-1.pdf", plot = p1, device = 'pdf', width = 26, height = 14, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-1.png" 
     alt="Neuron-1.png" 
     title="Neuron-1.png">

```R
p2 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "Biaoqian", split.by = "Group", pt.size=0.5, label = T, repel = TRUE, raster=FALSE, cols = cell_type_cols) + labs(x = "UMAP1", y = "UMAP2") + theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = "Figure 3A-2.pdf", plot = p2, device = 'pdf', width = 26, height = 14, units = 'cm')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-2.png" 
     alt="Neuron-2.png" 
     title="Neuron-2.png">


Export cell proportions
---
```R
Table1 <- table(testAB.integrated$Group, testAB.integrated$ranse)
write.table(Table1, file = "Cell counts-group.txt", sep ="\t")
```

Plot cells elevated compared to MCAO group
---
```R
tb <- data.frame(table(testAB.integrated$ranse,testAB.integrated$Sample, testAB.integrated$Group))
tb=tb[,c(1,3,4)]
```

Calculate Percentages
---
```R
tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var3 == x[2],3]))
tb<- tb %>% mutate(Percentage = round(Freq/Total,3) * 100)
tb=tb[,c(1,2,5)]
tb$Var1=as.factor(tb$Var1)
tb$Var3=as.factor(tb$Var3)
head(tb)
```



Perform t-Tests
---
```R
df= do.call(rbind,
            lapply(split(tb,tb$Var1), function(x){
              # x= split(tb,tb$Var1)[[1]]
              tmp = t.test(x$Percentage ~ x$Var3)
              return(c(tmp$p.value, tmp$estimate[1]-tmp$estimate[2]))
            }))
```

Add Threshold Labels
---
```R
colnames(df) = c("pval","Difference")
df = as.data.frame(df)
df$threshold = factor(ifelse(df$Difference > 0 ,'Down','Up'))
```

Visualization
---

```R
ggplot(df,aes(x=Difference,y=-log10(pval),color=threshold))+
  geom_point()+
  geom_text_repel(
    aes(label = rownames(df)),
    size = 4,
    segment.color = "black", show.legend = FALSE ) + #add cell name
  theme_bw()+# Modify plot background
  theme(legend.title = element_blank()) +  # Hide legend title
  ylab('-log10(pval)')+  # Update Y-axis label
  xlab('Difference')+  # Update X-axis label
  geom_vline(xintercept=c(0),lty=3,col="black",lwd=0.5)
ggsave("Figure 3C.pdf",width = 5,height = 3.8)
save(testAB.integrated,file = 'GSE232429 Neuron.Rdata')
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Neuron-3.png" 
     alt="Neuron-3.png" 
     title="Neuron-3.png">

