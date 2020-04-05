library(Seurat)
library(dplyr)
library(cowplot)

raw_data <- read.table(file = "D:/Sc-RNA Seq/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", header = T, sep="\t", as.is=T)
bc.data <- raw_data[,18:352]

rownames(bc.data) = raw_data[,1]
bc <- CreateSeuratObject(counts = bc.data, min.cells = 3, min.features  = 200, project = "SC_RNA", assay = "RNA")

bc <- FindVariableFeatures(object = bc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = bc))
mito.genes <- grep(pattern = "^MT-", x = rownames(bc@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(bc@assays[["RNA"]][mito.genes, ])/Matrix::colSums(bc@assays[["RNA"]])

bc <- AddMetaData(object = bc, metadata = percent.mito, col.name = "percent.mito") 
#in case the above function does not work simply do:
bc$percent.mito <- percent.mito

bc <- ScaleData(object = bc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
bc <- RunPCA(object = bc,  npcs = 30, verbose = FALSE)

DimPlot(object = bc, reduction = "pca")

#Clustering
bc <- FindNeighbors(bc, reduction = "pca", dims = 1:20)
bc <- FindClusters(bc, resolution = 0.5, algorithm = 1)

#Using TSNE
bc <- RunTSNE(object = bc, dims.use = 4:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = bc, reduction = "tsne")

# Running UMap
bc <- RunUMAP(bc, reduction = "pca", dims = 1:20)
DimPlot(bc, reduction = "umap", split.by = "seurat_clusters")
FeaturePlot(object = bc, features = c("ENSG00000107796.8","ENSG00000078098.9","ENSG00000150093.14","ENSG00000105974.7","ENSG00000041982.10","ENSG00000162493.12"), cols = c("green", "blue"))
