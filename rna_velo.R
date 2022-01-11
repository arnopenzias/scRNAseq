library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)

sars.list <- c(
  "SRR12878682.loom",
  "SRR12878683.loom",
  "SRR12878684.loom",
  "SRR12878685.loom",
  "SRR13960237.loom",
  "SRR13960238.loom",
  "SRR13960239.loom",
  "SRR13960240.loom"
)

setwd("/media/patrick/GERVAZIO/Bioinfo/neuroimuno_lab/scRNAseq-sarscov2_brain/velocyto/frontal_cortex/sars/")

ldat.sars5 <- ReadVelocity(file = "SRR13960237.loom")
bm.sars5 <- as.Seurat(x = ldat.sars5)
bm.sars5[["RNA"]] <- bm.sars5[["spliced"]]
bm.sars5 <- PercentageFeatureSet(bm.sars5, pattern = "^MT-", col.name = "percent.mt")
bm.sars5 <- subset(bm.sars5, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
bm.sars5 <- SCTransform(bm.sars5, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
bm.sars5$sars <- "SARS"




bm <- RunPCA(bm, npcs = 50, verbose = FALSE)
ElbowPlot(bm, ndims = 20, reduction = "pca") 
bm <- RunUMAP(bm, reduction = "pca", dims = 1:15, verbose = FALSE)
bm <- FindNeighbors(bm, reduction = "pca", dims = 1:15, verbose = FALSE)
bm <- FindClusters(bm, resolution = 0.5)
# save the Seurat object
saveRDS(bm, file = "./bm_control_A.rds")
# load back Seurat object
bm <- readRDS("bm_control_A.rds")
DimPlot(bm, reduction = "umap", label = TRUE)

DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "control_A_BM.h5Seurat")
Convert("control_A_BM.h5Seurat", dest = "h5ad")
