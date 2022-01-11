#################        Single-cell  RNAseq Seurat   #############################

# calling packages
library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(cowplot)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(Matrix)
library(celldex)
library(scater)
library(ensembldb)
library(car)
library(DESeq2)
library(DoubletFinder)

# Set working dir
setwd("/path/to/working/dir/")

# Load the dataset for each replicate
sarsA.data <- Read10X(data.dir = "/path/to/cellranger_count/dir/96hpi_infected_A/outs/filtered_feature_bc_matrix/")
ctrlA.data <- Read10X(data.dir = "/path/to/cellranger_count/dir/control_A/outs/filtered_feature_bc_matrix/")
sarsB.data <- Read10X(data.dir = "/path/to/cellranger_count/dir/96hpi_infected_B/outs/filtered_feature_bc_matrix/")
ctrlB.data <- Read10X(data.dir = "/path/to/cellranger_count/dir/control_B/outs/filtered_feature_bc_matrix/")
# Load the dataset for each condition (cellranger aggr group)
ctrl.data <- Read10X(data.dir = "/path/to/cellranger_aggr/dir/control/outs/count/raw_feature_bc_matrix/")
sars.data <- Read10X(data.dir = "/path/to/cellranger_aggr/dir/96hpi_infected/outs/count/raw_feature_bc_matrix/")

# Load the dataset for each condition (from cellranger count output)
setwd("/path/to/cellranger_count/dir/")

list1 <- c(
  "./sars1/",
  "./sars2/",
  "./sars3/",
  "./sars4/",
  "./sars5/"
)

list2 <- c(
  "./ctrl1/",
  "./ctrl2/",
  "./ctrl3/",
  "./ctrl4/",
  "./ctrl5/"
)


ctrl.data <- Read10X(data.dir = list2)
sars.data <- Read10X(data.dir = list1)


# Load a single aggr dataset
sars.data <- Read10X(data.dir = "/path/to/cellranger_aggr/dir/all/")


######################################     Regular pipe    ###########################################

# Initialize the Seurat object with the raw (non-normalized data).
# Set up Control:
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "Sars-cov_ctrl", min.cells = 3, min.features = 200)
ctrl$sars <- "CTRL"
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-") # Add mitochondrial percentage to the object metadata
#filter cells that have unique feature counts over 2,500 or less than 200 && cells that have >5% mitochondrial counts
ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
ctrl <- NormalizeData(ctrl, normalization.method = "LogNormalize", scale.factor = 10000) # Normalizing data -> Normalized values are stored in pbmc[["RNA"]]@data
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000) # Identification of highly variable features

# Set up Infected:
sars <- CreateSeuratObject(counts = sars.data, project = "Sars-cov_infec", min.cells = 3, min.features = 200)
sars$sars <- "SARS"
sars[["percent.mt"]] <- PercentageFeatureSet(sars, pattern = "^MT-") # Add mitochondrial percentage to the object metadata
#filter cells that have unique feature counts over 2,500 or less than 200 && cells that have >5% mitochondrial counts
sars <- subset(sars, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sars <- NormalizeData(sars, normalization.method = "LogNormalize", scale.factor = 10000) # Normalizing data -> Normalized values are stored in pbmc[["RNA"]]@data
sars <- FindVariableFeatures(sars, selection.method = "vst", nfeatures = 2000) # Identification of highly variable features


# Perform integration 
sars.anchors <- FindIntegrationAnchors(object.list = list(ctrl, sars), dims = 1:20)
sars.combined <- IntegrateData(anchorset = sars.anchors, dims = 1:20)

# Perform an integrated analysis
DefaultAssay(sars.combined) <- "integrated"

# Scaling the data -> shifts mean to 0 and variance to 1
# Gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(sars.combined)
sars.combined <- ScaleData(sars.combined, features = all.genes, verbose = FALSE)

# Perform linear dimensional reduction from the section right after ScTransform --------->>>>>>>>>




# More tips:
# Visualize QC metrics as a violin plot
VlnPlot(sars, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(sars, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sars, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sars), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(sars)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

# Examine and visualize PCA results a few different ways
print(sars[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sars, dims = 1:2, reduction = "pca")
DimPlot(sars, reduction = "pca")
DimHeatmap(sars, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sars, dims = 1:15, cells = 500, balanced = TRUE)




#################################       ScTransform Normalization method        #####################################

# This procedure omits the need for heuristic steps including pseudocount addition or log-transformation 
# and improves common downstream analytical tasks such as variable gene selection, dimensional reduction, and differential expression.
# Regress out heterogeneity associated with mitochondrial contamination


#### Introduction:

# Initialize the Seurat object with the raw (non-normalized data)
sars <- CreateSeuratObject(counts = sars.data)

# store mitochondrial percentage in object meta data
sars <- PercentageFeatureSet(sars, pattern = "^MT-", col.name = "percent.mt")

#filter cells that have unique feature counts over 2,500 or less than 200 && cells that have >5% mitochondrial counts
sars <- subset(sars, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# run sctransform
# this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
sars <- SCTransform(sars, vars.to.regress = "percent.mt", verbose = FALSE)
# or with a Gamma-Poisson Generalized Linear Model
sars <- SCTransform(sars, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# Perform dimensionality reduction by PCA and UMAP embedding
# These are now standard steps in the Seurat workflow for visualization and clustering
sars <- RunPCA(sars, verbose = FALSE)
sars <- RunUMAP(sars, dims = 1:30, verbose = FALSE)
sars <- FindNeighbors(sars, dims = 1:30, verbose = FALSE)
sars <- FindClusters(sars, verbose = FALSE)
DimPlot(sars, label = TRUE) + NoLegend()



##### Real pipeline:

# You can run the entire workflow at once after loading the data

### Pipeline for separate assays that will follow to integration: 
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "Sars-cov_ctrl", min.cells = 3, min.features = 200) %>% 
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt") #%>% 
  # RunPCA() %>% FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30) %>% FindClusters()
ctrl$sars <- "CTRL"

sars <- CreateSeuratObject(counts = sars.data, project = "Sars-cov_infec", min.cells = 3, min.features = 200) %>% 
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt") #%>% 
# RunPCA() %>% FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30) %>% FindClusters()
sars$sars <- "SARS"

# Perform integration
sars.list <- list(ctrl, sars)
features <- SelectIntegrationFeatures(object.list = sars.list, nfeatures = 3000)
sars.list <- PrepSCTIntegration(object.list = sars.list, anchor.features = features)
sars.anchors <- FindIntegrationAnchors(object.list = sars.list, normalization.method = "SCT", anchor.features = features)
sars.combined <- IntegrateData(anchorset = sars.anchors, normalization.method = "SCT")
# Perform an integrated analysis
DefaultAssay(sars.combined) <- "integrated"


### Pipeline from a single aggr dataset (integration performed through cellranger):
sars.combined <- CreateSeuratObject(counts = sars.data, min.cells = 3, min.features = 200) %>% 
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt") #%>% 
# RunPCA() %>% FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30) %>% FindClusters()

sars.combined$patients <- sub(".*-(.*)", "\\1", sars.combined@assays[["RNA"]]@data@Dimnames[[2]])
sars.combined$patients <- as.numeric(sars.combined$patients)
no_patho <- "No Pathology"
early_patho <- "Early Pathology"
late_patho <- "Late Pathology"
patho <- "Pathology"
no_patho.list <- scan("list1.txt")
early_patho.list <- scan("list3.txt")
late_patho.list <- scan("list2.txt")
patho.list <- scan("list4.txt")
sars.combined$patho <- recode(sars.combined$patients, "late_patho.list=late_patho; early_patho.list=early_patho; else=no_patho")
sars.combined$patho <- recode(sars.combined$patients, "patho.list=patho; else=no_patho")



# Perform linear dimensional reduction from next section ------------->>>>>>>>>>>




#########################   Continue both regular and ScTransform pipelines from here    ##################################


# Perform linear dimensional reduction
sars.combined <- RunPCA(sars.combined, npcs = 50, verbose = FALSE) # Perform linear dimensional reduction


# Optional:
# Determine the ‘dimensionality’ of the dataset -> how many components should we choose to include? 10? 20? 100?
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
sars.combined <- JackStraw(sars.combined, num.replicate = 100) #Cannot be performed with SCtransform normalized data
sars.combined <- ScoreJackStraw(sars.combined)#, dims = 1:30)
JackStrawPlot(sars.combined)#, dims = 1:15)
ElbowPlot(sars.combined, ndims = 20, reduction = "pca") #instead of jackstraw (for SCtransform normalized data)
DimHeatmap(sars.combined, dims = 1:20, cells = 500, balanced = TRUE)

# t-SNE or UMAP and Clustering
sars.combined <- RunTSNE(sars.combined, reduction = "pca", verbose = FALSE)
sars.combined <- RunUMAP(sars.combined, reduction = "pca", dims = 1:20, verbose = FALSE)
sars.combined <- FindNeighbors(sars.combined, reduction = "pca", dims = 1:20, verbose = FALSE)
sars.combined <- FindClusters(sars.combined, resolution = 0.5)

# save the Seurat object
saveRDS(sars.combined, file = "./sars.combined.rds")
# load back Seurat object
sars.combined <- readRDS("sars.combined.rds")

# Finding doublets
sweep.res <- paramSweep_v3(sars.combined, PCs = 1:20, sct = TRUE) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats) 
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
# define the expected number of doublet cells
nExp <- round(ncol(sars.combined) * 0.075)  # expect 7.5% doublets
sars.combined <- doubletFinder_v3(sars.combined, pN = 0.15, pK = 0.05, nExp = nExp, PCs = 1:20, sct = TRUE)

cowplot::plot_grid(ncol = 2, DimPlot(sars.combined, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(sars.combined, group.by = "DF.classifications_0.15_0.05_2077") + NoAxes()) ### the DF.classification always presents a unique name
sars.comb.no.doublet <- subset(sars.combined, DF.classifications_0.15_0.05_2077 == "Singlet")


# save the Seurat object
saveRDS(sars.comb.no.doublet, file = "./sars.comb.5000feats.nodoubs.rds")
# load back Seurat object
sars.combined <- readRDS("sars.comb.5000feats.nodoubs.rds")

# Cluster Visualization
p1 <- DimPlot(sars.combined, reduction = "umap", group.by = "patho") ## umap could be changed to tsne
p2 <- DimPlot(sars.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
# visualize the two conditions side-by-side
DimPlot(sars.combined, reduction = "umap", label = TRUE, split.by = "sars") #, raster=FALSE) #, split.by = "patho")



# scREAD Nervous System cells markers:
Astrocytes <- c("GFAP", "EAAT1", "AQP4", "LCN2", "GJA1", "SLC1A2", "FGFR3", "NKAIN4")
Endothelial_cells <- c("FLT1", "CLDN5", "VTN", "ITM2A", "VWF", "FAM167B", "BMX", "CLEC1B")
Excitatory_neurons <- c("SLC17A6",  "SLC17A7",  "NRGN", "CAMK2A", "SATB2", "COL5A1", "SDK2", "NEFM")
Inhibitory_neurons <- c("SLC32A1",  "GAD1", "GAD2", "TAC1", "PENK", "SST",  "NPY",  "MYBPC1", "PVALB", "GABBR2")
Microglia <- c("IBA-1", "P2RY12", "CSF1R",  "CD74", "C3", "CST3", "HEXB", "C1QA", "CX3CR1", "AIF-1")
Oligodendrocytes <- c("OLIG2",  "MBP",  "MOBP", "PLP1", "MOG",  "CLDN11", "MYRF", "GALC", "ERMN", "MAG")
Oligodendrocyte_precursor_cells <- c("VCAN", "CSPG4", "PDGFRA", "SOX10", "NEU4", "PCDG15", "GPR37L1", "C1QL1", "CDO1", "EPN2")
Pericytes <- c("AMBP",  "HIGD1B", "COX4I2", "AOC3", "PDE5A",  "PTH1R",  "P2RY14", "ABCC9", "KCNJ8", "CD248")
Macrophages <- c("CD14", "CD16", "CD64", "CD68", "CD71", "CCR5")


# Visualization of markers within the clusters. Obs: modify the targets as you wish
VlnPlot(sars.combined, features = c("CD14", "CD16"), pt.size = 0.01, ncol = 1, split.by = "sars")
# Visualize canonical marker genes on the sctransform embedding. Obs: modify the targets as you wish
FeaturePlot(sars.combined, features = "CD14", pt.size = 0.1, ncol = 2, reduction = "umap", split.by = "sars")


# find markers for every cluster compared to all remaining cells, report only positive ones
sars.markers <- FindAllMarkers(sars.combined, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25) %>% group_by(cluster)
# Follow All_markers.csv file to SCSA.py cellular identification
write.csv(as.data.frame(sars.markers), file="All_markers.csv")


# Rename clusters according to SCSA (most reliable cell typing method)
# The following line is just a example
sars.combined <- RenameIdents(sars.combined, `0` = "Oli_1", `1` = "Oli_2", `2` = "Ex_1", 
                              `3` = "Ast_1", `4` = "Opc", `5` = "Ast_2", `6` = "Oli_3", 
                              `7` = "Oli_4", `8` = "In_1", `9` = "Mic", 
                              `10` = "Ex_2", `11` = "In_2", `12` = "Ex_3",
                              `13` = "Ex_4", `14` = "In_3", `15` = "Oli_5",
                              `16` = "In_4", `17` = "Ex_5", `18` = "Ex_6", 
                              `19` = "Ex_7", `20` = "Peri", `21` = "Endo",
                              `22` = "Ex_8")
sars.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(sars.combined))
Idents(sars.combined) <- "celltypes"


DotPlot(sars.combined, features = c("CD14", "CD16"), cols = c("blue", "red"), col.max = 20, dot.scale = 15, cluster.idents = F,
split.by = "sars") + RotatedAxis()

DoHeatmap(sars.combined, features = GO_0008306, cells = 1:30000, size = 3)


# Differential expression analysis of two different conditions: Astrocytes infected and non-infected previously identified
sars.combined$celltype.sars <- paste(Idents(sars.combined), sars.combined$sars, sep = "_")
sars.combined$celltype <- Idents(sars.combined)
Idents(sars.combined) <- "seurat_clusters"
Idents(sars.combined) <- "celltype.sars"
sars.response.MAST2 <- FindMarkers(sars.combined, ident.1 = "Ast_SARS", 
                             ident.2 = "Ast_CTRL", 
                             #reduction = "umap",
                             test.use = "MAST", 
                             min.pct = 0.05, 
                             logfc.threshold = 0,
                             #slot = "counts",
                             verbose = FALSE
                             )
write.csv(as.data.frame(sars.response.MAST2), file="Markers_ctrl_x_sars_ast2.csv")

conv_markers <- read.delim("convergence_list.txt", header = F)
tconv_markers <- t(conv_markers)

top20 <- astro.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(astro.cells, features = top20$gene, cells = 1:4000, size = 2.5) + NoLegend()
