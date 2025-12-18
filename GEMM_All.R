####----------------------------------------------------------------------------
setwd("/Users/ytian2/Desktop")
## library packages 
suppressMessages({ 
  library(Seurat) 
  library(dplyr)
  library(SeuratObject) 
  library(viridis)
  library(ggplot2)
  library(pheatmap)
  library(tidyverse)
  library(harmony)
})
####----------------------------------------------------------------------------
obj <- readRDS("/Users/ytian2/Desktop/cluster.rds")
DimPlot(obj, reduction = "umap",label =T,group.by = "seurat_clusters")
DimPlot(obj, reduction = "umap",label =T,group.by = "orig.ident")
####---------include KRAS G12D mutation information-----------------------------
obj <- readRDS("/Users/ytian2/Desktop/KP_GEMM_KrasG12D_added.rds")
DimPlot(obj, reduction = "umap",label =T,group.by = "mutation_info")
DimPlot(obj, reduction = "umap",label =T,group.by = "KrasG12D_mutation")
####----------QC----------------------------------------------------------------
obj@active.assay <- "RNA"
obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
Idents(obj)  <- obj$orig.ident
Idents(obj)  <- paste0(Idents(obj)) 
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = F)
####----------Assign sample name------------------------------------------------
{
  metaData <- obj@meta.data
  metaData$orig.ident <- as.character(metaData$orig.ident)
  metaData$Sample[metaData$orig.ident %in% c("JZH_56524")] <- "GEMM_B0"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56525")] <- "GEMM_L0"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56514")] <- "GEMM_B1"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56515")] <- "GEMM_L1"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56516")] <- "GEMM_B2"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56517")] <- "GEMM_L2"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56518")] <- "GEMM_B3"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56519")] <- "GEMM_L3"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56520")] <- "GEMM_B4"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56521")] <- "GEMM_L4"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56522")] <- "GEMM_B5"
  metaData$Sample[metaData$orig.ident %in% c("JZH_56523")] <- "GEMM_L5"
  obj@meta.data <- metaData
}
####----------Assign sample_set name--------------------------------------------
{
  metaData <- obj@meta.data
  metaData$Sample <- as.character(metaData$Sample)
  metaData$Sample_set[metaData$Sample %in% c("GEMM_B0","GEMM_B1","GEMM_B2",
                                             "GEMM_B3","GEMM_B4","GEMM_B5")] <- "GEMM_B"
  metaData$Sample_set[metaData$Sample %in% c("GEMM_L0","GEMM_L1","GEMM_L2",
                                             "GEMM_L3","GEMM_L4","GEMM_L5")] <- "GEMM_L"
  obj@meta.data <- metaData
}
saveRDS(obj,file = "/Users/ytian2/Desktop/cluster_Sample.rds")
####----------------------------------------------------------------------------
obj <- readRDS("/Users/ytian2/Desktop/cluster_Sample.rds")
DimPlot(obj, reduction = "umap",label =T,raster=FALSE,group.by = "Sample")
DimPlot(obj, reduction = "umap",label =T,raster=FALSE,group.by = "Sample_set")
# Malignant/Epithelia cells
FeaturePlot(obj,features = c("Krt7","Krt8","Krt18","Clu"))
# T cells
FeaturePlot(obj,features = c("Ptprc","Cd2","Cd3d","Cd3g"))
# NK cells
FeaturePlot(obj,features = c("Ncr1","Nkg7","Klrd1","Gzma"))
# B cells
FeaturePlot(obj,features = c("Cd19","Cd79a","Cd79b","Ms4a1"))
# Mono/Macro/DC cells
FeaturePlot(obj,features = c("Csf1r","Fcgr2b","Lyz2","Cd68"))
# pDC (cluster20)
FeaturePlot(obj,features = c("Siglech","Ly6c2","Mpeg1","Cd209d"))
# Neutrophils
FeaturePlot(obj,features = c("Csf3r","S100a9","S100a8","G0s2"))
# Mast cells
FeaturePlot(obj,features = c("Ms4a2","Cpa3","Fcer1a","Cd200r3")) # Other mast cell markers: Tpsb2  Mcpt4
# Megakaryocytes
FeaturePlot(obj,features = c("Ppbp","Itga2b","Tubb1","Pf4"))
# Erythroid cells
FeaturePlot(obj,features = c("Gypa","Hbb-bs","Alas2","Bpgm","Snca","Mkrn1","Fech")) 
# Fibroblast
FeaturePlot(obj,features = c("Dcn","Col1a1","Pdgfra","Col1a2"))
# Endothelia
FeaturePlot(obj,features = c("Eng","Egfl7","Cdh5","Cldn5"))
####----------------------------------------------------------------------------
DEGs <- FindAllMarkers(obj, only.pos = T, min.pct = 0.1)
write_csv(DEGs,"Malignant_DEGs.csv")
####---------Assign cell types--------------------------------------------------
{
  metaData <- obj@meta.data
  metaData$seurat_clusters <- as.character(metaData$seurat_clusters)
  metaData$CellType[metaData$seurat_clusters %in% c(9,12,20)] <- "Epithelial/Malignant cells"
  metaData$CellType[metaData$seurat_clusters %in% c(0,6,7)]   <- "T cells"
  metaData$CellType[metaData$seurat_clusters %in% c(15)]      <- "NK cells"
  metaData$CellType[metaData$seurat_clusters %in% c(1,24)]    <- "B/Plasma cells"
  metaData$CellType[metaData$seurat_clusters %in% c(2,4,8,14,19)] <- "Mono/Macro/DC"
  metaData$CellType[metaData$seurat_clusters %in% c(10)]      <- "Neutrophils"
  metaData$CellType[metaData$seurat_clusters %in% c(3,16,21,25)] <- "Fibroblast" 
  metaData$CellType[metaData$seurat_clusters %in% c(18)]      <- "Mast cells"
  metaData$CellType[metaData$seurat_clusters %in% c(5,11,13,22)] <- "Endothelia" 
  metaData$CellType[metaData$seurat_clusters %in% c(17)]      <- "Megakaryocytes" 
  metaData$CellType[metaData$seurat_clusters %in% c(23,26)]   <- "Erythroid cells"
  obj@meta.data <- metaData
}
saveRDS(obj,file = "/Users/ytian2/Desktop/Celltype.rds")
####----------------------------------------------------------------------------
obj <- readRDS("/Users/ytian2/Desktop/Celltype.rds")
## UMAP
library(ArchR)
cols <- ArchR::paletteDiscrete(values = unique(obj$CellType),set = "stallion2") # find more color schemes at: https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
DimPlot(obj, reduction = "umap",label =T,group.by = "seurat_clusters")
DimPlot(obj, reduction = "umap",label =T,group.by = "CellType") + NoLegend()
DimPlot(obj, reduction = "umap",label =T,group.by = "CellType", cols=cols) 
DimPlot(obj, reduction = "umap",label =F,group.by = "CellType", cols=cols, raster=FALSE) + NoLegend()
DimPlot(obj, reduction = "umap",label =T,group.by = "Sample")  # stallion2, which included 20 colors,dose not work here, because there are 21 samples
##
col <- c("T cells"            = "#D51F26",  # Red
         "Mast cells"         = "#272E6A",  # Dark Blue
         "NK cells"           = "#208A42",  # Green
         "Mono/Macro/DC"      = "#89288F",  # Purple
         "B/Plasma cells"     = "#F47D2B",  # Orange
         "Megakaryocytes"     = "#FEE500",  # Yellow
         "Erythroid cells"    = "#90D5E4",  # Light Cyan
         "Neutrophils"        = "#89C75F",  # Lime Green
         "Fibroblast"         = "#F37B7D",  # Salmon Pink
         "Endothelia"         = "#6E4B9E",  # Indigo
         "Epithelial/Malignant cells" = "#0C727C")  # Teal

DimPlot(obj, reduction = "umap", label = TRUE, group.by = "CellType") +
        scale_color_manual(values = col) # + NoLegend()        # Figure 1d (6*6 inches)          
####----------------------------------------------------------------------------
Markers <- c("Krt7","Krt8","Krt18","Krt19",         # Epithelial/Malignant cells
             "Cd2","Cd3d","Cd3e","Cd3g",            # T cells
             "Ncr1","Nkg7","Klrd1","Gzma",          # NK cells
             "Cd19","Cd79a","Ms4a1","Cd79b",        # B/Plasma cells
             "Csf1r","Mpeg1","Lyz2","Cd68",         # Mono/Macro/DC
             "Csf3r","S100a9","S100a8","G0s2",      # Neutrophils
             "Ms4a2","Cpa3","Fcer1a","Cd200r3",     # Mast cells
             "Alas2","Fech","Hbb-bs","Hba-a1",      # Erythroid cells
             "Ppbp","Itga2b","Tubb1","Pf4",         # Mega
             "Dcn","Col1a1","Pdgfra","Col1a2",      # Fibroblast
             "Eng","Egfl7","Cdh5","Cldn5")          # Endothelia
iOrd <- c("Epithelial/Malignant cells","T cells","NK cells","B/Plasma cells","Mono/Macro/DC",
          "Neutrophils","Mast cells","Erythroid cells","Megakaryocytes","Fibroblast","Endothelia")
Idents(obj) <- factor(obj$CellType,levels=iOrd)
DotPlot(obj,features = Markers,) + 
  scale_color_gradientn(colors=rev(viridis(10))) + theme_bw()+
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_discrete(breaks=Markers,labels=Markers)      # Fig.S1a (4*15 inches)
####----------------------------------------------------------------------------
No <- as_tibble(obj@meta.data)%>%group_by(obj$Sample,CellType)%>%count
write_csv(No, "No.csv")

No <- as_tibble(obj@meta.data)%>%group_by(obj$seurat_clusters,Sample_set)%>%count
write_csv(No, "No.csv")
####---------KRAS G12D Mutation information-------------------------------------
obj <- readRDS("/Users/ytian2/Desktop/KP_GEMM_KrasG12D_added.rds")
library(ArchR)
cols <- ArchR::paletteDiscrete(values = unique(obj$CellType), set = "stallion2") # find more color schemes at: https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
DimPlot(obj, reduction = "umap",label =T,group.by = "seurat_clusters")
DimPlot(obj, reduction = "umap",label =T,group.by = "CellType",cols=cols)  #+ NoLegend() 
DimPlot(obj, reduction = "umap",label =F,group.by = "CellType",cols=cols,raster=FALSE) + NoLegend()
DimPlot(obj, reduction = "umap",label =T,group.by = "KrasG12D_mutation")
DimPlot(obj, reduction = "umap",label =T,group.by = "mutation_info")
##
No <- as_tibble(obj@meta.data)%>%group_by(obj$Sample,mutation_info) %>% count
write_csv(No, "No.csv")
####-----------------------------END--------------------------------------------

