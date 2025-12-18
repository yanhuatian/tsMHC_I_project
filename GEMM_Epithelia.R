## set working directory 
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
####------------------------------------------------------------------
obj <- readRDS("/Users/ytian2/Desktop/Celltype.rds")
####------------------------------------------------------------------
Epi <- subset(obj,subset = CellType == c("Epithelial/Malignant cells"))
DimPlot(Epi, reduction = "umap",label =T,group.by = "CellType") 
{
  aaa<- NormalizeData(Epi) 
  aaa<- FindVariableFeatures(object = aaa)
  aaa <- ScaleData(aaa)
  aaa <- RunPCA(object = aaa, verbose = FALSE,features = VariableFeatures(aaa))
  aaa<- FindNeighbors(aaa, dims=1:50)
  Epi <- FindClusters(object = aaa,resolution = 0.5)
} # Norm
Epi<- RunUMAP(Epi,reduction = "pca",seed.use = 800,dims=1:50,min.dist=0.5,repulsion.strength = 2,n.neighbors=50)
DimPlot(Epi, reduction = "umap",label =T)
DimPlot(Epi, reduction = "umap",label =T,group.by = "Sample")
DimPlot(Epi, reduction = "umap",label =T,group.by = "Sample_set")
DimPlot(Epi, reduction = "umap",label =T,group.by = "mutation_info")

saveRDS(Epi,file = "/Users/ytian2/Desktop/Epithelia.rds")
####-----------------------------------------------------------------
Epi <- readRDS("/Users/ytian2/Desktop/Epithelia.rds")
##
VlnPlot(Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = F)
##
No <- as_tibble(Epi@meta.data)%>%group_by(Epi$seurat_clusters,Sample)%>%count
write_csv(No, "No.csv")
##
DEGs <- FindAllMarkers(Epi,only.pos = T, min.pct = 0.1)
write_csv(DEGs,"Epi.csv")
##
# cluster 0 (n=6892) # Maligntnt cells (AT2-like)
# cluster 1 (n=6264) # Low quality 
# cluster 2 (n=6108) # Low quality 
# cluster 3 (n=6264) # Malignant cells (AT2/AT1 mixed)
# cluster 4 (n=5126) # Malignant cells (EMT)
# cluster 5 (n=6264) # Malignant cells  (similar with cluster 4)
# cluster 6 (n=6264) # AT1 cells (Low quality)
# cluster 7 (n=6264) # Low quality (similar with cluster 1)
# cluster 8 (n=1229) # Ciliated cells
# cluster 9 (n=1112) # Immune cell
# cluster 10 (n=1000)# Low quality (Removed) and expressed Myeloid and Neutrophil genes
# cluster 11 (n=839) # Endothelia  (Removed)
# cluster 12 (n=636) # Clara/Goblet
# cluster 13 (n=377) # Normal AT2 (not malignant cells)
# cluster 14 (n=347) # Basal cells
# cluster 15 (n=286) # T cells (Removed)
####-----------Epi_new --------------------------------------------------------------
Epi_new <- subset(Epi, subset = seurat_clusters %in% c("0","3","4","5","6","8","12","13","14")) # remove cluster 1,2,7,9,10,11,15
DimPlot(Epi_new, reduction = "umap",label =T)
DimPlot(Epi_new, reduction = "umap",label =F,group.by = "Sample")
{
  aaa<- NormalizeData(Epi_new) 
  aaa<- FindVariableFeatures(object = aaa)
  aaa <- ScaleData(aaa)
  aaa <- RunPCA(object = aaa, verbose = FALSE,features = VariableFeatures(aaa))
  aaa<- FindNeighbors(aaa, dims=1:50)
  Epi_new <- FindClusters(object = aaa,resolution = 0.5)
} # Norm
Epi_new<- RunUMAP(Epi_new,reduction = "pca",seed.use = 800,dims=1:50,min.dist=0.5,repulsion.strength = 2,n.neighbors=50)

DimPlot(Epi_new, reduction = "umap",label =T)
DimPlot(Epi_new, reduction = "umap",label =F,group.by = "Sample")
DimPlot(Epi_new, reduction = "umap",label =F,group.by = "mutation_info")
##
VlnPlot(Epi_new, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = F)
##
DEGs <- FindAllMarkers(Epi_new,only.pos = T, min.pct = 0.1)
write_csv(DEGs,"Epi_new.csv")
##
Markers <- c("Hopx","Ager","Emp2",     # AT1
             "Sftpb","Lyz2","Lamp3",   # AT2
             "Scgb3a2","Scgb1a1","Hp", # Club cells
             "Muc5ac","Agr2","Muc5b",  # Goblet cells
             "Foxj1","Tppp3","Ccdc153",# Ciliated cells
             "Krt5","Aqp3","Krt14")    # Basal
DotPlot(Epi_new, features = Markers)
##
{
  metaData <- Epi_new@meta.data
  metaData$seurat_clusters <- as.character(metaData$seurat_clusters)
  metaData$CellType_2[metaData$seurat_clusters %in% c(0,1,3,4,6,9,10)] <- "Malignant cells"
  metaData$CellType_2[metaData$seurat_clusters %in% c(2)] <- "AT1 cells"
  metaData$CellType_2[metaData$seurat_clusters %in% c(5)]<- "Cilated cells"
  metaData$CellType_2[metaData$seurat_clusters %in% c(12)]  <- "Club cells"
  metaData$CellType_2[metaData$seurat_clusters %in% c(7)]  <- "Goblet cells"
  metaData$CellType_2[metaData$seurat_clusters %in% c(8)]   <- "AT2 cells"
  metaData$CellType_2[metaData$seurat_clusters %in% c(11,13)] <- "Basal cells"
  Epi_new@meta.data <- metaData
}
DimPlot(Epi_new, reduction = "umap",label =T)                             + NoLegend()
DimPlot(Epi_new, reduction = "umap",label =F, group.by = "Sample")        + NoLegend()
DimPlot(Epi_new, reduction = "umap",label =T, group.by = "CellType_2")    + NoLegend()
DimPlot(Epi_new, reduction = "umap",label =T, group.by = "mutation_info") + NoLegend()
##
saveRDS(Epi_new,file = "/Users/ytian2/Desktop/GEMM_Epi_CellType_2.rds")
####--------------------------------------------------------------------------------------
Epi_new <- readRDS("/Users/ytian2/Desktop/GEMM_Epi_new.rds")
Epi_new <- readRDS("/Users/ytian2/Desktop/GEMM_Epi_CellType_2.rds")
DimPlot(Epi_new, reduction = "umap",label =T)                             + NoLegend()
DimPlot(Epi_new, reduction = "umap",label =F, group.by = "Sample")        + NoLegend()
DimPlot(Epi_new, reduction = "umap",label =T, group.by = "CellType_2")    + NoLegend()
DimPlot(Epi_new, reduction = "umap",label =T, group.by = "mutation_info") + NoLegend()
##
Normal_Epi <- subset(Epi_new, subset = CellType_2 %in% c("AT1 cells","Cilated cells","Club cells","Goblet cells","AT2 cells","Basal cells"))
DimPlot(Normal_Epi, reduction = "umap",label =F)                            #+NoLegend()
DimPlot(Normal_Epi, reduction = "umap",label =F,group.by = "Sample")        #+NoLegend()
DimPlot(Normal_Epi, reduction = "umap",label =F,group.by = "CellType_2")    #+NoLegend()
DimPlot(Normal_Epi, reduction = "umap",label =F,group.by = "mutation_info") #+NoLegend()
##
Idents(Epi_new) <- Epi_new$CellType_2
DEGs <- FindAllMarkers(Epi_new,only.pos = T, min.pct = 0.1)
write_csv(DEGs,"DEGs_Epi_new.csv")
##
Markers <- c("Hopx","Ager","Emp2",            # AT1
             "Sftpb","Lyz2","Lamp3",          # AT2
             "Scgb3a2","Scgb1a1","Hp",        # Club cells
             "Muc5ac","Agr2","Muc5b",         # Goblet cells
             "Foxj1","Tppp3","Ccdc153",       # Ciliated cells
             "Krt5","Aqp3","Krt14")           # Basal
iOrd <- c("Basal cells","Ciliated cells","Club cells","Goblet cells","AT2 cells","AT1 cells")
Idents(Normal_Epi) <- factor(Normal_Epi$CellType_2,levels=iOrd)
##
DotPlot(Normal_Epi,features = Markers,) + 
  scale_color_gradientn(colors=rev(viridis(10))) + theme_bw()+
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_discrete(breaks=Markers,labels=Markers)
##
groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}
mat <- GetAssayData(Normal_Epi,slot="data")[Markers,]
mat <- groupMeans(mat,paste0(Normal_Epi$CellType_2),sparse=T)
bk <- c(seq(-2,0,by=0.02),seq(0.02,2,by=0.02))
colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(86),
               colorRampPalette(c("#d1e5f0","#f7f7f7"))(16),
               colorRampPalette(c("#f7f7f7","#fddbc7"))(16),
               colorRampPalette(c("#fddbc7","#b2182b"))(86))
pheatmap(mat, cluster_rows=F, cluster_cols=F,show_colnames=T,show_rownames=T, 
         scale="row",
         breaks = bk,color = colour_bk,
         clustering_method = "ward.D",
         display_numbers = F,
         angle_col=0,
         cutree_cols = 2)
####----------------------------------------------------------------------------
Epi_new <- readRDS("/Users/ytian2/Desktop/GEMM_Epi_CellType_2.rds")
##
p1 <- DimPlot(Epi_new, reduction = "umap",label =T, group.by = "CellType_2")    + NoLegend()
p2 <- DimPlot(Epi_new, reduction = "umap",label =T, group.by = "mutation_info") + NoLegend()
p1+p2        # Fig. S1e
##
No <- as_tibble(Epi_new@meta.data)%>%group_by(Epi_new$CellType_2,mutation_info)%>%count
write_csv(No, "No.csv")
####-------------------------inferCNV-------------------------------------------
{
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(infercnv)
  library(readr)
}
##------------Creat geneAnno for mouse------------
{
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  mm10.genes <- genes(txdb)
  geneRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
  mm10.genes$symbol <- mapIds(org.Mm.eg.db,
                              keys = mm10.genes$gene_id,
                              column = "SYMBOL",
                              keytype = "ENTREZID",
                              multiVals = "first")
  df <- data.frame(mm10.genes)
  df <- df[!is.na(df$symbol),]
  df <- data.frame(V1=df$symbol,V2=df$seqnames,V3=df$start,V4=df$end,
                   row.names = df$symbol,stringsAsFactors = F)
  df <- df[order(df$V2),]
  #df$V2 <- paste0(df$V2)
  df <- df[df$V2 %in% levels(df$V2)[1:20],]
  df <- as_tibble(df) %>%
    arrange(V2, V3)
  df$V2 <- stringr::str_extract(string = df$V2, pattern = "(?<=chr).*")
  #df$V2 <- paste0(df$V2)
  write.table(df,
              file="~/Desktop/mm10_gene_pos.txt",
              sep="\t",quote = F,row.names = F, col.names = F)
}
####---------InferCNV for GEMM (Fibroblast as reference)------------------------
obj <- readRDS("/Users/ytian2/Desktop/Celltype.rds")
Idents(obj) = obj$CellType
data_infercnv <- subset(obj, idents = c("Epithelial/Malignant cells", "Fibroblast"))

obj_epi <- readRDS("/Users/ytian2/Desktop/GEMM_Epi_CellType_2.rds")

data_infercnv$CellType_2 <- "Fibroblast"
data_infercnv$CellType_2[data_infercnv$CellType!="Fibroblast"] <- "removed"
data_infercnv$CellType_2[data_infercnv$CellType!="Fibroblast" & Cells(data_infercnv) %in% Cells(obj_epi)] <- 
  obj_epi$CellType_2[match(Cells(data_infercnv)[data_infercnv$CellType!="Fibroblast" & Cells(data_infercnv) %in% Cells(obj_epi)],
                           Cells(obj_epi))]

Idents(data_infercnv) <- data_infercnv$CellType_2
data_infercnv <- subset(data_infercnv, idents = c("removed"), invert = T)
data_infercnv <- subset(data_infercnv, cells = Cells(data_infercnv)[! is.na(data_infercnv$CellType_2)])
##
md <- as_tibble(data_infercnv@meta.data) %>%
  mutate(Barcode = Cells(data_infercnv)) %>%
  dplyr::select(Barcode, CellType_2)
write_tsv(md, "/Users/ytian2/Desktop/annotation.txt", col_names = F)
##
infercnvobj = infercnv::CreateInfercnvObject(
  raw_counts_matrix=GetAssayData(data_infercnv, assay = "RNA", slot = "counts"),
  annotations_file="/Users/ytian2/Desktop/annotation.txt",
  delim="\t",
  gene_order_file="/Users/ytian2/Desktop/mm10_gene_pos.txt",
  ref_group_names="Fibroblast")

print('...run infercnv...')
infercnvobj = infercnv::run(infercnvobj,
                            cutoff = 0.1,
                            out_dir = "/Users/ytian2/Desktop/inferCNV",
                            num_threads = 4,
                            denoise = T,
                            HMM = T,
                            tumor_subcluster_partition_method = "qnorm",
                            cluster_by_groups = T)
####-------------------------Monocle2--------------------------------------------
obj_Monocle2 <- subset(Epi_new,subset = seurat_clusters %in% c(0,1,2,3,6,8,9,11))
DimPlot(obj_Monocle2, reduction = "umap",label =T)
{
  metaData <- obj_Monocle2@meta.data
  metaData$seurat_clusters <- as.character(metaData$seurat_clusters)
  metaData$CellType_2[metaData$seurat_clusters %in% c(0)] <- "Malignant_C0"
  metaData$CellType_2[metaData$seurat_clusters %in% c(1)] <- "Malignant_C1"
  metaData$CellType_2[metaData$seurat_clusters %in% c(2)] <- "Malignant_C2"
  metaData$CellType_2[metaData$seurat_clusters %in% c(3)] <- "Malignant_C3"
  metaData$CellType_2[metaData$seurat_clusters %in% c(6)] <- "Malignant_C6"
  metaData$CellType_2[metaData$seurat_clusters %in% c(8)] <- "Malignant_C8"
  metaData$CellType_2[metaData$seurat_clusters %in% c(11)] <- "Malignant_C11"
  metaData$CellType_2[metaData$seurat_clusters %in% c(9)]  <- "AT2 cells"
  obj_Monocle2@meta.data <- metaData
}
DimPlot(obj_Monocle2, reduction = "umap",label =F,group.by = "CellType_2") #+NoLegend()
DimPlot(obj_Monocle2, reduction = "umap",label =F) # +NoLegend ()
saveRDS(obj_Monocle2,file = "/Users/ytian2/Desktop/Monocle2/obj_Monocle2.rds")
####-------------------------------------------------------------------------------------
obj_Monocle2 <- readRDS("/Users/ytian2/Desktop/Monocle2/obj_Monocle2.rds")
obj_Monocle2 <- subset(obj_Monocle2,subset = seurat_clusters %in% c(0,1,2,3,6,9,11))
Markers <- c("Sftpc","Lyz2","Scd1","Sftpa1","Sftpd","Sftpb","Hp","Lamp3",  # AT2 marker
             "Napsa","Nkx2-1", # LUAD Marker
             "Lcn2","Cxcl15","Cd74","H2-Aa","C3","C6","Ly6c1") # immune molecules
groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}
mat <- GetAssayData(obj_Monocle2,slot="data")[Markers,]
mat <- groupMeans(mat,paste0(obj_Monocle2$CellType_2),sparse=T)
bk <- c(seq(-2,0,by=0.02),seq(0.02,2,by=0.02))
colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(86),
               colorRampPalette(c("#d1e5f0","#f7f7f7"))(16),
               colorRampPalette(c("#f7f7f7","#fddbc7"))(16),
               colorRampPalette(c("#fddbc7","#b2182b"))(86))
pheatmap(mat, cluster_rows=F, cluster_cols=T,show_colnames=T,show_rownames=T, 
         scale="row",
         breaks = bk,color = colour_bk,
         clustering_method = "ward.D",
         display_numbers = F,
         angle_col=0,
         cutree_cols = 2)
write.table(mat, file = "~/Desktop/mat.txt", sep =  "\t", row.names = T, col.names = T)
####--------------------END----------------------------------------------------------






