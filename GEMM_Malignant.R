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
####--------------Malignant-----------------------------------------------------
Malignant <- readRDS("/Users/ytian2/Desktop/GEMM_Malignant.rds")
DimPlot(Malignant, reduction = "umap",label =T, group.by = "CellType_3")    + NoLegend()

##
Markers <- c("H2-D1","H2-K1",
             "H2-Q7","H2-Q4","H2-T23","H2-T22","H2-T24",  
             "Mr1","Hfe","Cd1d1",
             "B2m","Tap1","Tap2","Tapbp","Canx","Calr","Pdia3","Pdia6","Hspa5",
             "Hsp90b1","Selenos","Ddost","Erp29","Crtap",
             "Psmb4","Psmb8","Psmb9","Psmb10","Psma1","Psma3","Psmb5",
             "Psme1","Psme2","Psme3",
             "Erap1","Lnpep","Dpp9","Dpp8","Ofd1",
             "Nlrc5","Irf1","Irf2","Irf8","Ciita","Rfx5","Nfkb1","Rela","Stat2",
             "Tnf","Ifng","Stat1","Stat3","Stat5a","Stat5b","Jak1","Jak2",
             "Socs1","Socs3","Usp18","Il6","Il1b","Cxcl9","Cxcl10","Pml",
             "Ulbp1","Raet1d","Raet1e") # All

Markers <- c("H2-K1","H2-D1","H2-Q7",  
             "H2-T23","H2-T22","H2-T24","H2-Q4",   
             "Mr1","Hfe","Cd1d1") # MHC-I Molecules

Markers <- c("B2m", "Tap1", "Tap2", "Tapbp", "Canx", "Calr", "Pdia3", "Pdia6", "Hspa5",
             "Hsp90b1", "Selenos", "Ddost", "Erp29", "Crtap",
             "Psmb4", "Psmb8", "Psmb9", "Psmb10", "Psma1", "Psma3", "Psmb5",
             "Psme1", "Psme2", "Psme3",
             "Erap1", "Lnpep", "Dpp9", "Dpp8", "Ofd1") # Antigen Processing Machinery (APM)

Markers <- c("Nlrc5", "Irf1", "Irf2", "Irf8", "Ciita", "Rfx5", "Nfkb1", "Rela", "Stat2",
             "Tnf", "Ifng", "Stat1", "Stat3", "Stat5a", "Stat5b", "Jak1", "Jak2",
             "Socs1", "Socs3", "Usp18", "Il6", "Il1b", "Cxcl9", "Cxcl10", "Pml") # Transcriptional Regulators

Markers <- c("Ulbp1","Raet1d","Raet1e") # Immune Surveilance Ligands
##
Malignant <- AddModuleScore(object = Malignant, 
                            features = list(Markers),
                            name = 'Antigen_presentation', 
                            ctrl = 10)
##
levels <- c("C0_AT2-like","C1_Str","C2_AT2/AT1","C3_HPCS","C4_hE/M","C5_Pro","C6_Rib","C7_EMT")

Malignant$CellType_3 <- factor(Malignant$CellType_3,levels = levels)

Idents(Malignant) <- Malignant$CellType_3
##
VlnPlot(Malignant, features = "Antigen_presentation1", pt.size = 0) +  # pt.size = 0 to hide points
        stat_summary(fun = mean, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(-0.25, 0.5) + NoLegend() # ✅ Fix y-axis limits

# Create the violin plot and store it in an object
a <- VlnPlot(object = Malignant,
          features = "Antigen_presentation1",
          pt.size  = 0) + 
          stat_summary(fun = mean,aes(group = 1),geom = "line",color = "black",size = 1) +
          stat_summary(fun = mean,geom = "point",color = "black",size = 3) +
          ylim(-0.25, 0.5) +
          NoLegend()
##
a+b+c+d
####---------------Correlation of tsMHC-I score with 8 clusters-----------------
# 1) Fetch the module score and cluster identity for each cell
df <- FetchData(Malignant, vars = c("Antigen_presentation1", "CellType_3"))

# 2) Compute the mean module score per cluster
cluster_means <- df %>%
  group_by(CellType_3) %>%
  summarise(mean_AP = mean(Antigen_presentation1, na.rm = TRUE)) %>%
  ungroup()

# 3) Assign numeric indices for each cluster (C0=0, C1=1, ..., C7=7)
cluster_means$cluster_index <- 0:(nrow(cluster_means)-1)

# 4) Perform a correlation test (Pearson or Spearman)
cor_test <- cor.test(cluster_means$cluster_index, cluster_means$mean_AP, method = "pearson")
# If you prefer rank-based correlation:
# cor_test <- cor.test(cluster_means$cluster_index, cluster_means$mean_AP, method = "spearman")

# Print the results
print(cluster_means)
print(cor_test)
####----------------------------------------------------------------------------
# Extract Antigen_presentation1 scores
scores <- Malignant@meta.data$Antigen_presentation1

# Estimate density
density_data <- density(scores)

# Identify local minima (valleys) as potential cutoff points
minima_indices <- which(diff(sign(diff(density_data$y))) == 2) + 1
cutoff_values <- density_data$x[minima_indices]

# Select the most relevant cutoff (first minimum)
cutoff_value <- ifelse(length(cutoff_values) > 0, cutoff_values[1], NA)

# Print identified cutoff value
print(paste("Selected Cutoff Value:", round(cutoff_value, 3)))

# Create density histogram plot
df <- data.frame(Antigen_presentation1 = scores)

ggplot(df, aes(x = Antigen_presentation1)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "brown", color = "black", alpha = 0.3) +
  geom_density(color = "black", size = 1) +
  geom_vline(xintercept = cutoff_value, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = cutoff_value, y = max(density_data$y) * 0.9, 
           label = paste0("Cutoff: ", round(cutoff_value, 3)), 
           color = "red", fontface = "bold", size = 5) +
  labs(x = "Antigen_presentation1 Score", y = "Density (%)") +
  theme_minimal()
##
Malignant$Antigen_presentation_group <- ifelse(
  Malignant@meta.data$Antigen_presentation1 < cutoff_value, "Low", "High"
)

# Check cell distribution
DimPlot(Malignant, group.by = "Antigen_presentation_group", label = TRUE) +
        ggtitle("Antigen Presentation Groups") 

# Create a dataframe from Seurat metadata
df <- data.frame(
  tsMHC_Score = Malignant@meta.data$Antigen_presentation1,  # Replace with your actual score column
  Antigen_presentation_group = Malignant$Antigen_presentation_group  # "High" and "Low" groups
)

# Generate the violin plot
ggplot(df, aes(x = Antigen_presentation_group, y = tsMHC_Score, fill = Antigen_presentation_group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +  # Violin plot with width scaling
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +  # Boxplot inside violins
  scale_fill_manual(values = c("Low" = "lightblue", "High" = "darkred")) +  # Custom colors
  labs(x = "tsMHC", y = "tsMHC Score") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))

saveRDS(Malignant,file = "/Users/ytian2/Desktop/GEMM_Malignant_tsMHC.rds")
####----------------------------------------------------------------------------
Markers <- c("Mki67","Birc5","Top2a","Ccna2","Ccnb1","Tk1")

##
Malignant <- AddModuleScore(object = Malignant, 
                            features = list(Markers),
                            name = 'Proliferating', 
                            ctrl = 10)
##
levels <- c("C0_AT2-like","C1_Str","C2_AT2/AT1","C3_HPCS","C4_hE/M","C5_Pro","C6_Rib","C7_EMT")

Malignant$CellType_3 <- factor(Malignant$CellType_3,levels = levels)

Idents(Malignant) <- Malignant$CellType_3
##
VlnPlot(Malignant, features = "Proliferating1", pt.size = 0) +  # pt.size = 0 to hide points
  stat_summary(fun = mean, geom = 'point', size = 10, colour = "black", shape = 95) +
  ylim(-0.25, 0.5) + NoLegend() # ✅ Fix y-axis limits
##

# 1. Extract tsMHC and proliferation scores, plus the existing tsMHC group from your Seurat object.
df <- FetchData(Malignant, vars = c("Antigen_presentation1", "Proliferating1", "tsMHC_group"))

# 2. For cells in the tsMHC-Low group, calculate the median proliferating score and classify them.
df <- df %>%
  group_by(tsMHC_group) %>%
  mutate(final_group = ifelse(tsMHC_group == "tsMHC_High",
                              "tsMHC_High",
                              ifelse(Proliferating1 >= median(Proliferating1, na.rm = TRUE),
                                     "tsMHC_Low_Prolif_High",
                                     "tsMHC_Low_Prolif_Low"))) %>%
  ungroup()

# 3. Add the new grouping back to the Seurat object's metadata.
Malignant$new_group <- df$final_group

# (Optional) Set the cell identities to the new grouping for downstream analysis/visualization.
Idents(Malignant) <- "new_group"

# 4. To check the distribution of cells among groups:
table(Malignant$new_group)

##
Idents(Malignant) <- Malignant$new_group
DEGs <- FindAllMarkers(Malignant,only.pos = F, min.pct = 0.1)  
DEGs$Gene=rownames(DEGs)
write_csv(DEGs,"DEGs_123.csv")
####-------------------END------------------------------------------------------
Malignant <- readRDS("/Users/ytian2/Desktop/GEMM_Malignant.rds")

Malignant_new <- subset(Malignant, subset = Sample %in% c("GEMM_L4","GEMM_L5"))
##
Markers <- c("H2-K1","H2-D1","H2-Q7",  
             "H2-T23","H2-T22","H2-T24","H2-Q4",   
             "Mr1","Hfe","Cd1d1") # MHC-I Molecules
##
Malignant_new <- Seurat::AddModuleScore(
  object = Malignant_new,                               # The Seurat object
  features = list(Antigen_presentation = Markers),  # Markers used for the module score
  name = "Antigen_presentation")                    # Prefix for the new score column
##
Idents(Malignant_new) <- Malignant_new$Sample
VlnPlot(Malignant_new, features ="Antigen_presentation1",pt.size = F) + 
  stat_summary(fun.y = mean, geom='point', size = 10, colour = "black", shape =95) +NoLegend()
##
ggplot(Malignant_new@meta.data, aes(x = Antigen_presentation1, fill = Sample)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Sample, scales = "free_y") +
  theme_bw() + NoLegend()

ggplot(Malignant_new@meta.data, aes(x = Antigen_presentation1)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  theme_bw()




