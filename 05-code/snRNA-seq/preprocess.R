library(Seurat)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(tidyverse)
library(DoubletFinder)

# Set R's global seed
set.seed(42)

# Load count matrix data from CellRanger ----
# Control
ctrl.11.mtx <- Read10X(data.dir = "02-processed/cellranger_Control_SP_11_GEX_FL-Z0041/filtered_feature_bc_matrix/")
ctrl.12.mtx <- Read10X(data.dir = "02-processed/cellranger_Control_SP_12_GEX_FL-Z0044/filtered_feature_bc_matrix/")
ctrl.14.mtx <- Read10X(data.dir = "02-processed/cellranger_Control_SP_14_GEX_FL-Z0047/filtered_feature_bc_matrix/")
ctrl.15.mtx <- Read10X(data.dir = "02-processed/cellranger_Control_SP_15_GEX_FL-Z0054/filtered_feature_bc_matrix/")

# MCT-Water
water.1.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Water_1_GEX_FL-Z0043/filtered_feature_bc_matrix/")
water.5.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Water_5_GEX_FL-Z0046/filtered_feature_bc_matrix/")
water.6.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Water_6_GEX_FL-Z0053/filtered_feature_bc_matrix/")
water.9.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Water_9_GEX_FL-Z0056/filtered_feature_bc_matrix/")

# MCT-Blumeria
blum.1.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Blumeria_1_GEX_FL-Z0042/filtered_feature_bc_matrix/")
blum.3.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Blumeria_3_GEX_FL-Z0045/filtered_feature_bc_matrix/")
blum.4.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Blumeria_4_GEX_FL-Z0048/filtered_feature_bc_matrix/")
blum.7.mtx <- Read10X(data.dir = "02-processed/cellranger_MCT_Blumeria_7_GEX_FL-Z0055/filtered_feature_bc_matrix/")

# Create Seurat Objects from Count Matrices ----
ctrl.11.seu <- CreateSeuratObject(counts = ctrl.11.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control11")
ctrl.12.seu <- CreateSeuratObject(counts = ctrl.12.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control12")
ctrl.14.seu <- CreateSeuratObject(counts = ctrl.14.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control14")
ctrl.15.seu <- CreateSeuratObject(counts = ctrl.15.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control15")
water.1.seu <- CreateSeuratObject(counts = water.1.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water1")
water.5.seu <- CreateSeuratObject(counts = water.5.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water5")
water.6.seu <- CreateSeuratObject(counts = water.6.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water6")
water.9.seu <- CreateSeuratObject(counts = water.9.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water9")
blum.1.seu <- CreateSeuratObject(counts = blum.1.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria1")
blum.3.seu <- CreateSeuratObject(counts = blum.3.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria3")
blum.4.seu <- CreateSeuratObject(counts = blum.4.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria4")
blum.7.seu <- CreateSeuratObject(counts = blum.7.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria7")

# Strict QC - Filtering ----
# Check mitochondrial DNA
ctrl.11.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.11.seu, pattern = "^Mt-")
ctrl.12.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.12.seu, pattern = "^Mt-")
ctrl.14.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.14.seu, pattern = "^Mt-")
ctrl.15.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.15.seu, pattern = "^Mt-")
water.1.seu[['percent.mt']] <- PercentageFeatureSet(water.1.seu, pattern = "^Mt-")
water.5.seu[['percent.mt']] <- PercentageFeatureSet(water.5.seu, pattern = "^Mt-")
water.6.seu[['percent.mt']] <- PercentageFeatureSet(water.6.seu, pattern = "^Mt-")
water.9.seu[['percent.mt']] <- PercentageFeatureSet(water.9.seu, pattern = "^Mt-")
blum.1.seu[['percent.mt']] <- PercentageFeatureSet(blum.1.seu, pattern = "^Mt-")
blum.3.seu[['percent.mt']] <- PercentageFeatureSet(blum.3.seu, pattern = "^Mt-")
blum.4.seu[['percent.mt']] <- PercentageFeatureSet(blum.4.seu, pattern = "^Mt-")
blum.7.seu[['percent.mt']] <- PercentageFeatureSet(blum.7.seu, pattern = "^Mt-")

# Visualize Mitochondrial Distribution per Sample
seu_list <- list(
  "Ctrl_11" = ctrl.11.seu,
  "Ctrl_12" = ctrl.12.seu,
  "Ctrl_14" = ctrl.14.seu,
  "Ctrl_15" = ctrl.15.seu,
  "Water_1" = water.1.seu,
  "Water_5" = water.5.seu,
  "Water_6" = water.6.seu,
  "Water_9" = water.9.seu,
  "Blum_1" = blum.1.seu,
  "Blum_3" = blum.3.seu,
  "Blum_4" = blum.4.seu,
  "Blum_7" = blum.7.seu
)

mt_qc_check <- map_df(seu_list, ~ .x@meta.data, .id = "Sample") %>%
  group_by(Sample) %>%
  summarize(
    Total_Nuclei = n(),
    Median_mt = round(median(percent.mt), 2),
    `Pass_<2%` = paste0(round(mean(percent.mt < 2) * 100, 1), "%"),
    `Pass_<3%` = paste0(round(mean(percent.mt < 3) * 100, 1), "%"),
    `Pass_<5%` = paste0(round(mean(percent.mt < 5) * 100, 1), "%")
  )

print(as.data.frame(mt_qc_check))


# 1. Extract just the metadata from all 12 objects into one data frame
meta_df <- map_df(seu_list, ~ .x@meta.data, .id = "List_Name")

# 2. Reshape the data so we can plot all 3 metrics side-by-side
meta_long <- meta_df %>%
  select(orig.ident, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(
    cols = c(nFeature_RNA, nCount_RNA, percent.mt),
    names_to = "Feature",
    values_to = "Value"
  )

# 3. Recreate the Seurat VlnPlot using standard ggplot2
ggplot(meta_long, aes(x = orig.ident, y = Value, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5) +
  facet_wrap(~Feature, scales = "free_y", ncol = 3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(x = "Sample", y = "Value", title = "Pre-QC Metrics per Sample")

# Function to filter seurat objects
filter_sample_qc <- function(seu_obj, mt_cutoff = 2) {
  # Calculate exact quantiles first
  f_min <- quantile(seu_obj$nFeature_RNA, 0.01)
  f_max <- quantile(seu_obj$nFeature_RNA, 0.99)
  c_min <- quantile(seu_obj$nCount_RNA, 0.01)
  c_max <- quantile(seu_obj$nCount_RNA, 0.99)
  
  # Apply filtering
  seu_filtered <- subset(
    seu_obj,
    subset = nFeature_RNA > f_min &
      nFeature_RNA < f_max &
      nCount_RNA > c_min &
      nCount_RNA < c_max &
      percent.mt < mt_cutoff
  )
  
  return(seu_filtered)
}

# Apply to one object:
ctrl.11.seu <- filter_sample_qc(ctrl.11.seu, mt_cutoff = 3)
ctrl.12.seu <- filter_sample_qc(ctrl.12.seu, mt_cutoff = 3)
ctrl.14.seu <- filter_sample_qc(ctrl.14.seu, mt_cutoff = 3)
ctrl.15.seu <- filter_sample_qc(ctrl.15.seu, mt_cutoff = 3)
water.1.seu <- filter_sample_qc(water.1.seu, mt_cutoff = 3)
water.5.seu <- filter_sample_qc(water.5.seu, mt_cutoff = 3)
water.6.seu <- filter_sample_qc(water.6.seu, mt_cutoff = 3)
water.9.seu <- filter_sample_qc(water.9.seu, mt_cutoff = 3)
blum.1.seu <- filter_sample_qc(blum.1.seu, mt_cutoff = 3)
blum.3.seu <- filter_sample_qc(blum.3.seu, mt_cutoff = 3)
blum.4.seu <- filter_sample_qc(blum.4.seu, mt_cutoff = 3)
blum.7.seu <- filter_sample_qc(blum.7.seu, mt_cutoff = 3)

# 1. Update the list with your newly filtered Seurat objects
seu_list_filtered <- list(
  "Ctrl_11" = ctrl.11.seu,
  "Ctrl_12" = ctrl.12.seu,
  "Ctrl_14" = ctrl.14.seu,
  "Ctrl_15" = ctrl.15.seu,
  "Water_1" = water.1.seu,
  "Water_5" = water.5.seu,
  "Water_6" = water.6.seu,
  "Water_9" = water.9.seu,
  "Blum_1"  = blum.1.seu,
  "Blum_3"  = blum.3.seu,
  "Blum_4"  = blum.4.seu,
  "Blum_7"  = blum.7.seu
)

# 2. Extract post-QC metadata
meta_df_post <- map_df(seu_list_filtered, ~ .x@meta.data, .id = "List_Name")

# 3. Reshape for plotting
meta_long_post <- meta_df_post %>%
  select(orig.ident, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(
    cols = c(nFeature_RNA, nCount_RNA, percent.mt),
    names_to = "Feature",
    values_to = "Value"
  )

# 4. Generate Post-QC Violin Plot
ggplot(meta_long_post, aes(x = orig.ident, y = Value, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5) +
  facet_wrap(~Feature, scales = "free_y", ncol = 3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(x = "Sample", y = "Value", title = "Post-QC Metrics per Sample")

qc_summary_post <- meta_df_post %>%
  group_by(List_Name) %>%
  summarize(
    Remaining_Nuclei = n(),
    Min_nFeature = min(nFeature_RNA),
    Max_nFeature = max(nFeature_RNA),
    Min_nCount   = min(nCount_RNA),
    Max_nCount   = max(nCount_RNA),
    Max_percent_mt = round(max(percent.mt), 2)
  )

print(as.data.frame(qc_summary_post))


#Violin Plots and subsetting
# VlnPlot(ctrl.11.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(ctrl.11.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(ctrl.11.seu$nCount_RNA, c(0.01, 0.99))
# quantile(ctrl.11.seu$percent.mt, 0.99)
# ctrl.11.seu <- subset(ctrl.11.seu,
#                       subset = nFeature_RNA > quantile(ctrl.11.seu$nFeature_RNA, 0.01) &
#                         nFeature_RNA < quantile(ctrl.11.seu$nFeature_RNA, 0.99) &
#                         nCount_RNA > quantile(ctrl.11.seu$nCount_RNA, 0.01) & 
#                         nCount_RNA < quantile(ctrl.11.seu$nCount_RNA, 0.99) &
#                         # percent.mt < quantile(ctrl.11.seu$percent.mt, 0.99))
#                         percent.mt < 2)

# VlnPlot(ctrl.12.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(ctrl.12.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(ctrl.12.seu$nCount_RNA, c(0.01, 0.99))
# quantile(ctrl.12.seu$percent.mt, 0.99)
# ctrl.12.seu <- subset(ctrl.12.seu,
                      # subset = nFeature_RNA > quantile(ctrl.12.seu$nFeature_RNA, 0.01) &
                        # nFeature_RNA < quantile(ctrl.12.seu$nFeature_RNA, 0.99) &
                        # nCount_RNA > quantile(ctrl.12.seu$nCount_RNA, 0.01) & 
                        # nCount_RNA < quantile(ctrl.12.seu$nCount_RNA, 0.99) &
                        # percent.mt < quantile(ctrl.12.seu$percent.mt, 0.99))
                        # percent.mt <2)

# VlnPlot(ctrl.14.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(ctrl.14.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(ctrl.14.seu$nCount_RNA, c(0.01, 0.99))
# quantile(ctrl.14.seu$percent.mt, 0.99)
# ctrl.14.seu <- subset(ctrl.14.seu,
                      # subset = nFeature_RNA > quantile(ctrl.14.seu$nFeature_RNA, 0.01) &
                        # nFeature_RNA < quantile(ctrl.14.seu$nFeature_RNA, 0.99) &
                        # nCount_RNA > quantile(ctrl.14.seu$nCount_RNA, 0.01) & 
                        # nCount_RNA < quantile(ctrl.14.seu$nCount_RNA, 0.99) &
                        # percent.mt < quantile(ctrl.14.seu$percent.mt, 0.99))
                        # percent.mt <2)

# VlnPlot(ctrl.15.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(ctrl.15.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(ctrl.15.seu$nCount_RNA, c(0.01, 0.99))
# quantile(ctrl.15.seu$percent.mt, 0.99)
# ctrl.15.seu <- subset(ctrl.15.seu,
                      # subset = nFeature_RNA > quantile(ctrl.15.seu$nFeature_RNA, 0.01) &
                        # nFeature_RNA < quantile(ctrl.15.seu$nFeature_RNA, 0.99) &
                        # nCount_RNA > quantile(ctrl.15.seu$nCount_RNA, 0.01) & 
                        # nCount_RNA < quantile(ctrl.15.seu$nCount_RNA, 0.99) &
                        # percent.mt < quantile(ctrl.15.seu$percent.mt, 0.99))
                        # percent.mt <2)

# VlnPlot(water.1.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(water.1.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(water.1.seu$nCount_RNA, c(0.01, 0.99))
# quantile(water.1.seu$percent.mt, 0.99)
# water.1.seu <- subset(water.1.seu,
#                       subset = nFeature_RNA > quantile(water.1.seu$nFeature_RNA, 0.01) &
#                         nFeature_RNA < quantile(water.1.seu$nFeature_RNA, 0.99) &
#                         nCount_RNA > quantile(water.1.seu$nCount_RNA, 0.01) & 
#                         nCount_RNA < quantile(water.1.seu$nCount_RNA, 0.99) &
#                         # percent.mt < quantile(water.1.seu$percent.mt, 0.99))
#                         percent.mt <2)

# VlnPlot(water.5.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(water.5.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(water.5.seu$nCount_RNA, c(0.01, 0.99))
# quantile(water.5.seu$percent.mt, 0.99)
# water.5.seu <- subset(water.5.seu,
#                       subset = nFeature_RNA > quantile(water.5.seu$nFeature_RNA, 0.01) &
#                         nFeature_RNA < quantile(water.5.seu$nFeature_RNA, 0.99) &
#                         nCount_RNA > quantile(water.5.seu$nCount_RNA, 0.01) &
#                         nCount_RNA < quantile(water.5.seu$nCount_RNA, 0.99) &
#                         # percent.mt < quantile(water.5.seu$percent.mt, 0.99))
#                         percent.mt<2)
# 
# VlnPlot(water.6.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(water.6.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(water.6.seu$nCount_RNA, c(0.01, 0.99))
# quantile(water.6.seu$percent.mt, 0.99)
# water.6.seu <- subset(water.6.seu,
#                       subset = nFeature_RNA > quantile(water.6.seu$nFeature_RNA, 0.01) &
#                         nFeature_RNA < quantile(water.6.seu$nFeature_RNA, 0.99) &
#                         nCount_RNA > quantile(water.6.seu$nCount_RNA, 0.01) & 
#                         nCount_RNA < quantile(water.6.seu$nCount_RNA, 0.99) &
#                         # percent.mt < quantile(water.6.seu$percent.mt, 0.99))
#                         percent.mt <2)
# 
# VlnPlot(water.9.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(water.9.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(water.9.seu$nCount_RNA, c(0.01, 0.99))
# quantile(water.9.seu$percent.mt, 0.99)
# water.9.seu <- subset(water.9.seu,
#                       subset = nFeature_RNA > quantile(water.9.seu$nFeature_RNA, 0.01) &
#                         nFeature_RNA < quantile(water.9.seu$nFeature_RNA, 0.99) &
#                         nCount_RNA > quantile(water.9.seu$nCount_RNA, 0.01) & 
#                         nCount_RNA < quantile(water.9.seu$nCount_RNA, 0.99) &
#                         # percent.mt < quantile(water.9.seu$percent.mt, 0.99))
#                         percent.mt <2)
# 
# VlnPlot(blum.1.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(blum.1.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(blum.1.seu$nCount_RNA, c(0.01, 0.99))
# quantile(blum.1.seu$percent.mt, 0.99)
# blum.1.seu <- subset(blum.1.seu,
#                      subset = nFeature_RNA > quantile(blum.1.seu$nFeature_RNA, 0.01) &
#                        nFeature_RNA < quantile(blum.1.seu$nFeature_RNA, 0.99) &
#                        nCount_RNA > quantile(blum.1.seu$nCount_RNA, 0.01) & 
#                        nCount_RNA < quantile(blum.1.seu$nCount_RNA, 0.99) &
#                        # percent.mt < quantile(blum.1.seu$percent.mt, 0.99))
#                        percent.mt <2)
# 
# VlnPlot(blum.3.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(blum.3.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(blum.3.seu$nCount_RNA, c(0.01, 0.99))
# quantile(blum.3.seu$percent.mt, 0.99)
# blum.3.seu <- subset(blum.3.seu,
#                      subset = nFeature_RNA > quantile(blum.3.seu$nFeature_RNA, 0.01) &
#                        nFeature_RNA < quantile(blum.3.seu$nFeature_RNA, 0.99) &
#                        nCount_RNA > quantile(blum.3.seu$nCount_RNA, 0.01) & 
#                        nCount_RNA < quantile(blum.3.seu$nCount_RNA, 0.99) &
#                        # percent.mt < quantile(blum.3.seu$percent.mt, 0.99))
#                        percent.mt <2)
# 
# VlnPlot(blum.4.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(blum.4.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(blum.4.seu$nCount_RNA, c(0.01, 0.99))
# quantile(blum.4.seu$percent.mt, 0.99)
# blum.4.seu <- subset(blum.4.seu,
#                      subset = nFeature_RNA > quantile(blum.4.seu$nFeature_RNA, 0.01) &
#                        nFeature_RNA < quantile(blum.4.seu$nFeature_RNA, 0.99) &
#                        nCount_RNA > quantile(blum.4.seu$nCount_RNA, 0.01) & 
#                        nCount_RNA < quantile(blum.4.seu$nCount_RNA, 0.99) &
#                        # percent.mt < quantile(blum.4.seu$percent.mt, 0.99))
#                        percent.mt <2)
# 
# VlnPlot(blum.7.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
# quantile(blum.7.seu$nFeature_RNA, c(0.01, 0.99))
# quantile(blum.7.seu$nCount_RNA, c(0.01, 0.99))
# quantile(blum.7.seu$percent.mt, 0.99)
# blum.7.seu <- subset(blum.7.seu,
#                      subset = nFeature_RNA > quantile(blum.7.seu$nFeature_RNA, 0.01) &
#                        nFeature_RNA < quantile(blum.7.seu$nFeature_RNA, 0.99) &
#                        nCount_RNA > quantile(blum.7.seu$nCount_RNA, 0.01) &
#                        nCount_RNA < quantile(blum.7.seu$nCount_RNA, 0.99) &
#                        # percent.mt < quantile(blum.7.seu$percent.mt, 0.99))
#                        percent.mt <2)


# Normalize, FindVariableFeatures, ScaleData, PCA (pK optimization needs this) ----
ctrl.11.seu <- NormalizeData(ctrl.11.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(ctrl.11.seu)
ctrl.12.seu <- NormalizeData(ctrl.12.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(ctrl.12.seu)
ctrl.14.seu <- NormalizeData(ctrl.14.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(ctrl.14.seu)
ctrl.15.seu <- NormalizeData(ctrl.15.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(ctrl.15.seu)
water.1.seu <- NormalizeData(water.1.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(water.1.seu)
water.5.seu <- NormalizeData(water.5.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(water.5.seu)
water.6.seu <- NormalizeData(water.6.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(water.6.seu)
water.9.seu <- NormalizeData(water.9.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(water.9.seu)
blum.1.seu <- NormalizeData(blum.1.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(blum.1.seu)
blum.3.seu <- NormalizeData(blum.3.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(blum.3.seu)
blum.4.seu <- NormalizeData(blum.4.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# ElbowPlot(blum.4.seu)
blum.7.seu <- NormalizeData(blum.7.seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(blum.7.seu)

# Cluster ----

ctrl.11.seu <- FindNeighbors(ctrl.11.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
ctrl.12.seu <- FindNeighbors(ctrl.12.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
ctrl.14.seu <- FindNeighbors(ctrl.14.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
ctrl.15.seu <- FindNeighbors(ctrl.15.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
water.1.seu <- FindNeighbors(water.1.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
water.5.seu <- FindNeighbors(water.5.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
water.6.seu <- FindNeighbors(water.6.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
water.9.seu <- FindNeighbors(water.9.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
blum.1.seu <- FindNeighbors(blum.1.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
blum.3.seu <- FindNeighbors(blum.3.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
blum.4.seu <- FindNeighbors(blum.4.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
blum.7.seu <- FindNeighbors(blum.7.seu, dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)

# Doublets ----

# Control11
sweep.res <- paramSweep(ctrl.11.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.11.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.15*length(Cells(ctrl.11.seu))) #Assuming 15% doublet formation rate based on .8% doublets per 1000 cells and 16622 cells recovered
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.11.seu <- doubletFinder(ctrl.11.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.11.seu@meta.data), value = TRUE)

DimPlot(ctrl.11.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.11.seu@meta.data[[df_col]])

#subset data to keep only singlets
# ctrl.11.seu <- subset(ctrl.11.seu, subset = ctrl.11.seu@meta.data[[df_col]] == 'Singlet')
ctrl.11.seu <- ctrl.11.seu[, ctrl.11.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.11.seu@meta.data[[df_col]])

# Control12
sweep.res <- paramSweep(ctrl.12.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.12.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.17*length(Cells(ctrl.12.seu))) #Assuming 17% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.12.seu <- doubletFinder(ctrl.12.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.12.seu@meta.data), value = TRUE)

DimPlot(ctrl.12.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.12.seu@meta.data[[df_col]])

#subset data to keep only singlets
ctrl.12.seu <- ctrl.12.seu[, ctrl.12.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.12.seu@meta.data[[df_col]])

# Control14
sweep.res <- paramSweep(ctrl.14.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.14.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.19*length(Cells(ctrl.14.seu))) #Assuming 19% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.14.seu <- doubletFinder(ctrl.14.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.14.seu@meta.data), value = TRUE)

DimPlot(ctrl.14.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.14.seu@meta.data[[df_col]])

#subset data to keep only singlets
ctrl.14.seu <- ctrl.14.seu[, ctrl.14.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.14.seu@meta.data[[df_col]])

# Control15
sweep.res <- paramSweep(ctrl.15.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.15.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.20*length(Cells(ctrl.15.seu))) #Assuming 20% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.15.seu <- doubletFinder(ctrl.15.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.15.seu@meta.data), value = TRUE)

DimPlot(ctrl.15.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.15.seu@meta.data[[df_col]])

#subset data to keep only singlets
ctrl.15.seu <- ctrl.15.seu[, ctrl.15.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.15.seu@meta.data[[df_col]])

# Water1
sweep.res <- paramSweep(water.1.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.1.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.15*length(Cells(water.1.seu))) #Assuming 15% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.1.seu <- doubletFinder(water.1.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.1.seu@meta.data), value = TRUE)

DimPlot(water.1.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.1.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.1.seu <- water.1.seu[, water.1.seu@meta.data[[df_col]] == "Singlet"]

table(water.1.seu@meta.data[[df_col]])

# Water5
sweep.res <- paramSweep(water.5.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.5.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.19*length(Cells(water.5.seu))) #Assuming 19% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.5.seu <- doubletFinder(water.5.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.5.seu@meta.data), value = TRUE)

DimPlot(water.5.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.5.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.5.seu <- water.5.seu[, water.5.seu@meta.data[[df_col]] == "Singlet"]

table(water.5.seu@meta.data[[df_col]])

# Water6
sweep.res <- paramSweep(water.6.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.6.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.18*length(Cells(water.6.seu))) #Assuming 18% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.6.seu <- doubletFinder(water.6.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.6.seu@meta.data), value = TRUE)

DimPlot(water.6.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.6.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.6.seu <- water.6.seu[, water.6.seu@meta.data[[df_col]] == "Singlet"]

table(water.6.seu@meta.data[[df_col]])

# Water9
sweep.res <- paramSweep(water.9.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.9.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.16*length(Cells(water.9.seu))) #Assuming 16% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.9.seu <- doubletFinder(water.9.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.9.seu@meta.data), value = TRUE)

DimPlot(water.9.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.9.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.9.seu <- water.9.seu[, water.9.seu@meta.data[[df_col]] == "Singlet"]

table(water.9.seu@meta.data[[df_col]])

# Blumeria1
sweep.res <- paramSweep(blum.1.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.1.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.17*length(Cells(blum.1.seu))) #Assuming 17% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.1.seu <- doubletFinder(blum.1.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.1.seu@meta.data), value = TRUE)

DimPlot(blum.1.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.1.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.1.seu <- blum.1.seu[, blum.1.seu@meta.data[[df_col]] == "Singlet"]

table(blum.1.seu@meta.data[[df_col]])

# Blumeria3
sweep.res <- paramSweep(blum.3.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.3.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.18*length(Cells(blum.3.seu))) #Assuming 18% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.3.seu <- doubletFinder(blum.3.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.3.seu@meta.data), value = TRUE)

DimPlot(blum.3.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.3.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.3.seu <- blum.3.seu[, blum.3.seu@meta.data[[df_col]] == "Singlet"]

table(blum.3.seu@meta.data[[df_col]])

# Blumeria4
sweep.res <- paramSweep(blum.4.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.4.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.16*length(Cells(blum.4.seu))) #Assuming 16% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.4.seu <- doubletFinder(blum.4.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.4.seu@meta.data), value = TRUE)

DimPlot(blum.4.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.4.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.4.seu <- blum.4.seu[, blum.4.seu@meta.data[[df_col]] == "Singlet"]

table(blum.4.seu@meta.data[[df_col]])

# Blumeria7
sweep.res <- paramSweep(blum.7.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.7.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.20*length(Cells(blum.7.seu))) #Assuming 20% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.7.seu <- doubletFinder(blum.7.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.7.seu@meta.data), value = TRUE)

DimPlot(blum.7.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.7.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.7.seu <- blum.7.seu[, blum.7.seu@meta.data[[df_col]] == "Singlet"]

table(blum.7.seu@meta.data[[df_col]])


# Add metadata ----

ctrl.11.seu$sample <- "Control11"
ctrl.11.seu$condition <- "Control"

ctrl.12.seu$sample <- "Control12"
ctrl.12.seu$condition <- "Control"

ctrl.14.seu$sample <- "Control14"
ctrl.14.seu$condition <- "Control"

ctrl.15.seu$sample <- "Control15"
ctrl.15.seu$condition <- "Control"

water.1.seu$sample <- "Water1"
water.1.seu$condition <- "MCT-Water"

water.5.seu$sample <- "Water5"
water.5.seu$condition <- "MCT-Water"

water.6.seu$sample <- "Water6"
water.6.seu$condition <- "MCT-Water"

water.9.seu$sample <- "Water9"
water.9.seu$condition <- "MCT-Water"

blum.1.seu$sample <- "Blumeria1"
blum.1.seu$condition <- "MCT-Blumeria"

blum.3.seu$sample <- "Blumeria3"
blum.3.seu$condition <- "MCT-Blumeria"

blum.4.seu$sample <- "Blumeria4"
blum.4.seu$condition <- "MCT-Blumeria"

blum.7.seu$sample <- "Blumeria7"
blum.7.seu$condition <- "MCT-Blumeria"

# Combine and prepare for Seurat v5 integration ----

seu_merge <- merge(
  x = ctrl.11.seu,
  y = list(
    ctrl.12.seu, ctrl.14.seu, ctrl.15.seu,
    water.1.seu, water.5.seu, water.6.seu, water.9.seu,
    blum.1.seu, blum.3.seu, blum.4.seu, blum.7.seu
  ),
  add.cell.ids = c(
    "Control11", "Control12", "Control14", "Control15",
    "Water1", "Water5", "Water6", "Water9",
    "Blumeria1", "Blumeria3", "Blumeria4", "Blumeria7"
  ),
  project = "All.Samples"
)

seu_merge$condition <- factor(
  seu_merge$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)

seu_merge$sample <- as.character(seu_merge$sample)

seu_merge <- NormalizeData(seu_merge) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

saveRDS(
  seu_merge,
  file = "03-analysis_scratch/seu_v5_preintegration.rds"
)




