# The purpose of this script is to run shared-nearest neighbor clustering on the 
# Harmony-corrected integrated cells
#-------------------------------------------------------------------------------#

# Clear environment
rm(list = ls())

# Load in libraries
library(imcRtools)
library(cytomapper)
library(RColorBrewer)
library(dittoSeq)
library(viridis)
library(gridGraphics)
library(cowplot)
library(tidyverse)
library(scater)
library(scran) 
library(patchwork)
library(harmony)
library(BiocSingular)
library(bluster)
library(scuttle)

# Prepare an output directory
opDir <- "Outputs/002_SNN_Clustering_Outputs"
if (!dir.exists(opDir)) {
  dir.create(opDir)
}

# Prepare a data subdirectory
dataDir <- "Data/Processed_Spe_Object/002_SNN_Clustering"
if (!dir.exists(dataDir)) {
  dir.create(dataDir)
}

################################################################################

# Read in the dimRedux spatialExperiment object
spe <- readRDS("Data/Processed_Spe_Object/DimRedux_Spe.Rds")

# Run SNN clustering
set.seed(03061999)
clustering <- clusterCells(spe[rowData(spe)$use_channel],
                           use.dimred = "harmony",
                           BLUSPARAM = SNNGraphParam(k = 20,
                                                     cluster.fun = "louvain",
                                                     type = "rank"))

# Append clusters to the spe object
spe$NN <- clustering

# Assign cluster colors
clusterCols <- setNames(brewer.pal(12, "Paired"), unique(spe$NN))

# Append cluster colors to metadata
metadata(spe)$color_vectors$colors <- clusterCols


dittoDimPlot(spe, var = "Indication", reduction.use = "UMAP_Harmony", size = 0.2) +
  scale_color_manual(values = metadata(spe)$color_vectors$Indication) +
  labs(color = "Indication",
       title = "UMAP: Indication") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 18, face = "bold"),
        legend.position = "bottom")

################################################################################

# Plot the SNN clusters on Harmony-corrected UMAP
SNN_UMAP <- dittoDimPlot(spe, 
                         var = "NN",
                         reduction.use = "UMAP_Harmony",
                         size = 0.2,
                         do.label = TRUE) +
  ggtitle("UMAP Harmony: SNN") +
  theme(axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 18, face = "bold"),
        legend.position = "right")
ggsave("Outputs/002_SNN_Clustering_Outputs/SNN_Clusters_Harmony_UMAP.tiff", SNN_UMAP, width = 6, height = 6, dpi = 300)

# Plot markers on UMAP
toPlot <- multi_dittoDimPlot(spe, var = rownames(spe)[rowData(spe)$use_channel], reduction.use = "UMAP_Harmony",
                             assay = "exprs", size = 0.2, list.out = TRUE)

plot_list <- lapply(toPlot, function(x) x + scale_color_viridis())
markers <- plot_grid(plotlist = plot_list)
ggsave("Outputs/002_SNN_Clustering_Outputs/Markers_UMAP.tiff", markers, width = 10, height = 10, dpi = 300)


################################################################################

# Aggregate across cells
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$NN,
                                   statistics = "median",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)

# Plot heatmap Z-scaled
ZScaled <- dittoHeatmap(image_mean, genes = rownames(image_mean)[rowData(image_mean)$use_channel],
                        assay = "exprs", 
                        cluster_cols = TRUE,
                        cluster_rows = TRUE, 
                        show_colnames = TRUE,
                        heatmap.colors = 
                        rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
                        annot.by = c("NN", "ncells"))
ggsave("Outputs/002_SNN_Clustering_Outputs/Median_ZScaled_SNN_Clusters.tiff", ZScaled, width = 6, height = 6, dpi = 300)

################################################################################

# Prepare an output directory to hold the plotCells outputs
cellsDir <- "Outputs/002_SNN_Clustering_Outputs/Plot_Cells"
if (!dir.exists(cellsDir)) {
  dir.create(cellsDir)
}

# Prepare an output directory to hold the plotPixels outputs
pixelsDir <- "Outputs/002_SNN_Clustering_Outputs/Plot_Pixels"
if (!dir.exists(pixelsDir)) {
  dir.create(pixelsDir)
}

# Prepare a data subdirectory to hold the normalized images
normImgDir <- "Data/Processed_Images/Normalized_Images"
if (!dir.exists(normImgDir)) {
  dir.create(normImgDir)
}

# Prepare an output directory for Lymphocyte plotPixels
lympho <- "Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Lymphocyte_Populations"
if (!dir.exists(lympho)) {
  dir.create(lympho)
}

# Prepare an output directory for structure plotPixels
structure <- "Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Structure_Populations"
if (!dir.exists(structure)) {
  dir.create(structure)
}

# Prepare an output directory for activated T-cells
cytoT <- "Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Cytotoxic_T"
if (!dir.exists(cytoT)) {
  dir.create(cytoT)
}

# Prepare an output directory for T-herlper/Tregs
helper <- "Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Th_Treg"
if (!dir.exists(helper)) {
  dir.create(helper)
}

# Visualize the clusters on the images
masks <- readRDS("Data/Processed_Images/WBH_Masks.Rds")
cur_ids <- names(masks)

# Plot cells belonging to each cluster one by one
for (i in cur_ids) {
  clusterVal <- plotCells(mask = masks[i],
                          object = spe,
                          cell_id = "ObjectNumber",
                          img_id = "sample_id",
                          colour_by = "NN",
                          legend = list(colour_by.labels.cex = 1.05,
                                        colour_by.legend.cex = 1.05,
                                        margin = 1),
                          return_plot = TRUE)
  outputPlot <- ggdraw(clusterVal$plot, clip = "on")
  fileName <- paste(i, "NN_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs/Plot_Cells", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

# Read in the images
images <- readRDS("Data/Processed_Images/WBH_Raw_Images.Rds")

# Normalize the images
img <- cytomapper::normalize(images, separateImages = TRUE)
images <- cytomapper::normalize(img, inputRange = c(0,0.1))

# Save the normalized images
saveRDS(images, file = "Data/Processed_Images/Normalized_Images/Normalized_WBH_Images.Rds")
images <- readRDS("Data/Processed_Images/Normalized_Images/Normalized_WBH_Images.Rds")

# Visualize some key lymphocyte markers (pixels)
for (i in cur_ids) {
  clusterVal <- plotPixels(images[i],
                           colour_by = c("Ecadherin", "CD3", "CD8a", "CD4", "FoxP3", "CD163"),
                           return_plot = TRUE)
  outputPlot <- ggdraw(clusterVal$plot, clip = "on")
  fileName <- paste(i, "NN_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Lymphocyte_Populations/", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


# Visualize some key structure markers (pixels)
for (i in cur_ids) {
  clusterVal <- plotPixels(images[i],
                           colour_by = c("Ecadherin", "aSMA", "Ki67"),
                           return_plot = TRUE)
  outputPlot <- ggdraw(clusterVal$plot, clip = "on")
  fileName <- paste(i, "NN_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Structure_Populations//", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


# Visualize cytototoxic T
for (i in cur_ids) {
  clusterVal <- plotPixels(images[i],
                           colour_by = c("Ecadherin", "CD3", "CD8a"),
                           return_plot = TRUE)
  outputPlot <- ggdraw(clusterVal$plot, clip = "on")
  fileName <- paste(i, "NN_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Cytotoxic_T", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

# Visualize T-helper/Tregs
for (i in cur_ids) {
  clusterVal <- plotPixels(images[i],
                           colour_by = c("Ecadherin", "CD8a", "CD4"),
                           colour = list(Ecadherin = c("black", "red"),
                                         CD8a = c("black", "green"),
                                         CD4 = c("black", "blue")),
                           return_plot = TRUE)
  outputPlot <- ggdraw(clusterVal$plot, clip = "on")
  fileName <- paste(i, "NN_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs/Plot_Pixels/Th_Treg", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

################################################################################

# Prepare an output directory for cell types on tissue
cellTypeDir <- "Outputs/002_SNN_Clustering_Outputs/CellTypes"
if (!dir.exists(cellTypeDir)) {
  dir.create(cellTypeDir)
}

# Based on the images, assign phenotypes
spe$CellType <- ifelse(spe$NN == 1, "Stroma",
                       ifelse(spe$NN == 2, "Granulocytes",
                              ifelse(spe$NN == 3, "Proliferating Epithelia",
                                     ifelse(spe$NN == 4, "Stroma",
                                            ifelse(spe$NN == 5, "Undefined",
                                                   ifelse(spe$NN == 6, "Stroma",
                                                          ifelse(spe$NN == 7, "Lymphocytes",
                                                                 ifelse(spe$NN == 8, "Macrophage",
                                                                        ifelse(spe$NN == 9, "Lymphocytes",
                                                                               ifelse(spe$NN == 10, "Tumor", 
                                                                                      ifelse(spe$NN == 11, "Proliferating Epithelia",
                                                                                             ifelse(spe$NN == 12, "Tumor", "Undefined"))))))))))))

# Assign color vector to CellTypes
celltype <- setNames(c("#3F1B03", "#F4AD31", "#894F36", "#1C750C", "#EF8ECC", 
                       "#6471E2", "#4DB23B", "grey"),
                     c("Stroma", "Granulocytes", "Proliferating Epithelia", "Lymphocytes", "Macrophage", 
                       "Tumor", "IELs", "Undefined"))
metadata(spe)$color_vectors$CellType <- celltype

# Visualize the clusters on the tissue
# Plot cells belonging to each cluster one by one
for (i in cur_ids) {
  clusterVal <- plotCells(mask = masks[i],
                          object = spe,
                          cell_id = "ObjectNumber",
                          img_id = "sample_id",
                          colour_by = "CellType",
                          colour = list(CellType = metadata(spe)$color_vectors$CellType),
                          legend = list(colour_by.labels.cex = 1,
                                        colour_by.legend.cex = 1,
                                        margin = 2),
                          return_plot = TRUE)
  outputPlot <- ggdraw(clusterVal$plot, clip = "on")
  fileName <- paste(i, "CellType_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs/CellTypes", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


for (i in cur_ids) {
OUT <- plotCells(masks[i],
          object = spe, 
          cell_id = "ObjectNumber", 
          img_id = "sample_id",
          colour_by = c("Ecadherin", "CD3", "CD8a", "CD4"),
          exprs_values = "exprs",
          return_plot = TRUE)
  outputPlot <- ggdraw(OUT$plot, clip = "on")
  fileName <- paste(i, "CellType_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}


for (i in cur_ids) {
  OUT <- plotPixels(image = images[i],
                   mask = masks[i],
                   object = spe, 
                   cell_id = "ObjectNumber", 
                   img_id = "sample_id",
                   colour_by = c("Ecadherin", "CD3", "CD4"),
                   outline_by = "CellType",
                   bcg = list(Ecadherin = c(0, 5, 1),
                              CD3 = c(0, 5, 1),
                              CD4 = c(0, 5, 1)),
                   colour = list(celltype = metadata(spe)$color_vectors$CellType),
                   thick = TRUE,
                   return_plot = TRUE)
  outputPlot <- ggdraw(OUT$plot, clip = "on")
  fileName <- paste(i, "CellType_Clusters.tiff", sep = "_")
  ggsave(paste("Outputs/002_SNN_Clustering_Outputs", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)
}

  

################################################################################

# Aggregate across cells
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$CellType,
                                   statistics = "mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)

# Plot heatmap Z-scaled
ZScaled <- dittoHeatmap(image_mean, genes = rownames(image_mean)[rowData(image_mean)$use_channel],
                        assay = "exprs", 
                        cluster_cols = TRUE,
                        cluster_rows = TRUE, 
                        show_colnames = TRUE,
                        heatmap.colors = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
                        annot.by = c("CellType", "ncells"),
                        annotation_colors = list(CellType = metadata(spe)$color_vectors$CellType,
                                                 ncells = rev(magma(100))))
ggsave("Outputs/002_SNN_Clustering_Outputs/ZScaled_CellTypes_Clusters.tiff", ZScaled, width = 6, height = 6, dpi = 300)


# Save Spe object
saveRDS(spe, file = "Data/Processed_Spe_Object/002_SNN_Clustering/SNN_Clustered_Spe.Rds")

################################################################################

# Differential marker expression
clusterNames <- unique(spe$CellType)

# Split spe into early and late
early <- spe[,spe$Indication == "EOCRC"]
late <- spe[,spe$Indication == "LOCRC"]

# Subset the early spatial experiment object
EOsamples <- unique(early$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(EOsamples), ncol = length(clusterNames))
propDataE <- as.data.frame(propData)
rownames(propDataE) <- EOsamples
colnames(propDataE) <- clusterNames

# Calculate the proportion of each cell type per patient
for(i in EOsamples) {
  patient_cells <- sum(early$sample_id == i)
  
  # Count the number of cells belonging to each cellType for each tissue
  for(j in clusterNames) {
    cluster_cells <- sum(early$sample_id == i & early$CellType == j)
    proportion <- cluster_cells / patient_cells 
    
    # Append to the results matrix
    propDataE[i,j] <- proportion
    
  }
}

# Add a group column
propDataE$Group <- "EOCRC"

# Pivot longer 
propDataE <- propDataE %>%
  pivot_longer(-c("Group"), names_to = "CellTypes", values_to = "values")


# Repeat this process for the late samples
LOsamples <- unique(late$sample_id)

# Initialize an empty matrix
propData <- matrix(0, nrow = length(LOsamples), ncol = length(clusterNames))
propDataL <- as.data.frame(propData)
rownames(propDataL) <- LOsamples
colnames(propDataL) <- clusterNames

# Iterate through the samples and cellTypes
for(i in LOsamples) {
  patient_cells <- sum(late$sample_id == i)
  
  for(j in clusterNames) {
    cluster_cells <- sum(late$sample_id == i & late$CellType == j)
    proportion <- cluster_cells / patient_cells 
    
    propDataL[i,j] <- proportion
  }
}

# Append a group column
propDataL$Group <- "LOCRC"

# Pivot longer
propDataL <- propDataL %>%
  pivot_longer(-c("Group"), names_to = "CellTypes", values_to = "values")

# Combine the two results dataframes
propData <- rbind(propDataE, propDataL)

# Rename the cell types for graphing
propData$CellTypes <- ifelse(propData$CellTypes == "Non-proliferating epithelia", "NPE", propData$CellTypes)
propData$CellTypes <- ifelse(propData$CellTypes == "Proliferating epithelia", "PE", propData$CellTypes)
propData$CellTypes <- ifelse(propData$CellTypes == "Cytotoxic immune cells", "Cytotoxic", propData$CellTypes)


# Set colors for early and late
custom_colors <- c("EOCRC" = "#1B9E77", "LOCRC" = "#7570B3")

# Plot
props <- ggplot(propData, aes(x = Group, y = values, fill = Group)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~CellTypes) +
  theme_bw() +
  labs(y = "Expression",
       title = "Proportion of Cell Types") +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(method = "t.test", label = "p.format") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 20),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        title = element_text(size = 18, face = "bold"))
ggsave("Outputs/004_Run_Clustering_Outputs/Phenotyped/Proportions_Per_Group.tiff", props, width = 10, height = 10, dpi = 300)

################################################################################

bp <- dittoBarPlot(spe, 
                   var = "CellType", 
                   group.by = "sample_id",
                   split.by = "Indication",
                   split.adjust = list(scales = "free")) +
  scale_fill_manual(values = metadata(spe)$color_vectors$CellType) +
  labs(y = "Percent of cells",
       title = "",
       fill = "Cell Type") +
  scale_y_continuous(labels = percent) +
  theme_classic() +
  theme(strip.text = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 12))

