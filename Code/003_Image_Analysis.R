# The purpose of this script is to perform image analysis on the WBH EOLO Samples

# Clear the environment
rm(list = ls())

# Load libraries
library(imcRtools)
library(immunoCluster)
library(flowCore)
library(cytomapper)

# Generate an output directory
opDir <- "Outputs/003_Image_Analysis_Outputs"
if (!dir.exists(opDir)){
  dir.create(opDir)
}

################################################################################

# Read in the clustered spe object
spe <- readRDS("Data/Processed_Spe_Object/002_SNN_Clustering/SNN_Clustered_Spe.Rds")

# Read in the normalized images and masks
images <- readRDS("Data/Processed_Images/WBH_Raw_Images.Rds")

# Normalize the images
img <- cytomapper::normalize(images, separateImages = TRUE)
images <- cytomapper::normalize(img, inputRange = c(0,0.1))

# Save the normalized images
saveRDS(images, file = "Data/Processed_Images/Normalized_Images/WBH_Normalized_Images.Rds")

# Read in the masks
masks <- readRDS("Data/Processed_Images/WBH_Masks.Rds")

################################################################################

# Differential marker abundance
image_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel], "SingleCellExperiment"), 
                                   ids = spe$sample_id,
                                   statistics = "mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean)/1)

dittoHeatmap(image_mean, genes = rownames(image_mean)[rowData(image_mean)$use_channel],
             assay = "exprs", 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             show_colnames = TRUE,
             heatmap.colors = 
               rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
             annot.by = c("sample_id", "Indication", "ncells"))
################################################################################

# Split the Spe object into EO and LO
eo <- spe[,spe$Indication == "EOCRC",]

# Aggregate
eo_mean <- aggregateAcrossCells(as(eo[rowData(eo)$use_channel], "SingleCellExperiment"), 
                                   ids = eo$sample_id,
                                   statistics = "mean",
                                   use.assay.type = "counts")
assay(eo_mean, "exprs") <- asinh(counts(eo_mean)/1)
EO <- as.data.frame(t(assay(eo_mean, "exprs")))
EO$Group <- "EOCRC"
EO <- EO %>%
  pivot_longer(-c("Group"), names_to = "Marker", values_to = "values")

lo <- spe[,spe$Indication == "LOCRC",]

# Aggregate
lo_mean <- aggregateAcrossCells(as(lo[rowData(lo)$use_channel], "SingleCellExperiment"), 
                                ids = lo$sample_id,
                                statistics = "mean",
                                use.assay.type = "counts")
assay(lo_mean, "exprs") <- asinh(counts(lo_mean)/1)
LO <- as.data.frame(t(assay(lo_mean, "exprs")))
LO$Group <- "LOCRC"
LO <- LO %>%
  pivot_longer(-c("Group"), names_to = "Marker", values_to = "values")


results <- rbind(EO, LO)

difExpr <- ggplot(results, aes(x = Group, y = values, fill = Group)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~Marker) +
  theme_bw() +
  labs(y = "Expression",
       title = "Mean Marker Expression") +
  stat_compare_means(method = "t.test", label = "p.format") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 20),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        title = element_text(size = 18, face = "bold"))
ggsave("Outputs/003_Image_Analysis_Outputs/Differential_Marker_Expression.tiff", difExpr, width = 14, height = 14, dpi = 300)

################################################################################

# Get a vector of early onset samples and late onset samples
eocrc <- unique(spe[,spe$Indication == "EOCRC"]$sample_id)
locrc <- unique(spe[,spe$Indication == "LOCRC"]$sample_id)

# Subset the images into two groups
eoImg <- images[names(images) %in% eocrc]
loImg <- images[names(images) %in% locrc]

# Generate an output subdirectory to hold the CD3+, CD4+ T cells
Th <- "Outputs/003_Image_Analysis_Outputs/CD3_CD4_Pos"
if (!dir.exists(Th)) {
  dir.create(Th)
}


# For early onset
clusterVal <- plotPixels(eoImg,
                         colour_by = c("Ecadherin", "CD3", "CD4"),
                         colour = list(Ecadherin = c("black", "burlywood"),
                                       CD3 = c("black", "cyan2"),
                                       CD4 = c("black", "firebrick1")),
                         return_plot = TRUE)
outputPlot <- ggdraw(clusterVal$plot, clip = "on")
fileName <- paste("EarlyOnset", "CD3_CD4.tiff", sep = "_")
ggsave(paste("Outputs/003_Image_Analysis_Outputs/CD3_CD4_Pos", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)

# For late onset
clusterVal <- plotPixels(loImg,
                         colour_by = c("Ecadherin", "CD3", "CD4"),
                         colour = list(Ecadherin = c("black", "burlywood"),
                                       CD3 = c("black", "cyan2"),
                                       CD4 = c("black", "firebrick1")),
                         return_plot = TRUE)
outputPlot <- ggdraw(clusterVal$plot, clip = "on")
fileName <- paste("LateOnset", "CD3_CD4.tiff", sep = "_")
ggsave(paste("Outputs/003_Image_Analysis_Outputs/CD3_CD4_Pos", fileName, sep = "/"), outputPlot, width = 15, height = 15, dpi = 300)





