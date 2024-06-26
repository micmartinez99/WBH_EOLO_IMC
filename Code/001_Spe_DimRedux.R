# The purpose of this script is to prepare a spatialExperiment object from the 
# combined Steinbock output and perform dimensionality reduction
# These samples are Early and Late onset tissues from Waterbury Hospital 
# preserved in LWRN (contrasted to JDH -- archived, not included in analysis)
#-------------------------------------------------------------------------------#

# Clear environment
rm(list = ls())

# Load libraries
library(imcRtools) # For spatialExperiment object generation
library(cytomapper) # For plotting pixels and cells
library(tidyverse) # Data manipulation
library(RColorBrewer) # Color palettes
library(viridis) # Color palettes
library(scater) # Single-cell analysis
library(harmony) # Batch correction
library(BiocSingular) # Batch correction
library(patchwork) # Plotting
library(cowplot) # Plotting
library(dittoSeq) # Single-cell analysis

# Prepare an output directory to hold the processed spe object
dataDir <- "Data/Processed_Spe_Object"
if (!dir.exists(dataDir)) {
  dir.create(dataDir)
}

# Prepare an output directory for dimensionality reduction 
dimReduxDir <- "Outputs/001_Spe_DimRedux_Outputs"
if (!dir.exists(dimReduxDir)) {
  dir.create(dimReduxDir)
}

################################################################################

# Read in the Steinbock folder
spe <- read_steinbock("Data/combined_steinbock/")

# Make sample_id more interpretable
spe$sample_id <- sub("Slide1-6_CRCTMA1_CRCTMA2_", "ROI_", spe$sample_id)
spe$sample_id <- sub("Slide2-5_CRCTMA3_CRCTMA4_", "ROI_", spe$sample_id)

# Check that there are 26 unique sample_id
length(unique(spe$sample_id))

# For any markers that have a hyphen in their name, replace the hypen with an underscore
rownames(spe) <- gsub("\\-", "_", rownames(spe))

# Specify a subset of the channels that will be used for all downstream analysis
rowData(spe)$use_channel <- !grepl("DNA1|DNA2|ICSK1|ICSK2|ICSK3|CD14|CD31|CD11c|CD45RA|Perforin|Pan_actin|CD45RA", rownames(spe))

################################################################################

# Metadata processing
metadata <- read.csv("Data/Updated_Metadata.csv")

# Add a column for Indication
spe$Indication <- metadata$Indication[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Patient ID
spe$Patient_ID <- metadata$Patient_ID[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for Source
spe$Source <- metadata$Source[match(spe$sample_id, metadata$Sample_ID)]

# Add a column for "Keep"
spe$Keep <- metadata$Keep[match(spe$sample_id, metadata$Sample_ID)]

# Subset the spe object based on the "keep" vector
spe <- spe[,spe$Source == "WBH"]
unique(spe$sample_id)

# Set cell names
colnames(spe) <- paste0('cell', seq_len(ncol(spe)))

################################################################################

# Initialize an empty list to hold the color vectors in the spe object
color_vectors <- list()

# Set colors for EOCRC and LOCRC
indication_colors <- c("EOCRC" = "#1B9E77",
                       "LOCRC" = "#7570B3")

cols <- brewer.pal(10, "Paired")
# Set colors for Patient_ID
cols <- c("grey", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
          "#FFFF99", "#B15928", "#FDCDAC", "#B3CDE3", "black")
patient_colors <- setNames(cols, unique(spe$Patient_ID)) 

# Append vectors to color list
color_vectors$Indication <- indication_colors
color_vectors$Patient_ID <- patient_colors

# Append color list to metadata slot in spe
metadata(spe)$color_vectors <- color_vectors

################################################################################

# Read in the Images
img <- loadImages("Data/combined_steinbock/img/")

# Adjust sample IDs for images/masks to match spe
names(img) <- sub("Slide1-6_CRCTMA1_CRCTMA2_", "ROI_", names(img))
names(img) <- sub("Slide2-5_CRCTMA3_CRCTMA4_", "ROI_", names(img))

# Set a vector of the samples that we are keeping (only WBH and only tumor)
keep <- metadata[metadata$Source == "WBH",]
keep <- keep[keep$Type == "Tumor",]
keep <- keep[keep$Keep == "Yes",]$Sample_ID

# Filter the images based on the vector of samples to keep
img <- img[names(img) %in% keep]

# Read in the masks
masks <- loadImages("Data/combined_steinbock/masks_deepcell/", as.is = TRUE)

# Adjust sample IDs for images/masks to match spe
names(masks) <- sub("Slide1-6_CRCTMA1_CRCTMA2_", "ROI_", names(masks))
names(masks) <- sub("Slide2-5_CRCTMA3_CRCTMA4_", "ROI_", names(masks))

# Filter the masks based on the vector of samples to keep
masks <- masks[names(masks) %in% keep]

# Ensure that channel names for images match the spe object
channelNames(img) <- rownames(spe)

# Check that all images/masks are there and sample_ids are the same as spe
length(names(img))
names(img)
length(names(masks))
names(masks)

# Set image names across images/masks and spe
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img))

# Generate vectors to append metadata to images
WBH <- metadata[metadata$Source == "WBH",]
EorL <- WBH[WBH$Type == "Tumor",]
EorL <- EorL[EorL$Keep == "Yes",]$Indication
length(EorL)

WBH <- metadata[metadata$Source == "WBH",]
Pat <- WBH[WBH$Type == "Tumor",]
Pat <- Pat[Pat$Keep == "Yes",]$Patient_ID
length(Pat)

# Assign the metadata to the images and masks as well
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img),
                                        Indication = EorL,
                                        Patient = Pat)

# Set a directory for pre-processed images
imageDir <- "Data/Processed_Images"
if (!dir.exists(imageDir)) {
  dir.create(imageDir)
}

# Save Rds objects of images and masks
saveRDS(img, file = "Data/Processed_Images/WBH_Raw_Images.Rds")
saveRDS(masks, file = "Data/Processed_Images/WBH_Masks.Rds")

################################################################################

# Transform the raw pixel counts with arcsine transformation (cofactor of 1)
spe <- spe[,spe$sample_id %in% keep]
assay(spe, "exprs") <- asinh(counts(spe)/1)

################################################################################

# Run the UMAP embedding on the expression values 
set.seed(030161999)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs")

################################################################################

# Run PCA on the data to prepare for Harmony batch correction
set.seed(03061999)
spe <- runPCA(spe,
              subset_row = rownames(spe)[rowData(spe)$use_channel],
              exprs_values = "exprs",
              ncomponents = 50,
              BSPARAM = ExactParam())

# Run Harmony correction
set.seed(03061999)
out <- RunHarmony(spe, group.by.vars = c("sample_id"))
reducedDim(spe, "harmony") <- reducedDim(out, "HARMONY")

# Run UMAP on the Harmony batch corrected files
set.seed(03061999)
spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_Harmony")

################################################################################

# Save spatialExperiment object
saveRDS(spe, file = "Data/Processed_Spe_Object/DimRedux_Spe.Rds")
