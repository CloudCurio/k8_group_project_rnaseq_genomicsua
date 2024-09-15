###############################################################################
#'This script takes annot_data_*resolution*.rds files, turns data into 
#'pseudobulk, and performs DEG analysis and Functional Gene Set Analysis on it
###############################################################################

library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(here)
library(dplyr)
library(glmGamPoi)
library(clustree)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

source("functions//fun.R")

# Install the Azimuth human cortex reference dataset
InstallData("humancortexref")

###############################################################################
#'Set up input values
###############################################################################

i_am("data_preparation//seurat_integration.R")
setwd(here())

in_dir <- "input_data\\GSE129308\\rds_data"
out_dir <- "outputs\\GSE129308\\pseudobulk"

dir.create(out_dir)

LoadData("humancortexref", "azimuth")
cell_types <- c("L2/3 IT", "L5 IT", "L5 ET", "L6 IT", "L6 CT", "Astro")

###############################################################################
#'Load data, filter it for desired cell subtypes
###############################################################################

#load seurat object with annotated cells
seurat <- readRDS("input_data\\GSE129308\\seurat\\annot_data_1.4.rds")

#add metadata column for sample name+condition
seurat$sample_id <- paste(seurat$orig.ident, seurat$condition, sep = ".")

#subset data to only include cell subtypes of interest
seurat <- seurat[, seurat@meta.data$predicted.subclass %in% cell_types]

###############################################################################
#'Load data, filter it for desired cell subtypes
###############################################################################

#create pseudobulks for each cell type (AD and Ctrl combined)
#pseudobulk_list <- list()

#extract count data and metadata to create SingleCellExperiment object
counts <- seurat[["RNA"]]$counts
metadata <- seurat@meta.data

#set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat$predicted.subclass)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

#Extract unique names of clusters (celltype+cond info)
cluster_names <- levels(colData(sce)$cluster_id)

#Extract sample id's
sample_names <- levels(as.factor(colData(sce)$sample_id))

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum")

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)


# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(condition, orig.ident, sample_id, )

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
# Rename rows
rownames(metadata) <- metadata$sample_id

# Number of cells per sample and cluster
t <- table(colData(sce)$sample_id,
           colData(sce)$cluster_id)


# Creating metadata list

## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

#'Important variables going forward:
#'`sce` - SingleCellExperiment object of our data
#'`counts_ls` - counts data by cell type
#'`metadata_ls` - metadata by cell type 
#'
#'Nice metrics to understand our data:
#'*t* - number of cells of each cell type in different samples
#'*aggr_counts* - aggregated counts in gene/sample format
#'*groups* - groups for aggregation per cell (celltype, sample)


###############################################################################
#'Prepare a DESeq2 object and perform PCA analysis and hclust
###############################################################################

#create a storage df for numbers of significantly up/downregulated genes for
#each cell type
up_down_reg_df <- data.frame(matrix(ncol = 2, nrow = 0))

for (celltype in cell_types){
  celltype_gene_counts <- DESeq2_workflow(counts_ls, metadata_ls, celltype, out_dir)
  up_down_reg_df <- rbind(up_down_reg_df, celltype_gene_counts)
  rownames(up_down_reg_df)[nrow(up_down_reg_df)] <- celltype
}

colnames(up_down_reg_df) <- c("upregulation >50%", "downregulation <30%")
write.csv(up_down_reg_df, paste(out_dir, "up_down_reg_genes.csv", sep = "\\"))
