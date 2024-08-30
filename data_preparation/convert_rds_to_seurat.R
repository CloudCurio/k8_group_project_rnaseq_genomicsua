###############################################################################
#'This script sets up Seurat object for every .rds file in the input folder, 
#'and stores them in the output folder.
###############################################################################

library(Matrix)
library(Seurat)
library(patchwork)

###############################################################################
#'Set up input values
###############################################################################

setwd("C:\\Users\\liber\\Desktop\\Study\\Genomics UA\\RNA-seq analysis with R\\Group Project")

in_dir <- "..\\..\\hard_data\\GSE129308\\in_data"
out_dir <- "..\\..\\hard_data\\GSE129308\\seurat_data"

###############################################################################
#'Iterate through files in @in_dir and assign metadata
###############################################################################

in_files <- list.files(in_dir, full.names = T)
seurat_list <- list()

for (filename in in_files){
  #read the file
  experiment <- as.matrix(readRDS(filename))
  
  #prepare metadata
  metadata_rownames <- colnames(experiment)
  ##check if control or sample for condition selection
  if (grepl("Control", filename)){
    condition <- rep("Control", length(metadata_rownames))
  } else {
    condition <- rep("MAP2", length(metadata_rownames))
  }
  
  #combine metadata into a data.frame
  metadata_df <- data.frame(row.names = metadata_rownames,
                            condition = condition)
  
  #make a nice name for the list element
  list_element_name <- gsub(".*/","", filename)
  list_element_name <- sub(".h5.rds", "", list_element_name)
  
  #'add the seurat file to @seurat_list
  seurat_object <- CreateSeuratObject(counts = experiment, meta.data = metadata_df)
  seurat_list[list_element_name] <- seurat_object
}

###############################################################################
#'Perform analysis for separate objects
###############################################################################

seurat_list <- lapply(seurat_list, NormalizeData)
seurat_list <- lapply(seurat_list, FindVariableFeatures)
seurat_list <- lapply(seurat_list, ScaleData)
seurat_list <- lapply(seurat_list, RunPCA)

seurat_list <- lapply(seurat_list, FindNeighbors, dims = 1:30, reduction = "pca")
seurat_list <- lapply(seurat_list, FindClusters, resolution = 2, 
                                                cluster.name = "unintegrated_clusters")

seurat_list <- lapply(seurat_list, RunUMAP, dims = 1:30, reduction = "pca", 
                                           reduction.name = "umap.unintegrated")
lapply(seurat_list, DimPlot, reduction = "umap.unintegrated", 
                            group.by = c("stim", "seurat_clusters"))
###############################################################################
#'Perform integration
###############################################################################

