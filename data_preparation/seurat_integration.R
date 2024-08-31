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

in_dir <- "input_data\\GSE129308\\rds_data"
out_dir <- "input_data\\GSE129308\\seurat"

dir.create(out_dir)

###############################################################################
#'Iterate through files in @in_dir and assign metadata
###############################################################################

in_files <- list.files(in_dir, full.names = T)[2:11]
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
  
  #extract experiment name
  #check https://hypebright.nl/index.php/en/2020/05/25/ultimate-cheatsheet-for-regex-in-r-2/
  #on how to work with regex
  #remove all characters before / (inclusive)
  experiment_name <- gsub(".*/","", filename)
  #remove all characters after _ (inclusive)
  experiment_name <- gsub("_.*", "", experiment_name)
  
  #combine metadata into a data.frame
  metadata_df <- data.frame(row.names = metadata_rownames,
                            condition = condition,
                            experiment = experiment_name)
  
  #'add the seurat file to @seurat_list
  seurat_object <- CreateSeuratObject(counts = experiment, meta.data = metadata_df)
  seurat_object$orig.ident <- seurat_object$experiment
  seurat_object$experiment <- NULL
  seurat_list[experiment_name] <- seurat_object

  # #save metadata as csv
  # write.csv(metadata_df, paste(out_dir, paste(experiment_name, "_metadata.csv", sep = ""), sep = "\\"))
}

###############################################################################
#'Perform integration
###############################################################################

#merge objects together
merged_seurat <- merge(seurat_list[[1]], seurat_list[-1])

#preprocessing
# run standard anlaysis workflow
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)

#integrate layers
merged_seurat <- IntegrateLayers(object = merged_seurat, method = CCAIntegration, 
                                 orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)

#save it because it takes A WHILE TO RUN
saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\integrated_data.rds")

#join layers after integration
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\joined_data.rds")
