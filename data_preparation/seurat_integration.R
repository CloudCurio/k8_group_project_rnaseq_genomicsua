###############################################################################
#'This script sets up Seurat object for every .rds file in the input folder, 
#'and stores them in the output folder.
###############################################################################

library(Matrix)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(here)
library(dplyr)
library(glmGamPoi)
library(clustree)
library(Azimuth)

# Install the Azimuth human cortex reference dataset
InstallData("humancortexref")

###############################################################################
#'Set up input values
###############################################################################

i_am("data_preparation//seurat_integration.R")
setwd(here())

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
    condition <- rep("AD", length(metadata_rownames))
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

#seurat_list$GSM3704357@assays$RNA$counts

###############################################################################
#'Perform integration
###############################################################################

#merge objects together
merged_seurat <- merge(seurat_list[[1]], seurat_list[-1])

merged_seurat

# preprocessing
# run standard anlaysis workflow
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)

#integrate layers
merged_seurat <- IntegrateLayers(object = merged_seurat, method = CCAIntegration, 
                                 orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)

# save it because it takes A WHILE TO RUN
saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\integrated_data.rds")

# join layers after integration
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\joined_data.rds")

###############################################################################
#'Shortcut of 70:97 uncomment this chunk to load the integrated file directly
#'TODO: just cut it into a separate script
###############################################################################

#merged_seurat <- readRDS("input_data\\GSE129308\\seurat\\joined_data.rds")


#'@NOTE: this chunk is obsolete, but there was a problem with it. Figure out
#'why controls switch to NA without @param .default specified, but specifying
#'@param .default breaks the RunUMAP() function
#'  
#'fix labeling (from MAP2 to AD for disease)
#'merged_seurat <- merged_seurat@meta.data %>% mutate(condition = case_when(
#'                         condition == "MAP2" ~ "AD",
#'                         .default = "Control"
#'))

###############################################################################
#'Integration cont'd: make a UMAP of the data, check that all clusters are
#'nicely integrated.
###############################################################################

#Make a UMAP of the merged object
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, reduction = "integrated.cca",
                         reduction.name = "umap.integrated")

DimPlot(merged_seurat, reduction = "umap.integrated", 
        group.by = c("orig.ident", "condition"))

###############################################################################
#' QC
###############################################################################
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

#save cell number for further measurements
orig_count <- ncol(merged_seurat)

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 13)

#check the retention percentage
print(paste("Percentage of cells left after MT-genes filtering: ", round(ncol(merged_seurat)*100/orig_count,3), "%", sep = ""))

#write merged_seurat into a file
saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\cleaned_data.rds")

###############################################################################
#'Perform clustering
###############################################################################

# reading 1: https://satijalab.org/seurat/articles/sctransform_vignette.html
# reading 2: https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html

#load obj
merged_seurat <- readRDS("input_data\\GSE129308\\seurat\\cleaned_data.rds")

#replaced with NormalizeData(). Reason: dataset exceeds max size (500 Mb) (local machine limitations)
#merged_seurat <- SCTransform(merged_seurat, verbose = FALSE)
merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, verbose = FALSE)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30, verbose = FALSE)

merged_seurat <- FindClusters(merged_seurat, verbose = FALSE, resolution = seq(0.2, 2.4, 0.2))

clustree(merged_seurat, prefix = "RNA_snn_res.")
#'conclusion: @param resolution = 1.4 is chosen, since a lot of crossovers 
#'happen at 1.4/1.6 with not much change in cluster size

merged_seurat <- FindClusters(merged_seurat, verbose = FALSE, resolution = 1.4)

saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\clustered_data.rds")

###############################################################################
#'Cell annotation
###############################################################################

# https://azimuth.hubmapconsortium.org/ - automatic label transfer
# https://cran.r-project.org/web/packages/scAnnotate/index.html - alternative

# https://satijalab.org/seurat/articles/pbmc3k_tutorial

merged_seurat <- readRDS("input_data\\GSE129308\\seurat\\clustered_data_1.4.rds")

#load reference dataset
LoadData("humancortexref", "azimuth")

# The RunAzimuth function can take a Seurat object as input
#for understanding the scores: https://github.com/satijalab/azimuth/issues/93
merged_seurat <- RunAzimuth(merged_seurat, reference = "humancortexref")

#plot clusters, condition and other useful parameters (such as metrics)
p1 <- DimPlot(merged_seurat, group.by = "predicted.subclass", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(merged_seurat, group.by = "condition")
p3 <- DimPlot(merged_seurat, group.by = "seurat_clusters", label = TRUE, label.size = 3) + NoLegend()
p4 <- FeaturePlot(merged_seurat, features = "mapping.score")
p5 <- FeaturePlot(merged_seurat, features = "predicted.subclass.score")

p1 + p2 #predicted subclass and condition
p1 + p3 #predicted subclass and seurat clustering comparison
p4 + p5 #mapping score and pred subclass score visualization

#add a subclass+condition metadata column
merged_seurat$subclass_cond <- paste(merged_seurat$predicted.subclass, merged_seurat$condition, sep = ".")


#save seurat obj for future use
saveRDS(merged_seurat, "input_data\\GSE129308\\seurat\\annot_data_1.4.rds")

#TODO: Selected neuron subtypes for further work: L2/3 (IT), L5 (IT+ET), L6(IT+CT), Astro

#merged_seurat <- FindAllMarkers(merged_seurat, only.pos = TRUE)

# p_val # do not use
# avg_log2FC # filter by > 0.5
# pct.1 # % of cells expressing this gene in the cluster (see cluster column)
# pct.2 # % of cells eexpressing this gene in the whole dataset
# pct.diff # do it yourself - mutate(pct.diff = pct.1 - pct.2)
# p_val_adj # < 0.05
# cluster # number of cluster
# gene # gene name

# cell annotation: check markers here: https://www.proteinatlas.org/
# do the cell annotation based on the markers
# and check the UMAP (whether it makes sense)
# Neurons_1_Gene1_Gene2

###############################################################################
#'Differential expression
###############################################################################

# Let's meet after you have annotated the cells


merged_seurat@assays$RNA$counts
