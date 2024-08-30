###############################################################################
#'This script converts h5 files into Seurat-suitable format
#'
#'Tutorial on how to set up the required packages: 
#'https://github.com/cellgeni/sceasy
#'
#'conda env (called GenomicsUA) start:
#'conda activate GenomicsUA
#'
#'conda env location:
#'C:\Users\liber\Desktop\Study\Genomics UA\RNA-seq analysis with R\Group Project
###############################################################################

###############################################################################
#'load required libraries
###############################################################################

library(Seurat)
library(SeuratDisk)

###############################################################################
#'initial settings
###############################################################################

setwd("C:\\Users\\liber\\Desktop\\Study\\Genomics UA\\RNA-seq analysis with R\\Group Project")

in_path <- "input_data//GSE129308//"
out_path <- "input_data//GSE129308//rds_data//"

###############################################################################
#'read input files
###############################################################################

file_list <- list.files(paste(in_path, "raw_data//", sep = ""), all.files = T, full.names = F)
file_list <- file_list[-c(1:2)] #drop /. and /..

dir.create(paste(in_path, "rds_data"))

###############################################################################
#'load h5 files as dgcMatrix objects and save them as .rds in the
#'rds_data directory
#'TODO: figure out why saving files to rds_data doesn't work
###############################################################################

for (filename in file_list){
  foo <- Read10X_h5(paste(in_path, "raw_data//", filename, sep = ""))
  saveRDS(foo, file = paste(out_path, filename, ".rds", sep = ""))
}
