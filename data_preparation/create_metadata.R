###############################################################################
#'This script makes a metadata df for every .rds file in the input folder, 
#'and stores them in the output folder.
###############################################################################

library(Matrix)
library(patchwork)

###############################################################################
#'Set up input values
###############################################################################

setwd("C:\\Users\\liber\\Desktop\\Study\\Genomics UA\\RNA-seq analysis with R\\Group Project")

in_dir <- "input_data\\GSE129308\\rds_data"
out_dir <- "input_data\\GSE129308\\metadata"

###############################################################################
#'Iterate through files in @in_dir and assign metadata
###############################################################################

in_files <- list.files(in_dir, full.names = T)[2:11]

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
  
  #save metadata as csv
  write.csv(metadata_df, paste(out_dir, paste(experiment_name, "_metadata.csv", sep = ""), sep = "\\"))
}
