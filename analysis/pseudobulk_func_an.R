###############################################################################
#'This script performs enrichment analysis for each cell type individually
#'99% of this script was copied from O.Petrenko's file:
#'https://github.com/GenomicsUA/2024-daad-rnaseq-course/blob/main/02_classes/10_functional_analysis/99_overrepresentation_filled.Rmd
###############################################################################

#load libraries (install missing ones)
pacman::p_load("here", "tidyverse", "enrichR", "gprofiler2")

###############################################################################
#'Set up initial values
###############################################################################

#set wd
i_am("analysis//pseudobulk_func_an.R")
setwd(here())
setwd("outputs//GSE129308//pseudobulk")

source("..//..//..//functions//fun.R")

#designate cell types of interest
#cell_types <- c("L2/3 IT", "L5 IT", "L5 ET", "L6 IT", "L6 CT", "Astro")

#load a table with numbers of up/downregulated genes per cell type
de_gene_numbers <- read.csv("up_down_reg_genes.csv", row.names = 1)

#iterate through cell types
for (celltype in rownames(de_gene_numbers)){
  print(paste("Starting analysis on", celltype))
  #check if strong log2 fold change is observed in any genes (>50% up/ <30% down)
  if (sum(de_gene_numbers[celltype,]) == 0){
    de_signif_flag <- FALSE
    warning("No genes with log2 fold change of >0.58 present") 
  } else {
    de_signif_flag <- TRUE
  }
  
  #load a table with significantly up/downregulated genes
  
  ##NOTE: L2/3 IT has to be loaded a bit differently due to Windows folder 
  ##naming rules
  if (celltype == "L2/3 IT"){
    file_celltype <- "L2_3 IT"
    de_results <- read.csv(paste(file_celltype, "DE_result_tables//AD_vs_Control_signif_genes.csv", 
                                 sep = "//"))
  } else {
    file_celltype <- celltype
    de_results <- read.csv(paste(file_celltype, "DE_result_tables//AD_vs_Control_signif_genes.csv", 
                                 sep = "//"))
  }
  
  #remove genes with padj>=0.05
  de_results <- de_results %>% filter(padj < 0.05)
  
  ###############################################################################
  #'Perform functional analysis
  ###############################################################################
  #select databases to use
  dbs <- listEnrichrDbs()
  dbs <- c("Human_Gene_Atlas", "GO_Molecular_Function_2023", 
           "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", 
           "Aging_Perturbations_from_GEO_down", "Aging_Perturbations_from_GEO_up", 
           "RNAseq_Automatic_GEO_Signatures_Human_Down", 
           "RNAseq_Automatic_GEO_Signatures_Human_Up",
           "Reactome_2016", "CellMarker_2024")
  
  #run enrichment analysis on all selected db's
  enriched <- enrichr(rownames(de_results[1:300,] %>% filter(log2FoldChange > 0) %>%
                                 arrange(desc(log2FoldChange))), dbs)
  
  #iterate through all databases and plot results
  for (db in names(enriched)){
    pdf(paste(file_celltype, "enrichment_results.pdf", sep = "//"))
    p <- plot_enrichr_results(enriched[[db]])
    print(p)
    dev.off()
  }
}