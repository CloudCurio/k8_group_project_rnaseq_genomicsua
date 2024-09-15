###############################################################################
#'This function subsets a Seurat object by celltype and creates a pseudobulk
#'SingleCellExperiment object as an output
#'@param seurat_obj - a Seurat object containing. Must have a metadata column 
#'"predicted.subclass"
#'@param celltype - character variable containing a cell type that exists in the
#'"predicted.subclass" metadata column
###############################################################################
# create_pseudobulk <- function(seurat_obj, celltype){
#   library(Seurat)
#   library(SingleCellExperiment)
#   
#   #parameter check
#   
#   #subset data
#   seurat <- seurat_obj[, seurat_obj@meta.data$predicted.subclass == celltype]
#   
#   
#   
#   return(sce)
# }

###############################################################################
#'@description
#'This function conducts full DE analysis, creates output plots and tables,
#'and returns a vector of two integer values (significantly up/downregulated
#'genes)
#'@param counts_ls list of dgCmatrices containing counts information for each
#'celltype
#'@param metadata_ls list of data.frames containing metadata for each celltype
#'@param celltype celltype to perform DE analysis on
#'@param out_dir path to the output directory, where results should be saved.
#'The function will create a folder for `celltype`, containing subfolders
#'`QC_plots`, `DE_plots`, and `DE_result_tables`
#'
#'@returns vector of form c(*number of genes upregulated by >50%*, 
#'                          *number of genes downregulated by <30%*)
###############################################################################

DESeq2_workflow <- function(counts_ls, metadata_ls, celltype, out_dir){
  #TODO: add param check if time allows
  
  #check compliance with Windows folder naming rules
  if (celltype == "L2/3 IT"){
    file_celltype <- "L2_3 IT"
  } else {
    file_celltype <- celltype
  }
  print(paste("Starting analysis on cell type", celltype))
  
  #set up output directories
  dir.create(paste(out_dir, file_celltype, sep = "\\"))
  dir.create(paste(out_dir, file_celltype, "QC_plots", sep = "\\"))
  dir.create(paste(out_dir, file_celltype, "DE_plots", sep = "\\"))
  dir.create(paste(out_dir, file_celltype, "DE_result_tables", sep = "\\"))
  # dir.create(paste(out_dir, file_celltype, "", sep = "\\"))
  # dir.create(paste(out_dir, file_celltype, "", sep = "\\"))
  
  #'get the indices of cells of `celltype`
  idx <- which(names(counts_ls) == celltype)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  # Check matching of matrix columns and metadata rows
  all(colnames(cluster_counts) == rownames(cluster_metadata))
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ condition)
  
  
  #sample level QC
  ##PCA
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  
  # Plot PCA
  p1 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "condition")
  p2 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
  p <- p1 + p2
  ggsave(paste(out_dir, file_celltype, "QC_plots", "PCA_cond_cell_count.png", sep = "\\"), p)
  
  ##hierarchical clustering
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Plot heatmap
  p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop=F])
  ggsave(paste(out_dir, file_celltype, "QC_plots", "hclust_cond.png", sep = "\\"), p)
  
  ###############################################################################
  #'Perform DE analysis
  ###############################################################################
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  # Plot dispersion estimates
  png(paste(out_dir, file_celltype, "DE_plots", "dispersion_estimates.png", sep = "\\"), width = 800, height = 800)
  plotDispEsts(dds)
  dev.off()
  
  # Check the coefficients for the comparison
  resultsNames(dds)
  
  # Generate results object
  res <- results(dds, 
                 name = "condition_Control_vs_AD",
                 alpha = 0.05)
  
  # Shrink the log2 fold changes to be more appropriate using the apeglm method - 
  # should cite [paper]() when using this method
  res <- lfcShrink(dds, 
                   coef = "condition_Control_vs_AD",
                   res=res,
                   type = "apeglm")
  
  ###############################################################################
  #'prepare results and visualize them
  ###############################################################################
  
  ##table of results for all genes
  #Turn the DESeq2 results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(padj)
  
  #Check results output
  res_tbl
  
  #Write all results to file
  write.csv(res_tbl,
            paste(out_dir, file_celltype, "DE_result_tables", "AD_vs_Control_all_genes.csv", sep = "\\"),
            quote = FALSE, 
            row.names = FALSE)
  
  ##table of results for significant genes
  #Set thresholds
  padj_cutoff <- 0.05
  
  #Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  #Check significant genes output
  sig_res
  
  #Write significant results to file
  write.csv(res_tbl,
            paste(out_dir, file_celltype, "DE_result_tables", "AD_vs_Control_signif_genes.csv", sep = "\\"),
            quote = FALSE, 
            row.names = FALSE)
  
  ##check for genes with log2 fold change >0.58 (~50% increase/~30% decrease expr)
  #Set thresholds
  log2fc_cutoff <- 0.58
  
  #Count significantly up/down genes above threshold
  n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>% 
    nrow()
  n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>% 
    nrow()
  
  cat(paste("genes with log2 fold change >0.58 for", celltype, "cells:\n",
            n_sig_up, "(>50% up);\n",
            n_sig_dn, "(<30% down);", sep = " "))
  
  print("CSV tables saved successfully")
  
  #if no significantly affected genes are present - terminate analysis
  if (n_sig_up+n_sig_dn == 0){
    message("No significantly altered genes found, terminating analysis")
    return(c(n_sig_up, n_sig_dn))
  }
  
  ##scatterplot of norm expr of top20 significant genes
  #Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  
  #Extract top 20 DEG from resLFC (make sure to order by padj)
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)
  
  #Extract matching normalized count values from matrix
  top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  
  #Convert wide matrix to long data frame for ggplot2
  top20_sig_df <- data.frame(top20_sig_counts)
  top20_sig_df$gene <- rownames(top20_sig_counts)
  
  top20_sig_df <- melt(setDT(top20_sig_df), 
                       id.vars = c("gene"),
                       variable.name = "cluster_sample_id") %>% 
    data.frame()
  
  ##NOTE: since L2/3 is messing formatting up, it has a separate chunk here -
  ##It's ugly, but it works and we don't have much time to do it properly
  if (celltype == "L2/3 IT"){
    #replace "." in "L2.3" with "/"
    top20_sig_df$cluster_sample_id <- sub("\\.", "/", top20_sig_df$cluster_sample_id)
  }
  
  #Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
  
  #extract metadata and turn cluster_sample_id into a factor
  metadata_for_merge <- as.data.frame(colData(dds))
  metadata_for_merge$cluster_sample_id <- as.factor(metadata_for_merge$cluster_sample_id)
  
  #Replace "." by " " again
  metadata_for_merge$cluster_sample_id <- gsub("\\.", " ", metadata_for_merge$cluster_sample_id)
  
  #Join counts data frame with metadata
  top20_sig_df <- plyr::join(top20_sig_df, metadata_for_merge,
                             by = "cluster_sample_id")
  top20_sig_df
  
  #Generate plot
  p <- ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
    geom_jitter(height = 0, width = 0.15) +
    scale_y_continuous(trans = 'log10') +
    ylab("log10 of normalized expression level") +
    xlab("condition") +
    ggtitle("Top 20 Significant DE Genes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ gene)
  ggsave(paste(out_dir, file_celltype, "DE_plots", "scatter_norm_expr_top20sig.png", sep = "\\"), p)
  print("Scatter plot saved successfully")
  
  ##heatmap of all significant genes
  #Extract normalized counts for significant genes only
  sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
  
  #Set a color-blind friendly palette
  heat_colors <- rev(brewer.pal(11, "PuOr"))
  
  #Run pheatmap using the metadata data frame for the annotation
  p <- pheatmap(sig_counts, 
           color = heat_colors, 
           cluster_rows = TRUE, 
           show_rownames = FALSE,
           annotation = cluster_metadata[, c("condition", "cluster_id")], 
           border_color = NA, 
           fontsize = 10, 
           scale = "row", 
           fontsize_row = 10, 
           height = 20)
  ggsave(paste(out_dir, file_celltype, "DE_plots", "heatmap_all_sig_genes.png", sep = "\\"), p)
  print("Heatmap saved successfully")
  
  #volcano plot
  res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
    mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
  min(log10(res_table_thres$padj))
  max(log10(res_table_thres$padj))
  
  ## Generate plot
  p <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    ggtitle("Volcano plot of stimulated B cells relative to control") +
    xlab("log2 fold change") +
    xlim(-4.5, 12) +
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0, 2)) +
    scale_color_manual(values = c("grey60", "red3")) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.3), hjust = 0.5),
          axis.title = element_text(size = rel(1.15)))
  ggsave(paste(out_dir, file_celltype, "DE_plots", "volcano_AD_vs_Ctrl.png", sep = "\\"), p)
  print("Volcano plot saved successfully")
  
  #return numbers of significantly up/downregulated genes for an overall table
  return(c(n_sig_up, n_sig_dn))
}