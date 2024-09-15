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
library(SeuratData)
library(Matrix)
library(patchwork)
library(here)
library(dplyr)
library(glmGamPoi)
library(clustree)
#library(Azimuth)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(biovizBase)
#InstallData("humancortexref")


###############################################################################
#'initial settings
###############################################################################


setwd("C:\\Users\\liber\\Desktop\\Study\\Genomics UA\\RNA-seq analysis with R\\Group Project")

in_path <- "input_data//GSE129308//"
out_path <- "input_data//GSE129308//rds_data//"

###############################################################################
#'read input files
###############################################################################

data_cells <- readRDS("input_data\\GSE129308\\seurat\\annot_data_1.4.rds")

dir.create(paste(in_path, "sinlgecell_data"))

###############################################################################
#'data
###############################################################################

VlnPlot(data_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#immunecells <- subset(immunecells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 13)


pbmc <- data_cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Add mitochondrial gene percentage to metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 13)


pbmc <- NormalizeData(pbmc)
# Visualize UMAP by condition (AD vs. control) to check for batch effects
# Visualize UMAP by condition (AD vs. control) to check for batch effects
DimPlot(pbmc, reduction = "umap", group.by = "condition")
DimPlot(pbmc, reduction = "umap", group.by = "predicted.cluster")

# Perform differential expression between AD and Control across all cells
markers_AD_vs_Control <- FindMarkers(pbmc, ident.1 = "AD", ident.2 = "Control", group.by = "condition")

# View the top differentially expressed genes
head(markers_AD_vs_Control)

# Perform DE analysis for a specific cell type 
markers_neurons <- FindMarkers(pbmc, ident.1 = "AD", ident.2 = "Control", group.by = "condition", subset.ident = "predicted.subclass")

# View the top differentially expressed genes for neurons
head(markers_neurons)

# Subset L5, L6, and L2/3 cells 
L5IT_cells <- subset(pbmc, subset = predicted.subclass == "L5 IT")
L6IT_cells <- subset(pbmc, subset = predicted.subclass == "L6 IT")
L2_3IT_cells <- subset(pbmc, subset = predicted.subclass == "L2/3 IT")
L5ET_cells <- subset(pbmc, subset = predicted.subclass == "L5 ET")
L6CT_cells <- subset(pbmc, subset = predicted.subclass == "L6 CT")
Astro_cells <- subset(pbmc, subset = predicted.subclass == "Astro")

# Differential expression between AD and control for each subset
markers_L5IT_cells <- FindMarkers(L5IT_cells, ident.1 = "AD", ident.2 = "Control", group.by = "condition")
markers_L6IT_cells <- FindMarkers(L6IT_cells, ident.1 = "AD", ident.2 = "Control", group.by = "condition")
markers_L2_3IT_cells <- FindMarkers(L2_3IT_cells, ident.1 = "AD", ident.2 = "Control", group.by = "condition")
markers_L5ET_cells <- FindMarkers(L5ET_cells, ident.1 = "AD", ident.2 = "Control", group.by = "condition")
markers_L6CT_cells <- FindMarkers(L6CT_cells, ident.1 = "AD", ident.2 = "Control", group.by = "condition")
markers_Astro_cells <- FindMarkers(Astro_cells, ident.1 = "AD", ident.2 = "Control", group.by = "condition")

head(markers_L5IT_cells)
head(markers_L6IT_cells)
head(markers_L2_3IT_cells)
head(markers_L5ET_cells)
head(markers_L6CT_cells)
head(markers_Astro_cells)

# Select top genes from L5
top_genes_L5IT_cells <- rownames(markers_L5IT_cells)[1:10]
top_genes_L6IT_cells <- rownames(markers_L6IT_cells)[1:10]
top_genes_L2_3IT_cells <- rownames(markers_L2_3IT_cells)[1:10]
top_genes_L5ET_cells <- rownames(markers_L5ET_cells)[1:10]
top_genes_L6CT_cells <- rownames(markers_L6CT_cells)[1:10]# Select top 10 genes
top_genes_Astro_cells <- rownames(markers_Astro_cells)[1:10]

###############################################################################
#'data plots
###############################################################################


# Save DE results to CSV files
write.csv(markers_L5IT_cells, "DE_markers_L5IT_cells_AD_vs_Control.csv")
write.csv(markers_L6IT_cells, "DE_markers_L6IT_cells_AD_vs_Control.csv")
write.csv(markers_L2_3IT_cells, "DE_markers_L2_3IT_cells_AD_vs_Control.csv")
write.csv(markers_L5ET_cells, "DE_markers_L5ET_cells_AD_vs_Control.csv")
write.csv(markers_L6CT_cells, "DE_markers_L6CT_cells_AD_vs_Control.csv")
write.csv(markers_Astro_cells, "DE_markers_Astro_cells_AD_vs_Control.csv")

# Create a heatmap for L5 cells
DoHeatmap(top_genes_L5IT_cells, features = top_genes_L5) + scale_fill_viridis()
DoHeatmap(top_genes_L6IT_cells, features = top_genes_L5) + scale_fill_viridis()
DoHeatmap(top_genes_L2_3IT_cells, features = top_genes_L5) + scale_fill_viridis()
DoHeatmap(top_genes_L5ET_cells, features = top_genes_L5) + scale_fill_viridis()
DoHeatmap(top_genes_L6CT_cells, features = top_genes_L5) + scale_fill_viridis()
DoHeatmap(top_genes_Astro_cells, features = top_genes_L5) + scale_fill_viridis()

# Create volcano plots for L5 IT cells
markers_L5IT_cells$gene <- rownames(markers_L5IT_cells)
volcano_L5IT <- ggplot(markers_L5IT_cells, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), size = 1.5) +
  labs(title = "Volcano Plot for L5 IT Cells (AD vs Control)", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot for L5 IT cells
ggsave("Volcano_L5IT_cells.png", plot = volcano_L5IT)

# Create volcano plots for L6 IT cells
markers_L6IT_cells$gene <- rownames(markers_L6IT_cells)
volcano_L6IT <- ggplot(markers_L6IT_cells, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), size = 1.5) +
  labs(title = "Volcano Plot for L6 IT Cells (AD vs Control)", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot for L6 IT cells
ggsave("Volcano_L6IT_cells.png", plot = volcano_L6IT)

# Create volcano plots for L2/3 IT cells
markers_L2_3IT_cells$gene <- rownames(markers_L2_3IT_cells)
volcano_L2_3IT <- ggplot(markers_L2_3IT_cells, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), size = 1.5) +
  labs(title = "Volcano Plot for L2/3 IT Cells (AD vs Control)", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot for L2/3 IT cells
ggsave("Volcano_L2_3IT_cells.png", plot = volcano_L2_3IT)

# Create volcano plots for L5 ET cells
markers_L5ET_cells$gene <- rownames(markers_L5ET_cells)
volcano_L5ET <- ggplot(markers_L5ET_cells, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), size = 1.5) +
  labs(title = "Volcano Plot for L5 ET Cells (AD vs Control)", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot for L5 ET cells
ggsave("Volcano_L5ET_cells.png", plot = volcano_L5ET)

# Create volcano plots for L6 CT cells
markers_L6CT_cells$gene <- rownames(markers_L6CT_cells)
volcano_L6CT <- ggplot(markers_L6CT_cells, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), size = 1.5) +
  labs(title = "Volcano Plot for L6 CT Cells (AD vs Control)", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot for L6 CT cells
ggsave("Volcano_L6CT_cells.png", plot = volcano_L6CT)

# Create volcano plots for astrocytes
markers_Astro_cells$gene <- rownames(markers_Astro_cells)
volcano_astrocytes <- ggplot(markers_Astro_cells, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), size = 1.5) +
  labs(title = "Volcano Plot for Astro Cells (AD vs Control)", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot for astrocytes
ggsave("Volcano_Astro_cells.png", plot = volcano_astrocytes)
# Create heatmap for top 10 DE genes in L5 IT cells


# Create heatmap for top 10 DE genes in L5 IT cells
top_genes_L5IT <- rownames(markers_L5IT_cells)[1:10]
png("Heatmap_L5IT_cells.png", width = 1000, height = 800)
DoHeatmap(L5IT_cells, features = top_genes_L5IT, group.by = "condition")
dev.off()

# Create heatmap for top 10 DE genes in L6 IT cells
top_genes_L6IT <- rownames(markers_L6IT_cells)[1:10]
png("Heatmap_L6IT_cells.png", width = 1000, height = 800)
DoHeatmap(L6IT_cells, features = top_genes_L6IT, group.by = "condition")
dev.off()

# Create heatmap for top 10 DE genes in L2/3 IT cells
top_genes_L2_3IT <- rownames(markers_L2_3IT_cells)[1:10]
png("Heatmap_L2_3IT_cells.png", width = 1000, height = 800)
DoHeatmap(L2_3IT_cells, features = top_genes_L2_3IT, group.by = "condition")
dev.off()

# Create heatmap for top 10 DE genes in L5 ET cells
top_genes_L5ET <- rownames(markers_L5ET_cells)[1:10]
png("Heatmap_L5ET_cells.png", width = 1000, height = 800)
DoHeatmap(L5ET_cells, features = top_genes_L5ET, group.by = "condition")
dev.off()

# Create heatmap for top 10 DE genes in L6 CT cells
top_genes_L6CT <- rownames(markers_L6CT_cells)[1:10]
png("Heatmap_L6CT_cells.png", width = 1000, height = 800)
DoHeatmap(L6CT_cells, features = top_genes_L6ET, group.by = "condition")
dev.off()

# Get the expression data for the top DE genes
top_genes <- rownames(markers_AD_vs_Control)[1:10]  # Adjust the number of genes as needed

# Extract normalized/scaled data from the Seurat object
expr_data <- GetAssayData(pbmc, assay = "RNA", slot = "scale.data")[top_genes, ]

# Optionally, transpose the matrix if needed for clustering by rows (genes) and columns (cells/clusters)
expr_data <- t(expr_data)


# Create heatmap for top 10 DE genes in Astro cells
top_genes_Astro <- rownames(markers_Astro_cells)[1:10]
png("Heatmap_Astro_cells.png", width = 1000, height = 800)
DoHeatmap(Astro_cells, features = top_genes_Astro, group.by = "condition")
dev.off()


###############################################################################
#'profiling
###############################################################################

library(pheatmap)

de_genesL5IT_cells <- rownames(markers_L5IT_cells)[markers_L5IT_cells$p_val_adj < 0.05]
de_genesL6IT_cells <- rownames(markers_L5IT_cells)[markers_L6IT_cells$p_val_adj < 0.05]
de_genesL2_3_IT_cells <- rownames(markers_L2_3IT_cells)[markers_L2_3IT_cells$p_val_adj < 0.05]
de_genesL5ET_cells <- rownames(markers_L5ET_cells)[markers_L5ET_cells$p_val_adj < 0.05]
de_genesL6CT_cells <- rownames(markers_L6CT_cells)[markers_L6CT_cells$p_val_adj < 0.05]
de_genesAstro_cells <- rownames(markers_Astro_cells)[markers_Astro_cells$p_val_adj < 0.05]

length(de_genesL5IT_cells)
length(de_genesL6IT_cells)
length(de_genesL2_3_IT_cells)
length(de_genesL5ET_cells)
length(de_genesL6CT_cells)
length(de_genesAstro_cells)

# List available databases and select those of interest
dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2016")

# Run enrichment analysis for each cell type

# L5 IT cells
enriched_L5IT <- enrichr(de_genesL5IT_cells, dbs)
print(enriched_L5IT)

# L6 IT cells
enriched_L6IT <- enrichr(de_genesL6IT_cells, dbs)
print(enriched_L6IT)

# L2/3 IT cells
enriched_L2_3IT <- enrichr(de_genesL2_3_IT_cells, dbs)
print(enriched_L2_3IT)

# L5 ET cells
enriched_L5ET <- enrichr(de_genesL5ET_cells, dbs)
print(enriched_L5ET)

# L6 CT cells
enriched_L6CT <- enrichr(de_genesL6CT_cells, dbs)
print(enriched_L6CT)

# Astrocytes
enriched_Astro <- enrichr(de_genesAstro_cells, dbs)
print(enriched_Astro)


plot_and_save_enrichr_results <- function(enrichr_results, db_name, cell_type_name, pval_threshold = 0.05, top_n = 10) {
  
  # Filter and select the top enriched terms based on adjusted p-value
  results_to_plot <- enrichr_results[[db_name]] %>%
    filter(Adjusted.P.value < pval_threshold) %>%
    head(top_n)
  
  # Create the bar plot
  p <- ggplot(results_to_plot, aes(x = reorder(Term, -Adjusted.P.value), 
                                   y = -log10(Adjusted.P.value), 
                                   fill = Combined.Score)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top Enriched Terms:", db_name, "in", cell_type_name),
         x = "Enrichment Term",
         y = "-log10(Adjusted P-value)") +
    theme_minimal()
  
  # Save the plot as PNG
  ggsave(filename = paste0("Enriched_", cell_type_name, "_", db_name, ".png"), plot = p, width = 10, height = 8)
}

# L5 IT cells
#enriched_L5IT <- enrichr(de_genesL5IT_cells, dbs)
plot_and_save_enrichr_results(enriched_L5IT, "GO_Biological_Process_2021", "L5IT")
plot_and_save_enrichr_results(enriched_L5IT, "KEGG_2021_Human", "L5IT")
plot_and_save_enrichr_results(enriched_L5IT, "Reactome_2016", "L5IT")

# L6 IT cells
enriched_L6IT <- enrichr(de_genesL6IT_cells, dbs)
plot_and_save_enrichr_results(enriched_L6IT, "GO_Biological_Process_2021", "L6IT")
plot_and_save_enrichr_results(enriched_L6IT, "KEGG_2021_Human", "L6IT")
plot_and_save_enrichr_results(enriched_L6IT, "Reactome_2016", "L6IT")

# L2/3 IT cells
enriched_L2_3IT <- enrichr(de_genesL2_3_IT_cells, dbs)
plot_and_save_enrichr_results(enriched_L2_3IT, "GO_Biological_Process_2021", "L2_3_IT")
plot_and_save_enrichr_results(enriched_L2_3IT, "KEGG_2021_Human", "L2_3_IT")
plot_and_save_enrichr_results(enriched_L2_3IT, "Reactome_2016", "L2_3_IT")

# L5 ET cells
enriched_L5ET <- enrichr(de_genesL5ET_cells, dbs)
plot_and_save_enrichr_results(enriched_L5ET, "GO_Biological_Process_2021", "L5ET")
plot_and_save_enrichr_results(enriched_L5ET, "KEGG_2021_Human", "L5ET")
plot_and_save_enrichr_results(enriched_L5ET, "Reactome_2016", "L5ET")

# L6 CT cells
enriched_L6CT <- enrichr(de_genesL6CT_cells, dbs)
plot_and_save_enrichr_results(enriched_L6CT, "GO_Biological_Process_2021", "L6CT")
plot_and_save_enrichr_results(enriched_L6CT, "KEGG_2021_Human", "L6CT")
plot_and_save_enrichr_results(enriched_L6CT, "Reactome_2016", "L6CT")

# Astrocytes
enriched_Astro <- enrichr(de_genesAstro_cells, dbs)
plot_and_save_enrichr_results(enriched_Astro, "GO_Biological_Process_2021", "Astro")
plot_and_save_enrichr_results(enriched_Astro, "KEGG_2021_Human", "Astro")
plot_and_save_enrichr_results(enriched_Astro, "Reactome_2016", "Astro")


# Save enrichment results for each cell type
write.csv(enriched_L5IT$GO_Biological_Process_2021, "enriched_L5IT_GO_BP.csv")
write.csv(enriched_L6IT$GO_Biological_Process_2021, "enriched_L6IT_GO_BP.csv")
write.csv(enriched_L2_3IT$GO_Biological_Process_2021, "enriched_L2_3IT_GO_BP.csv")
write.csv(enriched_L5ET$GO_Biological_Process_2021, "enriched_L5ET_GO_BP.csv")
write.csv(enriched_L6CT$GO_Biological_Process_2021, "enriched_L6CT_GO_BP.csv")
write.csv(enriched_Astro$GO_Biological_Process_2021, "enriched_Astro_GO_BP.csv")




