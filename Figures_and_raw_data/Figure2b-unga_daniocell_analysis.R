# This is an analysis pipeline for the Daniocell (https://daniocell.nichd.nih.gov) 
# scRNAseq dataset.

library(Seurat)
library(dplyr)
library(ggplot2)

# Load the Daniocell Seurat object and the latest annotations
path.to.seurat.object <- "daniocell_data/Daniocell2023_SeuratV4.rds"
path.to.annotations <- "daniocell_data/cluster_annotations.csv"

daniocell <- readRDS(path.to.seurat.object)
daniocell.annot <- read.csv(path.to.annotations)

FeaturePlot(daniocell, features = "unga") + ggtitle("Expression of unga across all cells")

# Violin plot of unga expression by tissue
VlnPlot(daniocell, features = "unga", group.by = "subset.full", pt.size = 0) + 
  ggtitle("unga expression by tissue")

# Violin plot of unga expression by developmental stage
VlnPlot(daniocell, features = "unga", group.by = "stage.group", pt.size = 0) + 
  ggtitle("unga expression by developmental stage")

# DotPlot by tissue
DotPlot(daniocell, features = c("unga", "ungb"), group.by = "subset.full") + 
  RotatedAxis() + ggtitle("unga and ungb expression by tissue")

ggsave("20250521-unga_ungb_danioCell_tissue.pdf", units = "cm", height = 15, width = 17.5)

# DotPlot by developmental stage
DotPlot(daniocell, features = c("unga", "ungb"), group.by = "stage.group") + 
  RotatedAxis() + ggtitle("unga and ungb expression by developmental stage")

ggsave("20250521-unga_ungb_danioCell_devStage.pdf", units = "cm", height = 10, width = 12.5)

expr_mat <- GetAssayData(daniocell, assay = "RNA", slot = "data")

# For example, select genes expressed in at least 10% of cells
expr_frac <- rowSums(expr_mat > 0) / ncol(expr_mat)
genes_to_use <- names(expr_frac[expr_frac > 0.1])

expr_mat_sub <- expr_mat[genes_to_use, ]

chunk_size <- 1000
cor_vals <- numeric(length(genes_to_use))

for (i in seq(1, length(genes_to_use), by = chunk_size)) {
  idx <- i:min(i + chunk_size - 1, length(genes_to_use))
  chunk <- expr_mat_sub[idx, ]
  
  cor_chunk <- sapply(1:nrow(chunk), function(j) {
    cor(unga_expr, chunk[j, ], method = "spearman")
  })
  
  cor_vals[idx] <- cor_chunk
  gc()
}

names(cor_vals) <- genes_to_use

# Select genes with correlation above a threshold (e.g., 0.3)
coexpressed_genes <- names(cor_vals[cor_vals > 0.3])

# Optionally remove unga itself
coexpressed_genes <- setdiff(coexpressed_genes, "unga")

library(biomaRt)

# Connect to Ensembl for zebrafish
mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Map gene symbols to Ensembl IDs
gene_map <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                  filters = "external_gene_name",
                  values = coexpressed_genes,
                  mart = mart)

# Use Ensembl IDs for GO enrichment
gene_list_ensembl <- gene_map$ensembl_gene_id

library(clusterProfiler)
library(org.Dr.eg.db)  # Zebrafish annotation package from Bioconductor

# Perform GO enrichment (Biological Process ontology)
go_results <- enrichGO(gene = gene_list_ensembl,
                       OrgDb = org.Dr.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",          # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)     # Converts IDs back to gene symbols

# View top enriched GO terms
head(go_results)

# Visualize results
dotplot(go_results, showCategory = 20) + ggtitle("GO Enrichment of Genes Coexpressed with unga")

ggsave("20250521-unga_danioCell_GO_analysis.pdf", units = "cm", height = 15, width = 20)

