# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Read the dataset
phylofish_unga_data <- read.csv("phylofish_unga_dataset-2025-05-27.csv", row.names = 1)

# Convert to matrix for heatmap
data_matrix <- as.matrix(phylofish_unga_data)

# Log2 transform the data (add 1 to avoid log(0))
data_log <- log2(data_matrix + 1)

# Create and save the heatmap
gt <- pheatmap(data_log,
         main = "Gene Expression Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",  # Scale by rows (genes)
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         #clustering_distance_cols = "euclidean",
         #clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10,
         fontsize_row = 8,
         fontsize_col = 10,
         angle_col = 45,
         border_color = NA,
         cellwidth = 30,
         cellheight = 15)$gtable

ggsave("gene_expression_heatmap.pdf", plot=gt)

