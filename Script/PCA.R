# Load the data into a data frame
gene_expression <- read.delim("/restricted/projectnb/ncrna/minty/samb_data4/DCC_DEseq2/Processing/CircRNACount_prepared.tsv", header = TRUE, row.names = 1)

# Transpose the data so that each row corresponds to a sample and each column corresponds to a gene
gene_expression_transposed <- t(gene_expression)

# Check that the data has been transposed correctly
head(gene_expression_transposed)

# Normalize the data
gene_expression_transposed_norm <- scale(gene_expression_transposed)

# Read the colData from a TSV file
colData <- read.delim("/restricted/projectnb/ncrna/minty/samb_data4/DCC_DEseq2/Processing/group_info.csv", stringsAsFactors = FALSE,sep=",")

# Set the row names of colData to match the sample names
rownames(colData) <- colData$sampleID

# Remove the sampleID column if not needed
colData$sampleID <- NULL

# Reorder colData based on colnames of gene_expression_transposed_norm
ordered_colData <- colData[rownames(gene_expression_transposed_norm), , drop = FALSE]

# Perform PCA on the data
pca_result <- prcomp(gene_expression_transposed_norm, scale = TRUE)

# View the summary of the PCA results
summary(pca_result)

# Plot the first two principal components
plot(pca_result$x[,1], pca_result$x[,2], xlab = "PC1", ylab = "PC2")

# Color the samples based on their grouping
groups <- ordered_colData$group
print(groups)
if (length(unique(groups)) > 1) {
  colors <- rainbow(length(unique(groups)))
  legend_labels <- levels(factor(groups, levels = unique(groups)))
  for (i in 1:length(unique(groups))) {
    points(pca_result$x[groups == unique(groups)[i],1], pca_result$x[groups == unique(groups)[i],2], col = colors[i], pch = 20)
  }
  legend("topright", legend = legend_labels, fill = colors)
}

# Save the plot
pdf("PCA_plot.pdf")
plot(pca_result$x[,1], pca_result$x[,2], xlab = "PC1", ylab = "PC2")
if (length(unique(groups)) > 1) {
  for (i in 1:length(unique(groups))) {
    points(pca_result$x[groups == unique(groups)[i],1], pca_result$x[groups == unique(groups)[i],2], col = colors[i], pch = 20)
  }
  legend("topright", legend = legend_labels, fill = colors)
}
dev.off()