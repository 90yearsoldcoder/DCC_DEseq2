library(DESeq2)


# Read in count data
countData <- read.delim("your_file.tsv", header = TRUE, row.names = 1)


# Read the colData from a TSV file
colData <- read.delim("colData.tsv", stringsAsFactors = FALSE)  # adjust file name and options as needed

# Set the row names of colData to match the sample names
rownames(colData) <- colData$sampleID

# Remove the sampleID column if not needed
colData$sampleID <- NULL


dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ Group)
# Run DESeq2 analysis
dds <- DESeq(dds)
# Replace 'Group1' and 'Group2' with the actual group names you want to compare
result <- results(dds, contrast = c("Group", "Group1", "Group2"))

# Filter for statistically significant differentially expressed genes
DEGs <- subset(result, padj < 0.05 & abs(log2FoldChange) > 1)


