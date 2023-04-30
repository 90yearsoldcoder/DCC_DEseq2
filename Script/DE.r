library(DESeq2)
library(tidyverse)
library(dplyr)


# Read in count data
countData <- read.delim("/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Processing/CircRNACount_prepared.tsv", header = TRUE)
countData <- countData[!duplicated(countData$gene), ]
row.names(countData) <- countData$gene
countData$gene <- NULL

# Read the colData from a TSV file
colData <- read.delim("/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Processing/group_info.csv", stringsAsFactors = FALSE,sep=",")  # adjust file name and options as needed

# Set the row names of colData to match the sample names
rownames(colData) <- colData$sampleID

# Remove the sampleID column if not needed
colData$sampleID <- NULL

# Reorder colData based on colnames of countData
ordered_colData <- colData[colnames(countData), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData, ordered_colData, design = ~ group)
# Run DESeq2 analysis
dds <- DESeq(dds)




generate_DEG_table <- function(dds, contrast, p_value_flag = FALSE, res_path, sig_path) {
  #res_path stores all gene DE result
  #sig_path stores the genes passing the filter
  
  # Generate results table for specified contrast
  result <- results(dds, contrast = contrast)
  
  # Filter for either p-value or adjusted p-value as specified
  if (p_value_flag) {
    DEGs <- subset(result, pvalue < 0.05 & abs(log2FoldChange) > 1)
  } else {
    DEGs <- subset(result, padj < 0.05 & abs(log2FoldChange) > 1)
  }
  
  # Convert results table to a tibble, add gene names as a column, and sort by p-value
  res_table <- as_tibble(result) %>%
    mutate(gene = rownames(result)) %>%
    relocate(gene) %>%
    arrange(pvalue)
  
  DEGs_table <- as_tibble(DEGs) %>%
    mutate(gene = rownames(DEGs)) %>%
    relocate(gene) %>%
    arrange(pvalue)
  
  # Save the table to output_path
  write_csv(res_table, res_path)
  write_csv(DEGs_table, sig_path)
  
  # Return the tibble
  return(res_table)
}

plot_volcano_pvalue <- function(res_table, title, fig_path){
  data <- mutate(
    res_table,
    `-log10(p_value)`=-log10(pvalue),
    `p_value<0.05`=pvalue<0.05
  ) 
  plot <- ggplot(data, aes(x=log2FoldChange,y=`-log10(p_value)`,color=`p_value<0.05`)) +
    ggtitle(title) +
    geom_point()
  
  #save the picture
  print(paste0(fig_path, "/", title, ".png"))
  ggsave(paste0(fig_path, "/", title, ".png"), plot)
  
}


# control1 vs otao1
table <- generate_DEG_table(dds, contrast =  c("group", "Control1", "oTau1"), p_value_flag = TRUE, 
                   res_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/control1_vs_oTau1.csv", 
                  sig_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/control1_vs_oTau1_sig.csv"            
                                      )
plot_volcano_pvalue(table, "control1_vs_oTau1", "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result")

# control2 vs otao2
table <- generate_DEG_table(dds, contrast =  c("group", "Control2", "oTau2"), p_value_flag = TRUE, 
                   res_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/control2_vs_oTau2.csv", 
                   sig_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/control2_vs_oTau2_sig.csv"            
)
plot_volcano_pvalue(table, "control2_vs_oTau2", "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result")

# control1 vs control2
table <- generate_DEG_table(dds, contrast =  c("group", "Control1", "Control2"), p_value_flag = TRUE, 
                   res_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/control1_vs_control2.csv", 
                   sig_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/control1_vs_control2_sig.csv"            
)
plot_volcano_pvalue(table, "control1_vs_control2", "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result")

# otao1 vs otao2
table <- generate_DEG_table(dds, contrast =  c("group", "oTau1", "oTau2"), p_value_flag = TRUE, 
                   res_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/oTau1_vs_oTau2.csv", 
                   sig_path = "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result/oTau1_vs_oTau2_sig.csv"            
)
plot_volcano_pvalue(table, "oTau1_vs_oTau2", "/restricted/projectnb/ncrna/minty/samb_data4_removeOL/DCC_DEseq2/Result")










