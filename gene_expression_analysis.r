# Gene Expression Analysis in R
# Bioinformatics Pipeline: Differential Gene Expression (RNA-seq) + GO Enrichment

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load count data and sample metadata
count_data <- read.csv("counts_matrix.csv", row.names = 1)
sample_info <- read.csv("sample_metadata.csv", row.names = 1)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~ condition)

# Pre-filtering low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res_ordered), file = "deseq2_results.csv")

# Plot MA plot
plotMA(res, main = "DESeq2 MA Plot")

# Volcano Plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'Differential Expression',
                subtitle = 'Treated vs Control')

# Heatmap of top 30 genes
top_genes <- head(order(res$padj), 30)
norm_counts <- assay(rlog(dds))
pheatmap(norm_counts[top_genes,], cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = TRUE, annotation_col = sample_info)

# PCA Plot
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = "condition")

# GO Enrichment Analysis
sig_genes <- rownames(res[which(res$padj < 0.05), ])
# Convert gene IDs (assuming Gene Symbols)
entrez_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Run GO enrichment
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

# Save and plot GO results
write.csv(as.data.frame(go_results), file = "GO_enrichment_results.csv")
barplot(go_results, showCategory = 15, title = "Top GO Biological Processes")
