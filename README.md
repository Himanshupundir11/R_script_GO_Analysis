
# 🧬 Gene Expression Analysis with DESeq2 & GO Enrichment

This repository demonstrates a full RNA-seq differential expression analysis workflow in R using DESeq2, along with functional enrichment using Gene Ontology (GO).

---

## 📁 Project Structure

```
.
├── gene_expression_analysis.R       # Main analysis script
├── counts_matrix.csv                # Simulated gene expression count matrix
├── sample_metadata.csv              # Metadata with sample condition info
├── deseq2_results.csv               # DE results (auto-generated)
├── GO_enrichment_results.csv        # Enriched GO terms (auto-generated)
```

---

## 🔬 Workflow Overview

### Step 1: Load and Preprocess
- Reads in raw count matrix and sample metadata
- Filters low-expression genes

### Step 2: Differential Expression (DE)
- Conducted using **DESeq2**
- Results filtered by adjusted p-value (padj < 0.05)
- Outputs:
  - MA Plot
  - Volcano Plot
  - Top 30 genes heatmap
  - PCA plot for sample clustering

### Step 3: Functional Enrichment
- Converts DE gene symbols to Entrez IDs
- Performs GO Biological Process enrichment via `clusterProfiler`
- Plots top 15 enriched GO terms

---

## 📦 Dependencies

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "pheatmap", "clusterProfiler", "org.Hs.eg.db"))
```

---

## 📈 Sample Output

You can view the differential expression result in:
- `deseq2_results.csv`
- `GO_enrichment_results.csv`

Plots will be generated automatically by the R script.

---

## 🚀 Run the Script

```r
source("gene_expression_analysis.R")

Prepared for educational and demonstration purposes for bioinformatics pipelines using R.
