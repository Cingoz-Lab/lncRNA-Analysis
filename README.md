# lncRNA Analysis

End-to-end downstream analysis pipeline for RNA-seq differential expression data, with a focus on identifying and characterising **long non-coding RNAs (lncRNAs)**.

**Author:** Furkan Emre Bora — [Cingoz Lab](https://github.com/Cingoz-Lab)

---

## What this pipeline does

Starting from a DESeq2 results table and a GENCODE GTF annotation file, the script produces:

| # | Output | Description |
|---|--------|-------------|
| 01 | log2FC distribution | QC — should be roughly symmetric around 0 |
| 02 | -log10(pvalue) distribution | QC — enrichment at high values expected |
| 03 | -log10(padj) distribution | QC — FDR landscape after BH correction |
| 04 | baseMean violin by class | lncRNA vs mRNA expression levels |
| 05 | Exon count histogram | lncRNAs typically have fewer exons |
| 06 | Gene length histogram | lncRNA length distribution vs mRNA |
| 07 | MA plot | log fold change vs mean expression, all genes |
| 08 | Volcano plot (lncRNA) | key lncRNAs labelled |
| 09 | Top lncRNA bar plot | top hits ranked by padj |
| 10 | Top lncRNA lollipop plot | alternative view of top hits |
| 11 | Biotype composition bar | biotypes among significant genes |
| 12 | Expression vs \|log2FC\| scatter | detect LFC inflation at low counts |
| 13 | Expression vs -log10(padj) scatter | power increases with expression |
| 14 | GSEA Hallmark dotplot | top enriched pathways |
| 15 | GSEA enrichment curve | top-ranked pathway ES curve |
| 16 | GSEA ridgeplot | per-gene rank distribution per pathway |
| 17 | GSEA emapplot | pathway similarity network |
| 18-21 | PCA / heatmaps | sample-level QC (requires DESeq2 `dds` object) |

---

## Requirements

### R version
R >= 4.2 recommended.

### Packages
All packages are **installed automatically** the first time the script is run.

**CRAN:** `data.table`, `dplyr`, `stringr`, `tibble`, `ggplot2`, `openxlsx`, `ggrepel`, `forcats`, `scales`

**Bioconductor:** `clusterProfiler`, `enrichplot`, `msigdbr`, `fgsea`

---

## Input files

| File | Description |
|------|-------------|
| `deg_results.xlsx` | Excel workbook with a sheet called **`deg_results`** containing DESeq2 output. Required columns: `Column1` (ENSEMBL gene ID), `baseMean`, `log2FoldChange`, `pvalue`, `padj`, `Gene` (symbol). |
| `gencode.annotation.gtf` | GENCODE GTF file (tested with v36–v46, hg38). Any recent build works. Download from [gencodegenes.org](https://www.gencodegenes.org/human/). |

---

## How to run

1. Clone the repository and place your input files in the same folder:
   ```bash
   git clone https://github.com/Cingoz-Lab/lncRNA-Analysis.git
   cd lncRNA-Analysis
   ```

2. Open `All_In.R` and edit **Section 0** (the only section that needs changes):
   ```r
   base_dir  <- "."                        # folder containing your input files
   xlsx_file <- file.path(base_dir, "deg_results.xlsx")
   gtf_file  <- file.path(base_dir, "gencode.annotation.gtf")
   outdir    <- file.path(base_dir, "lncRNA_Analysis_Output")
   ```

3. Run from the terminal:
   ```bash
   Rscript All_In.R
   ```
   Or source inside an R session:
   ```r
   source("All_In.R")
   ```

4. All outputs are written to `lncRNA_Analysis_Output/plots/`, `tables/`, and `gsea/`.

---

## Significance thresholds

Three tiers are pre-defined to balance stringency and power:

| Set | padj | pvalue | \|log2FC\| | Use |
|-----|------|--------|-----------|-----|
| **strict** | < 0.05 | < 0.05 | ≥ 1.0 | Publication-level reporting |
| **broad** | < 0.05 | — | ≥ 1.0 | Default filter |
| **moderate** | < 0.10 | < 0.05 | ≥ 0.58 | Pathway/network enrichment |

> **Why so few significant lncRNAs?**
> DESeq2's independent filtering removes low-count genes before multiple testing correction, leaving ~77 % of genes with `padj = NA`. This is expected behaviour. Using raw p-value for the volcano plot and relaxed thresholds for exploratory analysis is standard practice in this situation.

---

## Sample-level plots (optional)

Plots 18–21 (PCA, correlation heatmap, top-variable-gene heatmap, lncRNA heatmap) require a DESeq2 `dds` object in your R session. Run your DESeq2 pipeline first and then source this script:

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)
dds <- DESeq(dds)
source("All_In.R")   # vsd is created automatically from dds
```

---

## Output structure

```
lncRNA_Analysis_Output/
├── plots/
│   ├── 01_log2FC_distribution.png
│   ├── 02_neglog10pvalue_distribution.png
│   ├── ...
│   └── 17_GSEA_Hallmark_emapplot.png
├── tables/
│   ├── DEG_annotated_all_genes.csv
│   ├── DEG_annotated_lncRNA_only.csv
│   ├── SIG_lncRNA_STRICT.csv
│   ├── SIG_lncRNA_MODERATE.csv
│   ├── SIG_lncRNA_padj0.05_lfc1.csv
│   └── QC_summary.csv
└── gsea/
    └── GSEA_Hallmark_results.csv
```
