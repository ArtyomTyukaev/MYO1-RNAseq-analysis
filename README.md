# **MYO1-RNAseq-analysis**

This repository contains all scripts, environments, and processed data used in the RNA-seq analysis of the **MYO1 gene family (MYO1A–MYO1H)** across a large panel of human cell lines (GSE240542).  
It provides a fully reproducible workflow — from raw FASTQ preprocessing to gene expression quantification, TPM calculation, and figure generation.

---

## **Directory structure**

```
MYO1-RNAseq-analysis/
│
├── data/                # Processed input data used in all R / Python scripts
│   ├── expression_tables_fraction/
│   ├── expression_tables_unique/
│   ├── cell_line_annotation_update.csv
│   ├── cell_line_annotation.tsv
│   ├── cell_lines_into_groups.tsv
│   └── GSE240542_SraRunTable.csv
│
├── envs/                # Conda environments & R session info
│   ├── define_chain_specificity.yml
│   ├── python_analysis_environment.yml
│   ├── RNAseq_preprocessing.yml
│   └── R_sessionInfo.txt
│
├── scripts/             # All analysis scripts
│   ├── 1_fasterq_SSD.sh                 # fastq generation from .sra (parallel, SSD tmp)
│   ├── 2_run_fastp.sh                   # adapter trimming + QC
│   ├── 3_run_STAR.sh                    # genome alignment (STAR)
│   ├── 4_chain_specificity.sh           # strand specificity inference (RSeQC)
│   ├── 5_run_featureCounts.sh           # gene-level quantification
│   ├── 6_get_genes_counts.sh            # extraction of counts for each SRR
│   ├── calculate_TPM_and_merge_SRR_to_cell_lines.ipynb
│   └── R_plots.Rmd                      # all final visualization code (boxplots, heatmaps, UMAP, violins)
│
└── README.md
```

---

## **Project overview**

This repository implements a complete RNA-seq workflow focused on the **MYO1 gene family**:

### **1. Raw data preprocessing**
- Conversion of `.sra` to `.fastq.gz` (via `fasterq-dump`)
- Adapter trimming and QC (`fastp`)
- Strand specificity detection (`RSeQC`)

### **2. Alignment & quantification**
- Genome alignment with **STAR** (GENCODE v48, GRCh38.p14)
- Gene-level counts from **featureCounts**
- Unique-only and multimapping-fraction modes

### **3. Expression metrics**
- TPM and log2(TPM+1) calculation per SRR
- Aggregation into per–cell-line mean expression

### **4. Annotation**
- Metadata from SraRunTable
- Cell-line grouping into tumor types and categories

### **5. Figures and downstream analysis**
All figures used in the manuscript can be reproduced using:

- `R_plots.Rmd`  
- `calculate_TPM_and_merge_SRR_to_cell_lines.ipynb`

These include heatmaps, boxplots, violins, UMAPs, and more.

## Precomputed expression matrices

All final expression matrices used in downstream analysis are provided in:
```
/data/expression_tables_unique/
/data/expression_tables_fraction/
```
### Unique vs Fraction (what they mean)

**unique**  
- Only uniquely mapping read pairs are counted  
- Multimapping reads are *discarded*  

**fraction**  
- Multimapping reads are included  
- Each read is assigned fractionally across all compatible loci using the  
  `-M --fraction` mode of **featureCounts**  

Both versions are provided for transparency and reproducibility. We used unique.

---
## **Reproducible environments**

All environments required to run the pipeline are in `envs/`:

- `RNAseq_preprocessing.yml`  
- `define_chain_specificity.yml`  
- `python_analysis_environment.yml`  
- `R_sessionInfo.txt` (exact R setup used to generate the figures)

Create environments with:

```bash
conda env create -f envs/RNAseq_preprocessing.yml
conda env create -f envs/define_chain_specificity.yml
conda env create -f envs/python_analysis_environment.yml
```

---
## **How to reproduce the pipeline**
1. Run preprocessing scripts sequentially
```
scripts/1_fasterq_SSD.sh
scripts/2_run_fastp.sh
scripts/3_run_STAR.sh
scripts/4_chain_specificity.sh
scripts/5_run_featureCounts.sh
scripts/6_get_genes_counts.sh
```
2. Calculate TPM and other expression metrics and merge replicates
``scripts/calculate_TPM_and_merge_SRR_to_cell_lines.ipynb``
3. Generate all figures
```scripts/R_plots.Rmd```
