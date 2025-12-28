# Tumor scRNA-seq Analysis Toolkit (Scanpy)

## Problem
Single-cell RNA-seq analysis of a tumor microenvironment to (1) identify major cell populations, (2) interpret cluster markers, and (3) compare malignant vs immune compartments with pathway-level interpretation.

## Dataset
- Public melanoma tumor scRNA-seq dataset (GEO: GSE72056; processed expression matrix).
- Analysis starts from the provided processed matrix (not FASTQs).

## Methods (high level)
1. **QC & filtering**: genes/cell, counts/cell, mitochondrial fraction; basic filtering to remove low-quality cells.
2. **Preprocessing**: normalize, log-transform, highly-variable gene selection, scaling.
3. **Embedding & clustering**: PCA → neighbors → UMAP → Leiden clustering.
4. **Markers & annotation**: Wilcoxon marker ranking per cluster + marker-based cell type labeling (immune vs malignant).
5. **Compartment DE (no pseudoreplication)**: differential expression computed using **cluster-level pseudobulk** (clusters as replicates) rather than treating cells as independent replicates.
6. **GO enrichment**: GO Biological Process enrichment on significant gene lists (malignant-up vs immune-up).

## Results
### Clustering + annotation
- Leiden clustering produced **18 clusters**.
- Marker genes recapitulate canonical tumor microenvironment biology (e.g., cytotoxic lymphocytes show CD8A/NKG7/PRF1 signatures).

### Malignant vs immune differential expression (cluster pseudobulk)
- Output table: `results/tables/de_malignant_vs_immune_cluster_pseudobulk.tsv`
- Volcano plot: `figures/de/volcano_malignant_vs_immune.png`

Top signals include melanoma lineage/tumor-associated genes (e.g., **PMEL**, **GPNMB**) enriched in malignant clusters and immune-associated programs enriched in immune clusters.

### Pathway enrichment (GO Biological Process)
- Tables:
  - `results/enrichment/go_bp_malignant_up.csv`
  - `results/enrichment/go_bp_immune_up.csv`
- Plots:
  - `figures/enrichment/go_bp_malignant_up_dotplot.png`
  - `figures/enrichment/go_bp_immune_up_dotplot.png`

## How to reproduce
### Environment
Create and activate the conda environment:
```bash
conda env create -f envs/scanpy.yaml
conda activate scrna-tumor

python scripts/01_qc.py
python scripts/02_preprocess.py
python scripts/03_cluster.py
python scripts/04_markers.py
python scripts/05_annotation.py
python scripts/06_attach_geo_metadata.py
python scripts/07_de_malignant_vs_immune_cluster_pseudobulk.py
python scripts/08_go_enrichment.py
python scripts/09_volcano.py
