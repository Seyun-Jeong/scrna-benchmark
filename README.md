# PBMC3k Single-Cell RNA-seq Analysis with Scanpy
<<<<<<< HEAD
=======

Reproducible implementation of a standard PBMC3k single-cell RNA-seq analysis workflow using Scanpy.

This repository documents an end-to-end **CPU-only Scanpy pipeline** on the canonical 10x Genomics PBMC3k dataset. It does **not** yet include a CPU–GPU benchmark, RAPIDS-singlecell, Seurat/ScaleSC comparisons, runtime measurements, or ARI/module-profile validation.

## Dataset

10x Genomics PBMC3k (`filtered_gene_bc_matrices/hg19/`):

| Stage | Cells | Genes |
|---|---|---|
| Raw load | 2,700 | 32,738 |
| After initial filter (`min_genes=100`, `min_cells=3`) | 2,700 | 13,714 |
| After QC filter (200–2,500 genes, MT < 5%) | 2,638 | 13,714 |

## Workflow

`notebooks/pbmc3k/scanpy.ipynb`:

1. Load 10x matrix
2. Quality-control plots
3. Cell and gene filtering
4. Total-count normalization
5. `log1p`
6. Selection of 2,000 highly variable genes
7. PCA
8. Neighbor graph and UMAP
9. Leiden clustering
10. Wilcoxon marker-gene analysis

## Validation (seeds 0–9)

`notebooks/pbmc3k/validation.ipynb` and `scripts/seed_sweep.py` run the pipeline across random seeds 0–9 and compare Leiden cluster marker genes against published PBMC3k cell-type marker sets (Jaccard overlap). Results are saved under `results/`.

## Repository layout

```
data/                          # PBMC3k 10x matrix
notebooks/pbmc3k/
  scanpy.ipynb                   # Main Scanpy workflow
  validation.ipynb               # Seed sweep + marker validation
scripts/
  seed_sweep.py                  # Batch seed runner
  ground_truth.py                # PBMC3k marker gene reference
results/                         # Seed sweep outputs
```

## Requirements

- Python 3.10+
- scanpy
- pandas
- numpy

## Run

```bash
jupyter notebook notebooks/pbmc3k/scanpy.ipynb
```

Or run the seed sweep:

```bash
python scripts/seed_sweep.py
```

## Scope and related work

**What this repo demonstrates:** a reproducible PBMC3k Scanpy analysis workflow with seed-based marker validation.

**What is not in this repository:** CPU–GPU benchmarking, RAPIDS-singlecell, Seurat/ScaleSC comparisons, runtime tables, ARI comparisons, module-profile correlation, or mouse-brain (1.3M-cell) experiments from my SC Workshops 2025 paper and M.S. thesis. Those experiments were executed on SLU BioHPC infrastructure; release is subject to thesis review and storage constraints.

**Next step:** add a RAPIDS-singlecell implementation for side-by-side comparison with this Scanpy pipeline.

## Citation

If you use this workflow, please cite the SC Workshops 2025 paper for the broader CPU/GPU benchmarking context:

> Jeong, S., et al. *Evaluating Accuracy and Performance Tradeoffs in GPU Accelerated Single Cell RNA-seq Analysis.* DRBSD-11, SC Workshops 2025. https://doi.org/10.1145/3731599.3767378

## License

See [LICENSE](LICENSE).
>>>>>>> dev
