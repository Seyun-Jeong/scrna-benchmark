import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ground_truth import PBMC3K_MARKERS, best_cell_type_match

ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = ROOT / "data/filtered_gene_bc_matrices/hg19"
RESULTS_DIR = ROOT / "results"
SEEDS = range(10)

os.makedirs(RESULTS_DIR, exist_ok=True)

def build_adata():
    adata = sc.read_10x_mtx(DATA_PATH, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None)
    adata = adata[
        (adata.obs.n_genes_by_counts < 2500) &
        (adata.obs.n_genes_by_counts > 200) &
        (adata.obs.pct_counts_mt < 5)
    ].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata)
    return adata

print("Building base AnnData (deterministic up to PCA)...")
adata_base = build_adata()

summary_rows = []

for seed in SEEDS:
    print(f"\n--- Seed {seed} ---")
    adata = adata_base.copy()

    sc.pp.neighbors(adata, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, random_state=seed)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    n_clusters = adata.obs["leiden"].nunique()
    print(f"  Clusters found: {n_clusters}")

    result = sc.get.rank_genes_groups_df(adata, group=None)

    rows = []
    for cluster in adata.obs["leiden"].unique():
        top_genes = (
            result[result["group"] == cluster]
            .sort_values("scores", ascending=False)
            .head(10)["names"]
            .tolist()
        )
        cell_type, score, all_scores = best_cell_type_match(top_genes)
        rows.append({
            "seed": seed,
            "cluster": cluster,
            "top_genes": ",".join(top_genes),
            "best_cell_type": cell_type,
            "jaccard": round(score, 4),
        })
        print(f"  Cluster {cluster} → {cell_type} (J={score:.3f})")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / f"seed_{seed}_markers.csv", index=False)
    summary_rows.extend(rows)

summary = pd.DataFrame(summary_rows)
summary.to_csv(RESULTS_DIR / "summary_all_seeds.csv", index=False)
print(f"\nDone. Results in {RESULTS_DIR}/")
