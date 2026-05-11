#!/usr/bin/env python
#Load each count matrix from the previous alignment step into Scanpy (Python) or Seurat (R).
#This is to filter real cells from empty droplets.

import scanpy as sc
adata = sc.read_10x_mtx("Tumor1/Solo.out/GeneFull/raw/")

# Standard QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Filter thresholds (tune to your data's knee plots)
adata = adata[adata.obs.n_genes_by_counts > 200]
adata = adata[adata.obs.n_genes_by_counts < 6000]
adata = adata[adata.obs.pct_counts_mt < 20]
