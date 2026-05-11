#!/usr/bin/env python

#Merge all three AnnData objects and correct for sample-of-origin batch effects

import scanpy.external as sce

# Merge
adata = sc.concat(
    {"Tumor1": adata_t1, "Tumor2": adata_t2, "Normal": adata_n},
    label="sample"
)

# Normalise and find HVGs
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample")
sc.pp.pca(adata, use_highly_variable=True)

# Harmony integration
sce.pp.harmony_integrate(adata, "sample")

##Clustering & cell type annotation

# Neighbourhood graph and UMAP on Harmony embedding
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
sc.tl.umap(adata)

# Leiden clustering (tune resolution 0.3–1.0)
sc.tl.leiden(adata, resolution=0.5)

# Visualise clusters coloured by sample and condition
sc.pl.umap(adata, color=["leiden", "sample", "condition"])
