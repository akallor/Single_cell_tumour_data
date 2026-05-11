#Perform differential analysis using the DESeq2 alogrithm to analyze tumour vs normal

library(DESeq2)

# Aggregate counts per sample per cluster
pseudo <- aggregateAcrossCells(sce, ids = colData(sce)[,c("cluster","sample")])

# For each cluster, run DESeq2 comparing tumour vs normal
dds <- DESeqDataSet(pseudo, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumour", "normal"))
