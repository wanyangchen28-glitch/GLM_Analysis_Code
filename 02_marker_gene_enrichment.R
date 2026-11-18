library(Seurat)
library(Matrix)
library(presto)
library(clusterProfiler)
library(org.Hs.eg.db)

obj <- readRDS("annotated_single_cell.rds")
expr <- LayerData(obj, assay = "RNA", layer = "data")
expr <- expr[rowSums(expr > 0) > 0, ]
markers <- wilcoxauc(expr, y = obj$cell_type)
linear_expr <- expr
linear_expr@x <- expm1(linear_expr@x)
total <- rowSums(linear_expr)
for (type in unique(obj$cell_type)) {
  take <- obj$cell_type == type
  inside <- rowSums(linear_expr[, take])
  log2fc <- log2((inside + 1) / sum(take)) - log2((total - inside + 1) / sum(!take))
  rows <- markers$group == type
  markers$avg_log2FC[rows] <- log2fc[match(markers$feature[rows], rownames(expr))]
}
markers <- subset(markers, padj < 0.05 & avg_log2FC > 0.5 & pct_in >= 25 & auc > 0.5)
go <- lapply(split(markers$feature, markers$group), function(genes) {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  enrichGO(ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
})
write.csv(markers, "cell_type_markers.csv", row.names = FALSE)
saveRDS(go, "marker_go_enrichment.rds")
