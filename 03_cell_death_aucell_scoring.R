library(Seurat)
library(AUCell)
library(clusterProfiler)

obj <- readRDS("annotated_single_cell.rds")
expr <- LayerData(obj, assay = "RNA", layer = "data")
gmt <- read.gmt("gene_sets/cell_death_programs.gmt")
gene_sets <- lapply(split(gmt$gene, gmt$term), intersect, rownames(expr))
set.seed(123)
rankings <- AUCell_buildRankings(expr, plotStats = FALSE, keepZeroesAsNA = FALSE)
auc <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(nrow(expr) * 0.10), normAUC = TRUE)
scores <- scale(t(getAUC(auc)))
summary <- do.call(rbind, lapply(colnames(scores), function(program) {
  do.call(rbind, lapply(unique(obj$cell_type), function(type) {
    x <- scores[obj$cell_type == type, program]
    data.frame(program, cell_type = type, mean_z = mean(x), fraction_z_gt0 = mean(x > 0))
  }))
}))
write.csv(summary, "cell_death_score_summary.csv", row.names = FALSE)
