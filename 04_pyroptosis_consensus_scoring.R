library(Seurat)
library(AUCell)
library(UCell)
library(singscore)
library(GSVA)

score_signature <- function(obj, genes) {
  expr <- LayerData(obj, assay = "RNA", layer = "data")
  genes <- intersect(genes, rownames(expr))
  gene_sets <- list(Signature = genes)
  set.seed(123)
  rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
  auc <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(nrow(expr) * 0.10))
  ucell <- ScoreSignatures_UCell(expr, features = gene_sets, maxRank = 1500, missing_genes = "skip")
  sing <- simpleScore(rankGenes(as.matrix(expr)), upSet = genes, centerScore = TRUE, knownDirection = TRUE)
  ssgsea <- gsva(ssgseaParam(as.matrix(expr), gene_sets, normalize = FALSE))
  obj <- AddModuleScore(obj, features = gene_sets, name = "Signature", nbin = 24, ctrl = 100, seed = 123)
  scores <- cbind(
    AUCell = as.numeric(getAUC(auc)[1, colnames(obj)]),
    UCell = ucell[colnames(obj), 1],
    singscore = sing[colnames(obj), "TotalScore"],
    ssGSEA = as.numeric(ssgsea[1, colnames(obj)]),
    AddModuleScore = obj$Signature1
  )
  setNames(rowMeans(scale(scores)), colnames(obj))
}

obj <- readRDS("annotated_single_cell.rds")
genes <- readLines("gene_sets/pyroptosis_40_genes.txt")
obj$pyroptosis_score <- score_signature(obj, genes)[colnames(obj)]
obj$score_group <- ifelse(obj$pyroptosis_score > median(obj$pyroptosis_score), "High", "Low")
meta <- obj[[]]
groups <- c(list(All = meta), split(meta, meta$cell_type))
tests <- do.call(rbind, lapply(names(groups), function(type) {
  d <- groups[[type]]
  data.frame(cell_type = type, pvalue = wilcox.test(pyroptosis_score ~ condition_paper, data = d)$p.value)
}))
tests$p_adjust <- p.adjust(tests$pvalue, method = "BH")
write.csv(tests, "pyroptosis_group_comparisons.csv", row.names = FALSE)
saveRDS(obj, "pyroptosis_scored_single_cell.rds")
