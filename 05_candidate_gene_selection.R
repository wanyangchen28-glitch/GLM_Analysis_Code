library(Seurat)

obj <- readRDS("pyroptosis_scored_single_cell.rds")
glm <- subset(obj, subset = condition_paper == "GLM")
expr <- LayerData(glm, assay = "RNA", layer = "data")
rho <- apply(expr, 1, function(x) cor(x, glm$pyroptosis_score, method = "spearman"))
rho <- sort(rho[is.finite(rho) & rho > 0], decreasing = TRUE)
top100 <- names(head(rho, 100))
deg <- FindMarkers(obj, group.by = "score_group", ident.1 = "High", ident.2 = "Low", test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
up <- rownames(deg)[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.5]
pyroptosis <- readLines("gene_sets/pyroptosis_40_genes.txt")
candidates <- Reduce(intersect, list(top100, up, pyroptosis))
writeLines(candidates, "string_candidate_genes.txt")

edges <- read.delim("string_interactions.tsv")
edges <- subset(edges, combined_score >= 0.4)
nodes <- unique(c(edges$preferredName_A, edges$preferredName_B))
write.table(edges, "ppi_network.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(intersect(candidates, nodes), "ppi_retained_nodes.txt")
rankings <- lapply(c("MCC", "EPC", "MNC", "Degree"), function(method) {
  d <- read.csv(paste0("cytoHubba_", method, ".csv"))
  head(d$Name[order(d$Rank)], 2)
})
writeLines(Reduce(intersect, rankings), "key_candidate_genes.txt")
