library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(igraph)

mac <- readRDS("final_macrophages.rds")
rna <- FindMarkers(mac, group.by = "condition_paper", ident.1 = "GLM", ident.2 = "HC", test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1)
rna$gene <- rownames(rna)
rna_up <- subset(rna, p_val_adj < 0.05 & avg_log2FC > 0)$gene
protein <- read.csv("protein_differential_statistics.csv")
mapping <- read.csv("human_rat_ortholog_mapping.csv")
mapping <- subset(mapping, unambiguous & ensembl_supported & support_n >= 3)
protein <- merge(protein, mapping, by.x = "GeneName", by.y = "symbol")
protein <- subset(protein, representative & mapping_ok)
protein_up <- subset(protein, pvalue < 0.05 & FC >= 1.5)
common <- intersect(rna_up, protein_up$human_symbol)
background <- intersect(rna$gene, protein$human_symbol)
writeLines(common, "shared_upregulated_genes.txt")

enrichment <- function(genes, background, database, organism) {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = database)$ENTREZID
  bg <- if (is.null(background)) NULL else bitr(background, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = database)$ENTREZID
  go <- enrichGO(ids, universe = bg, OrgDb = database, ont = "ALL", pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)
  kegg <- enrichKEGG(ids, universe = bg, organism = organism, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)
  go <- setReadable(go, OrgDb = database, keyType = "ENTREZID")@result
  kegg <- setReadable(kegg, OrgDb = database, keyType = "ENTREZID")@result
  go$category <- go$ONTOLOGY
  kegg$category <- "KEGG"
  columns <- c("ID", "Description", "pvalue", "p.adjust", "geneID", "Count", "category")
  rbind(go[, columns], kegg[, columns])
}

terms <- enrichment(common, background, org.Hs.eg.db, "hsa")
terms <- subset(terms, pvalue < 0.05 & category %in% c("BP", "KEGG"))
sets <- strsplit(terms$geneID, "/", fixed = TRUE)
pairs <- combn(seq_len(nrow(terms)), 2)
edges <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  shared <- length(intersect(sets[[i]], sets[[j]]))
  jaccard <- shared / length(union(sets[[i]], sets[[j]]))
  if (shared >= 3 & jaccard >= 0.25) data.frame(from = terms$ID[i], to = terms$ID[j], weight = jaccard)
}))
network <- graph_from_data_frame(edges, directed = FALSE, vertices = terms$ID)
set.seed(42)
modules <- membership(cluster_louvain(network, weights = E(network)$weight))
terms$module <- modules[terms$ID]
representatives <- do.call(rbind, lapply(split(terms, terms$module), function(x) x[which.min(x$pvalue), ]))
write.csv(terms, "functional_network_nodes.csv", row.names = FALSE)
write.csv(edges, "functional_network_edges.csv", row.names = FALSE)
write.csv(representatives, "functional_module_representatives.csv", row.names = FALSE)

rna_terms <- enrichment(rna_up, rna$gene, org.Hs.eg.db, "hsa")
rat <- read.csv("protein_differential_statistics.csv")
rat_up <- subset(rat, pvalue < 0.05 & FC >= 1.5)$GeneName
rat_terms <- enrichment(rat_up, NULL, org.Rn.eg.db, "rno")
write.csv(rna_terms, "macrophage_upregulated_gene_enrichment.csv", row.names = FALSE)
write.csv(rat_terms, "rat_upregulated_protein_enrichment.csv", row.names = FALSE)
rna_terms$key <- ifelse(rna_terms$category == "KEGG", sub("^hsa", "", rna_terms$ID), rna_terms$ID)
rat_terms$key <- ifelse(rat_terms$category == "KEGG", sub("^rno", "", rat_terms$ID), rat_terms$ID)
shared_terms <- merge(subset(rna_terms, pvalue < 0.05), subset(rat_terms, pvalue < 0.05), by = c("key", "category"), suffixes = c("_RNA", "_protein"))
shared_terms$shared_genes <- mapply(function(human, rat) {
  rat <- strsplit(rat, "/", fixed = TRUE)[[1]]
  human <- strsplit(human, "/", fixed = TRUE)[[1]]
  paste(intersect(human, mapping$human_symbol[match(rat, mapping$symbol)]), collapse = "/")
}, shared_terms$geneID_RNA, shared_terms$geneID_protein)
write.csv(shared_terms, "shared_enriched_pathways.csv", row.names = FALSE)
