library(Seurat)
library(Matrix)
library(scTenifoldKnk)
library(clusterProfiler)
library(org.Hs.eg.db)

mac <- readRDS("final_macrophages.rds")
counts <- LayerData(mac, assay = "RNA", layer = "counts")
set.seed(1234)
counts <- counts[, sample(colnames(counts), min(3000, ncol(counts)))]
counts <- counts[rowSums(counts > 0) >= 10, ]
x <- NormalizeData(CreateSeuratObject(counts))
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
markers <- c("CASP1", "NLRP3", "PYCARD", "CASP4", "CASP5", "GSDMD", "IL1B", "IL18", "AIM2", "NINJ1")
genes <- head(unique(c(intersect(markers, rownames(counts)), VariableFeatures(x))), 3000)
counts <- counts[genes, ]
ko <- scTenifoldKnk(
  countMatrix = counts, gKO = "CASP1", qc = TRUE,
  qc_mtThreshold = 0.20, qc_minLSize = 500, qc_minCells = 10,
  nc_lambda = 0, nc_nNet = 10, nc_nCells = min(500, floor(ncol(counts) * 0.8)),
  nc_nComp = 3, nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9,
  td_K = 3, td_maxIter = 1000, td_maxError = 1e-5, td_nDecimal = 3, ma_nDim = 2
)
dr <- ko$diffRegulation
dr$gene <- rownames(dr)
significant <- subset(dr, p.adj < 0.05)
ids <- bitr(significant$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go <- enrichGO(ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
go <- subset(go@result, p.adjust < 0.05)
write.csv(dr, "casp1_knockout_perturbations.csv", row.names = FALSE)
write.csv(go, "casp1_knockout_go_enrichment.csv", row.names = FALSE)
