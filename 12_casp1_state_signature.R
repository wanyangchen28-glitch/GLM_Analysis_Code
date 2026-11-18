library(Seurat)
library(Matrix)

mac <- readRDS("final_macrophages.rds")
mac <- subset(mac, subset = condition_paper == "GLM")
expr <- LayerData(mac, assay = "RNA", layer = "data")
positive <- expr["CASP1", ] > 0
mean_positive <- rowMeans(expr[, positive])
mean_negative <- rowMeans(expr[, !positive])
var_positive <- pmax(rowMeans(expr[, positive]^2) - mean_positive^2, 0)
var_negative <- pmax(rowMeans(expr[, !positive]^2) - mean_negative^2, 0)
variance_sum <- var_positive + var_negative
v_score <- (mean_negative - mean_positive) / (sqrt(variance_sum) + (variance_sum == 0))
write.csv(data.frame(Gene = rownames(expr), v_score), "casp1_state_transition_signature.csv", row.names = FALSE)
