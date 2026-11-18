library(Seurat)
library(CellChat)

obj <- readRDS("annotated_single_cell.rds")
mac <- readRDS("final_macrophages.rds")
obj <- subset(obj, cells = colnames(obj)[obj$cell_type != "Macrophage" | colnames(obj) %in% colnames(mac)])

run_cellchat <- function(obj, split_casp1 = FALSE) {
  expr <- LayerData(obj, assay = "RNA", layer = "data")
  meta <- obj[[]]
  meta$group <- as.character(meta$cell_type)
  if (split_casp1) {
    take <- meta$group == "Macrophage"
    meta$group[take] <- ifelse(expr["CASP1", take] > 0, "CASP1+ Macrophage", "CASP1- Macrophage")
  }
  keep <- meta$group %in% names(which(table(meta$group) >= 10))
  cc <- createCellChat(expr[, keep], meta = meta[keep, ], group.by = "group")
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE, population.size = FALSE, nboot = 100)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc, thresh = 0.05)
  cc <- aggregateNet(cc, thresh = 0.05)
  netAnalysis_computeCentrality(cc, slot.name = "netP")
}

hc <- run_cellchat(subset(obj, subset = condition_paper == "HC"))
glm <- run_cellchat(subset(obj, subset = condition_paper == "GLM"))
casp1 <- run_cellchat(subset(obj, subset = condition_paper == "GLM"), split_casp1 = TRUE)
comparison <- data.frame(group = c("HC", "GLM"), count = c(sum(hc@net$count), sum(glm@net$count)), strength = c(sum(hc@net$weight), sum(glm@net$weight)))
write.csv(comparison, "communication_counts_and_strength.csv", row.names = FALSE)
write.csv(subsetCommunication(glm, sources.use = "Macrophage", thresh = 0.05), "glm_macrophage_outgoing_communication.csv", row.names = FALSE)
write.csv(subsetCommunication(glm, targets.use = "Macrophage", thresh = 0.05), "glm_macrophage_incoming_communication.csv", row.names = FALSE)
write.csv(subsetCommunication(casp1, thresh = 0.05), "casp1_positive_negative_communication.csv", row.names = FALSE)
saveRDS(list(HC = hc, GLM = glm, CASP1 = casp1), "cellchat_results.rds")
