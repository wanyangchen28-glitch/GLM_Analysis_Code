library(Seurat)
library(dplyr)

obj <- readRDS("pyroptosis_scored_single_cell.rds")
d <- FetchData(obj, vars = c("GBP2", "CASP1", "cell_type", "condition_paper", "pyroptosis_score", "score_group"))
expression <- d %>% group_by(cell_type, condition_paper) %>%
  summarise(across(c(GBP2, CASP1), list(mean = mean, positive_fraction = ~mean(.x > 0))), .groups = "drop")
mac <- subset(d, cell_type == "Macrophage")
glm <- subset(mac, condition_paper == "GLM")
correlation <- do.call(rbind, lapply(c("GBP2", "CASP1"), function(gene) {
  test <- cor.test(glm[[gene]], glm$pyroptosis_score, method = "spearman", exact = FALSE)
  data.frame(gene, rho = unname(test$estimate), pvalue = test$p.value)
}))
comparisons <- data.frame(
  comparison = c("GLM_vs_HC", "High_vs_Low"),
  pvalue = c(wilcox.test(CASP1 ~ condition_paper, data = mac)$p.value,
             wilcox.test(CASP1 ~ score_group, data = mac)$p.value)
)
write.csv(expression, "candidate_expression_summary.csv", row.names = FALSE)
write.csv(correlation, "candidate_pyroptosis_correlations.csv", row.names = FALSE)
write.csv(comparisons, "casp1_expression_comparisons.csv", row.names = FALSE)
