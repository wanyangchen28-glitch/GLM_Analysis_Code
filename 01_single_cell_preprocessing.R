library(Seurat)
library(harmony)

samples <- read.csv("sample_metadata.csv")
objects <- lapply(seq_len(nrow(samples)), function(i) {
  x <- CreateSeuratObject(Read10X(samples$directory[i]))
  x$sample_id <- samples$sample_id[i]
  x$condition_paper <- samples$condition[i]
  x
})
obj <- merge(objects[[1]], y = objects[-1], add.cell.ids = samples$sample_id)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)
obj <- RunHarmony(obj, group.by.vars = "sample_id")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.7)
obj <- JoinLayers(obj, assay = "RNA")
cluster_markers <- FindAllMarkers(obj, only.pos = TRUE, test.use = "wilcox")
write.csv(cluster_markers, "cluster_markers.csv", row.names = FALSE)

annotation <- read.csv("cluster_annotations.csv")
obj$cell_type <- annotation$cell_type[match(obj$seurat_clusters, annotation$cluster)]
obj <- subset(obj, subset = seurat_clusters != "24")
saveRDS(obj, "annotated_single_cell.rds")

composition <- as.data.frame(table(obj$sample_id, obj$cell_type))
names(composition) <- c("sample_id", "cell_type", "n_cells")
composition$proportion <- composition$n_cells / ave(composition$n_cells, composition$sample_id, FUN = sum)
write.csv(composition, "cell_composition_by_sample.csv", row.names = FALSE)
composition <- as.data.frame(table(obj$condition_paper, obj$cell_type))
names(composition) <- c("condition", "cell_type", "n_cells")
composition$proportion <- composition$n_cells / ave(composition$n_cells, composition$condition, FUN = sum)
write.csv(composition, "cell_composition_by_condition.csv", row.names = FALSE)
