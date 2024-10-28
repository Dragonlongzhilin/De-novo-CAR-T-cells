#' @description: clustering and annotation

library(Seurat)
library(harmony)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(clustree)
library(dplyr)
library(ggpubr)
set.seed(101)
library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #50G
setwd("/data/scRNAseq-CARTcell-20220501") # 91 server
set.resolutions <- seq(0.1, 1.2, by = 0.1)

CART.scRNA <- readRDS("CART.scRNA.rds")

# remove MT and RPS gene
genes <- rownames(CART.scRNA)[grep("^MT-|RPL|RPS", rownames(CART.scRNA))]
CART.scRNA <- subset(CART.scRNA, features = setdiff(rownames(CART.scRNA), genes))

#####################2.correct batch effect########################
DefaultAssay(CART.scRNA) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
CART.scRNA <- NormalizeData(CART.scRNA, verbose = FALSE)
CART.scRNA <- CellCycleScoring(CART.scRNA, s.features = s.genes, g2m.features = g2m.genes)
CART.scRNA@meta.data$CC.Difference <- CART.scRNA@meta.data$S.Score - CART.scRNA@meta.data$G2M.Score
## select the top variable genes
DefaultAssay(CART.scRNA) <- "RNA"
CART.scRNA.list <- SplitObject(CART.scRNA, split.by = "orig.ident")
CART.scRNA.list.Standard <- lapply(CART.scRNA.list, function(x){
	DefaultAssay(x) <- "RNA"
	x <- NormalizeData(x, verbose = FALSE)
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
	return(x)
})

source(file = "/home/scRNA.Integrate.multipleSample.R")
pdf("2.Cluster/Seurat.Standard.PC30.feature3000.pdf")
Seurat.Standard <- Seurat.integration.reduceDimension(seurat.lists = CART.scRNA.list.Standard, assay = "RNA", set.resolutions = set.resolutions, vars.to.regress = c("mt_ratio", "nCount_RNA"),
                                                      adjusted.cellCycle = T, PC = 30, nfeatures = 3000, npcs = 50)
dev.off()
Seurat.Standard$seurat_clusters <- Seurat.Standard$integrated_snn_res.0.3
saveRDS(Seurat.Standard, file = "2.Cluster/Seurat.Standard.rds")

#####################3.Annotated cell cluster########################
DefaultAssay(Seurat.Standard) <- "RNA"
pdf("2.Cluster/cluster.FeaturePlot.pdf")
FeaturePlot(Seurat.Standard, reduction = "umap", cols = c("lightgrey", "red"), features = c("CD3D", "CD68", "CD4", "CD8A"), ncol = 2)
FeaturePlot(Seurat.Standard, reduction = "tsne", cols = c("lightgrey", "red"), features = c("CD3D", "CD68", "CD4", "CD8A"), ncol = 2)
dev.off()

featureMarker <- c("CD3D", "CD3E", "CD4", "IL7R", "CD40LG", "FOXP3", "IL2RA",
                   "CD8A", "CD8B", "MKI67", "TOP2A", "STMN1", "CXCL10", "CCL3", "CCND1",
                   "SELL", "CCR7", "TCF7", "CD69", "LAMP1", "ZNF683", "ITGAE", "PRDM1",
                   "GZMA", "GZMB", "GZMK", "GZMH", "IFNG", "GNLY", "KLRD1", "KLRB1",
                   "PDCD1", "CTLA4", "HAVCR2", "LAG3", "VSIR", "TOX")
pdf("2.Cluster/cluster.Dotplot.pdf", height = 5, width = 6)
DotPlot(Seurat.Standard, features = featureMarker, assay = "RNA", group.by = "seurat_clusters", cols = c("#1e90ff", "#ff5a36"), dot.scale = 4) + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)) + NoLegend()
dev.off()

## differential analysis
Idents(Seurat.Standard) <- Seurat.Standard$seurat_clusters
cluster.all.markers <- FindAllMarkers(Seurat.Standard, only.pos = TRUE,
                                      group.by = "seurat_clusters",
                                      test.use = "MAST")
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05 & cluster.all.markers$avg_log2FC>0.25),]
saveRDS(cluster.sig.markers, file = "2.Cluster/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "2.Cluster/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)

## assign cell population
Seurat.Standard$seurat_clusters <- Seurat.Standard$integrated_snn_res.0.3
clusters <- plyr::mapvalues(Seurat.Standard$seurat_clusters, 
                            from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
                            to = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
Seurat.Standard$cellType <- paste0("C", clusters)
Seurat.Standard$cellType <- factor(Seurat.Standard$cellType, levels = c("C5", "C12", "C3", "C11", "C4", "C1", "C2", "C7", "C9", "C6", "C8", "C10"))

pdf("2.Cluster/AnnotatedCellType/celltype.pdf")
DimPlot(object = Seurat.Standard, reduction = 'umap', label = TRUE, repel = T, group.by = "cellType")+NoLegend()
DimPlot(object = Seurat.Standard, reduction = 'umap', label = F, group.by = "orig.ident")
DimPlot(object = Seurat.Standard, reduction = 'umap', label = T, repel = T, group.by = "cellType", split.by = "orig.ident", ncol=2)+NoLegend()
dev.off()

DefaultAssay(Seurat.Standard) <- "RNA"
pdf("2.Cluster/AnnotatedCellType/phenotype.markerExpression.pdf")
FeaturePlot(Seurat.Standard, reduction = "umap", pt.size = 0.8, cols = c("lightgrey", "red"), , features = "GZMB", split.by = "orig.ident")+patchwork::plot_layout(ncol = 2, nrow = 2)
FeaturePlot(Seurat.Standard, reduction = "umap", pt.size = 0.8, cols = c("lightgrey", "red"), features = "MKI67", split.by = "orig.ident")+patchwork::plot_layout(ncol = 2, nrow = 2)
dev.off()

#### cell ratio
cellNum <- as.data.frame(table(Seurat.Standard@meta.data[,c("orig.ident")]))
cellRatio <- as.data.frame(table(Seurat.Standard@meta.data[,c("orig.ident", "cellType")]))
cellRatio$cellNum <- rep(cellNum$Freq, times = length(unique(Seurat.Standard@meta.data$cellType)))
cellRatio$ratio <- round(cellRatio$Freq/cellRatio$cellNum,4)

GSC.EGFR_Binder <- cellRatio[which(cellRatio$orig.ident == "EGFR Binder + GSC"),]
GSC.EGFR_scFv <- cellRatio[which(cellRatio$orig.ident == "EGFR scFv + GSC"),]
Control.EGFR_Binder <- cellRatio[which(cellRatio$orig.ident == "EGFR Binder"),]
Control.EGFR_scFv <- cellRatio[which(cellRatio$orig.ident == "EGFR scFv"),]

## Control
EGFR_Binder.vs.scFv <- Control.EGFR_Binder$ratio/Control.EGFR_scFv$ratio
Control.change1 <- data.frame(type = rep("EGFR Binder", length(EGFR_Binder.vs.scFv)), group = Control.EGFR_Binder$cellType, changes = EGFR_Binder.vs.scFv)
EGFR_scFv.vs.Binder <- Control.EGFR_scFv$ratio/Control.EGFR_Binder$ratio
Control.change2 <- data.frame(type = rep("EGFR scFv", length(EGFR_scFv.vs.Binder)), group = Control.EGFR_scFv$cellType, changes = EGFR_scFv.vs.Binder)
Control.change <- rbind(Control.change1, Control.change2)
Control.change$type <- factor(Control.change$type, levels = c("EGFR scFv", "EGFR Binder"))
Control.change$group <- factor(Control.change$group, levels = rev(names(sort(table(Seurat.Standard$cellType)))))

## GSC
GSC.EGFR_Binder.vs.scFv <- GSC.EGFR_Binder$ratio/GSC.EGFR_scFv$ratio
GSC.change1 <- data.frame(type = rep("EGFR Binder", length(GSC.EGFR_Binder.vs.scFv)), group = GSC.EGFR_Binder$cellType, changes = GSC.EGFR_Binder.vs.scFv)
GSC.EGFR_scFv.vs.Binder <- GSC.EGFR_scFv$ratio/GSC.EGFR_Binder$ratio
GSC.change2 <- data.frame(type = rep("EGFR scFv", length(GSC.EGFR_scFv.vs.Binder)), group = GSC.EGFR_scFv$cellType, changes = GSC.EGFR_scFv.vs.Binder)
GSC.change <- rbind(GSC.change1, GSC.change2)
GSC.change$type <- factor(GSC.change$type, levels = c("EGFR scFv", "EGFR Binder"))
GSC.change$group <- factor(GSC.change$group, levels = rev(names(sort(table(Seurat.Standard$cellType)))))

## GSC VS control
EGFR_Bind.GSC.vs.Control <- GSC.EGFR_Binder$ratio/Control.EGFR_Binder$ratio
change1 <- data.frame(type = rep("EGFR Binder", length(EGFR_Bind.GSC.vs.Control)), group = GSC.EGFR_Binder$cellType, changes = EGFR_Bind.GSC.vs.Control)
change1 <- change1[order(change1$changes),]
EGFR_scFv.GSC.vs.Control <- GSC.EGFR_scFv$ratio/Control.EGFR_scFv$ratio
change2 <- data.frame(type = rep("EGFR scFv", length(EGFR_scFv.GSC.vs.Control)), group = GSC.EGFR_scFv$cellType, changes = EGFR_scFv.GSC.vs.Control)
GSC.vs.Control.change <- rbind(change1, change2)
GSC.vs.Control.change$type <- factor(GSC.vs.Control.change$type, levels = c("EGFR scFv", "EGFR Binder"))
GSC.vs.Control.change$changes <- log2(GSC.vs.Control.change$changes)
GSC.vs.Control.change$group <- factor(GSC.vs.Control.change$group, levels = change1$group)

library(ggpubr)
pdf("2.Cluster/AnnotatedCellType/cellRatio.change.pdf")
ggbarplot(GSC.vs.Control.change, x = "group", y = "changes", fill = "type", color = "type",
          palette = c("#214DA9", "#FF1205"), position = position_dodge(0.8), ylab = "Log2(Fold change stimulation vs unstimulation)", xlab = "")
dev.off()
