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
setwd("/data/ExtraDisk/sdf/longzhilin/ExtraWork/xiezhen/scRNAseq-CARTcell-20220501") # 91 server
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
CART.scRNA.list.Standard < lapply(CART.scRNA.list, function(x){
	DefaultAssay(x) <- "RNA"
	x <- NormalizeData(x, verbose = FALSE)
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
	return(x)
})

source(file = "/home/longzhilin/Analysis_Code/SingleCell/scRNA.Integrate.multipleSample.R")
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

#### phenotype marker expression
DefaultAssay(Seurat.Standard) <- "RNA"
library(ComplexHeatmap)
library(circlize)
library(ggsci)
genes <- read.table("2.Cluster/AnnotatedCellType/markers.txt", header = T, stringsAsFactors = F, sep = "\t")
genes.exp <- AverageExpression(Seurat.Standard, features = genes$Gene, group.by = "cellType")
genes.exp <- genes.exp$RNA
genes.exp <- scale(t(genes.exp))

idx <- match(colnames(genes.exp), genes$Gene)
row_split <- genes$Type[idx]
genes.exp <- genes.exp[as.character(change1$group),]
pdf("2.Cluster/AnnotatedCellType/celltype.phenotypeMarker.Heatmap.pdf")
p <- Heatmap(t(genes.exp), cluster_rows = T, cluster_columns = F, show_row_dend = F, row_split = row_split, name = "Z-score",
        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), height = unit(14, "cm"), width = unit(8, "cm"))
print(p)
dev.off()

#### load RNA-seq result
library(readxl)
source("/home/longzhilin/Analysis_Code/Combined.P.FC.R")
GSC_EGFR_Binder_scFv.up <- read_xlsx("/data/ExtraDisk/sdd/longzhilin/Project/ExtraWork/xiezhen/RNAseq-CARTcell-20220501/Differential.result.xlsx", sheet = "GSC_EGFR_Binder_scFv.up")
GSC_EGFR_Binder_scFv.up <- as.data.frame(GSC_EGFR_Binder_scFv.up)
pi <- Combined.P.FC(GSC_EGFR_Binder_scFv.up[,c(2,6)], log2FC = T, log10P = F)
GSC_EGFR_Binder_scFv.up$pi <- pi$pi
GSC_EGFR_Binder_scFv.up <- arrange(GSC_EGFR_Binder_scFv.up, desc(pi))
GSC_EGFR_Binder_scFv.down <- read_xlsx("/data/ExtraDisk/sdd/longzhilin/Project/ExtraWork/xiezhen/RNAseq-CARTcell-20220501/Differential.result.xlsx", sheet = "GSC_EGFR_Binder_scFv.down")
GSC_EGFR_Binder_scFv.down <- as.data.frame(GSC_EGFR_Binder_scFv.down)
pi <- Combined.P.FC(GSC_EGFR_Binder_scFv.down[,c(2,6)], log2FC = T, log10P = F)
GSC_EGFR_Binder_scFv.down$pi <- pi$pi
GSC_EGFR_Binder_scFv.down <- arrange(GSC_EGFR_Binder_scFv.down, desc(pi))
EGFR_Binder_scFv.up <- read_xlsx("/data/ExtraDisk/sdd/longzhilin/Project/ExtraWork/xiezhen/RNAseq-CARTcell-20220501/Differential.result.xlsx", sheet = "EGFR_Binder_scFv.up")
EGFR_Binder_scFv.up <- as.data.frame(EGFR_Binder_scFv.up)
pi <- Combined.P.FC(EGFR_Binder_scFv.up[,c(2,6)], log2FC = T, log10P = F)
EGFR_Binder_scFv.up$pi <- pi$pi
EGFR_Binder_scFv.up <- arrange(EGFR_Binder_scFv.up, desc(pi))
EGFR_Binder_scFv.down <- read_xlsx("/data/ExtraDisk/sdd/longzhilin/Project/ExtraWork/xiezhen/RNAseq-CARTcell-20220501/Differential.result.xlsx", sheet = "EGFR_Binder_scFv.down")
EGFR_Binder_scFv.down <- as.data.frame(EGFR_Binder_scFv.down)
pi <- Combined.P.FC(EGFR_Binder_scFv.down[,c(2,6)], log2FC = T, log10P = F)
EGFR_Binder_scFv.down$pi <- pi$pi
EGFR_Binder_scFv.down <- arrange(EGFR_Binder_scFv.down, desc(pi))

library(ggheatmap)
#### cluster C5
cluster4 <- cluster.sig.markers[which(cluster.sig.markers$cluster=="4"),]
cluster4 <- arrange(cluster4, desc(avg_log2FC))
Control.up.cluster4 <- intersect(cluster4$gene, EGFR_Binder_scFv.up$gene)
GSC.up.cluster4 <- intersect(cluster4$gene, GSC_EGFR_Binder_scFv.up$gene)
# top30
cluster4.top30 <- cluster4[match(GSC.up.cluster4[1:30], cluster4$gene),]
GSC_EGFR_Binder_scFv.up.top30 <- GSC_EGFR_Binder_scFv.up[match(GSC.up.cluster4[1:30], GSC_EGFR_Binder_scFv.up$gene),]
rownames(GSC_EGFR_Binder_scFv.up.top30) <- GSC_EGFR_Binder_scFv.up.top30$gene
rownames(cluster4.top30) <- cluster4.top30$gene
GSC_EGFR_Binder_scFv.up.top30 <- GSC_EGFR_Binder_scFv.up.top30[,2,drop = F]
cluster4.top30 <- cluster4.top30[,2,drop = F]
cluster4.p1 <- ggheatmap(GSC_EGFR_Binder_scFv.up.top30, color = colorRampPalette(c(rgb(255/255,160/255,122/255), "red"))(50), legendName = "log2(FC)", levels_rows = rev(rownames(GSC_EGFR_Binder_scFv.up.top30))) + ggtitle(label = "RNA-seq(GSC.UP)\nEGFR_Binder vs scFv") + theme(axis.text.y = element_text(size = 8))
cluster4.p2 <- ggheatmap(cluster4.top30, color = colorRampPalette(c(rgb(255/255,160/255,122/255), "red"))(50), legendName = "log2(FC)", levels_rows = rev(rownames(cluster4.top30))) + ggtitle(label = "scRNA-seq\nC4") + theme(axis.text.y = element_text(size = 8))

#### cluster C6
cluster5 <- cluster.sig.markers[which(cluster.sig.markers$cluster=="5"),]
cluster5 <- arrange(cluster5, desc(avg_log2FC))
Control.up.cluster5 <- intersect(cluster5$gene, EGFR_Binder_scFv.up$gene)
GSC.up.cluster5 <- intersect(cluster5$gene, GSC_EGFR_Binder_scFv.up$gene)
# top30
cluster5.top30 <- cluster5[match(GSC.up.cluster5[1:30], cluster5$gene),]
GSC_EGFR_Binder_scFv.up.top30 <- GSC_EGFR_Binder_scFv.up[match(GSC.up.cluster5[1:30], GSC_EGFR_Binder_scFv.up$gene),]
rownames(GSC_EGFR_Binder_scFv.up.top30) <- GSC_EGFR_Binder_scFv.up.top30$gene
rownames(cluster5.top30) <- cluster5.top30$gene
GSC_EGFR_Binder_scFv.up.top30 <- GSC_EGFR_Binder_scFv.up.top30[,2,drop = F]
cluster5.top30 <- cluster5.top30[,2,drop = F]
cluster5.p1 <- ggheatmap(GSC_EGFR_Binder_scFv.up.top30, color = colorRampPalette(c(rgb(255/255,160/255,122/255), "red"))(50), legendName = "log2(FC)", levels_rows = rev(rownames(GSC_EGFR_Binder_scFv.up.top30))) + ggtitle(label = "RNA-seq(GSC.UP)\nEGFR_Binder vs scFv") + theme(axis.text.y = element_text(size = 8))
cluster5.p2 <- ggheatmap(cluster5.top30, color = colorRampPalette(c(rgb(255/255,160/255,122/255), "red"))(50), legendName = "log2(FC)", levels_rows = rev(rownames(cluster5.top30))) + ggtitle(label = "scRNA-seq\nC5") + theme(axis.text.y = element_text(size = 8))

pdf("2.Cluster/AnnotatedCellType/Combined.RNAseq.pdf")
ggarrange(cluster4.p1, cluster4.p2, ncol = 3, nrow =2)
ggarrange(cluster5.p1, cluster5.p2, ncol = 3, nrow =2)
dev.off()

#### T state signature
DefaultAssay(Seurat.Standard) <- "RNA"
T.state1 <- read.table("2.Cluster/AnnotatedCellType/markers.txt", header = T, stringsAsFactors = F, sep = "\t")
T.state.list1 <- lapply(unique(T.state1$Type), function(x){
  a <- which(T.state1$Type == x)
  return(T.state1$Gene[a])
})
names(T.state.list1) <- unique(T.state1$Type)

T.state2 <- read.table("2.Cluster/AnnotatedCellType/T.signature.txt", header = T, sep = "\t", stringsAsFactors = F)
T.state.list2 <- lapply(unique(T.state2$Type), function(x){
  a <- which(T.state2$Type == x)
  return(T.state2$Gene[a])
})
names(T.state.list2) <- unique(T.state2$Type)

# merge T state
T.state.list <- c(T.state.list1, T.state.list2[-3])

# Method1: Signature score
SignatureScore <- AddModuleScore(Seurat.Standard, features = T.state.list[c(2,4,8,9,10)])
T.state.score <- SignatureScore@meta.data[,37:42]
colnames(T.state.score) <- c("CellType", names(T.state.list)[c(2,4,8,9,10)])
T.state.score.avg <- apply(T.state.score[,-1], 2, function(x){
    avg.score <- tapply(x, T.state.score$CellType, mean)
    return(avg.score)
})

# ggradar plot
library(ggsci)
library(ggradar)
library(ggsci)
library(dplyr)
library(scales)
library(tibble)

cols <- pal_npg("nrc")(3)
cols <- gsub("FF$", "", cols)
names(cols) <- c("C6", "C5", "C8")

# Normalization:0-1 for each clusters
T.state.score.avg.norm <- apply(T.state.score.avg, 1, function(x){
    return((x-min(x))/(max(x)-min(x)))
})
T.state.score.avg.norm <- t(T.state.score.avg.norm)
T.state.score.avg.norm <- cbind(data.frame(group = rownames(T.state.score.avg.norm)), T.state.score.avg.norm)
T.state.score.avg.norm$group <- factor(T.state.score.avg.norm$group, levels = c("C5", "C6", "C8"))
T.state.score.avg.norm <- T.state.score.avg.norm[,c("group", "Naive marker", "Activation markers", "Cytotoxic", "Exhaustion", "Cell cycle")]
colnames(T.state.score.avg.norm) <- c("group", "Naive", "Activation", "Cytotoxic", "Exhaustion", "Cell proliferation")
ggradar(T.state.score.avg.norm, values.radar = c("0", "0.5", "1"), gridline.mid.colour = "grey", background.circle.colour = "white", grid.mid = 0.5, grid.max = 1, group.point.size = 3, group.line.width = 0.75, group.colours = cols)
