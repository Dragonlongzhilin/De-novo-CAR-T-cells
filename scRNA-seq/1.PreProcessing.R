#' @description: processing the data

library(Seurat)
library(harmony)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(clustree)
library(dplyr)
library(ggpubr)
set.seed(101)
setwd("/data/scRNAseq-CARTcell-20220501")
set.resolutions <- seq(0.1, 1.2, by = 0.1)

CART.scRNA <- Read10X(data.dir = "/data/scRNAseq-CARTcell-20220501/0.CellRange_result/filtered_feature_bc_matrix")
CART.scRNA <- CreateSeuratObject(counts = CART.scRNA, project = "CART", min.cells = 3, min.features = 200)
CART.scRNA[["mt_ratio"]] <- PercentageFeatureSet(CART.scRNA, pattern = "^MT-")
CART.scRNA[["rp_ratio"]] <- PercentageFeatureSet(CART.scRNA, pattern = "^RPL|^RPS")

# Add the patient information of the Seurat object
groups <- rownames(CART.scRNA@meta.data)
groups <- gsub(".*-1", "EGFR Binder + GSC", groups)
groups <- gsub(".*-2", "EGFR Binder", groups)
groups <- gsub(".*-3", "EGFR scFv + GSC", groups)
groups <- gsub(".*-4", "EGFR scFv", groups)
CART.scRNA <- AddMetaData(object = CART.scRNA, metadata = groups, col.name = "orig.ident")

pdf("1.QualityControl/PreFiltering.QC.pdf")
VlnPlot(
  object = CART.scRNA,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(
  object = CART.scRNA,
  features = c("mt_ratio", "rp_ratio"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# log10
VlnPlot(
  object = CART.scRNA,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & yscale("log10", .format = TRUE) & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

qc.info <- CART.scRNA@meta.data
ggdensity(qc.info, x = "nCount_RNA", title = "nCount_RNA", xlab = "") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nCount_RNA", title = "nCount_RNA", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nFeature_RNA", title = "nFeature_RNA", xlab = "") & theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nFeature_RNA", title = "nFeature_RNA", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "mt_ratio", title = "mt_ratio", xlab = "") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "mt_ratio", title = "mt_ratio", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
dev.off()

## 1.filtering low quality cell
CART.scRNA.pro <- subset(CART.scRNA, subset = nCount_RNA > 1500 & nCount_RNA < 25000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & mt_ratio < 10)

pdf("1.QualityControl/AfterFiltering.QC.pdf")
VlnPlot(
  object = CART.scRNA.pro,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(
  object = CART.scRNA.pro,
  features = c("mt_ratio", "rp_ratio"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# log10
VlnPlot(
  object = CART.scRNA.pro,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & yscale("log10", .format = TRUE) & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# check the quality metrics
qc.info <- CART.scRNA.pro@meta.data
ggdensity(qc.info, x = "nCount_RNA", title = "nCount_RNA", xlab = "") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nCount_RNA", title = "nCount_RNA", xlab = "", color = "orig.ident", fill = "orig.ident", xscale = "log10", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nFeature_RNA", title = "nFeature_RNA", xlab = "", add_density = TRUE) & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nFeature_RNA", title = "nFeature_RNA", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "mt_ratio", title = "mt_ratio", xlab = "", add_density = TRUE) & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "mt_ratio", title = "mt_ratio", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
dev.off()

## 2.filtering doublets
source(file = "/data/doubletDetect.R")
doubletRate <- read.table("/data/activate_data/longzhilin/Analysis_code/SingleCell/doubletDetect.doubletRate.txt", header = T, sep = "\t", stringsAsFactors = F)
scRNA.QC <- function(sample_id, seurat_obj, PCs = 30, doubletRate = doubletRate){

    cat(sample_id, ":remove doublets ...\n")
    seurat_obj <- subset(seurat_obj, subset = orig.ident == sample_id) # single patient
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = c("nCount_RNA", "mt_ratio"), verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = PCs, verbose = FALSE)
    seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:PCs, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:PCs, verbose = FALSE)
    p <- DimPlot(seurat_obj, group.by = "seurat_clusters")
    print(p)

    # method1: DoubleFinder
    cell.ranges <- round(nrow(seurat_obj@meta.data)/1000, 0)*1000
    doublet.rate <- doubletRate[which(doubletRate[, 2] == cell.ranges), 3]/100
    seurat_obj <- doubletDetect(Seurat.object = seurat_obj, PCs = 1:PCs, doublet.rate = doublet.rate, sct = T)
    p <- DimPlot(seurat_obj, group.by = "Doublet")
    print(p)

    return(seurat_obj)
}
# plot doublets cell distribution
pdf("1.QualityControl/DoubleFinder.pdf")
CART.scRNA.list <- lapply(as.character(unique(CART.scRNA.pro$orig.ident)), function(x){
  a <- scRNA.QC(sample_id = x, seurat_obj = CART.scRNA.pro, doubletRate = doubletRate)
  return(a)
})
names(CART.scRNA.list) <- unique(CART.scRNA.pro$orig.ident)
dev.off()

# remove doublet
pdf("1.QualityControl/QC.after.doubletRemoval.pdf")
CART.scRNA.list.pro <- lapply(CART.scRNA.list, function(seurat_obj){
  seurat_obj <- subset(seurat_obj, subset = Doublet == "Singlet")
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  p <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = T) + NoLegend()
  print(p)

  ## cluster distribution: doublets or multiples likely form distinct clusters with hybrid expression features and exhibit an aberrantly high gene count
  p <- FeaturePlot(seurat_obj, cols = c("lightgrey", "red"), features = c("nCount_RNA"))
  print(p)
  p <- FeaturePlot(seurat_obj, cols = c("lightgrey", "red"), features = c("mt_ratio"))
  print(p)
  p <- FeaturePlot(seurat_obj, cols = c("lightgrey", "red"), features = c("nFeature_RNA"))
  print(p)
  return(seurat_obj)
})
dev.off()

CART.scRNA <- merge(CART.scRNA.list.pro[[1]], y = CART.scRNA.list.pro[2:length(CART.scRNA.list.pro)], project = "CART")

pdf("1.QualityControl/PassFiltering.QC.pdf")
VlnPlot(
  object = CART.scRNA,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(
  object = CART.scRNA,
  features = c("mt_ratio", "rp_ratio"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# log10
VlnPlot(
  object = CART.scRNA,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & yscale("log10", .format = TRUE) & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# check the quality metrics
qc.info <- CART.scRNA@meta.data
ggdensity(qc.info, x = "nCount_RNA", title = "nCount_RNA", xlab = "") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nCount_RNA", title = "nCount_RNA", xlab = "", color = "orig.ident", fill = "orig.ident", xscale = "log10", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nFeature_RNA", title = "nFeature_RNA", xlab = "", add_density = TRUE) & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "nFeature_RNA", title = "nFeature_RNA", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "mt_ratio", title = "mt_ratio", xlab = "", add_density = TRUE) & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
ggdensity(qc.info, x = "mt_ratio", title = "mt_ratio", xlab = "", color = "orig.ident", fill = "orig.ident", facet.by = "orig.ident", scales = "free_y") & theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position="none")
dev.off()
saveRDS(CART.scRNA, file = "CART.scRNA.rds")

# filtering the cell without expression of CD3D and CD3E or expressed CD8A and CD4
# CD8A & CD4
T.exp <- CART.scRNA@assays$RNA@data[c("CD3D", "CD3E","CD8A", "CD4"),]
T.exp <- as.data.frame(t(as.matrix(T.exp)))
CD3D_E_label <- rep("No", nrow(T.exp))
idx <- which(T.exp$CD3D >0 | T.exp$CD3E>0)
CD3D_E_label[idx] <- "Yes"

CD8A.CD4.coexpressed <- rep("No", nrow(T.exp))
idx <- which(T.exp$CD8A >0 & T.exp$CD4>0)
CD8A.CD4.coexpressed[idx] <- "Yes"
CART.scRNA$CD3DorCD3E.exressed <- CD3D_E_label
CART.scRNA$CD8A.CD4.coexpressed <- CD8A.CD4.coexpressed

CART.scRNA <- subset(CART.scRNA,, subset = CD3DorCD3E.exressed == "Yes" & CD8A.CD4.coexpressed == "No")
saveRDS(CART.scRNA, file = "CART.scRNA.rds")

# cells with in each patient
patient.cells <- as.data.frame(table(CART.scRNA$orig.ident))
pdf("1.QualityControl/patients.cellNumber.pdf")
ggbarplot(patient.cells, x = "Var1", y = "Freq", xlab = "", ylab = "Cell number", fill = "Var1", label = T) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="none")
dev.off()

pdf("1.QualityControl/FinalFiltering.QC.pdf")
VlnPlot(
  object = CART.scRNA,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(
  object = CART.scRNA,
  features = c("mt_ratio", "rp_ratio"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# log10
VlnPlot(
  object = CART.scRNA,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & yscale("log10", .format = TRUE) & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
