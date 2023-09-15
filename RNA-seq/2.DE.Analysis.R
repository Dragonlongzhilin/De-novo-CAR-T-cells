library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(DESeq2)
library(writexl)
setwd("/data/sde/longzhilin/ExtraWork/activate_data/Xiazhen/RNAseq-CARTcell-20220501/5.Analysis") # 231 server

quantFiles <- list.files("/data/sde/longzhilin/ExtraWork/activate_data/Xiazhen/RNAseq-CARTcell-20220501/4.salmon_quant", 
                         pattern=paste0("^quant.genes.sf$"), full.names=TRUE, recursive=T)
quantFiles <- quantFiles[1:8]
groups <- c("Stimulated-EGFR_Binder-rep1", "Stimulated-EGFR_Binder-rep2", "Unstimulated-EGFR_Binder-rep1", "Unstimulated-EGFR_Binder-rep2",
            "Stimulated-EGFR_scFv-rep1", "Stimulated-EGFR_scFv-rep2", "Unstimulated-EGFR_scFv-rep1", "Unstimulated-EGFR_scFv-rep2")

#####################DESeq analysis########################
gene.matrix <- lapply(quantFiles, function(x){
    sf <- read.table(x, stringsAsFactors = F, sep = "\t", header = T)
    sf$Name <- gsub("\\.\\d+", "", sf$Name)
    return(sf[,c("Name", "NumReads")])
})
gene.matrix <- Reduce(function(x,y) merge(x,y,by="Name"), gene.matrix)
rownames(gene.matrix) <- gene.matrix$Name
gene.matrix <- gene.matrix[,-1]
colnames(gene.matrix) <- groups
gene.matrix <- gene.matrix[which(rowSums(gene.matrix)>0),]

# ID convert
source(file = "/data/activate_data/longzhilin/Analysis_code/IDConvert.R")
genes <- rownames(gene.matrix)
genes.convert <- IDConvert(genes = genes, method = "clusterProfiler", fromType = "ENSEMBL", toType = c("SYMBOL"))
index <- match(genes.convert$ENSEMBL, rownames(gene.matrix))
gene.matrix.pro <- gene.matrix[index,]
rownames(gene.matrix.pro) <- genes.convert$SYMBOL #18835
saveRDS(gene.matrix.pro, file = "gene.matrix.pro.rds")

#### DEseq2 analysis
gene.matrix.pro <- round(gene.matrix.pro)
condition <- factor(c("Stimulated-EGFR_Binder", "Stimulated-EGFR_Binder", "Unstimulated-EGFR_Binder", "Unstimulated-EGFR_Binder",
                      "Stimulated-EGFR_scFv", "Stimulated-EGFR_scFv", "Unstimulated-EGFR_scFv", "Unstimulated-EGFR_scFv"))
coldata <- data.frame(row.names = colnames(gene.matrix.pro), condition)
dds <- DESeqDataSetFromMatrix(countData=gene.matrix.pro, colData=coldata, design=~condition)
dds <- DESeq(dds)
EGFR_Binder <- results(dds, contrast=c("condition","Stimulated-EGFR_Binder","Unstimulated-EGFR_Binder"))
EGFR_scFv <- results(dds, contrast=c("condition","Stimulated-EGFR_scFv","Unstimulated-EGFR_scFv"))
EGFR_Binder_scFv <- results(dds, contrast=c("condition","Unstimulated-EGFR_Binder","Unstimulated-EGFR_scFv"))
GSC_EGFR_Binder_scFv <- results(dds, contrast=c("condition","Stimulated-EGFR_Binder","Stimulated-EGFR_scFv"))

DESeq2.res <- list(EGFR_Binder = as.data.frame(EGFR_Binder), EGFR_scFv = as.data.frame(EGFR_scFv), 
                   EGFR_Binder_scFv = as.data.frame(EGFR_Binder_scFv), GSC_EGFR_Binder_scFv = as.data.frame(GSC_EGFR_Binder_scFv))
DESeq2.res$EGFR_Binder$gene <- rownames(DESeq2.res$EGFR_Binder)
DESeq2.res$EGFR_scFv$gene <- rownames(DESeq2.res$EGFR_scFv)
DESeq2.res$EGFR_Binder_scFv$gene <- rownames(DESeq2.res$EGFR_Binder_scFv)
DESeq2.res$GSC_EGFR_Binder_scFv$gene <- rownames(DESeq2.res$GSC_EGFR_Binder_scFv)
saveRDS(DESeq2.res, file = "DESeq2.result.rds")
write_xlsx(DESeq2.res, "DESeq2.result.xlsx")

# construct the GSEA input files
EGFR_Binder.GSEA <- data.frame(gene = rownames(DESeq2.res$EGFR_Binder), log2FoldChange = DESeq2.res$EGFR_Binder$log2FoldChange)
EGFR_Binder.GSEA <- EGFR_Binder.GSEA[order(EGFR_Binder.GSEA[,2], decreasing=T),]
EGFR_scFv.GSEA <- data.frame(gene = rownames(DESeq2.res$EGFR_scFv), log2FoldChange = DESeq2.res$EGFR_scFv$log2FoldChange)
EGFR_scFv.GSEA <- EGFR_scFv.GSEA[order(EGFR_scFv.GSEA[,2], decreasing=T),]
EGFR_Binder_scFv.GSEA <- data.frame(gene = rownames(DESeq2.res$EGFR_Binder_scFv), log2FoldChange = DESeq2.res$EGFR_Binder_scFv$log2FoldChange)
EGFR_Binder_scFv.GSEA <- EGFR_Binder_scFv.GSEA[order(EGFR_Binder_scFv.GSEA[,2], decreasing=T),]
GSC_EGFR_Binder_scFv.GSEA <- data.frame(gene = rownames(DESeq2.res$GSC_EGFR_Binder_scFv), log2FoldChange = DESeq2.res$GSC_EGFR_Binder_scFv$log2FoldChange)
GSC_EGFR_Binder_scFv.GSEA <- GSC_EGFR_Binder_scFv.GSEA[order(GSC_EGFR_Binder_scFv.GSEA[,2], decreasing=T),]
write.table(EGFR_Binder.GSEA, file = "EGFR_Binder-Stimulated.vs.Unstimulated.rnk", quote = F, sep = "\t", row.names=F, col.names=F)
write.table(EGFR_scFv.GSEA, file = "EGFR_scFv-Stimulated.vs.Unstimulated.rnk", quote = F, sep = "\t", row.names=F, col.names=F)
write.table(EGFR_Binder_scFv.GSEA, file = "Unstimulated-EGFR_Binder.vs.scFv.rnk", quote = F, sep = "\t", row.names=F, col.names=F)
write.table(GSC_EGFR_Binder_scFv.GSEA, file = "Stimulated-EGFR_Binder.vs.scFv.rnk", quote = F, sep = "\t", row.names=F, col.names=F)

#### Volcano plot
interest.genes <- c("GZMA", "GZMB", "MKI67", "CDK4", "DUSP14")
DESeq2.res$GSC_EGFR_Binder_scFv$log10FDR <- -log10(DESeq2.res$GSC_EGFR_Binder_scFv$padj)
DESeq2.res$GSC_EGFR_Binder_scFv$Status <- "None"
DESeq2.res$GSC_EGFR_Binder_scFv$Status[which(DESeq2.res$GSC_EGFR_Binder_scFv$log2FoldChange > 1 & DESeq2.res$GSC_EGFR_Binder_scFv$padj < 0.05)] <- "Up"
DESeq2.res$GSC_EGFR_Binder_scFv$Status[which(DESeq2.res$GSC_EGFR_Binder_scFv$log2FoldChange < -1 & DESeq2.res$GSC_EGFR_Binder_scFv$padj < 0.05)] <- "Down"
DESeq2.res$GSC_EGFR_Binder_scFv$Status <- factor(DESeq2.res$GSC_EGFR_Binder_scFv$Status, levels = c("None", "Down", "Up"))
p1 <- ggscatter(DESeq2.res$GSC_EGFR_Binder_scFv, x = "log2FoldChange", y = "log10FDR", color = "Status", fill = "Status", size = 3, 
                palette = c("grey30", "royalblue", "red2"), xlab = "Log2(fold change)", ylab = "-Log10(FDR)", alpha = 0.5,
                label = rownames(DESeq2.res$GSC_EGFR_Binder_scFv), label.select = interest.genes, repel = T)
p1 + geom_hline(yintercept = -log10(0.05), lty = "dashed", lwd = 1) + geom_vline(xintercept = c(-1, 1), lty = "dashed", lwd = 1)

#### gene signature
# add selected pathways
library(readxl)
Tcell.pathways <- read_excel("immune.pathways.xlsx")
Tcell.pathways <- as.data.frame(Tcell.pathways)
Tcell.pathways <- apply(Tcell.pathways, 1, function(x){
  name <- x[1]
  x <- as.character(na.omit(x[-c(1,2)]))
  names(x) <- NULL
  res <- list(a = x)
  names(res) <- name
  return(res)
})
Tcell.pathways <- unlist(Tcell.pathways, recursive = F)
Tcell.pathways <- Tcell.pathways[-7]

library(GSVA)
exp.matrix <- readRDS("gene.matrix.pro.rds")
# ref: ImmuneSignatureGeneSet.T.Chung.2017.NatureComm
T.state <- read.table("T.signature.txt", header = T, sep = "\t", stringsAsFactors = F)
T.state.list <- lapply(unique(T.state$Type), function(x){
  a <- which(T.state$Type == x)
  return(T.state$Gene[a])
})
names(T.state.list) <- unique(T.state$Type)
T.state.list <- T.state.list[-2]
T.state.list <- c(Tcell.pathways, T.state.list)
#T.state.list <- T.state.list[-c(1,3,4,5,7,9,11,13)]
gsva.score <- gsva(expr = as.matrix(exp.matrix), gset.idx.list = T.state.list, method = "gsva", kcdf = "Poisson")
column_ha <- HeatmapAnnotation(Condition = c(rep("Stimulate", 2), rep("Unstimulate", 2),rep("Stimulate", 2), rep("Unstimulate", 2)),
                               col = list(Condition = c("Unstimulate" = "#619CFF", "Stimulate" = "#F8766D")))
p <- Heatmap(gsva.score, show_row_dend = F, show_column_dend = F, name = "GSVA score", top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), height = unit(5, "cm"), width = unit(6, "cm"))
print(p)
