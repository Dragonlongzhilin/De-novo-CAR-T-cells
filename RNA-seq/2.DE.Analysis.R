library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(DESeq2)
library(writexl)
setwd("/data/RNAseq-CARTcell-20220501/5.Analysis") # 231 server

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

