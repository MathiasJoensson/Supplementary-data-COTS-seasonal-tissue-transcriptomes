### This script is used to perform the differential expression analysis of the RNA-seq data from the COTS tissues. The script is based on the DESeq2 tutorial by Stephen Turner (https://gist.github.com/aleccstansell/9fafa21d31cb5143fd2eba37ce59ec63).

## RNA-seq analysis with DESeq2
## Stephen Turner, @genetics_blog

# RNA-seq data from GSE52202
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse52202. All patients with
# ALS, 4 with C9 expansion ("exp"), 4 controls without expansion ("ctl")

# Import & pre-process ----------------------------------------------------

countdata <- read.table("[input raw reads]", header=TRUE, row.names=1)



# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first 7 are female samples of radial nerve collected during summer, second 6 are samples of radial nerve collected in winter)
(condition <- factor(c(rep("RNf", 7), rep("RNw", 6))))
#adding sex as another factor for pcas
#(sex <- factor(c(rep("female", 7), rep("male", 6), rep("female", 7), rep("male", 6), rep("female", 7), rep("male", 6), rep("female", 7), rep("male", 5), rep("female", 7), rep("male", 6), rep("female", 7), rep("male", 6), rep("female", 7), rep("male", 6), rep("female", 7), rep("male", 6))))
  #season <- factor(c(rep("#e41a1c", 13), rep("#1f78b4", 7)))
#group.colors <- c(Summer = "#e41a1c", Winter = "#1f78b4")
# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)
library(colorRamps)
hmcol <- colorRampPalette(c("blue", "white", "firebrick2"))(n=100)
hmcol = blue2red(400)
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) >= 1 ]
nrow(dds)


dds <- DESeq(dds)

# Plot dispersions
pdf("qc-dispersions.pdf", 7, 5, pointsize=10)
plotDispEsts(dds, main="Dispersion plot")
dev.off()


# Regularized log transformation for clustering/heatmaps, etc You can use varianceStabilizingTransformation or rlogTransformation
rld <- varianceStabilizingTransformation(dds,blind = FALSE)
head(assay(rld))
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
#library(gplots)
#pdf("qc-heatmap-samples.pdf", w=7, h=5, pointsize=10)
#heatmap.2(as.matrix(sampleDists), key=F, trace="none",
col=colorpanel(100, "black", "white"),
ColSideColors=mycols[condition], RowSideColors=mycols[condition],
margin=c(10, 10), main="Sample Distance Matrix")
#dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
library(ggplot2)
options(ggplot2.discrete.color = list(c("red", "#e41a1c"), "blue", "#1f78b4"))
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
names = colnames(counts)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,  sep="-" )
#colnames(sampleDistMatrix) <- NULL

ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(aes(size=3), show.legend = F) +
  #scale_shape_manual(values = c(17, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 15, 15, 15, 15, 15, 15, 15)) +
  geom_text(aes(label=colnames(sampleDistMatrix)), size=3, hjust=0.5, vjust=-1, check_overlap = FALSE , show.legend = F) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  stat_ellipse(aes(mapping=NULL, data=NULL), type = "norm", level=0.95) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values = c("#ED1C24", "#F4FF00", "#C49A00", "#53B400", "#4CBD94", "#004D99", "#764DFF", "#FB61D7"))

ggsave("qc_pca_names.pdf", plot = last_plot(), device = "pdf")

#no names pca
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(aes(size=5), show.legend = T) +
  #geom_text(aes(label=colnames(sampleDistMatrix)), size=3, hjust=0.5, vjust=-1, check_overlap = FALSE , show.legend = F) +
  #scale_shape_manual(values = c(17, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 15, 15, 15, 15, 15, 15, 15)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  stat_ellipse(aes(mapping=NULL, data=NULL), type = "norm", level=0.95) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values = c("#ED1C24", "#C49A00", "#53B400", "#4CBD94", "#004D99", "#764DFF", "#FB61D7"))

ggsave("qc_pca_nonames.pdf", plot = last_plot(), device = "pdf")
  
  #no names pca
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(aes(size=3,shape=sex), show.legend = F) +
  scale_shape +
  geom_text(aes(label=colnames(sampleDistMatrix)), size=3, hjust=0.5, vjust=-1, check_overlap = FALSE , show.legend = F) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #stat_ellipse(aes(mapping=NULL, data=NULL), type = "norm", level=0.95) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype="dotted") +
    geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values = c("#ED1C24", "#F4FF00", "#C49A00", "#53B400", "#4CBD94", "#004D99", "#764DFF", "#FB61D7"))

ggsave("qc_pca_noellipse.pdf", plot = last_plot(), device = "pdf")
  
  
  #no names pca
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(aes(size=3,shape = sex), show.legend = F) +
  #scale_shape_manual(values = c(17, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 15, 15, 15, 15, 15, 15, 15)) +
  #geom_text(aes(label=colnames(sampleDistMatrix)), size=3, hjust=0.5, vjust=-1, check_overlap = FALSE , show.legend = F) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #stat_ellipse(aes(mapping=NULL, data=NULL), type = "norm", level=0.95) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype="dotted") +
    geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values = c("#ED1C24", "#F4FF00", "#C49A00", "#53B400", "#4CBD94", "#004D99", "#764DFF", "#FB61D7"))

ggsave("qc_pca_nonamesnoellipse.pdf", plot = last_plot(), device = "pdf")

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
metadata(res)$filterThreshold
pdf("independent_filtering.pdf", 7, 5, pointsize=10)
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()


## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
pdf("diffexpr-maplot.pdf", 7, 5, pointsize=10)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=0, sigthresh=0.05, main="Volcano Plot", legendpos="top" , labelsig=T, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="black", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    #with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  #legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
pdf("diffexpr-volcanoplot.pdf", 7, 5, pointsize=10)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-5, 5))
dev.off()

#Enhanced volcanoplot
library(EnhancedVolcano)
pdf("test.pdf", 7, 5, pointsize=1)
EnhancedVolcano(resdata,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Summer versus Winter",
                selectLab = c("gbr.614.1.t1", "gbr.58.209.t1", "gbr.11.124.t1", "gbr.17.96.t1", "gbr.187.17.t1","gbr.227.3.t1_gbr.227.4.t1", "gbr.183.69.t1", "gbr.4.78.t1", "gbr.38.74.t1"),
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 0.5,
                labSize = 2.0,
                labCol = "black",
                labFace = "bold",
                boxedLabels = T,
                col=c("black", "black", "black", "#53B400"),
                colAlpha = 4/5,
                drawConnectors = T,
                cutoffLineType = "blank",
                gridlines.major = F,
                gridlines.minor = F,
                #legendLabels = c("Not sig.", "Log (base 2) FC", "p-value", "p-value & Log (base 2) FC"),
                legendPosition = "top")
dev.off()

x11()
library(gplots)
cbPalette <- c("#9e0142", "#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#2D5270")
heatmap.2(sampleDists, key=TRUE, symkey = TRUE, density.info = "none", scale="none", trace="n",
          col=cbPalette,
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margins =c(10,10), main="Sample Distance Matrix")
dev.off()

library(pheatmap)
pheatmap(mat = sampleDists, 
         color = rev(hmcol),
         #color = colorRampPalette(c(brewer.pal(n=7, name = "RdYlBu")))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_colnames = F,
         show_rownames = F,
         border_color=T,
         annotation_row = as.data.frame(coldata[,"condition"], row.names=rownames(coldata)),
         filename = "heatmap3.pdf") 
