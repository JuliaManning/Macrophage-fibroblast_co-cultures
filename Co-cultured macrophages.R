library(tidyverse)

# Read in all data
load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Preprocessing_of_all_Mo.RData")

#Select the macrophages 
# all macrophage samples includimg sample 16 (despite poor sequencing)
macrocountdata <- allcountdata[,c(1,4,6, 9, 12, 17, 22, 25, 28, 31, 32, 35, 38, 
                                  41, 44, 46, 49, 52, 54, 57, 60, 61, 64, 67, 
                                  70, 73,76, 79, 81 )]
macrosampleinfo <- allsampleinfo[c(1,4,6, 9, 12, 17, 22, 25, 28, 31, 32, 35, 38, 
                                   41, 44, 46, 49, 52, 54, 57, 60, 61, 64, 67, 
                                   70, 73,76, 79, 81 ),]

# load packages
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)
library(cowplot)
library(ggplot2)

# Initial analysis with the monocultures still in ----

# Rename macrophage to mono
macrosampleinfo$Diagnosis<-
  replace(macrosampleinfo$Diagnosis,macrosampleinfo$Diagnosis=="Macrophage","Mono")

vstcounts <- vst(macrocountdata)
vst_pcDat <- prcomp(t(vstcounts))

# Initial PCA
autoplot(vst_pcDat, 
         data = macrosampleinfo, 
         shape = 'CellType',
         colour = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)

# design, with batch and diagnosis
macro_design <- as.formula(~ Batch + Diagnosis)

# change y intercept to macrophage
macrosampleinfo$Diagnosis <- factor(macrosampleinfo$Diagnosis, levels = c("Mono", "Norm", "Res", "VeRA", "EstRA", "Jrep"))
macro_modelMatrix <- model.matrix(macro_design, data = macrosampleinfo)
macro_modelMatrix

# DESeq
#DESeq without filtering 
macro_dds <- DESeqDataSetFromMatrix(countData = macrocountdata,
                                    colData = macrosampleinfo,
                                    design = macro_design)

dds <- DESeq(macro_dds)

#To normalise the data and look at normalised counts
vsd <- vst(dds)

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, macro_design, macrocountdata, macrosampleinfo, macro_dds, dds, vsd,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_with_mono_dds.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_with_mono_dds.RData")

# Looking at the variably expressed genes

#Analysis of the most variable genes
# PCA
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Macrophage_PCA_initial_with_mono.tiff", 
     width = 500, height = 400)
plotPCA(vsd, intgroup = "Diagnosis") 
dev.off()

# Analysis of just the co-cultured cells ----

# Remove the monocultured cells 
macrocountdata <- macrocountdata[,c(-10, -21, -29 )]
macrosampleinfo <- macrosampleinfo[c(-10, -21, -29),]

vstcounts <- vst(macrocountdata)
vst_pcDat <- prcomp(t(vstcounts))

# Initial PCA
autoplot(vst_pcDat, 
         data = macrosampleinfo, 
         shape = 'CellType',
         colour = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)

# design, with batch and diagnosis
macro_design <- as.formula(~ Batch + Diagnosis)

# what is the y intercept as standard? 
macro_modelMatrix <- model.matrix(macro_design, data = macrosampleinfo)
macro_modelMatrix

# change y intercept to normal
macrosampleinfo$Diagnosis <- factor(macrosampleinfo$Diagnosis, levels = c("Norm", "Res", "VeRA", "EstRA", "Jrep"))
macro_modelMatrix <- model.matrix(macro_design, data = macrosampleinfo)
macro_modelMatrix

#DESeq without filtering 
macro_dds <- DESeqDataSetFromMatrix(countData = macrocountdata,
                                    colData = macrosampleinfo,
                                    design = macro_design)

dds <- DESeq(macro_dds)

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/DispersionEst.tiff", 
     width = 500, height = 300)
plotDispEsts(dds)
dev.off()

# To normalise the data and look at normalised counts
vsd <- vst(dds)

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, macro_design, macrocountdata, macrosampleinfo, macro_dds, dds, vsd,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds.RData")

# Looking at the variably expressed genes ----

# Explore with PCA explorer
library("pcaExplorer")
pcaExplorer(dds = dds)

# For PCAs, with biplots, loadings etc. 
vst_macro <- assay(vst(dds))

# To colour the graphs accoridng to Diagnosis, Joint etc. 
my_metadata <- data.frame(row.names = colnames(dds))
my_metadata$Diagnosis <- dds$Diagnosis
my_metadata$Joint <- dds$Joint
my_metadata$Batch <- dds$Batch
my_metadata$Age <- dds$Age
my_metadata$Culture <- dds$Co_vs_mono
my_metadata$Fibroblast_ID <- dds$SampleName
my_metadata$RIN <- dds$RIN

# To remove the lower 10% of variables based on variance
p <- pca(vst_macro, metadata = my_metadata, removeVar = 0.1)

# To get gene names rather than ensembl gene names
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

PCAtools::screeplot(p)
# goes to PC25
PCAtools::plotloadings(p, labSize = 5)

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Diagnosis biplot.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Diagnosis biplot with loadings.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Joint biplot.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Batch biplot.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Batch', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Batch', legendPosition = 'right')
dev.off()

# Pairsplots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Diagnosis pairsplot.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Batch pairsplot.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Batch', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Joint pairsplot.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=2)
dev.off()

# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples
vst_macro <- assay(vst(dds))

library(EnsDb.Hsapiens.v79)
topVarGenes <- head( order( rowVars(vst_macro), decreasing=TRUE ), 100 )
subset <- vst_macro[ topVarGenes, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(subset) <- newnames_subset

diagnosis_col <- as.data.frame(dds$Diagnosis)
Joint_col <- as.data.frame(dds$Joint)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(4, 12, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Macrophages 100 variable genes heatmap.tiff",
     width = 500, height = 800)
Heatmap(subset, column_labels = dds$FeatureCounts, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, Joint = dds$Joint, 
                                           Batch = dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))
dev.off()

# Correcting for batch with limma ----
mat <- assay(vsd)
mm <- model.matrix(~Diagnosis, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "Diagnosis")

batch_vst <- assay(vsd)

topVarGenes <- head( order( rowVars(batch_vst), decreasing=TRUE ), 100 )
subset <- batch_vst[ topVarGenes, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(subset) <- newnames_subset

diagnosis_col <- as.data.frame(dds$Diagnosis)
Joint_col <- as.data.frame(dds$Joint)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(2, 12, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Macrophages 100 variable genes heatmap only batch.tiff",
     width = 500, height = 800)
Heatmap(subset, column_labels = dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, Joint = dds$Joint, 
                                           Batch = dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))
dev.off()

# For PCAs, with biplots, loadings etc. 

# To remove the lower 10% of variables based on variance
p <- pca(batch_vst, metadata = my_metadata, removeVar = 0.1)

# To get gene names rather than ensembl gene names
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

PCAtools::screeplot(p)

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Loadings plot only batch corrected.tiff",
     width = 500, height = 400)
PCAtools::plotloadings(p, labSize = 5)
dev.off()

# Biplot of macrophages 
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Diagnosis biplot only batch loadings.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Diagnosis biplot only batch.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Batch biplot batch.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Batch', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Batch', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Joint biplot batch.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
dev.off()

# Pairsplots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Diagnosis pairsplot batch.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Batch pairsplot batch.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Batch', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Joint pairsplot batch.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=2)
dev.off()

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, macro_design, macrocountdata, macrosampleinfo, macro_dds, dds, vsd,
     my_metadata, batch_vst, culture_conditions, diagnosis_col, Joint_col, mat, mm, p, subset, vst_macro, vst_pcDat, 
     dists, macro_design, newnames, newnames_subset, o, rv, topVarGenes, col_fun,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_PCAs.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_PCAs.RData")

# Then make heatmap or counts plots of the genes from biplots

# Heatmap
Batch_PCA <- c("ENSG00000140678",
"ENSG00000139505",
"ENSG00000167118",
"ENSG00000148803",
"ENSG00000112394",
"ENSG00000278311",
"ENSG00000015475",
"ENSG00000138668",
"ENSG00000079277",
"ENSG00000150779")

subset <- vst_macro[ Batch_PCA, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(subset) <- newnames_subset

diagnosis_col <- as.data.frame(dds$Diagnosis)
Joint_col <- as.data.frame(dds$Joint)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

Heatmap(subset, column_labels = dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, Joint = dds$Joint, 
                                           Batch = dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))

Initial_PCA <- c("ENSG00000140678",
"ENSG00000139505",
"ENSG00000167118",
"ENSG00000148803",
"ENSG00000108468",
"ENSG00000107929",
"ENSG00000120963",
"ENSG00000106609",
"ENSG00000130304",
"ENSG00000109790")

subset <- vst_macro[ Initial_PCA, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(subset) <- newnames_subset

diagnosis_col <- as.data.frame(dds$Diagnosis)
Joint_col <- as.data.frame(dds$Joint)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

Heatmap(subset, column_labels = dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, Joint = dds$Joint, 
                                           Batch = dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))





# Counts plots 
# Those in both 
ITGAX <- plotCounts(dds, gene='ENSG00000140678', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(ITGAX, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ITGAX") +
  ylab("Normalised count")
MTMR6 <- plotCounts(dds, gene='ENSG00000139505', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(MTMR6, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MTMR6") +
  ylab("Normalised count")
URM1 <- plotCounts(dds, gene='ENSG00000167118', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(URM1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("URM1")+
  ylab("Normalised count")
FUOM <- plotCounts(dds, gene='ENSG00000148803', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(FUOM, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("FUOM")+
  ylab("Normalised count")

# From the initial PCA
SLC16A10 <- plotCounts(dds, gene='ENSG00000112394', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(SLC16A10, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SLC16A10")+
  ylab("Normalised count")
GGNBP2  <- plotCounts(dds, gene='ENSG00000278311', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(GGNBP2 , aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("GGNBP2 ")+
  ylab("Normalised count")
BID <- plotCounts(dds, gene='ENSG00000015475', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(BID, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("BID")+
  ylab("Normalised count")
HNRNPD <- plotCounts(dds, gene='ENSG00000138668', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(HNRNPD, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HNRNPD")+
  ylab("Normalised count")
MKNK1 <- plotCounts(dds, gene='ENSG00000079277', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(MKNK1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MKNK1")+
  ylab("Normalised count")
TIMM8B <- plotCounts(dds, gene='ENSG00000150779', intgroup="Diagnosis", returnData = TRUE)
j <- ggplot(TIMM8B, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TIMM8B")+
  ylab("Normalised count")

# Counts plots from the batch corrected biplot
CBX1 <- plotCounts(dds, gene='ENSG00000108468', intgroup="Diagnosis", returnData = TRUE)
k <- ggplot(CBX1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CBX1")+
  ylab("Normalised count")
LARP4B <- plotCounts(dds, gene='ENSG00000107929', intgroup="Diagnosis", returnData = TRUE)
l <- ggplot(LARP4B, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("LARP4B")+
  ylab("Normalised count")
ZNF706 <- plotCounts(dds, gene='ENSG00000120963', intgroup="Diagnosis", returnData = TRUE)
m <- ggplot(ZNF706, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ZNF706")+
  ylab("Normalised count")
TMEM248 <- plotCounts(dds, gene='ENSG00000106609', intgroup="Diagnosis", returnData = TRUE)
n <- ggplot(TMEM248, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TMEM248")+
  ylab("Normalised count")
SLC27A1 <- plotCounts(dds, gene='ENSG00000130304', intgroup="Diagnosis", returnData = TRUE)
o <- ggplot(SLC27A1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SLC27A1")+
  ylab("Normalised count")
KLHL5 <- plotCounts(dds, gene='ENSG00000109790', intgroup="Diagnosis", returnData = TRUE)
p <- ggplot(KLHL5, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("KLHL5")+
  ylab("Normalised count")

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/PCA plotcounts.tiff",
     width = 700, height = 600)
cowplot::plot_grid(a, b, c,d,  e, f, g, h, i, j,k, l,m, n, o, p,
                   labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
                              "K", "L", "M", "N", "O", "P", "Q"),
                   ncol = 4, nrow = 4)
dev.off()


# Countsplots seem bimodal, so with labels, see if they are the same outliers
aa <- ggplot(ITGAX, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ITGAX") +geom_text(hjust=0.5, vjust=1)
bb <- ggplot(MTMR6, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MTMR6") +geom_text(hjust=0.5, vjust=1)
cc <- ggplot(URM1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("URM1") +geom_text(hjust=0.5, vjust=1)
dd <- ggplot(FUOM, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("FUOM") +geom_text(hjust=0.5, vjust=1)

ee <- ggplot(SLC16A10, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SLC16A10") +geom_text(hjust=0.5, vjust=1)
ff <- ggplot(GGNBL2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("GGNBL2") +geom_text(hjust=0.5, vjust=1)
gg <- ggplot(BID, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("BID") +geom_text(hjust=0.5, vjust=1)
hh <- ggplot(HHRNPD, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HHRNPD") +geom_text(hjust=0.5, vjust=1)
ii <- ggplot(MKNK1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MKNK1") +geom_text(hjust=0.5, vjust=1)
jj <- ggplot(TIMM8B, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TIMM8D") +geom_text(hjust=0.5, vjust=1)

kk <- ggplot(CBX1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CBX1") +geom_text(hjust=0.5, vjust=1)
ll <- ggplot(LARP4B, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("LARP4B") +geom_text(hjust=0.5, vjust=1)
mm <- ggplot(ZNF706, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ZNF706") +geom_text(hjust=0.5, vjust=1)
nn <- ggplot(TMEM248, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TMEM248") +geom_text(hjust=0.5, vjust=1)
oo <- ggplot(SLC27A1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SLC27A1") +geom_text(hjust=0.5, vjust=1)
pp <- ggplot(KLHL5, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("KLHL5") +geom_text(hjust=0.5, vjust=1)


cowplot::plot_grid(aa, bb, cc, dd,  ee, ff, gg, hh, ii,jj, kk, ll, mm, nn, oo, pp,
                   labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
                   "K","L", "M", "N", "O", "P"),
                   ncol = 4, nrow = 4)
# Don't see (m)any of the same outliers across the board

# Differential expression analysis  ----

# VeRA vs Res
VeRA_vs_Res <- results(dds, contrast=c("Diagnosis", "VeRA", "Res"))

# Sum of significant and DEGs
sum(VeRA_vs_Res$padj < 0.05, na.rm = TRUE)
sum(VeRA_vs_Res$padj < 0.05 & VeRA_vs_Res$log2FoldChange > 1, na.rm = TRUE)
sum(VeRA_vs_Res$padj < 0.05 & VeRA_vs_Res$log2FoldChange < -1, na.rm = TRUE)

# To get the entrez ID and gene symbol names
# load
library('org.Hs.eg.db')

# To see what they're called
columns(org.Hs.eg.db)

# Get ensembl names as rownames
# Although this says VeRA_vs_Res, it's the same gene list for all
ensembl <- as.vector(rownames(VeRA_vs_Res))

entrez <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL'))
symbol <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'SYMBOL', 'ENSEMBL'))

# VeRA VS Res
# Get DEGs as a data frame
annotVeRA_vs_Res <- as.data.frame(VeRA_vs_Res) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_Res <- cbind (annotVeRA_vs_Res, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_Res), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Res.csv")

# export the significant results
VeRA_vs_Res_sig <- as.data.frame(annotVeRA_vs_Res)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig[order(VeRA_vs_Res_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Res_sig.csv")

# VeRA vs Jrep 
VeRA_vs_Jrep <- results(dds, contrast=c("Diagnosis", "VeRA", "Jrep"))

sum(VeRA_vs_Jrep$padj < 0.05, na.rm = TRUE)
sum(VeRA_vs_Jrep$padj < 0.05 & VeRA_vs_Jrep$log2FoldChange > 1, na.rm = TRUE)
sum(VeRA_vs_Jrep$padj < 0.05 & VeRA_vs_Jrep$log2FoldChange < -1, na.rm = TRUE)

# Get DEGs as a data frame
annotVeRA_vs_Jrep <- as.data.frame(VeRA_vs_Jrep) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_Jrep <- cbind (annotVeRA_vs_Jrep, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_Jrep), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Jrep.csv")
# export the significant results
VeRA_vs_Jrep_sig <- as.data.frame(annotVeRA_vs_Jrep)
VeRA_vs_Jrep_sig <- VeRA_vs_Jrep_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Jrep_sig <- VeRA_vs_Jrep_sig[order(VeRA_vs_Jrep_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Jrep_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Jrep_sig.csv")

# VeRA vs EstRA 
VeRA_vs_EstRA <- results(dds, contrast=c("Diagnosis", "VeRA", "EstRA"))

sum(VeRA_vs_EstRA$padj < 0.05, na.rm = TRUE)
sum(VeRA_vs_EstRA$padj < 0.05 & VeRA_vs_EstRA$log2FoldChange > 1, na.rm = TRUE)
sum(VeRA_vs_EstRA$padj < 0.05 & VeRA_vs_EstRA$log2FoldChange < -1, na.rm = TRUE)

# Get DEGs as a data frame
annotVeRA_vs_EstRA <- as.data.frame(VeRA_vs_EstRA) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_EstRA <- cbind (annotVeRA_vs_EstRA, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_EstRA), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_EstRA.csv")
# export the significant results
VeRA_vs_EstRA_sig <- as.data.frame(annotVeRA_vs_EstRA)
VeRA_vs_EstRA_sig <- VeRA_vs_EstRA_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_EstRA_sig <- VeRA_vs_EstRA_sig[order(VeRA_vs_EstRA_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_EstRA_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_EstRA_sig.csv")

# VeRA vs Norm 
VeRA_vs_Norm <- results(dds, contrast=c("Diagnosis", "VeRA", "Norm"))
sum(VeRA_vs_Norm$padj < 0.05, na.rm = TRUE)
sum(VeRA_vs_Norm$padj < 0.05 & VeRA_vs_Norm$log2FoldChange > 1, na.rm = TRUE)
sum(VeRA_vs_Norm$padj < 0.05 & VeRA_vs_Norm$log2FoldChange < -1, na.rm = TRUE)

# Get DEGs as a data frame
annotVeRA_vs_Norm <- as.data.frame(VeRA_vs_Norm) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_Norm <- cbind (annotVeRA_vs_Norm, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_Norm), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Norm.csv")

# export the significant results
VeRA_vs_Norm_sig <- as.data.frame(annotVeRA_vs_Norm)
VeRA_vs_Norm_sig <- VeRA_vs_Norm_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Norm_sig <- VeRA_vs_Norm_sig[order(VeRA_vs_Norm_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Norm_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Norm_sig.csv")

# Res vs Norm 
Res_vs_Norm <- results(dds, contrast=c("Diagnosis", "Res", "Norm"))
sum(Res_vs_Norm$padj < 0.05, na.rm = TRUE)
sum(Res_vs_Norm$padj < 0.05 & Res_vs_Norm$log2FoldChange > 1, na.rm = TRUE)
sum(Res_vs_Norm$padj < 0.05 & Res_vs_Norm$log2FoldChange < -1, na.rm = TRUE)

# Get DEGs as a data frame
annotRes_vs_Norm <- as.data.frame(Res_vs_Norm) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotRes_vs_Norm <- cbind (annotRes_vs_Norm, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotRes_vs_Norm), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Res_vs_Norm.csv")

# export the significant results
Res_vs_Norm_sig <- as.data.frame(annotRes_vs_Norm)
Res_vs_Norm_sig <- Res_vs_Norm_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
Res_vs_Norm_sig <- Res_vs_Norm_sig[order(Res_vs_Norm_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "Res", "Norm"), alpha = 0.05)
write.csv(as.data.frame(Res_vs_Norm_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Res_vs_Norm_sig.csv")

# Histogram of p values
hist(VeRA_vs_Res$pvalue)
hist(Res_vs_Norm$pvalue)
hist(VeRA_vs_Jrep$pvalue)
hist(VeRA_vs_EstRA$pvalue)
hist(VeRA_vs_Norm$pvalue)

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, macro_design, macrocountdata, macrosampleinfo, macro_dds, dds, vsd,
     my_metadata, batch_vst, culture_conditions, diagnosis_col, Joint_col, mat, mm, p, subset, vst_macro, vst_pcDat, 
     dists, macro_design, newnames, newnames_subset, o, rv, topVarGenes, col_fun,VeRA_vs_Res, annotVeRA_vs_Res, VeRA_vs_Res_sig, 
     VeRA_vs_Jrep, VeRA_vs_Jrep_sig, annotVeRA_vs_Jrep, VeRA_vs_EstRA, VeRA_vs_EstRA_sig, annotVeRA_vs_EstRA, VeRA_vs_Norm, annotVeRA_vs_Norm,
     VeRA_vs_Norm_sig, Res_vs_Norm, Res_vs_Norm_sig, annotRes_vs_Norm,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_DESeq2_run.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_DESeq2_run.RData")

# For number of DEGs in other comparisons
Res_vs_EstRA <- results(dds, contrast=c("Diagnosis", "Res", "EstRA"))
sum(Res_vs_EstRA$padj < 0.05 & (Res_vs_EstRA$log2FoldChange > 1|Res_vs_EstRA$log2FoldChange < -1), na.rm = TRUE)
Res_vs_Jrep <- results(dds, contrast=c("Diagnosis", "Res", "Jrep"))
sum(Res_vs_Jrep$padj < 0.05 & (Res_vs_Jrep$log2FoldChange > 1|Res_vs_Jrep$log2FoldChange < -1), na.rm = TRUE)
Norm_vs_EstRA <- results(dds, contrast=c("Diagnosis", "Norm", "EstRA"))
sum(Norm_vs_EstRA$padj < 0.05 & (Norm_vs_EstRA$log2FoldChange > 1|Norm_vs_EstRA$log2FoldChange < -1), na.rm = TRUE)
Norm_vs_Jrep <- results(dds, contrast=c("Diagnosis", "Norm", "Jrep"))
sum(Norm_vs_Jrep$padj < 0.05 & (Norm_vs_Jrep$log2FoldChange > 1|Norm_vs_Jrep$log2FoldChange < -1), na.rm = TRUE)
Jrep_vs_EstRA <- results(dds, contrast=c("Diagnosis", "Jrep", "EstRA"))
sum(Jrep_vs_EstRA$padj < 0.05 & (Jrep_vs_EstRA$log2FoldChange > 1|Jrep_vs_EstRA$log2FoldChange < -1), na.rm = TRUE)

# To get volcano plots ----

# VeRA_vs_Res
# Shrink LF change for volcano plot
VeRA_vs_Res_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','Res'), 
                                res=VeRA_vs_Res, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnslDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_Res_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_Res_shrunk) <- newnames

# VeRA_vs_Jrep
VeRA_vs_Jrep_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','Jrep'), 
                                 res=VeRA_vs_Jrep, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_Jrep_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_Jrep_shrunk) <- newnames

# VeRA_vs_Jrep
VeRA_vs_EstRA_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','EstRA'), 
                                  res=VeRA_vs_EstRA, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_EstRA_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_EstRA_shrunk) <- newnames

# VeRA_vs_Norm
VeRA_vs_Norm_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','Norm'), 
                                 res=VeRA_vs_Norm, type = 'ashr')
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_Norm_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_Norm_shrunk) <- newnames

# Res_vs_Norm
Res_vs_Norm_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','Res','Norm'), 
                                res=Res_vs_Norm, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(Res_vs_Norm_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(Res_vs_Norm_shrunk) <- newnames

library(EnhancedVolcano)

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, macro_design, macrocountdata, macrosampleinfo, macro_dds, dds, vsd,
     my_metadata, batch_vst, culture_conditions, diagnosis_col, Joint_col, mat, mm, p, subset, vst_macro, vst_pcDat, 
     dists, macro_design, newnames, newnames_subset, o, rv, topVarGenes, col_fun,VeRA_vs_Res, annotVeRA_vs_Res, VeRA_vs_Res_sig, 
     VeRA_vs_Jrep, VeRA_vs_Jrep_sig, annotVeRA_vs_Jrep, VeRA_vs_EstRA, VeRA_vs_EstRA_sig, annotVeRA_vs_EstRA, VeRA_vs_Norm, annotVeRA_vs_Norm,
     VeRA_vs_Norm_sig, Res_vs_Norm, Res_vs_Norm_sig, annotRes_vs_Norm, VeRA_vs_Res_shrunk, VeRA_vs_Jrep_shrunk, VeRA_vs_EstRA_shrunk,
     VeRA_vs_Norm_shrunk, Res_vs_Norm_shrunk,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_DESeq2_shrinkage.RData")

load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_DESeq2_shrinkage.RData")

# Not sure on the best way to present this data (heatmap, volcano, or MA plot?) ----

# MA plots
DESeq2::plotMA(VeRA_vs_Res)
DESeq2::plotMA(VeRA_vs_Res_shrunk)
DESeq2::plotMA(VeRA_vs_Norm)
DESeq2::plotMA(VeRA_vs_Norm_shrunk)
DESeq2::plotMA(VeRA_vs_EstRA)
DESeq2::plotMA(VeRA_vs_EstRA_shrunk)
DESeq2::plotMA(VeRA_vs_Jrep)
DESeq2::plotMA(VeRA_vs_Jrep_shrunk)
DESeq2::plotMA(Res_vs_Norm)
DESeq2::plotMA(Res_vs_Norm_shrunk)

# Heatmaps for top differentially expresssed genes
VeRA_vs_Res_TEST <- as.data.frame(annotVeRA_vs_Res)
VeRA_vs_Res_TEST <- VeRA_vs_Res_TEST %>% dplyr::filter (padj < 0.1 & log2FoldChange > 6 | padj < 0.1 & log2FoldChange < -6)
VeRA_vs_Res_TEST <- VeRA_vs_Res_TEST[order(VeRA_vs_Res_TEST$log2FoldChange),]

VeRA_vs_Res_HM_test <- c(VeRA_vs_Res_TEST$GeneID)

VeRA_vs_Res_HM <- c(VeRA_vs_Res_sig$GeneID) # VeRA_vs_Norm_sig$GeneID, VeRA_vs_Jrep$GeneID, VeRA_vs_EstRA$GeneID, Res_vs_Norm$GeneID)
VeRA_vs_Norm_HM <- c(VeRA_vs_Norm_sig$GeneID)
VeRA_vs_Jrep_HM <- c(VeRA_vs_Jrep_sig$GeneID)
VeRA_vs_EstRA_HM <- c(VeRA_vs_EstRA_sig$GeneID)
Res_vs_Norm_HM <- c(Res_vs_Norm_sig$GeneID)

#If using >1 comparison
Heatmap_data <- Heatmap_data[!duplicated(Heatmap_data)]

VeRA_vs_Res_HM_test <- vstcounts[VeRA_vs_Res_HM_test,]

VeRA_vs_Res_HM <- vstcounts[VeRA_vs_Res_HM,]
VeRA_vs_Norm_HM <- vstcounts[VeRA_vs_Norm_HM,]
VeRA_vs_Jrep_HM <- vstcounts[VeRA_vs_Jrep_HM,]
VeRA_vs_EstRA_HM <- vstcounts[VeRA_vs_EstRA_HM,]
Res_vs_Norm_HM <- vstcounts[Res_vs_Norm_HM,]

# Select only the samples of the comparison
VeRA_vs_Res_HM_test <- VeRA_vs_Res_HM_test[,c(1,5,6,7,9,13,15,17,18,19,20,24,26)]

VeRA_vs_Res_HM <- VeRA_vs_Res_HM[,c(1,5,6,7,9,13,15,17,18,19,20,24,26)]
VeRA_vs_Norm_HM <- VeRA_vs_Norm_HM[,c(1,6,7,17,18,19,20,4,21,10,12,16)]
VeRA_vs_EstRA_HM <- VeRA_vs_EstRA_HM[,c(1,6,7,17,18,19,20,2,11)]
VeRA_vs_Jrep_HM <- VeRA_vs_Jrep_HM[,c(1,6,7,17,18,19,20,8,14,22,23,25)]
Res_vs_Norm_HM <- Res_vs_Norm_HM[,c(5,9,13,15,24,26,4,21,10,12,16)]

library(ComplexHeatmap)

# To get correct labels etc
VeRA_vs_Res_dds <- dds[,c(1,5,6,7,9,13,15,17,18,19,20,24,26)]
VeRA_vs_Norm_dds <- dds[,c(1,6,7,17,18,19,20,4,21,10,12,16)]
VeRA_vs_EstRA_dds <- dds[,c(1,6,7,17,18,19,20,2,11)]
VeRA_vs_Jrep_dds <- dds[,c(1,6,7,17,18,19,20,8,14,22,23,25)]
Res_vs_Norm_dds <- dds[,c(5,9,13,15,24,26,4,21,10,12,16)]

#dend = cluster_between_groups(Heatmap_data, Heatmap_data_dds$Diagnosis)

col_fun <- colorRamp2(breaks = seq(4, 10, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

library(EnsDb.Hsapiens.v79)
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(VeRA_vs_Res_HM_test),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(VeRA_vs_Res_HM_test) <- newnames_subset

Heatmap(VeRA_vs_Res_HM_test, column_labels = VeRA_vs_Res_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
         col = col_fun, cluster_rows = FALSE, #show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Res_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Res_HM.tiff",
width = 280, height = 400)
Heatmap(VeRA_vs_Res_HM, column_labels = VeRA_vs_Res_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Res_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Norm_HM.tiff",
     width = 275, height = 400)
Heatmap(VeRA_vs_Norm_HM, column_labels = VeRA_vs_Norm_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE,  col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Norm_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_EstRA_HM.tiff",
     width = 250, height = 400)
Heatmap(VeRA_vs_EstRA_HM, column_labels = VeRA_vs_EstRA_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_EstRA_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Jrep_HM.tiff",
     width = 275, height = 400)
Heatmap(VeRA_vs_Jrep_HM, column_labels = VeRA_vs_Jrep_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE,  col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Jrep_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Res_vs_Norm_HM.tiff",
     width = 270, height = 400)
Heatmap(Res_vs_Norm_HM, column_labels = Res_vs_Norm_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=Res_vs_Norm_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

# Without genes clustered

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Res_HM_gene_not_clutsered.tiff",
     width = 270, height = 400)
Heatmap(VeRA_vs_Res_HM, column_labels = VeRA_vs_Res_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, cluster_rows = FALSE,col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Res_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Norm_HM_gene_not_clutsered.tiff",
     width = 265, height = 400)
Heatmap(VeRA_vs_Norm_HM, column_labels = VeRA_vs_Norm_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, cluster_rows = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Norm_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_EstRA_HM_gene_not_clutsered.tiff",
     width = 235, height = 400)
Heatmap(VeRA_vs_EstRA_HM, column_labels = VeRA_vs_EstRA_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, cluster_rows = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_EstRA_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Jrep_HM_gene_not_clutsered.tiff",
     width = 265, height = 400)
Heatmap(VeRA_vs_Jrep_HM, column_labels = VeRA_vs_Jrep_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, cluster_rows = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=VeRA_vs_Jrep_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Res_vs_Norm_HM_gene_not_clutsered.tiff",
     width = 260, height = 400)
Heatmap(Res_vs_Norm_HM, column_labels = Res_vs_Norm_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, cluster_rows = FALSE, col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=Res_vs_Norm_dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))
dev.off()


Res_vs_Norm_HM - rowMeans((Res_vs_Norm_HM)) -> Res_vs_Norm_HM_meanSubtract
Res_vs_Norm_HM_meanSubtract/rowSds(as.matrix(Res_vs_Norm_HM)) ->
  heatmap_data_zscores

Heatmap(heatmap_data_zscores, column_labels = Res_vs_Norm_dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, #col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=Res_vs_Norm_dds$Diagnosis, Joint = Res_vs_Norm_dds$Joint, 
                                           Batch = Res_vs_Norm_dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))


Heatmap(scaled, column_labels = dds$SampleName, #cluster_columns = dend,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))


Heatmap(heatmap_data_zscores, column_labels = dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, #col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, Joint = dds$Joint, 
                                           Batch = dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))

Heatmap(Heatmap_data, column_labels = dds$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        show_row_names = FALSE, #col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, Joint = dds$Joint, 
                                           Batch = dds$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))


# Volcano plots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes VeRA_vs_JRep.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_Jrep_shrunk,
                lab = rownames(VeRA_vs_Jrep_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,11),
                title = 'VeRA vs Jrep',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes VeRA_vs_EstRA.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_EstRA_shrunk,
                lab = rownames(VeRA_vs_EstRA_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,21), 
                title = 'VeRA vs EstRA',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes VeRA_vs_Res.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_Res_shrunk,
                lab = rownames(VeRA_vs_Res_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'VeRA vs Res',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes VeRA_vs_Norm.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_Norm_shrunk,
                lab = rownames(VeRA_vs_Norm_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'VeRA vs Norm',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/
     Final/Macrophages/Volcanoes Res_vs_Norm.tiff",
     width = 500, height = 500)
EnhancedVolcano(Res_vs_Norm_shrunk,
                lab = rownames(Res_vs_Norm_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'Res vs Norm',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

# Or don't use the shrunk versions?? 
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes_VeRA_vs_Res_notshrunk.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_Res,
                lab = rownames(VeRA_vs_Res_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'VeRA vs Res',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes_VeRA_vs_Norm_notshrunk.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_Norm,
                lab = rownames(VeRA_vs_Norm_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'VeRA vs EstRA',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes_VeRA_vs_EstRA_notshrunk.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_EstRA,
                lab = rownames(VeRA_vs_EstRA_shrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,21),
                title = 'VeRA vs EstRA',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes_VeRA_VS_Jrep_notshrunk.tiff",
     width = 500, height = 500)
EnhancedVolcano(VeRA_vs_Jrep,
                lab = rownames(VeRA_vs_Jrep_shrunk),
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,6),
                title = 'VeRA vs Jrep',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Volcanoes Res_vs_Norm_notshrunk.tiff",
     width = 500, height = 500)
EnhancedVolcano(Res_vs_Norm,
                lab = rownames(Res_vs_Norm_shrunk),
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,8),
                title = 'Res vs Norm',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)
dev.off()


# DESeq2 for histograms ----
# DESeq with filtering 
macro_dds <- DESeqDataSetFromMatrix(countData = macrocountdata,
                                    colData = macrosampleinfo,
                                    design = macro_design)

filtered_dds <- estimateSizeFactors(macro_dds)
idx <- rowSums( counts(filtered_dds, normalized=TRUE) >= 10 ) >= 3
#This would say, e.g. filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5.

filtered_dds <- filtered_dds[idx,]
filtered_dds <- DESeq(filtered_dds)

filtered_VeRA_vs_Res <- results(filtered_dds, contrast=c("Diagnosis", "VeRA", "Res"))
filtered_VeRA_vs_Norm <- results(filtered_dds, contrast=c("Diagnosis", "VeRA", "Norm"))
filtered_VeRA_vs_EstRA <- results(filtered_dds, contrast=c("Diagnosis", "VeRA", "EstRA"))
filtered_VeRA_vs_Jrep <- results(filtered_dds, contrast=c("Diagnosis", "VeRA", "Jrep"))
filtered_Norm_vs_Res <- results(filtered_dds, contrast=c("Diagnosis", "Norm", "Res"))
hist(filtered_VeRA_vs_Res$pvalue)
hist(filtered_VeRA_vs_Norm$pvalue)
hist(filtered_VeRA_vs_EstRA$pvalue)
hist(filtered_VeRA_vs_Jrep$pvalue)
hist(filtered_Norm_vs_Res$pvalue)

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, macro_design, macrocountdata, macrosampleinfo, macro_dds, dds, vsd,
     my_metadata, batch_vst, culture_conditions, diagnosis_col, Joint_col, mat, mm, p, subset, vst_macro, vst_pcDat, 
     dists, macro_design, newnames, newnames_subset, o, rv, topVarGenes, col_fun,VeRA_vs_Res, annotVeRA_vs_Res, VeRA_vs_Res_sig, 
     VeRA_vs_Jrep, VeRA_vs_Jrep_sig, annotVeRA_vs_Jrep, VeRA_vs_EstRA, VeRA_vs_EstRA_sig, annotVeRA_vs_EstRA, VeRA_vs_Norm, annotVeRA_vs_Norm,
     VeRA_vs_Norm_sig, Res_vs_Norm, Res_vs_Norm_sig, annotRes_vs_Norm, VeRA_vs_Res_shrunk, VeRA_vs_Jrep_shrunk, VeRA_vs_EstRA_shrunk,
     VeRA_vs_Norm_shrunk, Res_vs_Norm_shrunk, filtered_Norm_vs_Res, filtered_VeRA_vs_EstRA, filtered_VeRA_vs_Norm, filtered_VeRA_vs_Res,
     filtered_VeRA_vs_Jrep,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_variable_filtered.RData")

load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_dds_variable_filtered.RData")

# To compare between the filtered and not filtered results
df_VeRA_vs_Res <- as.data.frame(VeRA_vs_Res)
df_f_VeRA_vs_Res <- as.data.frame(filtered_VeRA_vs_Res)
View(df_f_VeRA_vs_Res)
View(df_VeRA_vs_Res)

write.csv(df_f_VeRA_vs_Res, 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Res_filtered.csv")

write.csv(df_VeRA_vs_Res, 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/VeRA_vs_Res_not_filtered.csv")

EnhancedVolcano(filtered_VeRA_vs_Res,
                lab = rownames(filtered_VeRA_vs_Res),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'VeRA vs Res filtered p value',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)

EnhancedVolcano(filtered_VeRA_vs_Res,
                lab = rownames(filtered_VeRA_vs_Res),
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,13),
                title = 'VeRA vs Res filtered padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)

EnhancedVolcano(VeRA_vs_Res,
                lab = rownames(VeRA_vs_Res),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'VeRA_vs_Res pvalue',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)

EnhancedVolcano(VeRA_vs_Res,
                lab = rownames(VeRA_vs_Res),
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,13),
                title = 'VeRA_vs_Res padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)

EnhancedVolcano(VeRA_vs_Res,
                lab = rownames(VeRA_vs_Res),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,13),
                title = 'Res vs Norm',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 5)

sum(VeRA_vs_Res$padj < 0.05 & (VeRA_vs_Res$log2FoldChange > 1|VeRA_vs_Res$log2FoldChange < -1), na.rm = TRUE)
sum(filtered_VeRA_vs_Res$padj < 0.05 & (filtered_VeRA_vs_Res$log2FoldChange > 1|filtered_VeRA_vs_Res$log2FoldChange < -1), na.rm = TRUE)

# For GSEA ----

# So that any pathview stuff goes to the right location 
setwd("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages")

# Select only the significant genes, for pathway analysis using loose p values
sigGenesVeRA_vs_Res <- annotVeRA_vs_Res %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(padj < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

sigGenesVeRA_vs_Norm <- annotVeRA_vs_Norm %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(pvalue < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

sigGenesVeRA_vs_EstRA <- annotVeRA_vs_EstRA %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(pvalue < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

sigGenesVeRA_vs_Jrep <- annotVeRA_vs_Jrep %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(pvalue < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

sigGenesRes_vs_Norm <- annotRes_vs_Norm %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(pvalue < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

# load cluster profiler
library(clusterProfiler)

go_VeRA_vs_Res <- enrichGO(gene = sigGenesVeRA_vs_Res, OrgDb = 'org.Hs.eg.db', ont = "CC", 
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
dotplot(go_VeRA_vs_Res)

# If you get an issue with KEGG then use this
#install.packages("R.utils")
#R.utils::setOption("clusterProfiler.download.method","auto")

# GO analysis
go_VeRA_vs_Res <- enrichGO(gene = sigGenesVeRA_vs_Res, OrgDb = 'org.Hs.eg.db')
go_VeRA_vs_Norm <- enrichGO(gene = sigGenesVeRA_vs_Norm, OrgDb = 'org.Hs.eg.db')
go_VeRA_vs_EstRA <- enrichGO(gene = sigGenesVeRA_vs_EstRA, OrgDb = 'org.Hs.eg.db')
go_VeRA_vs_JRep <- enrichGO(gene = sigGenesVeRA_vs_Jrep, OrgDb = 'org.Hs.eg.db')
go_Res_vs_Norm <- enrichGO(gene = sigGenesRes_vs_Norm, OrgDb = 'org.Hs.eg.db')

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/GO dotplot macrophages.tiff",
     width = 800, height = 1000)
a <- dotplot(go_VeRA_vs_Res)
b <- dotplot(go_VeRA_vs_Norm)
c <- dotplot(go_VeRA_vs_EstRA)
d <- dotplot(go_VeRA_vs_JRep)
e <- dotplot(go_Res_vs_Norm)
cowplot::plot_grid(a, b, c, d, e,
                   labels = c("A - VeRA vs Res", "B - VeRA vs Norm", "C - VeRA VS EstRA", "D - VeRA vs JRep", "E - Res vs Norm"),
                   ncol = 2, nrow = 3)
dev.off()

# For some reason KEGG sometimes can't connect, use this
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

# KEGG stuff
kk_VeRA_vs_Res <- enrichKEGG(gene = sigGenesVeRA_vs_Res, organism = 'hsa')
kk_VeRA_vs_Norm <- enrichKEGG(gene = sigGenesVeRA_vs_Norm, organism = 'hsa')
kk_VeRA_vs_JRep <- enrichKEGG(gene = sigGenesVeRA_vs_Jrep, organism = 'hsa')
kk_VeRA_vs_EstRA <- enrichKEGG(gene = sigGenesVeRA_vs_EstRA, organism = 'hsa')
kk_Res_vs_Norm <- enrichKEGG(gene = sigGenesRes_vs_Norm, organism = 'hsa')

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/KEGG dotplot macrophages.tiff",
     width = 800, height = 1000)
a <- dotplot(kk_VeRA_vs_Res)
b <- dotplot(kk_Res_vs_Norm)
c <- dotplot(kk_VeRA_vs_EstRA)
d <- dotplot(kk_VeRA_vs_JRep)
e <- dotplot(kk_VeRA_vs_Norm)
cowplot::plot_grid(a, b, c, d, e,
                   labels = c("A - VeRA vs Res", "B - VeRA vs Norm", "C - VeRA VS EstRA", "D - VeRA vs JRep", "E - Res vs Norm"),
                   ncol = 2, nrow = 3)
dev.off()

# To make file for GSEA with MSigDB
counts <- counts(dds)
write.csv(counts,file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/macrophage_dds_counts_table_all.csv")

write.csv(macrosampleinfo$Diagnosis,file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/macro_sample_info_all.csv")

# Make plots of the counts 

# Counts plots for genes of interest from pathway analysis ----

BIN1 <- plotCounts(dds, gene='ENSG00000136717', intgroup="Diagnosis", returnData = TRUE)
ggplot(BIN1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("BIN1")+
  ylab("Normalised count")

IRF7 <- plotCounts(dds, gene='ENSG00000185507', intgroup="Diagnosis", returnData = TRUE)
ggplot(IRF7, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IRF7")+
  ylab("Normalised count")

ggplot(MAP2K3, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MAP2K3")+
  ylab("Normalised count")


#IL-6-JAK-STAT from hallmark
TNFRSF1A <- plotCounts(dds, gene='ENSG00000067182', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(TNFRSF1A, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TNFRSF1A")+
  ylab("Normalised count")
ggplot(TNFRSF1A, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TNFRSF1A") +geom_text(hjust=0.5, vjust=1)
MAP3K8 <- plotCounts(dds, gene='ENSG00000107968', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(MAP3K8, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MAP3K8")+
  ylab("Normalised count")
IRF9 <- plotCounts(dds, gene='ENSG00000213928', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(IRF9, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IRF9")+
  ylab("Normalised count")
ggplot(IRF9, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IRF9") +geom_text(hjust=0.5, vjust=1)
PTPN1 <- plotCounts(dds, gene='ENSG00000196396', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(PTPN1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("PTPN1")+
  ylab("Normalised count")
HMOX1 <- plotCounts(dds, gene='ENSG00000100292', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(HMOX1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HMOX1")+
  ylab("Normalised count")
IFNGR1 <- plotCounts(dds, gene='ENSG00000027697', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(IFNGR1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IFNGR1")+
  ylab("Normalised count")
CXCL1 <- plotCounts(dds, gene='ENSG00000163739', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(CXCL1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CXCL1")+
  ylab("Normalised count")
IL2RG <- plotCounts(dds, gene='ENSG00000147168', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(IL2RG, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL2RG")+
  ylab("Normalised count")
IL1B<- plotCounts(dds, gene='ENSG00000125538', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(IL1B, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL1B")+
  ylab("Normalised count")

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/IL6-JAK-STAT genes macrophages.tiff",
     width = 400, height = 650)
cowplot::plot_grid(a,b,c,d,e,f,g,h,i,
                   ncol = 2, nrow =5)
dev.off()

# TNF from hallmark
ID2 <- plotCounts(dds, gene='ENSG00000115738', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(ID2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ID2")+
  ylab("Normalised count")
EHD1 <- plotCounts(dds, gene='ENSG00000110047', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(EHD1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("EHD1")+
  ylab("Normalised count")
CXCL2 <- plotCounts(dds, gene='ENSG00000081041', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(CXCL2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CXCL2")+
  ylab("Normalised count")
TRAF1 <- plotCounts(dds, gene='ENSG00000056558', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(TRAF1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TRAF1")+
  ylab("Normalised count")
NFE2L2 <- plotCounts(dds, gene='ENSG00000116044', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(NFE2L2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("NFE2L2")+
  ylab("Normalised count")
MAP2K3 <- plotCounts(dds, gene='ENSG00000034152', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(MAP2K3, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MAP2K3")+
  ylab("Normalised count")
ggplot(MAP2K3, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MAP2K3") +geom_text(hjust=0.5, vjust=1)
IL1A <- plotCounts(dds, gene='ENSG00000115008', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(IL1A, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL1A")+
  ylab("Normalised count")
BIRC3 <- plotCounts(dds, gene='ENSG00000023445', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(BIRC3, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("BIRC3")+
  ylab("Normalised count")
CEBPD<- plotCounts(dds, gene='ENSG00000221869', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(CEBPD, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CEBPD")+
  ylab("Normalised count")
SOD2<- plotCounts(dds, gene='ENSG00000112096', intgroup="Diagnosis", returnData = TRUE)
j <- ggplot(SOD2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SOD2")+
  ylab("Normalised count")
IRS2 <- plotCounts(dds, gene='ENSG00000185950', intgroup="Diagnosis", returnData = TRUE)
K <- ggplot(IRS2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IRS2")+
  ylab("Normalised count")

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/TNF genes macrophages.tiff",
     width = 400, height = 800)
cowplot::plot_grid(a,b,c,d,e,f,g,h,i,j,K,
                   ncol = 2, nrow =6)
dev.off()



#IFN-y from hallmark
SP110 <- plotCounts(dds, gene='ENSG00000135899', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(SP110, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SP110")+
  ylab("Normalised count")
SOD2 <- plotCounts(dds, gene='ENSG00000112096', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(SOD2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SOD2")+
  ylab("Normalised count")
CASP8 <- plotCounts(dds, gene='ENSG00000196396', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(CASP8, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CASP8")+
  ylab("Normalised count")
PTPN1 <- plotCounts(dds, gene='ENSG00000196396', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(PTPN1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("PTPN1")+
  ylab("Normalised count")
MARCHF1 <- plotCounts(dds, gene='ENSG00000145416', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(MARCHF1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MARCHF1")+
  ylab("Normalised count")
SAMHD1 <- plotCounts(dds, gene='ENSG00000101347', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(SAMHD1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SAMHD1")+
  ylab("Normalised count")
STAT1 <- plotCounts(dds, gene='ENSG00000115415', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(STAT1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("STAT1")+
  ylab("Normalised count")
ggplot(STAT1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis, label = macrosampleinfo$SampleName)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("STAT1") +geom_text(hjust=0.5, vjust=1)
TRAFD1 <- plotCounts(dds, gene='ENSG00000135148', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(TRAFD1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TRAFD1")+
  ylab("Normalised count")
OAS2 <- plotCounts(dds, gene='ENSG00000111335', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(OAS2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("OAS2")+
  ylab("Normalised count")
CMKLR1 <- plotCounts(dds, gene='ENSG00000174600', intgroup="Diagnosis", returnData = TRUE)
j <- ggplot(CMKLR1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CMKLR1")+
  ylab("Normalised count")
HLAA <- plotCounts(dds, gene='ENSG00000206503', intgroup="Diagnosis", returnData = TRUE)
k <- ggplot(HLAA, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HLAA")+
  ylab("Normalised count")
HIF1A <- plotCounts(dds, gene='ENSG00000100644', intgroup="Diagnosis", returnData = TRUE)
l <- ggplot(HIF1A, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HIF1A")+
  ylab("Normalised count")
SAMD9L <- plotCounts(dds, gene='ENSG00000177409', intgroup="Diagnosis", returnData = TRUE)
m <- ggplot(SAMD9L, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SAMD9L")+
  ylab("Normalised count")
HLAB <- plotCounts(dds, gene='ENSG00000234745', intgroup="Diagnosis", returnData = TRUE)
n <- ggplot(HLAB, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HLAB")+
  ylab("Normalised count")
MVP <- plotCounts(dds, gene='ENSG00000013364', intgroup="Diagnosis", returnData = TRUE)
o <- ggplot(MVP, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MVP")+
  ylab("Normalised count")

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/IFNy genes macrophages.tiff",
     width = 550, height = 600)
cowplot::plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,
                   ncol = 3, nrow =5)
dev.off()

# Make the IFN gene heatmap
IFN_all <- "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/IFN_A_B_genes.txt"
IFN_all <-readLines(IFN_all)

IFN_all_vst <- vstcounts[IFN_all,]

library(EnsDb.Hsapiens.v79)
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(IFN_all_vst),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(IFN_all_vst) <- newnames_subset

library(ComplexHeatmap)
dend = cluster_between_groups(IFN_all_vst, dds$Diagnosis)

IFN_all_vst - rowMeans((IFN_all_vst)) -> IFN_all_vst_meanSubtract
IFN_all_vst_meanSubtract/rowSds(as.matrix(IFN_all_vst)) ->
  IFN_all_vst_zscores

Heatmap(IFN_all_vst_zscores, column_labels = dds$SampleName, #cluster_columns = dend,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"))))

# Mariola paper
MerTK <- plotCounts(dds, gene='ENSG00000153208', intgroup="Diagnosis", returnData = TRUE) # MerTK
CD206 <- plotCounts(dds, gene='ENSG00000260314', intgroup="Diagnosis", returnData = TRUE) # CD206
CD163 <- plotCounts(dds, gene='ENSG00000177575', intgroup="Diagnosis", returnData = TRUE) # CD163
CD11b <- plotCounts(dds, gene='ENSG00000169896', intgroup="Diagnosis", returnData = TRUE) # CD11b - out of interest? 
# From the ICAM1 cluster - mariola
IL1B <- plotCounts(dds, gene='ENSG00000125538', intgroup="Diagnosis", returnData = TRUE) # IL-1B
CCL4 <- plotCounts(dds, gene='ENSG00000275302', intgroup="Diagnosis", returnData = TRUE) # CCL4
ICAM1 <- plotCounts(dds, gene='ENSG00000090339', intgroup="Diagnosis", returnData = TRUE) # ICAM1
# From the ISG15 cluster - mariola
ISG15 <- plotCounts(dds, gene='ENSG00000187608', intgroup="Diagnosis", returnData = TRUE) # ISG15
GBP1 <- plotCounts(dds, gene='ENSG00000117228', intgroup="Diagnosis", returnData = TRUE) # GBP1
IFI6  <- plotCounts(dds, gene='ENSG00000126709', intgroup="Diagnosis", returnData = TRUE) # IFI6 

a <- ggplot(MerTK, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MerTK") +
  ylab("Normalised count")
b <- ggplot(CD206, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CD206")+
  ylab("Normalised count")
c <- ggplot(CD163, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CD163")+
  ylab("Normalised count")
d <- ggplot(CD11b, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CD11b")+
  ylab("Normalised count")
e <- ggplot(CCL4, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CLL4")+
  ylab("Normalised count")
f <- ggplot(ICAM1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ICAM1")+
  ylab("Normalised count")
g <- ggplot(IL1B, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL1B")+
  ylab("Normalised count")
h <- ggplot(ISG15, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ISG15")+
  ylab("Normalised count")
i <- ggplot(GBP1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("GBP1")+
  ylab("Normalised count")
j <- ggplot(IFI6, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IFI6")+
  ylab("Normalised count")


# Other genes of interest from mariola paper
STAB1 <- plotCounts(dds, gene='ENSG00000010327', intgroup="Diagnosis", returnData = TRUE) #STAB1
k <- ggplot(STAB1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("STAB1")+
  ylab("Normalised count")
LYVE1 <- plotCounts(dds, gene='ENSG00000133800', intgroup="Diagnosis", returnData = TRUE) #STAB1
l <- ggplot(LYVE1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("LYVE1")+
  ylab("Normalised count")

S100A9 <- plotCounts(dds, gene='ENSG00000163220', intgroup="Diagnosis", returnData = TRUE) 
m <- ggplot(S100A9, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("S100A9")+
  ylab("Normalised count")

TREM2 <- plotCounts(dds, gene='ENSG00000095970', intgroup="Diagnosis", returnData = TRUE) 
n <- ggplot(TREM2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TREM2")+
  ylab("Normalised count")

ID2 <- plotCounts(dds, gene='ENSG00000115738', intgroup="Diagnosis", returnData = TRUE) 
o <- ggplot(ID2, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ID2")+
  ylab("Normalised count")
LY6E <- plotCounts(dds, gene='ENSG00000160932', intgroup="Diagnosis", returnData = TRUE)
p <- ggplot(LY6E, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("LY6E")+
  ylab("Normalised count")

# IFN genes from original donlin paper 
CXCL9 <- plotCounts(dds, gene='ENSG00000138755', intgroup="Diagnosis", returnData = TRUE)
q <- ggplot(CXCL9, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CXCL9") +
  ylab("Normalised count")
CXCL10 <- plotCounts(dds, gene='ENSG00000169245', intgroup="Diagnosis", returnData = TRUE)
r <- ggplot(CXCL10, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CXCL10")+
  ylab("Normalised count")
IFIT1 <- plotCounts(dds, gene='ENSG00000185745', intgroup="Diagnosis", returnData = TRUE)
s <- ggplot(IFIT1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IFIT1")+
  ylab("Normalised count")
MX1 <- plotCounts(dds, gene='ENSG00000157601', intgroup="Diagnosis", returnData = TRUE)
t <- ggplot(MX1, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MX1")+
  ylab("Normalised count")

NKG7 <- plotCounts(dds, gene='ENSG00000105374', intgroup="Diagnosis", returnData = TRUE)
u <- ggplot(NKG7, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("NKG7")+
  ylab("Normalised count")
CXCL5 <- plotCounts(dds, gene='ENSG00000163735', intgroup="Diagnosis", returnData = TRUE)
v <- ggplot(CXCL5, aes(x=Diagnosis, y=count, colour = macrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CXCL5")+
  ylab("Normalised count")

cowplot::plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,
                   labels = c("A", "B", "C", "D", "E", 
                              "F", "G", "H", "I", "J", 
                              "K", "L", "M", "N", "O", 
                              "P", "Q", "R", "S", "T", 
                              "U", "V"),
                   ncol = 6, nrow =4)

# To build a 3D volcano ----
# For this question leave only Res, Norm and VeRA 
RAcountdata <- macrocountdata[,c(4,5,6,7,9,10,12,13,15,16,17,18,19,20,21,24,26)]
RAsampleinfo <- macrosampleinfo[c(4,5,6,7,9,10,12,13,15,16,17,18,19,20,21,24,26),]

library(DESeq2)
library(volcano3D)

#DESeq without filtering 
ra_dds <- DESeqDataSetFromMatrix(countData = RAcountdata,
                                 colData = RAsampleinfo,
                                 design = ~ Diagnosis + Batch)

ra_dds_DE <- DESeq(ra_dds)

library('org.Hs.eg.db')

# To see what they're called
columns(org.Hs.eg.db)

# Get ensembl names as rownames
# Although this says VeRA_vs_Res, it's the same gene list for all
ensembl <- as.vector(rownames(RES))

entrez <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL'))
symbol <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'SYMBOL', 'ENSEMBL'))

# VeRA VS Res
# Get DEGs as a data frame
annotVeRA_vs_Res <- as.data.frame(RES) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_Res <- cbind (annotVeRA_vs_Res, entrez, symbol)

# export the results as csv
res <- as.data.frame(annotVeRA_vs_Res)
          

RES <- results(ra_dds_DE, contrast=c("Diagnosis", "VeRA", "Res"))
res_df <- as.data.frame(RES)

# likelihood ratio test on 'Diagnosis'
dds_LRT <- DESeq(ra_dds, test = "LRT", reduced = ~ Batch, parallel = TRUE) 

# To get gene names
library(EnsDb.Hsapiens.v79)
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(ra_dds_DE),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(ra_dds_DE) <- newnames_subset

# create 'volc3d' class object for plotting
res <- deseq_polar(ra_dds_DE, dds_LRT, "Diagnosis")

# plot 3d volcano plot
volcano3D(res)

save(ra_dds, ra_dds_DE, RAcountdata, RAsampleinfo,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Volcano_plot_withoutBX027.RData")

load(ra_dds, ra_dds_DE, RAcountdata, RAsampleinfo,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Volcano_plot_withoutBX027.RData")


save(ra_dds, ra_dds_DE, RAcountdata, RAsampleinfo,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Volcano_plot.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Volcano_plot.RData")

vst_macro <- assay(vst(ra_dds_DE))
topVarGenes <- head( order( rowVars(vst_macro), decreasing=TRUE ), 100 )
subset <- vst_macro[ topVarGenes, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(subset) <- newnames_subset


diagnosis_col <- as.data.frame(dds$Diagnosis)
Joint_col <- as.data.frame(dds$Joint)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(4, 12, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Macrophages/Macrophages 100 variable genes heatmap VeRA_Norm_Res.tiff",
     width = 400, height = 800)
Heatmap(subset, column_labels = ra_dds_DE$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)),col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=ra_dds_DE$Diagnosis, Joint = ra_dds_DE$Joint, 
                                           Batch = ra_dds_DE$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))
dev.off()


subset - rowMeans((subset)) -> subset_meanSubtract
subset_meanSubtract/rowSds(as.matrix(subset)) ->
  heatmap_data_zscores

Heatmap(heatmap_data_zscores, column_labels = ra_dds_DE$FeatureCounts, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)),# col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=ra_dds_DE$Diagnosis, Joint = ra_dds_DE$Joint, 
                                           Batch = ra_dds_DE$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))


vst_macro <- vst(ra_dds_DE)
plotPCA(vst_macro, intgroup = "Batch")
plotPCA(vst_macro, intgroup = "Diagnosis")


vsd <- vst(ra_dds_DE)
mat <- assay(vsd)
mm <- model.matrix(~Diagnosis, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "Diagnosis")

batch_vst <- assay(vsd)

# To colour the graphs accoridng to Diagnosis, Joint etc. 
my_metadata <- data.frame(row.names = colnames(ra_dds_DE))
my_metadata$Diagnosis <- ra_dds_DE$Diagnosis
my_metadata$Joint <- ra_dds_DE$Joint
my_metadata$Batch <- ra_dds_DE$Batch
my_metadata$Age <- ra_dds_DE$Age
my_metadata$Culture <- ra_dds_DE$Co_vs_mono
my_metadata$Fibroblast_ID <- ra_dds_DE$SampleName
my_metadata$RIN <- ra_dds_DE$RIN

# To remove the lower 10% of variables based on variance
p <- pca(batch_vst, metadata = my_metadata, removeVar = 0.1)

# To get gene names rather than ensembl gene names
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

PCAtools::screeplot(p)
# goes to PC25
PCAtools::plotloadings(p, labSize = 5)
PCAtools::biplot((p), colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')
PCAtools::biplot((p), colby = 'Batch', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')

Batch_PCA <- c("ENSG00000158517",
"ENSG00000048740",
"ENSG00000198839",
"ENSG00000105393",
"ENSG00000164091",
"ENSG00000138777",
"ENSG00000130119",
"ENSG00000169139",
"ENSG00000129255",
"ENSG00000132153")

vst_macro <- assay(vst(ra_dds_DE))
subset <- vst_macro[Batch_PCA, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames_subset)
rownames(subset) <- newnames_subset

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

Heatmap(subset, column_labels = ra_dds_DE$SampleName, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=ra_dds_DE$Diagnosis, Joint = ra_dds_DE$Joint, 
                                           Batch = ra_dds_DE$Batch, col = list
                                           (Diagnosis = c("Macrophage" = "grey", "Res" = "chartreuse4", "VeRA" = "dodgerblue3","Jrep" = "firebrick",
                                                          "Norm" = "magenta4", "EstRA" = "darkorange"),
                                             Joint = c("Macrophage" = "grey", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Batch = c("Batch_1" = "lightseagreen", "Batch_2" = "blueviolet", "Batch_3" = "coral1"))))


