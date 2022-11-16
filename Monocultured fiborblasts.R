library(tidyverse)

# Read in all data
load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Preprocessing_of_all_Mo.RData")

#Select the macrophages 
# all macrophage samples includimg sample 16 (despite poor sequencing)
fibrocountdata <- allcountdata[,-c(1,4,6, 9, 12, 17, 22, 25, 28, 31, 32, 35, 38, 
                                  41, 44, 46, 49, 52, 54, 57, 60, 61, 64, 67, 
                                  70, 73,76, 79, 81 )]
fibrosampleinfo <- allsampleinfo[-c(1,4,6, 9, 12, 17, 22, 25, 28, 31, 32, 35, 38, 
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

vstcounts_fib <- vst(fibrocountdata)
vst_pcDat_fib <- prcomp(t(vstcounts_fib))

# Initial PCA
autoplot(vst_pcDat_fib, 
         data = fibrosampleinfo, 
         shape = 'Co_vs_mono',
         colour = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)

# For this question, remove the co-cultured cells 
fibrocountdata_mono <- fibrocountdata[,-c(1,3,4,6,8,10,12,14,16,18,20,22,24,26,28,31,33,35,36,38,40,42,44,46,48,50,52)]
fibrosampleinfo_mono <- fibrosampleinfo[-c(1,3,4,6,8,10,12,14,16,18,20,22,24,26,28,31,33,35,36,38,40,42,44,46,48,50,52),]

vstcounts <- vst(fibrocountdata_mono)
vst_pcDat <- prcomp(t(vstcounts))

# Initial PCA
autoplot(vst_pcDat, 
         data = fibrosampleinfo_mono, 
         shape = 'Co_vs_mono',
         colour = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)

autoplot(vst_pcDat, 
         data = fibrosampleinfo_mono, 
         shape = 'Co_vs_mono',
         colour = "Batch",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)

# design, with batch and diagnosis
fibro_design <- as.formula(~ Batch + Diagnosis)

# what is the y intercept as standard? 
fibro_modelMatrix <- model.matrix(fibro_design, data = fibrosampleinfo_mono)
fibro_modelMatrix

# change y intercept to resolving
fibrosampleinfo_mono$Diagnosis <- factor(fibrosampleinfo_mono$Diagnosis, levels = c("Norm", "Res", "VeRA", "EstRA", "Jrep"))
fibro_modelMatrix <- model.matrix(fibro_design, data = fibrosampleinfo)
fibro_modelMatrix

# DESeq

#DESeq without filtering 
fibro_dds <- DESeqDataSetFromMatrix(countData = fibrocountdata_mono,
                                    colData = fibrosampleinfo_mono,
                                    design = fibro_design)

dds <- DESeq(fibro_dds)

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/DispersionEst mono.tiff", 
     width = 500, height = 300)
plotDispEsts(dds)
dev.off()

# To normalise the data and look at normalised counts
vsd <- vst(dds)

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, fibro_design, fibrocountdata, fibrosampleinfo, fibro_dds, dds, vsd,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Fibroblast_dds_mono.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Fibroblast_dds_mono.RData")

# Looking at the variably expressed genes ----

# Explore with PCA explorer
library("pcaExplorer")
pcaExplorer(dds = dds)

#Analysis of the most variable genes
# PCA
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Fibroblast_PCA_initial_mono.tiff", 
     width = 700, height = 400)
plotPCA(vsd, intgroup = "Diagnosis") +
  geom_text_repel(aes(x=PC1, y=PC2, label=fibro_dds$FeatureCounts), box.padding = 0.8)
dev.off()

# Dendogram
rv <- rowVars(assay(vsd))
o <- order(rv,decreasing=TRUE)
dists <- dist(t(assay(vsd)[head(o,500),]))

hc <- hclust(dists)
plot(hc, labels=vsd$SampleName)
plot(hc, labels=vsd$Diagnosis)

# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples

vst_fibro <- assay(vst(dds))

library(EnsDb.Hsapiens.v79)

topVarGenes <- head( order( rowVars(vst_fibro), decreasing=TRUE ), 100 )
subset <- vst_fibro[ topVarGenes, ]
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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Fibroblast 100 variable genes heatmap mono.tiff",
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
vst_fibro <- assay(vst(dds))

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
p <- pca(vst_fibro, metadata = my_metadata, removeVar = 0.1)

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

# Biplot of macrophages 

# Biplot of factors

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot with loadings mono.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot with labells mono.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, lab = dds$FeatureCounts, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right', lab = NULL)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot labelled mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right', lab = dds$SampleName)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Joint biplot mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Batch biplot mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Batch', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Batch', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Batch biplot PC2 and PC3 mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, "PC3", "PC2", colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Batch', legendPosition = 'right')
dev.off()

# Pairsplots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis pairsplot mono.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Batch pairsplot mono.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Batch', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Joint pairsplot mono.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=2)
dev.off()

# Correcting for batch effects -----

# Correcting for batch only
mat <- assay(vsd)
mm <- model.matrix(~Diagnosis, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "Diagnosis")
plotPCA(vsd, intgroup = "Joint")

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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Fibroblasts 100 variable genes heatmap batch only batch mono.tiff",
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
PCAtools::plotloadings(p, labSize = 5)

# Biplot of macrophages 
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot batch loadings only batch mono.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot batch loadings only batch labelled mono.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, lab = dds$FeatureCounts, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis biplot batch only batch mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right' )
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Batch biplot batch only batch mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Batch', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Batch', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Joint biplot batch only batch mono.tiff",
     width = 400, height = 300)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
dev.off()

# Pairsplots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Diagnosis pairsplot only batch mono.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Batch pairsplot only batch mono.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Batch', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Joint pairsplot only batch mono.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=2)
dev.off()

save(allcountdata, allsampleinfo, batch_vst, culture_conditions, dds, diagnosis_col, fibro_dds, fibro_modelMatrix, 
     fibrocountdata,  fibrocountdata_mono, fibrosampleinfo,  fibrosampleinfo_mono, 
      fibrosampleinfo_mono, fibrosampleinfo, fibrosampleinfo_mono, Joint_col, 
     mat, mm, my_metadata, p, subset, vsd, vst_fibro, vst_pcDat_fib,  vstcounts, vstcounts_fib, 
     dists, fibro_design, newnames, newnames_subset, o, rv, topVarGenes, col_fun, 
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Firboblast_dds_variable_mono.RData")

# Differential expression analysis

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
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_Res_mono.csv")

# export the significant results
VeRA_vs_Res_sig <- as.data.frame(annotVeRA_vs_Res)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig[order(VeRA_vs_Res_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_Res_sig_mono.csv")

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
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_Jrep_mono.csv")
# export the significant results
VeRA_vs_Jrep_sig <- as.data.frame(annotVeRA_vs_Jrep)
VeRA_vs_Jrep_sig <- VeRA_vs_Jrep_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Jrep_sig <- VeRA_vs_Jrep_sig[order(VeRA_vs_Jrep_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Jrep_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_Jrep_sig_mono.csv")

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
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_EstRA_mono.csv")
# export the significant results
VeRA_vs_EstRA_sig <- as.data.frame(annotVeRA_vs_EstRA)
VeRA_vs_EstRA_sig <- VeRA_vs_EstRA_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_EstRA_sig <- VeRA_vs_EstRA_sig[order(VeRA_vs_EstRA_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_EstRA_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_EstRA_sig_mono.csv")

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
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_Norm_mono.csv")

# export the significant results
VeRA_vs_Norm_sig <- as.data.frame(annotVeRA_vs_Norm)
VeRA_vs_Norm_sig <- VeRA_vs_Norm_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Norm_sig <- VeRA_vs_Norm_sig[order(VeRA_vs_Norm_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Norm_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/VeRA_vs_Norm_sig_mono.csv")

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
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Res_vs_Norm_mono.csv")

# export the significant results
Res_vs_Norm_sig <- as.data.frame(annotRes_vs_Norm)
Res_vs_Norm_sig <- Res_vs_Norm_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
Res_vs_Norm_sig <- Res_vs_Norm_sig[order(Res_vs_Norm_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "Res", "Norm"), alpha = 0.05)
write.csv(as.data.frame(Res_vs_Norm_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Res_vs_Norm_sig_mono.csv")

# Histogram of p values
hist(VeRA_vs_Res$pvalue)
hist(Res_vs_Norm$pvalue)
hist(VeRA_vs_Jrep$pvalue)
hist(VeRA_vs_EstRA$pvalue)
hist(VeRA_vs_Norm$pvalue)

# For number of DEGs in other comparisons
sum(VeRA_vs_Res$padj < 0.05 & (VeRA_vs_Res$log2FoldChange > 1|VeRA_vs_Res$log2FoldChange < -1), na.rm = TRUE)
sum(VeRA_vs_Norm$padj < 0.05 & (VeRA_vs_Norm$log2FoldChange > 1|VeRA_vs_Norm$log2FoldChange < -1), na.rm = TRUE)
sum(VeRA_vs_Jrep$padj < 0.05 & (VeRA_vs_Jrep$log2FoldChange > 1|VeRA_vs_Jrep$log2FoldChange < -1), na.rm = TRUE)
sum(VeRA_vs_EstRA$padj < 0.05 & (VeRA_vs_EstRA$log2FoldChange > 1|VeRA_vs_EstRA$log2FoldChange < -1), na.rm = TRUE)
sum(Res_vs_Norm$padj < 0.05 & (Res_vs_Norm$log2FoldChange > 1|Res_vs_Norm$log2FoldChange < -1), na.rm = TRUE)


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

# VeRA_vs_Res
# Shrink LF change for volcano plot
VeRA_vs_Res_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','Res'), 
                                res=VeRA_vs_Res, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
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

save(allcountdata, allsampleinfo, batch_vst, culture_conditions, dds, diagnosis_col, fibro_dds, fibro_modelMatrix, 
     fibrocountdata,  fibrocountdata_mono, fibrosampleinfo,  fibrosampleinfo_mono, 
     fibrosampleinfo_mono, fibrosampleinfo, fibrosampleinfo_mono, Joint_col, 
     mat, mm, my_metadata, p, subset, vsd, vst_fibro, vst_pcDat_fib, vstcounts, vstcounts_fib, 
     dists, fibro_design, newnames, newnames_subset, o, rv, topVarGenes, col_fun, VeRA_vs_Res, annotVeRA_vs_Res, VeRA_vs_Res_sig, 
     VeRA_vs_Jrep, VeRA_vs_Jrep_sig, annotVeRA_vs_Jrep, VeRA_vs_EstRA, VeRA_vs_EstRA_sig, annotVeRA_vs_EstRA, VeRA_vs_Norm, annotVeRA_vs_Norm,
     VeRA_vs_Norm_sig, Res_vs_Norm, Res_vs_Norm_sig, annotRes_vs_Norm, VeRA_vs_Res_shrunk, VeRA_vs_Jrep_shrunk, VeRA_vs_EstRA_shrunk,
     VeRA_vs_Norm_shrunk, Res_vs_Norm_shrunk,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Fibroblast_dds_diff_mono.RData")

load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Fibroblast_dds_diff_mono.RData")

# Do you see the same batch genes in mono-cultured cells? 

Batch_PCA <- c("ENSG00000140443",
               "ENSG00000164932",
               "ENSG00000119559",
               "ENSG00000042445", 
               "ENSG00000115649",
               "ENSG00000267534",
               "ENSG00000122863",
               "ENSG00000164548",
               "ENSG00000131475",
               "ENSG00000035687")
subset <- vst_fibro[ Batch_PCA, ]
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



# Batch from this analysis

ENSG00000105472
ENSG00000174442
ENSG00000148773
ENSG00000122203
ENSG00000204386
ENSG00000158859
ENSG00000115762
ENSG00000166016
ENSG00000151632
ENSG00000177105

Batch_PCA <- c("ENSG00000105472",
               "ENSG00000174442",
               "ENSG00000148773",
               "ENSG00000122203", 
               "ENSG00000204386",
               "ENSG00000158859",
               "ENSG00000115762",
               "ENSG00000166016",
               "ENSG00000151632",
               "ENSG00000177105")
subset <- vst_fibro[ Batch_PCA, ]
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





# Volcano plot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Volcanoes VeRA_vs_JRep_mono.tiff",
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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Volcanoes VeRA_vs_EstRA_mono.tiff",
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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Volcanoes VeRA_vs_Res_mono.tiff",
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


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Volcanoes VeRA_vs_Norm_mono.tiff",
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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/Volcanoes Res_vs_Norm_mono.tiff",
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


# For filtered analsyis and histograms-----
# Re-do filtering, but only use two samples due to EstRA dropping to 2 after outlier removal
fibro_dds_filter <- DESeqDataSetFromMatrix(countData = fibrocountdata,
                                           colData = fibrosampleinfo,
                                           design = fibro_design)

fibro_dds_filter <- estimateSizeFactors(fibro_dds_filter)

# Filter out the lower expressed genes 
# e.g. idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
# would filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5
# This is a fairly low cut off and seems to show good histograms
idx <- rowSums(counts(fibro_dds_filter, normalized=TRUE) >= 10) >= 3
fibro_dds_filter <- fibro_dds_filter[idx,]
fibro_dds_filter <- DESeq(fibro_dds_filter)

plotDispEsts(fibro_dds_filter)

filtered_VeRA_vs_Res <- results(fibro_dds_filter, contrast=c("Diagnosis", "VeRA", "Res"))
filtered_VeRA_vs_Norm <- results(fibro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"))
filtered_VeRA_vs_EstRA <- results(fibro_dds_filter, contrast=c("Diagnosis", "VeRA", "EstRA"))
filtered_VeRA_vs_Jrep <- results(fibro_dds_filter, contrast=c("Diagnosis", "VeRA", "Jrep"))
filtered_Norm_vs_Res <- results(fibro_dds_filter, contrast=c("Diagnosis", "Norm", "Res"))
hist(filtered_VeRA_vs_Res$pvalue)
hist(filtered_VeRA_vs_Norm$pvalue)
hist(filtered_VeRA_vs_EstRA$pvalue)
hist(filtered_VeRA_vs_Jrep$pvalue)
hist(filtered_Norm_vs_Res$pvalue)

# Pathway analysis ----

# To make file for GSEA
counts <- counts(dds)
write.csv(counts,file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/fibroblasts_from_macrophage_dds_counts_table_all.csv")

write.csv(fibrosampleinfo_co$Diagnosis,file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/fibroblast_from_macro_sample_info_all.csv")

# For KEGG and GO ----
library (clusterProfiler)
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

# GO analysis
go_VeRA_vs_Res <- enrichGO(gene = sigGenesVeRA_vs_Res, OrgDb = 'org.Hs.eg.db')
go_VeRA_vs_Norm <- enrichGO(gene = sigGenesVeRA_vs_Norm, OrgDb = 'org.Hs.eg.db')
go_VeRA_vs_EstRA <- enrichGO(gene = sigGenesVeRA_vs_EstRA, OrgDb = 'org.Hs.eg.db')
go_VeRA_vs_JRep <- enrichGO(gene = sigGenesVeRA_vs_Jrep, OrgDb = 'org.Hs.eg.db')
go_Res_vs_Norm <- enrichGO(gene = sigGenesRes_vs_Norm, OrgDb = 'org.Hs.eg.db')

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/GO dotplot fibroblasts.tiff",
     width = 1100, height = 700)
a <- dotplot(go_VeRA_vs_Res)
b <- dotplot(go_VeRA_vs_Norm)
c <- dotplot(go_VeRA_vs_EstRA)
d <- dotplot(go_VeRA_vs_JRep)
e <- dotplot(go_Res_vs_Norm)
cowplot::plot_grid(a, b, c, d, e,
                   labels = c("A - VeRA vs Res", "B - VeRA vs Norm", "C - VeRA VS EstRA", "D - VeRA vs JRep", "E - Res vs Norm"),
                   ncol = 3, nrow = 2)
dev.off()


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/GO dotplot fibroblasts long.tiff",
     width = 800, height = 1000)
cowplot::plot_grid(a, b, c, d, e,
                   labels = c("A - VeRA vs Res", "B - VeRA vs Norm", "C - VeRA VS EstRA", "D - VeRA vs JRep", "E - Res vs Norm"),
                   ncol = 2, nrow = 3)
dev.off()


# KEGG stuff
kk_VeRA_vs_Res <- enrichKEGG(gene = sigGenesVeRA_vs_Res, organism = 'hsa')
kk_VeRA_vs_Norm <- enrichKEGG(gene = sigGenesVeRA_vs_Norm, organism = 'hsa')
kk_VeRA_vs_JRep <- enrichKEGG(gene = sigGenesVeRA_vs_Jrep, organism = 'hsa')
kk_VeRA_vs_EstRA <- enrichKEGG(gene = sigGenesVeRA_vs_EstRA, organism = 'hsa')
kk_Res_vs_Norm <- enrichKEGG(gene = sigGenesRes_vs_Norm, organism = 'hsa')

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/KEGG dotplot fibroblasts.tiff",
     width = 800, height = 1000)
dotplot(kk_VeRA_vs_Res)
dotplot(kk_VeRA_vs_EstRA)
dotplot(kk_VeRA_vs_JRep)
dotplot(kk_VeRA_vs_Norm)
cowplot::plot_grid(a, b, c, d, e,
                   labels = c("A - VeRA vs Res", "B - VeRA vs Norm", "C - VeRA VS EstRA", "D - VeRA vs JRep", "E - Res vs Norm"),
                   ncol = 2, nrow = 3)
dev.off()


# Genes of interest ----

# From the biplot, corrected
CTHRC1 <- plotCounts(dds, gene='ENSG00000164932', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(CTHRC1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CTHRC1")
C19orf25 <- plotCounts(dds, gene='ENSG00000119559', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(C19orf25, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("C19orf25")
IGFR1 <- plotCounts(dds, gene='ENSG00000140443', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(IGFR1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IGFR1")

CNPPD1 <- plotCounts(dds, gene='ENSG00000115649', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(CNPPD1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CNPPD1")
RETSAT <- plotCounts(dds, gene='ENSG00000042445', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(RETSAT, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("RETSAT")

S1PR2 <- plotCounts(dds, gene='ENSG00000267534', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(S1PR2, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("S1PR2")
CHST3 <- plotCounts(dds, gene='ENSG00000122863', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(CHST3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CHST3")
# DON'T HAVE ADSS?? 
ADSS <- plotCounts(dds, gene='ENSG00000267312', intgroup="Diagnosis", returnData = TRUE)
ggplot(ADSS, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ADSS")

TRA2A <- plotCounts(dds, gene='ENSG00000164548', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(TRA2A, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TRA2A")
VPS25 <- plotCounts(dds, gene='ENSG00000131475', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(VPS25, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("VPS25")

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts.tiff",
     width = 700, height = 500)
cowplot::plot_grid(a, b, c, d,  e, f, g, h, i,
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                   ncol = 3, nrow = 3)
dev.off()

# From the GSEA
SMAD3 <- plotCounts(dds, gene='', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(SMAD3, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SMAD3")
ATP2B1 <- plotCounts(dds, gene='', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(ATP2B1, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ATP2B1")
PLAUR <- plotCounts(dds, gene='', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(PLAUR, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("PLAUR")

# For TGFB
XIAP <- plotCounts(dds, gene='ENSG00000101966', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(XIAP, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("XIAP")
SMAD3 <- plotCounts(dds, gene='ENSG00000166949', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(SMAD3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SMAD3")
THBS1 <- plotCounts(dds, gene='ENSG00000137801', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(THBS1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("THBS1")
FNTA <- plotCounts(dds, gene='ENSG00000168522', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(FNTA, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("FNTA")
SKIL <- plotCounts(dds, gene='ENSG00000136603', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(SKIL, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SKIL")
WWTR1 <- plotCounts(dds, gene='ENSG00000018408', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(WWTR1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("WWTR1")
ID1 <- plotCounts(dds, gene='ENSG00000125968', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(ID1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ID1")
BMPR1A <- plotCounts(dds, gene='ENSG00000107779', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(BMPR1A, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("BMPR1A")
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts TGF-B.tiff",
     width = 600, height = 400)
cowplot::plot_grid(a, b, c, d,  e, f, g, h, 
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                   ncol = 2, nrow = 4)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts TGF-B side.tiff",
     width = 600, height = 350)
cowplot::plot_grid(a, b, c, d,  e, f, 
                   #labels = c("A", "B", "C", "D", "E", "F"),
                   ncol = 3, nrow = 2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts TGF-B below.tiff",
     width = 1000, height = 150)
cowplot::plot_grid(a, b, c, d,  e, f, 
                   #labels = c("A", "B", "C", "D", "E", "F"),
                   ncol = 6, nrow = 1)
dev.off()

# For TNF
ATP2B1 <- plotCounts(dds, gene='ENSG00000070961', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(ATP2B1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ATP2B1")
CD83 <- plotCounts(dds, gene='ENSG00000112149', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(CD83, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CD83")
SMAD3 <- plotCounts(dds, gene='ENSG00000166949', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(SMAD3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SMAD3")
KLF9 <- plotCounts(dds, gene='ENSG00000119138', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(KLF9, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("KLF9")
PLAUR <- plotCounts(dds, gene='ENSG00000011422', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(PLAUR, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("PLAUR")
PNRC1 <- plotCounts(dds, gene='ENSG00000146278', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(PNRC1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("PNRC1")
NFKB1 <- plotCounts(dds, gene='ENSG00000109320', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(NFKB1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("NFKB1")
NINJ1 <- plotCounts(dds, gene='ENSG00000131669', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(NINJ1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("NINJ1")
CXCL6 <- plotCounts(dds, gene='ENSG00000124875', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(CXCL6, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CXCL6")
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts TNF.tiff",
     width = 600, height = 400)
cowplot::plot_grid(a, b, c, d,  e, f, g, h, i,
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                   ncol = 3, nrow = 3)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts TNF side.tiff",
     width = 400, height = 800)
cowplot::plot_grid(a, b, c, d,  e, f, g, h, i,
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                   ncol = 2, nrow = 5)
dev.off()

# For IFN-A
IRF7 <- plotCounts(dds, gene='ENSG00000185507', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(IRF7, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IRF7")
IFI27 <- plotCounts(dds, gene='ENSG00000165949', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(IFI27, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IFI27")
MOV10 <- plotCounts(dds, gene='ENSG00000155363', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(MOV10, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MOV10")
HERC6 <- plotCounts(dds, gene='ENSG00000138642', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(HERC6, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HERC6")
CSF1 <- plotCounts(dds, gene='ENSG00000184371', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(CSF1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CSF1")
IFIT3 <- plotCounts(dds, gene='ENSG00000119917', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(IFIT3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IFIT3")
MVB12A <- plotCounts(dds, gene='ENSG00000141971', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(MVB12A, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MVB12A")
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts IFN.tiff",
     width = 600, height = 400)
cowplot::plot_grid(a, b, c, d,  e, f, g, 
                   #labels = c("A", "B", "C", "D", "E", "F", "G"),
                   ncol = 3, nrow = 3)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts IFN side.tiff",
     width = 400, height = 550)
cowplot::plot_grid(a, b, c, d,  e, f, g, 
                   #labels = c("A", "B", "C", "D", "E", "F", "G"),
                   ncol = 2, nrow = 4)
dev.off()


# For G2M
SLC7A1 <- plotCounts(dds, gene='ENSG00000139514', intgroup="Diagnosis", returnData = TRUE)
a <- ggplot(SLC7A1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SLC7A1")
PML <- plotCounts(dds, gene='ENSG00000140464', intgroup="Diagnosis", returnData = TRUE)
b <- ggplot(PML, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("PML")
SMAD3 <- plotCounts(dds, gene='ENSG00000166949', intgroup="Diagnosis", returnData = TRUE)
c <- ggplot(SMAD3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SMAD3")
HMGA1 <- plotCounts(dds, gene='ENSG00000137309', intgroup="Diagnosis", returnData = TRUE)
d <- ggplot(HMGA1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("HMGA1")
CCND1<- plotCounts(dds, gene='ENSG00000110092', intgroup="Diagnosis", returnData = TRUE)
e <- ggplot(CCND1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CCND1")
CDK1 <- plotCounts(dds, gene='ENSG00000170312', intgroup="Diagnosis", returnData = TRUE)
f <- ggplot(CDK1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CDK1")
TACC3 <- plotCounts(dds, gene='ENSG00000013810', intgroup="Diagnosis", returnData = TRUE)
g <- ggplot(TACC3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TACC3")
EZH2 <- plotCounts(dds, gene='ENSG00000106462', intgroup="Diagnosis", returnData = TRUE)
h <- ggplot(EZH2, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("EZH2")
CENPF <- plotCounts(dds, gene='ENSG00000117724', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(CENPF, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CENPF")
NASP <- plotCounts(dds, gene='ENSG00000132780', intgroup="Diagnosis", returnData = TRUE)
j <- ggplot(NASP, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("NASP")
ODF2<- plotCounts(dds, gene='ENSG00000136811', intgroup="Diagnosis", returnData = TRUE)
k <- ggplot(ODF2, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ODF2")
ODC1<- plotCounts(dds, gene='ENSG00000115758', intgroup="Diagnosis", returnData = TRUE)
l <- ggplot(ODC1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ODC1")
SFPQ<- plotCounts(dds, gene='ENSG00000116560', intgroup="Diagnosis", returnData = TRUE)
i <- ggplot(SFPQ, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SFPQ")
NUP50<- plotCounts(dds, gene='ENSG00000093000', intgroup="Diagnosis", returnData = TRUE)
m <- ggplot(NUP50, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("NUP50")
CDC6<- plotCounts(dds, gene='ENSG00000094804', intgroup="Diagnosis", returnData = TRUE)
n <- ggplot(CDC6, aes(x=Diagnosis, y=count, colour = fibrosampleinfo_co$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("CDC6")
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts G2M.tiff",
     width = 750, height = 500)
cowplot::plot_grid(a, b, c, d,  e, f, g, h, i, j, k,l, m, n,
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N"),
                   ncol = 4, nrow = 4)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/Fibroblasts/PCA plotcounts G2M side.tiff",
     width = 400, height = 1000)
cowplot::plot_grid(a, b, c, d,  e, f, g, h, i, j, k,l, m, n,
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N"),
                   ncol = 2, nrow = 7)
dev.off()

# From the volcano plots
IRF7 <- plotCounts(dds, gene='ENSG00000185507', intgroup="Diagnosis", returnData = TRUE)
ggplot(IRF7, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IRF7")
MAST4 <- plotCounts(dds, gene='ENSG00000069020', intgroup="Diagnosis", returnData = TRUE)
ggplot(MAST4, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MAST4")
MACROD1 <- plotCounts(dds, gene='ENSG00000133315', intgroup="Diagnosis", returnData = TRUE)
ggplot(MACROD1, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MACROD1")
ZNF200 <- plotCounts(dds, gene='ENSG00000010539', intgroup="Diagnosis", returnData = TRUE)
ggplot(ZNF200, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ZNF200")
STIM2 <- plotCounts(dds, gene='ENSG00000109689', intgroup="Diagnosis", returnData = TRUE)
ggplot(STIM2, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("STIM2")
MASP1 <- plotCounts(dds, gene='ENSG00000127241', intgroup="Diagnosis", returnData = TRUE)
ggplot(MASP1, aes(x=Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("MASP1")
TNFAIP6 <- plotCounts(dds, gene='ENSG00000123610', intgroup="Diagnosis", returnData = TRUE)
ggplot(TNFAIP6, aes(x=dds$Diagnosis, y=count)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TNFAIP6")
