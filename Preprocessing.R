# Loading data in ----
library (tidyverse)

# Read in sample information
allsampleinfo <- read_tsv("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Macrophage-fibroblast samples/Sample_info/Macrophage-fibroblast-sample-info.txt")

# To load in counts table
seqdata <- read_tsv("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Macrophage-fibroblast samples/Analysis/Mo_count_f.txt")

# Transform countdata into matrix
allcountdata <- seqdata %>%
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  rename_all(str_remove, ".bam") %>% # remove the ".bam" from the column names
  select(allsampleinfo$FeatureCounts) %>%
  as.matrix()
# Get a lot of 0s

# load packages
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)
library(cowplot)
library(ggplot2)

# Remove samples 6, 16 and 89
allcountdata <- allcountdata[,c (-6, -16, -83)]
allsampleinfo <- allsampleinfo[c (-6, -16, -83),]

# Look at boxplot and mean vs SD
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Boxplot_all_samples.tiff", 
     width = 2000, height = 600)
boxplot(allcountdata, main='Raw counts distribution', las=2, names = allsampleinfo$SampleName)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Mean_vs_SD_all_samples.tiff", 
     width = 450, height = 450)
plot(rowMeans(allcountdata), rowSds(allcountdata), 
     main='Raw Counts : Mean Vs SD', 
     xlim=c(0,10000), 
     ylim=c(0,5000))
dev.off()

#Vectors for colour
CellCol <- match(allsampleinfo$CellType, c("Macrophage", "Fibroblast")) + 1
DiagnosisCol <- match(allsampleinfo$Diagnosis, c("Macrophage", "Norm", "Res", "VeRA", "EstRA", "JRep")) + 1

# Transform data with variance stabalisation (VST)
vstcounts <- vst(allcountdata)

# Plot the counts distribution
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/VST_boxplots.tiff", 
     width = 1000, height = 450)
boxplot(vstcounts, 
        xlab="", 
        ylab="VST counts",
        las=2,
        col = CellCol) 
dev.off()

# Plot the mean vs Sd
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/VST_Mean_vs_SD.tiff", 
     width = 450, height = 450)
plot(rowMeans(vstcounts), rowSds(vstcounts), 
     main='VST Counts : Mean Vs SD')
dev.off()

# PCA
vst_pcDat <- prcomp(t(vstcounts))

# plot PCA
# with sample names
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/PCA_all.tiff", 
     width = 600, height = 350)
autoplot(vst_pcDat, 
         data = allsampleinfo, 
         shape = 'CellType',
         colour = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)
dev.off()

save(allcountdata, allsampleinfo, vstcounts, vst_pcDat, file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Preprocessing_of_all_Mo.RData")

# Based on this and the MultiQC plot, sample 1, 6, 16 are removed (1 and 6 macrophages, 16 fibroblast)

# Macrophages only ----

# To pick up from this point 
load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Preprocessing_of_all_Mo.RData")

# all macrophage samples includimg sample 16 (despite poor sequencing)
macrocountdata <- allcountdata[,c(1,4,6, 9, 12, 17, 22, 25, 28, 31, 32, 35, 38, 
                                  41, 44, 46, 49, 52, 54, 57, 60, 61, 64, 67, 
                                  70, 73,76, 79, 81 )]
macrosampleinfo <- allsampleinfo[c(1,4,6, 9, 12, 17, 22, 25, 28, 31, 32, 35, 38, 
                                  41, 44, 46, 49, 52, 54, 57, 60, 61, 64, 67, 
                                  70, 73,76, 79, 81 ),]

# Re load packages if needed
# Quick sanity check everything is working
macro_vstcounts <- vst(macrocountdata)
macro_vst_pcDat <- prcomp(t(macro_vstcounts))

autoplot(macro_vst_pcDat, 
         data = macrosampleinfo, 
         colour = 'Diagnosis',
         main = 'Macrophage VST transform PCA',
         max.overlaps = 10,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=FeatureCounts), box.padding = 0.8)

autoplot(macro_vst_pcDat, 
         data = macrosampleinfo, 
         colour = 'Batch',
         main = 'Macrophage VST transform PCA',
         max.overlaps = 10,
         size = 3) 

autoplot(macro_vst_pcDat, 
         data = macrosampleinfo, 
         colour = 'Joint',
         main = 'Macrophage VST transform PCA',
         max.overlaps = 10,
         size = 3) 

# Then run DESeq2 ----

# design suggesting diagnosis causes differences in groups
macro_design <- as.formula(~ Diagnosis)

# what is the y intercept as standard? 
macro_modelMatrix <- model.matrix(macro_design, data = macrosampleinfo)
macro_modelMatrix

# change y intercept to resolving
macrosampleinfo$Diagnosis <- factor(macrosampleinfo$Diagnosis, levels = c("Macrophage", "Norm", "Res", "VeRA", "EstRA", "Jrep"))
macro_modelMatrix <- model.matrix(macro_design, data = macrosampleinfo)
macro_modelMatrix

# DESeq

#DESeq without filtering 
macro_dds <- DESeqDataSetFromMatrix(countData = macrocountdata,
                                    colData = macrosampleinfo,
                                    design = macro_design)

macro_DEdds <- DESeq(macro_dds)

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/DispersionEstimate_macrophage_pre-filter_Mo.tiff", 
     width = 500, height = 300)
plotDispEsts(macro_DEdds)
dev.off()

save(macro_dds, macro_DEdds, macro_modelMatrix, macro_design, macrosampleinfo, macrocountdata,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_all_DESeq_no_filter_Mo.RData")

load("~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_all_DESeq_no_filter_Mo.RData")

# This dispersion estimate seems ok, so no need at this stage to filter
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/DispersionEstimate_macrophage_with-filter-Mo.tiff", 
     width = 500, height = 300)
plotDispEsts(macro_DEdds)
dev.off()

# pca 
vsd <- vst(macro_DEdds)
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophage_PCA_initial_Mo.tiff", 
     width = 700, height = 400)
plotPCA(vsd, intgroup = "Diagnosis") +
  geom_text_repel(aes(x=PC1, y=PC2, label=macro_dds$FeatureCounts), box.padding = 0.8)
dev.off()

plotPCA(vsd, intgroup = "Batch") +
  geom_text_repel(aes(x=PC1, y=PC2, label=macro_dds$FeatureCounts), box.padding = 0.8)
plotPCA(vsd, intgroup = "Joint") +
  geom_text_repel(aes(x=PC1, y=PC2, label=macro_dds$FeatureCounts), box.padding = 0.8)

save(macro_dds, macro_DEdds, macro_modelMatrix, macro_design, macrosampleinfo, macrocountdata, vsd,
     file="~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Preprocessing/Macrophage_all_DESeq_with_filter_pcas.RData")

# mean vs sd plot of the filtered genes
library("vsn")
ntd <- normTransform(macro_DEdds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

# To look at PCA, boxplot and counts use PCAExplorer, has tendency to crash R
library("pcaExplorer")
pcaExplorer(dds = macro_DEdds)

# To get the normalised counts
normalisedCounts_filter <- counts (macro_DEdds, normalized = TRUE)
logNormalizedCounts_filter <- log2(normalisedCounts_filter + 1)
limma::plotMA(logNormalizedCounts_filter, "5")
abline(h=0, col="red")

# mean vs sd plot of the filtered genes
library("vsn")
ntd <- normTransform(macro_dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

# More PCAs, with biplots, loadings etc. 
vst_macro <- assay(vst(macro_DEdds))

# To colour the graphs accoridng to Diagnosis, Joint etc. 
my_metadata <- data.frame(row.names = colnames(macro_DEdds))
my_metadata$Diagnosis <- macro_DEdds$Diagnosis
my_metadata$Joint <- macro_DEdds$Joint
my_metadata$Batch <- macro_DEdds$Batch
my_metadata$Age <- macro_DEdds$Age
my_metadata$Culture <- macro_DEdds$Co_vs_mono
my_metadata$Fibroblast_ID <- macro_DEdds$SampleName

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

# Biplot of macrophages 
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Batch biplot.tiff",
     width = 700, height = 500)
PCAtools::biplot(p, colby = 'Batch', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Batch', legendPosition = 'right',
                 encircle = TRUE,   encircleAlpha = 1/6)
# Or can draw ellipse in - these are statistically sig, the encirle ones aren't, but for this analysis not sig 
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Joint biplot.tiff",
     width = 700, height = 500)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right',
                 encircle = TRUE,   encircleAlpha = 1/6)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Diagnosis biplot.tiff",
     width = 700, height = 500)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right',
                 encircle = TRUE,   encircleAlpha = 1/6)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Culture biplot.tiff",
     width = 700, height = 500)
PCAtools::biplot(p, colby = 'Culture', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Culture', legendPosition = 'right',
                 encircle = TRUE,   encircleAlpha = 1/6)
dev.off()

# Pairsplots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Batch pairsplot.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Batch', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Joint pairsplot.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Diagnosis pairsplot.tiff",
     width = 1000, height = 700)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=2)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Macrophage RNA seq/Outputs/Final/QC/Macrophages/Culture pairsplot.tiff",
     width = 800, height = 400)
PCAtools::pairsplot(p, colby = 'Culture', pointSize=2)
dev.off()


