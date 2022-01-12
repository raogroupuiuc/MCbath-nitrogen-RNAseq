getwd()

setwd("D:/Box Sync/02 - RNA seq work/2002-Methylococcus_bath_rnaseq/2021/05_R_analysis/")

options(stringsAsFactors = FALSE)

# BiocManager::install("PCAtools")


# Load libraries ####

library(edgeR)
library(statmod)
library(affycoretools)
library(limma)
library(gplots)
library(rgl)
library(Glimma)
library(rtracklayer)
library(sva)
library(ggplot2)
library(PCAtools)
# library(WGCNA)

ls()
rm(list = ls())

# [ Load this to add all data to the environment ] #### 

# load("2021-10-22-05-after-OutConts.Rdata")

# FileName <- dir(path = "../featCounts/", pattern = "featCounts")
# write.csv(FileName, file = "MCbath-CountFileNames.csv")


targets <- readTargets("MC_bath_ReadFates.txt")

targets$Group <- paste(targets$Treatment, targets$Time, sep = ".")
targets$Label <- paste(targets$Treatment, targets$Time, targets$Rep, sep = ".")
targets$Group %>% unique()
#  "AMS.37" "AMS.42" "NMS.37" "NMS.42"

targets$GpF <- factor(targets$Group, levels = c("AMS.37", "AMS.42", "NMS.37", "NMS.42"))
targets$Col <- as.numeric(targets$GpF)
head(targets)

# read in counts #### 

d <- readDGE(targets, path = "../02_featureCounts/", columns = c(1,7), labels = targets$Label, 
             group = targets$Group, comment.char = "#", header = TRUE)

dim(d)
# 3095 12

head(d$samples)

apply(d$counts,2, max) / d$samples$lib.size*100

jpeg("maximum count percentage of each library_210122.jpeg", width = 6, height = 6, units = "in", res = 300, quality = 100)
apply(d$counts,2,function(x){max(x)/sum(x)}) %>%
  barplot(las=2,col=targets$col,main="maximum count percentage of each library", cex.axis = 0.8)
abline(h=mean(apply(d$counts,2,function(x){max(x)/sum(x)})))
dev.off()


# add annotation #### 

gtf0 <- import("../00_Genome/GCF_000008325.1_ASM832v1_genomic.gtf")

table(gtf0$type)

head(d)
head(gtf0)

gtf1 <- gtf0[gtf0$type == "gene"]
head(gtf1)
table(gtf1$type)

gtf1$gene_id %>% unique() %>% length()
# 3095

gtf1$gene %>% unique() %>% length()
# 698 

gtf2 <- gtf0[gtf0$type == "CDS"]
head(gtf2)

gtf2$gene_id %>% unique() %>% length()
# 3039
gtf2$gene %>% unique() %>% length()
# 693 

gtf2$product %>% unique() %>%  length()
# 1872 


all.equal(rownames(d$counts),gtf1$gene_id)
#TRUE 

d$genes <- mcols(gtf1)[,c("gene_id", "gene")]


annot.gene <- unique(mcols(gtf1[,c("gene_id", "gene")]))
annot.CDS <- unique(mcols(gtf2[,c("gene_id", "product", "protein_id")]))

dim(annot.gene)
# 3095 2 
dim(annot.CDS)
# 3039 3

annot.gene$gene_id %>% unique() %>%  length()
# 3095

annot.CDS$gene_id %>% unique() %>% length()
# 3039 

annot <- merge(annot.gene, annot.CDS, by = "gene_id",  
               all.x = TRUE, all.y = FALSE, sort = F)

dim(annot)
# 3095 4


all.equal(d$genes$gene_id,annot$gene_id)
# FALSE 


rm("gtf0", "gtf1", "gtf2", "annot.CDS", "annot.gene")

# Filtering zero count reads #### 
x11()
plotDensities( d$counts, group = d$samples$GpF, col = 1:4, main = "raw counts" )

jpeg("Raw_counts_log_scale.jpeg", width=6, height=6, units="in", res=300, quality=100)
#x11()
plotDensities( log2( d$counts + 0.1 ), group = d$samples$GpF, col = 1:7, main = "Raw Counts (Log scale)" )
dev.off()

jpeg("librarysizes.jpeg", width=6, height=6, units="in", res=300, quality=100)
barplot( d$samples$lib.size / 1e6, ylab = "total number of reads (million)", las = 2, ylim = c(0,50),
         col = d$samples$Col, main = "Library sizes", cex.names = 0.6, names.arg = targets$Label)
dev.off()


jpeg("sequencingsize.jpeg", width=6, height=6, units="in", res=300, quality=100)
barplot( targets$NumReads / 1e6, ylab = "total number of reads (million)", las = 2, ylim = c(0,50),
         col = d$samples$Col, main = "Library sizes", cex.names = 0.6, names.arg = targets$Label)
dev.off()


x11()
hist(log2(rowSums(d$counts + 0.1)), 1000, main = "Total counts per gene, log2 scale")

d <- calcNormFactors(d)

cpm.values <- cpm(d$counts)
head(cpm.values)

above1cpm <- rowSums(cpm.values  >= 1)
head(above1cpm)

x11()
hist(above1cpm, xlab = "Number of samples with > 1 CPM", ylab = "Number of genes")

sum(above1cpm >= 3)
# This is number of genes that passed the filtering & should be kept
# 3033

mean(above1cpm >= 3) 
# 0.9799677 

x11()
plotDensities( log2( d$counts[above1cpm >= 3 , ] + 0.1 ), group = d$samples$GpF, 
               col = 1:11, main = "Counts (log scale) with CPM>3" )

jpeg("Above_1_CPM.jpeg", width=12, height=6, units="in", res=300, quality=100)
hist(above1cpm, xlab = "Samples with CPM > 1", ylab = "Number of genes")
dev.off()

d.filt <- d[above1cpm >= 3, , keep.lib.sizes = FALSE]

d.filt$samples$lib.size / d$samples$lib.size
#  [1] 0.9999974 0.9999801 0.9999951 0.9999943 0.9999321 0.9999903 0.9999948 0.9999721 0.9999938 0.9999935 0.9999743 0.9999943

nrow(d.filt) / nrow(d)
# [1] 0.9799677

x11()
hist(log2(rowSums(d.filt$counts + 0.1)), 1000, main = "Total counts per gene, log2 scale", xlim = c(2,30))

save.image("2021-10-22-01-post_filtering.Rdata")

d.filt <- calcNormFactors(d.filt)

d.filt$samples

jpeg("Tmm-normalization.jpeg", width=12, height=6, units="in", res=300, quality=100)
barplot( d.filt$samples$norm.factors, col = targets$Col, las = 2, legend.text = FALSE, cex.name = 0.7,
         ylab = "TMM norm values", main = "TMM normalization factors", names.arg = targets$Label )
# add line:
abline(h = 1)
dev.off()

all.equal(d.filt$samples$norm.factors,d$samples$norm.factors)
# [1] "Mean relative difference: 0.000665963"

x11()
plot(d.filt$samples$norm.factors,d$samples$norm.factors, col=targets$Col)
abline(0,1)
cor(d.filt$samples$norm.factors,d$samples$norm.factors)
#  0.9999442 corelation b/w d and d-filt.

# clustering ####

x11()
plotMDS(d.filt, main = "MDS plot", col = d.filt$samples$Col, top = 1000 )

glMDSPlot(d.filt, top = 1000, labels = d.filt$samples$Label,
          groups = d.filt$samples[,c("Treatment","Time","Group")],
          folder = "interactive_MDS", launch = FALSE)

### CPM values #### 

logCPM <- cpm(d.filt, log = TRUE)
class(logCPM)
dim(logCPM)
# 3033 12
head(logCPM)

jpeg("number of unique genes.jpeg", width = 10, height = 6, units = "in", res = 300, quality = 100)
unique_index <- rowSums(logCPM>log2(1))==1
barplot(colSums((logCPM>log2(1))[unique_index,]), col = targets$Col, names.arg = targets$Label, 
        las=2, cex.names=0.8, main=paste0("number of unique genes out of ", nrow(logCPM)))
dev.off()



x11()
plotDensities( logCPM, group = d$samples$GpF, col = 1:11, main = "edgeR logCPM" )

hc.edgeR <- hclust(dist(t(logCPM)), method = "average")


jpeg("Dendogram_vertical.jpeg", width=8, height=6, units="in", res=300, quality=100)
plot(hc.edgeR, hang = -1, main = "MCbath - cluster dendogram", sub = "", 
     xlab = "", cex = 0.9)
dev.off()

jpeg("Dendogram_horizontal.jpeg", width=8, height=6, units="in", res=300, quality=100)
nodePar <- list(lab.cex = 0.7, pch = c(NA, 19), cex = 0, col = "black")

plot(as.dendrogram(hc.edgeR), horiz = TRUE, nodePar = nodePar, xlab = "Height", 
     main = "MCbath - cluster dendogram")

dev.off()

# PCA ####


jpeg("PCA_variance.jpeg", width=6, height=6, units="in", res=300, quality=100)
plotPCA(logCPM, screeplot= TRUE) #variance 
dev.off()



jpeg("PCA-axis1-2.jpeg", width=6, height=6, units="in", res=300, quality=100)
plotPCA(logCPM, pch = 16, col = d.filt$samples$Col, groupnames = levels(d.filt$samples$GpF), 
        main = "MCbath - PCA", outside = TRUE)
dev.off()

jpeg("PCA-axis2-3.jpeg", width=6, height=6, units="in", res=300, quality=100)
plotPCA(logCPM, pch = 16, col = d.filt$samples$Col, groupnames = levels(d.filt$samples$GpF), plot3d = FALSE, 
        pcs = 2:3, outside = TRUE)
dev.off()


jpeg("PCA-nicer.jpeg", width=5, height=6, units="in", res=300, quality=100)
test_new <- pca(logCPM)
setEPS()
postscript("PCA_nicer.eps", width = 5, height = 6)
biplot(test_new, #xlim = c(-50,50), ylim = c(-40,50), 
       lab = d.filt$samples$group,
       pointSize = 4, gridlines.major = FALSE, gridlines.minor = FALSE)
dev.off()

x11()
screeplot(test_new)

x11(w=10,h = 10)
pairsplot(test_new, pointSize = 2, gridlines.major = FALSE, gridlines.minor = FALSE)

save.image("2021-10-22-02-after-PCA.Rdata")

# Design matrix #### 

design <- model.matrix(~0 + d.filt$samples$group)
colnames(design) <- levels(d.filt$samples$group)
rownames(design) <- d.filt$samples$Label  


cont.matrix <- makeContrasts(temp.A.42vsA.37 = AMS.42 - AMS.37, 
                             temp.N.42vsN.37 = NMS.42 - NMS.37,
                             nitr.A.42vsN.42 = AMS.42 - NMS.42, 
                             nitr.A.37vsN.37 = AMS.37 - NMS.37, 
                             both.N.42vsA.37 = NMS.42 - AMS.37, 
                             both.A.42vsN.37 = AMS.42 - NMS.37, levels=design)

fit <- lmFit(logCPM, design = design)

summary(decideTests(fit))
# AMS.37 AMS.42 NMS.37 NMS.42
# Down        8      0      6      3
# NotSig    156    134    143    140
# Up       2869   2899   2884   2890

fit2 <- contrasts.fit(fit, contrasts = cont.matrix)

fit2 <- eBayes(fit2, trend = T)

summary(decideTests(fit2, p.value = 0.05))
#        temp.A.42vsA.37 temp.N.42vsN.37 nitr.A.42vsN.42 nitr.A.37vsN.37 both.N.42vsA.37 both.A.42vsN.37
# Down                 5             571              13             140              70             725
# NotSig            3027            1954            2991            2674            2910            1681
# Up                   1             508              29             219              53             627


jpeg("Venn_nitrogen-comparision.jpeg", width=12, height=6, units="in", res=300, quality=100)
layout(matrix(2:1,1,2)) #this divides the graph into 3 parts
vennDiagram(decideTests(fit2[,c(1,2)]), include = "down", main = "Downregulated genes", circle.col = c("red", "blue", "green3"))
vennDiagram(decideTests(fit2[,c(1,2)]), include = "up", main = "Upregulated genes", circle.col = c("red", "blue", "green3"))
dev.off()

jpeg("Venn_temperature-comparision.jpeg", width=12, height=6, units="in", res=300, quality=100)
layout(matrix(2:1,1,2)) #this divides the graph into 3 parts
vennDiagram(decideTests(fit2[,c(3,4)]), include = "down", main = "Downregulated genes", circle.col = c("red", "blue", "green3"))
vennDiagram(decideTests(fit2[,c(3,4)]), include = "up", main = "Upregulated genes", circle.col = c("red", "blue", "green3"))
dev.off()


# jpeg("All.jpeg", width=12, height=6, units="in", res=300, quality=100)
setEPS()
postscript("venn_all_compared.eps", width = 12, height = 6)
layout(matrix(2:1,1,2)) #this divides the graph into 2 parts
vennDiagram(decideTests(fit2[,c(1:4)]), include = "down", main = "Downregulated genes", circle.col = c("red", "blue", "green3", "yellow"))
vennDiagram(decideTests(fit2[,c(1:4)]), include = "up", main = "Upregulated genes", circle.col = c("red", "blue", "green3", "yellow"))
dev.off()


# jpeg("All_combined.jpeg", width=6, height=6, units="in", res=300, quality=100)
setEPS()
postscript("venn_all_combined_compared.eps", width = 6, height = 6)
vennDiagram(decideTests(fit2[,c(1:4)]), circle.col = c("red", "blue", "green3", "yellow"))

dev.off()


# 1  way test #### 

res.anova <- topTable(fit2, coef = 1:6, num = Inf, sort.by = "none" )

sum(res.anova$adj.P.Val < 0.05)
#1190

Sig_index <- topTable(fit2, coef = 1:6, n=Inf, sort.by = "none")$adj.P.Val < 0.05
table(Sig_index)
# Sig_index
# FALSE  TRUE 
# 1843  1190


source("summarizeFit.R")

out.conts <- summarizeFit(fit2, addAnova = TRUE)
names(out.conts)
# [1] "Amean"                "FC.temp.A.42vsA.37"   "rawP.temp.A.42vsA.37" "FDR.temp.A.42vsA.37"  "FC.temp.N.42vsN.37"  
# [6] "rawP.temp.N.42vsN.37" "FDR.temp.N.42vsN.37"  "FC.nitr.A.42vsN.42"   "rawP.nitr.A.42vsN.42" "FDR.nitr.A.42vsN.42" 
# [11] "FC.nitr.A.37vsN.37"   "rawP.nitr.A.37vsN.37" "FDR.nitr.A.37vsN.37"  "FC.both.N.42vsA.37"   "rawP.both.N.42vsA.37"
# [16] "FDR.both.N.42vsA.37"  "FC.both.A.42vsN.37"   "rawP.both.A.42vsN.37" "FDR.both.A.42vsN.37"  "F"                   
# [21] "F.rawP"               "F.FDR"  

sum(out.conts$F.FDR < 0.05)
# 1190 - same as before 

colnames(out.conts)[20:22] <- c("Fstat.ANOVA","rawP.ANOVA","FDR.ANOVA")
rownames(out.conts) <- (rownames(fit2$coefficients))
out.conts$gene <- rownames(out.conts) 

save.image("2021-10-22-03-after-anovatesting.Rdata")


#### heatmap ####

getwd()

setwd("./heatmaps/")

N.42 = c("NMS.42.1", "NMS.42.2", "NMS.42.3")
A.42 = c("AMS.42.1", "AMS.42.2", "AMS.42.3")
N.37 = c("NMS.37.1", "NMS.37.2", "NMS.37.3")
A.37 = c("AMS.37.1", "AMS.37.2", "AMS.37.3")

logCPM.heat.A.42vsN.42 <- t(scale(t(logCPM[Sig_index,c(A.42, N.42)])))
logCPM.heat.A.37vsN.37 <- t(scale(t(logCPM[Sig_index,c(A.37, N.37)])))
logCPM.heat.N.42vsN.37 <- t(scale(t(logCPM[Sig_index,c(N.42, N.37)])))
logCPM.heat.A.42vsA.37 <- t(scale(t(logCPM[Sig_index,c(A.42, A.37)])))
logCPM.heat <- t(scale(t(logCPM[Sig_index,])))

dim(logCPM.heat)

col.pan <- colorpanel(100, "blue", "white", "red")

logCPM.heat.reorg <- logCPM.heat[,c(7:9, 1:3, 10:12, 4:6)]

# x11(10,15)
jpeg("Heatmap_MC_bath_overall.jpeg", width = 6, height = 10, units = "in", res = 300, quality = 100)
heatmap.2(logCPM.heat.reorg, col=col.pan, Rowv=TRUE,Colv=FALSE, scale="none",
          trace="none", dendrogram="none", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3,6,9), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat), " genes" ))
dev.off()



#Do blocks on heatmap

x.cluster.h1 <- hclust(dist(logCPM.heat))

x11(10,15)
heatmap.2(logCPM.heat, col=col.pan, Rowv=as.dendrogram(x.cluster.h1),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3,6,9), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat), " genes" ))

x11(12,6)
plot(x.cluster.h1, labels = FALSE)
abline(h = 4.5 )

rowcols.h1 <- cutree(x.cluster.h1, h = 4.5)
table(rowcols.h1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14 
# 168 258  43  48 174 210  70 109  31  35   7  31   4   3 

# jpeg("Heatmap_overall_clustered.jpeg", width = 6, height = 10, units = "in", res = 300, quality = 100)
setEPS()
postscript("Heatmap_overall_clustered.eps", width = 6, height = 10)

heatmap.2(logCPM.heat.reorg, col=col.pan, Rowv=as.dendrogram(x.cluster.h1),Colv=FALSE, scale="none",
          trace="none", dendrogram="none", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3,6,9), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat), " genes" ))
dev.off()

## individual heatmaps #### 

# logCPM.heat.A.42vsN.42 <- t(scale(t(logCPM[Sig_index,c(A.42, N.42)])))
# logCPM.heat.A.37vsN.37 <- t(scale(t(logCPM[Sig_index,c(A.37, N.37)])))
# logCPM.heat.N.42vsN.37 <- t(scale(t(logCPM[Sig_index,c(N.42, N.37)])))
# logCPM.heat.A.42vsA.37 <- t(scale(t(logCPM[Sig_index,c(A.42, A.37)])))

jpeg("Heatmap_A42_vs_N42.jpeg", width = 6, height = 10, units = "in", res = 300, quality = 100)
heatmap.2(logCPM.heat.A.42vsN.42, col=col.pan, Rowv=as.dendrogram(hclust(dist(logCPM.heat.A.42vsN.42))),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat.A.42vsN.42), " genes" ))
dev.off()

jpeg("Heatmap_A37_vs_N37.jpeg", width = 6, height = 10, units = "in", res = 300, quality = 100)
heatmap.2(logCPM.heat.A.37vsN.37, col=col.pan, Rowv=as.dendrogram(hclust(dist(logCPM.heat.A.37vsN.37))),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat.A.37vsN.37), " genes" ))
dev.off()

jpeg("Heatmap_N42_vs_N37.jpeg", width = 6, height = 10, units = "in", res = 300, quality = 100)
heatmap.2(logCPM.heat.N.42vsN.37, col=col.pan, Rowv=as.dendrogram(hclust(dist(logCPM.heat.N.42vsN.37))),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat.N.42vsN.37), " genes" ))
dev.off()

jpeg("Heatmap_A42_vs_A37.jpeg", width = 6, height = 10, units = "in", res = 300, quality = 100)
heatmap.2(logCPM.heat.A.42vsA.37, col=col.pan, Rowv=as.dendrogram(hclust(dist(logCPM.heat.A.42vsA.37))),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3), sepcolor = "black",
          margins = c(7,2),density.info = "none",
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat.A.42vsA.37), " genes" ))
dev.off()


names(out.conts)
# [1] "Amean"                "FC.temp.A.42vsA.37"   "rawP.temp.A.42vsA.37" "FDR.temp.A.42vsA.37"  "FC.temp.N.42vsN.37"  
# [6] "rawP.temp.N.42vsN.37" "FDR.temp.N.42vsN.37"  "FC.nitr.A.42vsN.42"   "rawP.nitr.A.42vsN.42" "FDR.nitr.A.42vsN.42" 
# [11] "FC.nitr.A.37vsN.37"   "rawP.nitr.A.37vsN.37" "FDR.nitr.A.37vsN.37"  "FC.both.N.42vsA.37"   "rawP.both.N.42vsA.37"
# [16] "FDR.both.N.42vsA.37"  "FC.both.A.42vsN.37"   "rawP.both.A.42vsN.37" "FDR.both.A.42vsN.37"  "Fstat.ANOVA"         
# [21] "rawP.ANOVA"           "FDR.ANOVA"            "gene"     

all.equal(rownames(out.conts), rownames(logCPM))
# TRUE 

out.conts <- cbind(out.conts, logCPM)


save.image("../2021-10-22-04-after-heatmaps.RData")


# add annotation ##### 

dim(out.conts)
# 3033 35

out.conts2 <- merge(out.conts, annot, by.x = "gene", by.y = "gene_id", sort = FALSE)

dim(out.conts2)
# 3033 38

names(out.conts2)
# [1] "gene"                 "Amean"                "FC.temp.A.42vsA.37"   "rawP.temp.A.42vsA.37" "FDR.temp.A.42vsA.37" 
# [6] "FC.temp.N.42vsN.37"   "rawP.temp.N.42vsN.37" "FDR.temp.N.42vsN.37"  "FC.nitr.A.42vsN.42"   "rawP.nitr.A.42vsN.42"
# [11] "FDR.nitr.A.42vsN.42"  "FC.nitr.A.37vsN.37"   "rawP.nitr.A.37vsN.37" "FDR.nitr.A.37vsN.37"  "FC.both.N.42vsA.37"  
# [16] "rawP.both.N.42vsA.37" "FDR.both.N.42vsA.37"  "FC.both.A.42vsN.37"   "rawP.both.A.42vsN.37" "FDR.both.A.42vsN.37" 
# [21] "Fstat.ANOVA"          "rawP.ANOVA"           "FDR.ANOVA"            "AMS.37.1"             "AMS.37.2"            
# [26] "AMS.37.3"             "AMS.42.1"             "AMS.42.2"             "AMS.42.3"             "NMS.37.1"            
# [31] "NMS.37.2"             "NMS.37.3"             "NMS.42.1"             "NMS.42.2"             "NMS.42.3"            
# [36] "gene.y"               "product"              "protein_id" 




write.table(out.conts2[,c(1,38,36,37,2,21:23,3:20,24:35)], file = "../Results_MCbath_RNAseq_2021-10-22.txt", 
            row.names = FALSE, sep = "\t")

save.image("../2021-10-22-05-after-OutConts.Rdata")
sessionInfo()
# 
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] PCAtools_2.4.0       ggrepel_0.9.1        ggplot2_3.3.5        sva_3.40.0           BiocParallel_1.26.2 
# [6] genefilter_1.74.1    mgcv_1.8-36          nlme_3.1-152         rtracklayer_1.52.1   GenomicRanges_1.44.0
# [11] GenomeInfoDb_1.28.4  IRanges_2.26.0       S4Vectors_0.30.2     Glimma_2.2.0         rgl_0.107.14        
# [16] gplots_3.1.1         affycoretools_1.64.0 Biobase_2.52.0       BiocGenerics_0.38.0  statmod_1.4.36      
# [21] edgeR_3.34.1         limma_3.48.3        
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2                  R.utils_2.11.0              tidyselect_1.1.1            RSQLite_2.2.8              
# [5] AnnotationDbi_1.54.1        htmlwidgets_1.5.4           grid_4.1.1                  munsell_0.5.0              
# [9] ScaledMatrix_1.0.0          codetools_0.2-18            preprocessCore_1.54.0       withr_2.4.2                
# [13] colorspace_2.0-2            Category_2.58.0             filelock_1.0.2              OrganismDbi_1.34.0         
# [17] knitr_1.36                  rstudioapi_0.13             labeling_0.4.2              MatrixGenerics_1.4.3       
# [21] GenomeInfoDbData_1.2.6      hwriter_1.3.2               farver_2.1.0                bit64_4.0.5                
# [25] vctrs_0.3.8                 generics_0.1.0              xfun_0.26                   biovizBase_1.40.0          
# [29] BiocFileCache_2.0.0         R6_2.5.1                    rsvd_1.0.5                  locfit_1.5-9.4             
# [33] AnnotationFilter_1.16.0     bitops_1.0-7                cachem_1.0.6                reshape_0.8.8              
# [37] DelayedArray_0.18.0         assertthat_0.2.1            BiocIO_1.2.0                scales_1.1.1               
# [41] nnet_7.3-16                 gtable_0.3.0                beachmat_2.8.1              affy_1.70.0                
# [45] ggbio_1.40.0                ensembldb_2.16.4            rlang_0.4.11                splines_4.1.1              
# [49] lazyeval_0.2.2              dichromat_2.0-0             checkmate_2.0.0             BiocManager_1.30.16        
# [53] yaml_2.2.1                  reshape2_1.4.4              GenomicFeatures_1.44.2      backports_1.2.1            
# [57] Hmisc_4.6-0                 RBGL_1.68.0                 tools_4.1.1                 affyio_1.62.0              
# [61] ellipsis_0.3.2              ff_4.0.4                    RColorBrewer_1.1-2          Rcpp_1.0.7                 
# [65] plyr_1.8.6                  sparseMatrixStats_1.4.2     base64enc_0.1-3             progress_1.2.2             
# [69] zlibbioc_1.38.0             purrr_0.3.4                 RCurl_1.98-1.5              prettyunits_1.1.1          
# [73] rpart_4.1-15                cowplot_1.1.1               SummarizedExperiment_1.22.0 cluster_2.1.2              
# [77] magrittr_2.0.1              data.table_1.14.2           ProtGenerics_1.24.0         matrixStats_0.61.0         
# [81] hms_1.1.1                   xtable_1.8-4                XML_3.99-0.8                jpeg_0.1-9                 
# [85] gcrma_2.64.0                gridExtra_2.3               compiler_4.1.1              biomaRt_2.48.3             
# [89] tibble_3.1.5                KernSmooth_2.23-20          crayon_1.4.1                ReportingTools_2.32.1      
# [93] R.oo_1.24.0                 htmltools_0.5.2             GOstats_2.58.0              Formula_1.2-4              
# [97] geneplotter_1.70.0          DBI_1.1.1                   dbplyr_2.1.1                rappdirs_0.3.3             
# [101] Matrix_1.3-4                R.methodsS3_1.8.1           pkgconfig_2.0.3             GenomicAlignments_1.28.0   
# [105] foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1               annotate_1.70.0            
# [109] dqrng_0.3.0                 XVector_0.32.0              AnnotationForge_1.34.1      stringr_1.4.0              
# [113] VariantAnnotation_1.38.0    digest_0.6.28               graph_1.70.0                Biostrings_2.60.2          
# [117] htmlTable_2.3.0             DelayedMatrixStats_1.14.3   GSEABase_1.54.0             restfulr_0.0.13            
# [121] curl_4.3.2                  Rsamtools_2.8.0             gtools_3.9.2                rjson_0.2.20               
# [125] lifecycle_1.0.1             jsonlite_1.7.2              PFAM.db_3.13.0              BSgenome_1.60.0            
# [129] fansi_0.5.0                 pillar_1.6.4                lattice_0.20-44             GGally_2.1.2               
# [133] KEGGREST_1.32.0             fastmap_1.1.0               httr_1.4.2                  survival_3.2-11            
# [137] GO.db_3.13.0                glue_1.4.2                  png_0.1-7                   iterators_1.0.13           
# [141] bit_4.0.4                   Rgraphviz_2.36.0            stringi_1.7.5               blob_1.2.2                 
# [145] oligoClasses_1.54.0         DESeq2_1.32.0               BiocSingular_1.8.1          latticeExtra_0.6-29        
# [149] caTools_1.18.2              memoise_2.0.0               dplyr_1.0.7                 irlba_2.3.3