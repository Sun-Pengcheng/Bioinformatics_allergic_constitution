# prepare -------------------------------------------------------------------

rm(list=ls())
getwd()
options(stringsAsFactors = F)
library(stringr)  
library(data.table)

fpkm1 <- read.table("./rdata/GSE75011/mat_fpkm_norm.txt",header = T,sep = '\t')
colnames(fpkm1)
max(fpkm1)
min(fpkm1)
boxplot(fpkm1)
sum(is.na(fpkm1))


fpkm2 <- read.table("./rdata/GSE125916/mat_fpkm_norm.txt",header = T,sep = '\t')
colnames(fpkm2)
max(fpkm2)
min(fpkm2)
boxplot(fpkm2)
sum(is.na(fpkm2))
which(is.na(fpkm2))

fpkm3 <- read.table("./rdata/GSE184237/mat_fpkm_norm.txt",header = T,sep = '\t')
colnames(fpkm3)
max(fpkm3)
min(fpkm3)
boxplot(fpkm3)
sum(is.na(fpkm3))


all_gene <- intersect(rownames(fpkm1),intersect(rownames(fpkm3),rownames(fpkm2)))
length(all_gene)
fpkm1 <- fpkm1[all_gene,]
fpkm2 <- fpkm2[all_gene,]
fpkm3 <- fpkm3[all_gene,]
identical(rownames(fpkm1),rownames(fpkm3))
identical(rownames(fpkm1),rownames(fpkm2))

fpkm <- cbind(fpkm1,fpkm2,fpkm3)
sum(is.na(fpkm))
boxplot(fpkm)

group <- read.table("./rdata/WGCNA/group.txt", sep = '\t', header = T)

identical(group$id,colnames(fpkm))
fpkm <- fpkm[,group$id]
identical(group$id,colnames(fpkm))

colnames(fpkm) <- group$geo_accession


library(sva)
group$batch <- factor(group$gse)
batch = group$batch
# exp[1:5, 1:5]
boxplot(fpkm,range=0,col="blue")

sum(is.na(fpkm))

# combat校正
combat_edata = ComBat(dat=fpkm, batch=batch, mod=NULL, par.prior=TRUE)
# combat_edata[1:5, 1:5]

save(combat_edata,group,file = "./rdata/WGCNA/input.Rda")



# WGCNA -------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(WGCNA)

enableWGCNAThreads()

load("./rdata/WGCNA/input.Rda")

data <- combat_edata
head(data)
dim(data)
max(data)
min(data)

datExpr0 <- t(data)
head(datExpr0[,1:5])
dim(datExpr0)
tmp1 <- goodGenes(datExpr0)
sum(tmp1)==ncol(datExpr0)
tmp2 <- goodSamples(datExpr0)
sum(tmp2)==nrow(datExpr0)
rm(tmp1, tmp2)

data[1:8,1:8]

datExpr0[1:8,1:8]
dim(data)

sampleTree <- hclust(dist(datExpr0, method = "manhattan"), method = "average")
sampleTree

pdf("./figure/1.Sample_clustering.pdf", width = 16, height = 8)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1)
abline(h=2500, col="red")
dev.off()


# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2500, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
table(keepSamples)
datExpr0 = datExpr0[keepSamples, ]
data=data[,keepSamples]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

pheno_filter <- group[group$geo_accession%in%colnames(data),]
table(pheno_filter$group)

sampleTree2 <- hclust(dist(datExpr0, method = "manhattan"), method = "average")
sampleTree2

pdf("./figure/1.Sample_clustering_cut.pdf", width = 16, height = 8)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1)
dev.off()

save(datExpr0,pheno_filter,file = "./rdata/WGCNA/ready.Rda")


# input expression matrix, calculate Soft Threshold -------------------------------------------------------------------
rm(list = ls())
load("./data/WGCNA/ready.Rda")
dim(datExpr0)
# power
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
powers

sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, RsquaredCut = 0.9)
sft
sft$powerEstimate

pdf("./figure/2.power_0.9.pdf", width = 7, height = 4)
par(mfrow = c(1,2))
plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit, signed R^2", type="n",
  main=paste("Scale independence")
)
text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,
  cex=0.85,
  col="red"
)
abline(h=0.9, col="red")
cex1 = 0.85

plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab="Soft Threshold (power)",
  ylab="Mean Connectivity",
  type='n',
  main=paste("Mean connectivity")
)
text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

# module build --------------------------------------------------------------------

dim(datExpr0)

cor = WGCNA::cor
net = blockwiseModules(datExpr0,
                       power = sft$powerEstimate, 
                       maxBlockSize = 30000,
                       TOMType = "unsigned", 
                       minModuleSize = 50, 
                       reassignThreshold = 0,
                       mergeCutHeight = 0.15, 
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE, 
                       verbose = 3)
cor = stats::cor
save(net, sft, file = "./rdata/WGCNA/WGCNA_net.Rda")


# reload ------------------------------------------------------------------
library(WGCNA)
load(file = "./data/WGCNA/WGCNA_net.Rda")

mergedColors = labels2colors(net$colors)
table(mergedColors)
length(table(mergedColors))

gene_and_color <- data.frame(gene_id=colnames(datExpr0), colors=mergedColors)
head(gene_and_color)  
write.table(gene_and_color, file = "./rdata/gene_and_color.txt", sep = '\t',
            row.names = F, col.names = T, quote = F) 

pdf("./figure/3.Dendrogram.pdf", width = 14, height = 12)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modules", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.colorLabels = 3, marAll = c(1, 13, 3, 1),
                    cex.axis = 3, cex.main=3,cex.lab=3, mgp=c(5,1,0))
dev.off()

# calculate the adjacency of genes
adjacency = adjacency(datExpr0, power = sft$powerEstimate)
All_degrees=intramodularConnectivity(adjacency, mergedColors)
All_degrees <- All_degrees %>% rownames_to_column(.,var = "GeneID")
All_degrees <- merge(All_degrees,gene_and_color,by.x="GeneID",by.y="gene_id")

all_color <- unique(gene_and_color$colors)
gene_rank <- data.frame()
for (i in all_color) {
  # i <- "greenyellow"
  temp <- All_degrees[All_degrees$colors==i,]
  temp$Rank <- rank(-temp$kWithin)
  temp$GeneNumber <- rep(nrow(temp),nrow(temp))
  temp$Rank_ratio <- temp$Rank/temp$GeneNumber*100
  temp$Rank_ratio <- round(temp$Rank_ratio,1)
  gene_rank <- rbind(gene_rank,temp)
}

write.table(gene_rank,file =  "./output/WGCNA_gene_connectivity.txt", sep = '\t',
            row.names = F, col.names = T, quote = F)
connectivity <- read.table("./output/WGCNA_gene_connectivity.txt",sep = '\t',header = T)
colnames(connectivity)
colnames(gene_and_color)

gene_and_color <- read.table(file = "./rdata/gene_and_color.txt", sep = '\t',
                            header = T) 
colnames(connectivity)
connectivity <- connectivity[,-6]
gene_and_color <- merge(gene_and_color,connectivity,by.x="gene_id",by.y="GeneID")
colnames(gene_and_color)
write.table(gene_and_color,file =  "./output/WGCNA_gene_color_connectivity.txt", sep = '\t',
            row.names = F, col.names = T, quote = F)

# reload: gene color connectivity ---------------------------------------------


gene_and_color <- read.table( "./output/WGCNA_gene_color_connectivity.txt",sep = '\t',header = T)

# calculate the moduleEigengenes
MEList = moduleEigengenes(datExpr0, colors = mergedColors)
MEs = MEList$eigengenes
head(MEs)

MET = orderMEs(MEs)
rownames(MET) = rownames(datExpr0)
head(MET)[,1:5]

pdf(file = "./figure/eigengeneNetwork.pdf",width = 6,height = 6)
plotEigengeneNetworks(MET, '', marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2),
                      cex.lab = 0.8, xLabelsAngle = 90, printAdjacency = T, cex.adjacency = 0.8)
dev.off()


# modules and samples -------------------------------------------------------------------

moduleColors <- labels2colors(net$colors)

modules <- mergedColors
allModules <- unique(moduleColors)
allModules

Genes <- colnames(datExpr0)
inModule <- is.finite(match(moduleColors, modules))
modGenes <- Genes[inModule]

nSamples <- nrow(datExpr0)
sampleMatrix <- as.data.frame(diag(x=1, nrow = nSamples))

rownames(sampleMatrix) <- rownames(datExpr0)
colnames(sampleMatrix) <- rownames(datExpr0)

sampleMatrix[1:5,1:5]

MET[1:5,1:5]
moduleSampleCor <- cor(MET, sampleMatrix, use = "p")
moduleSampleCor[1:5,1:5]

moduleSamplePvalue <- corPvalueStudent(moduleSampleCor, nSamples)

textMatrix <- paste(signif(moduleSampleCor,2), "\n(", signif(moduleSamplePvalue,1),")", sep = "")
head(textMatrix)

pdf("./figure/module_sample.pdf", width = 30, height = 10)
par(mar = c(4,8,4,4))
labeledHeatmap(Matrix = moduleSampleCor, xLabels = names(sampleMatrix), yLabels = names(MET),
               ySymbols = names(MET), xLabelsAngle = 90, cex.lab = 0.6,
               colorLabels = FALSE, colors = blueWhiteRed(100),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1,1),
               yColorWidth = 0.03, xColorWidth = 0.05, main = paste("Module-Sample relationship"))
dev.off()



# modules and traits -------------------------------------------------------------------

nSamples <- nrow(datExpr0)
sampleMatrix <- as.data.frame(diag(x=1, nrow = nSamples))

dataTrait <- data.frame(trait=pheno_filter$geo_accession,
                        group=pheno_filter$group)
dataTrait <- dataTrait %>% column_to_rownames(.,var="trait")
dataTrait$Allergy <- ifelse(dataTrait$group=="Control",0,1)

dataTrait$AR <- ifelse(pheno_filter$class=="AR" & pheno_filter$group=="Patients",1,ifelse(pheno_filter$group=="Control",0,NA))
dataTrait$AA <- ifelse(pheno_filter$class=="AA" & pheno_filter$group=="Patients",1,ifelse(pheno_filter$group=="Control",0,NA))
dataTrait$AD <- ifelse(pheno_filter$class=="AD" & pheno_filter$group=="Patients",1,ifelse(pheno_filter$group=="Control",0,NA))

colnames(dataTrait)
dataTrait <- dataTrait[,-1]

moduleTraitCor = cor(MET, dataTrait, use = "p")
moduleTraitCor

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1),")", sep = "")


module_trait_pvalue <- data.frame("color"=rownames(moduleTraitPvalue), moduleTraitPvalue)
# module_trait_pvalue
module_trait_pvalue[module_trait_pvalue$allergy<0.05,]$color

pdf("./figure/4.trait.pdf", width = 5, height = 4)
par(mar = c(4,10,4,4))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(dataTrait), yLabels = names(MET), 
               ySymbols = names(MET),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
table(mergedColors)
table(pheno_filter$group)

# Modules and genes ----------------------------------------------------------------

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
head(MEs)[,1:8]
head(datExpr0)[,1:8]
head(geneModuleMembership)[,1:8]

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
head(MMPvalue)[,1:5]
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr0, dataTrait$Allergy, use = "p"));
head(geneTraitSignificance)

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", "Allergy", sep="");
names(GSPvalue) = paste("p.GS.", "Allergy", sep="");

# plot for sig modules

modules <- "blue" # first time, plot for blue module
modules <- "turquoise" # second time, plot for turquoise module

column = match(modules, modNames);
moduleGenes = mergedColors==modules;
gene_choose <- colnames(datExpr0)[moduleGenes]

pdf(file = paste0("./figure/MSvsGS_", modules, ".pdf"), width = 5, height = 5)
# sizeGrWindow(5, 5);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", modules, "module"),
                   ylab = "Gene significance for allergy",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = modules)
abline(v=0.7, h=0.2, col = "red", lty = 2, lwd=2)
dev.off()


# blue --------------------------------------------------------------------

table(mergedColors)
modules <- "blue"

column = match(modules, modNames);
moduleGenes = mergedColors==modules;
gene_choose <- colnames(datExpr0)[moduleGenes]

blue <- data.frame(gene=gene_choose)

dim(geneModuleMembership)
head(geneModuleMembership)
MM <- geneModuleMembership %>% rownames_to_column(.,var = "gene")
colnames(MM)
MM <- MM[,c("gene","MMblue")]
blue <- merge(blue,MM,by="gene")

dim(MMPvalue)
head(MMPvalue)
MM_P <- MMPvalue %>% rownames_to_column(.,var = "gene")
colnames(MM_P)
MM_P <- MM_P[,c("gene","p.MMblue")]
blue <- merge(blue,MM_P,by="gene")
colnames(blue) <- c("gene","MM","p.MM")

dim(geneTraitSignificance)
head(geneTraitSignificance)
GS <- geneTraitSignificance %>% rownames_to_column(.,var = "gene")
blue <- merge(blue,GS,by="gene")

GS_P <- GSPvalue %>% rownames_to_column(.,var = "gene")
blue <- merge(blue,GS_P,by="gene")

colnames(gene_and_color)
blue <- merge(blue,gene_and_color,by.x="gene",by.y = "gene_id")

colnames(blue)

# turquoise -------------------------------------------------------------------

modules <- "turquoise"

column = match(modules, modNames);
moduleGenes = mergedColors==modules;
gene_choose <- colnames(datExpr0)[moduleGenes]

turquoise <- data.frame(gene=gene_choose)

dim(geneModuleMembership)
head(geneModuleMembership)
MM <- geneModuleMembership %>% rownames_to_column(.,var = "gene")
colnames(MM)
MM <- MM[,c("gene","MMturquoise")]
turquoise <- merge(turquoise,MM,by="gene")

dim(MMPvalue)
head(MMPvalue)
MM_P <- MMPvalue %>% rownames_to_column(.,var = "gene")
colnames(MM_P)
MM_P <- MM_P[,c("gene","p.MMturquoise")]
turquoise <- merge(turquoise,MM_P,by="gene")
colnames(turquoise) <- c("gene","MM","p.MM")

dim(geneTraitSignificance)
head(geneTraitSignificance)
GS <- geneTraitSignificance %>% rownames_to_column(.,var = "gene")
turquoise <- merge(turquoise,GS,by="gene")

GS_P <- GSPvalue %>% rownames_to_column(.,var = "gene")
turquoise <- merge(turquoise,GS_P,by="gene")

colnames(gene_and_color)
turquoise <- merge(turquoise,gene_and_color,by.x="gene",by.y = "gene_id")

colnames(turquoise)


# filter ------------------------------------------------------------------


blue_filter <- blue[blue$kWithin>0.3,]
blue_filter <- blue_filter[abs(blue_filter$MM)>0.7 & abs(blue_filter$GS.Allergy)>0.2,]

turquoise_filter <- turquoise[turquoise$kWithin>0.3,]
turquoise_filter <- turquoise_filter[abs(turquoise_filter$MM)>0.7 & abs(turquoise_filter$GS.Allergy)>0.2,]

write.table(blue_filter,file = "./output/WGCNA_blue_filter.txt",sep = '\t',col.names = T,row.names = F,quote = F)
write.table(turquoise_filter,file = "./output/WGCNA_turquoise_filter.txt",sep = '\t',col.names = T,row.names = F,quote = F)


# recalculate the relation of modules and filtered genes -------------------------------------------------------

all <- rbind(blue_filter,turquoise_filter)

head(mergedColors)

datExpr0[1:5,1:5]
sum(colnames(datExpr0)%in%all$gene)

mergedColors <- mergedColors[colnames(datExpr0)%in%all$gene]
table(mergedColors)

datExpr0 <- datExpr0[,colnames(datExpr0)%in%all$gene]

MEList = moduleEigengenes(datExpr0, colors = mergedColors)
MEs = MEList$eigengenes
head(MEs)

MET = orderMEs(MEs)
rownames(MET) = rownames(datExpr0)
head(MET)[,1:2]

nSamples <- nrow(datExpr0)
sampleMatrix <- as.data.frame(diag(x=1, nrow = nSamples))

dataTrait <- data.frame(trait=pheno_filter$geo_accession,
                        group=pheno_filter$group)
dataTrait <- dataTrait %>% column_to_rownames(.,var="trait")
dataTrait$Allergy <- ifelse(dataTrait$group=="Control",0,1)

dataTrait$AR <- ifelse(pheno_filter$class=="AR" & pheno_filter$group=="Patients",1,ifelse(pheno_filter$group=="Control",0,NA))
dataTrait$AA <- ifelse(pheno_filter$class=="AA" & pheno_filter$group=="Patients",1,ifelse(pheno_filter$group=="Control",0,NA))
dataTrait$AD <- ifelse(pheno_filter$class=="AD" & pheno_filter$group=="Patients",1,ifelse(pheno_filter$group=="Control",0,NA))

colnames(dataTrait)
dataTrait <- dataTrait[,-1]

moduleTraitCor = cor(MET, dataTrait, use = "p")
moduleTraitCor

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1),")", sep = "")

module_trait_pvalue <- data.frame("color"=rownames(moduleTraitPvalue), moduleTraitPvalue)
module_trait_pvalue
