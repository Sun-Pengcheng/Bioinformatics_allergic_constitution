# RRA algorithm
rm(list = ls())
library(tidyverse)

path_AR1 <- "./rdata/GSE75011/DEG_DESeq2_AR1.Rdata"
path_AA1 <- "./rdata/GSE75011/DEG_DESeq2_AA1.Rdata"
path_AD1 <- "./rdata/GSE125916/DEG_DESeq2_AD1.Rdata"
path_AD2 <- "./rdata/DEG_DESeq2_AD2.Rdata"

# 1.1 load data AR1----
load(file = path_AR1)
head(DEG_DESeq2)

DEG_limma_voom = DEG_DESeq2
names(DEG_limma_voom)
colnames(DEG_limma_voom)[2]='logFC'
colnames(DEG_limma_voom)[5]='P.Value'
DEG_limma_voom=DEG_limma_voom[order(DEG_limma_voom$logFC),]
head(DEG_limma_voom)

## set threshold to filter up and down genes
## in general, we choose 1 as the logFC threshold, 0.05 as the pvalue threshold
cut_logFC <- with(DEG_limma_voom,mean(abs(DEG_limma_voom$logFC)) + 2*sd(abs(DEG_limma_voom$logFC)) )
cut_logFC

if (cut_logFC > 1) {
  logFC_t <- 1
}

if (cut_logFC < 1) {
  logFC_t <- cut_logFC
}
logFC_t

pvalue_t <- 0.05

deg1_up <- DEG_limma_voom[with(DEG_limma_voom,logFC > logFC_t & P.Value < pvalue_t),]
deg1_up <- deg1_up[order(deg1_up$logFC,decreasing = T),]
deg1_up1 <- rownames(deg1_up)

deg1_down <- DEG_limma_voom[with(DEG_limma_voom,logFC < -logFC_t & P.Value < pvalue_t),]
deg1_down <- deg1_down[order(deg1_down$logFC,decreasing = F),]
deg1_down1 <- rownames(deg1_down)

# 1.2 load data AA1----
load(file = path_AA1)
head(DEG_DESeq2)

DEG_limma_voom = DEG_DESeq2
names(DEG_limma_voom)
colnames(DEG_limma_voom)[2]='logFC'
colnames(DEG_limma_voom)[5]='P.Value'
DEG_limma_voom=DEG_limma_voom[order(DEG_limma_voom$logFC),]
head(DEG_limma_voom)

## set threshold to filter up and down genes
## in general, we choose 1 as the logFC threshold, 0.05 as the pvalue threshold
cut_logFC <- with(DEG_limma_voom,mean(abs(DEG_limma_voom$logFC)) + 2*sd(abs(DEG_limma_voom$logFC)) )
cut_logFC

if (cut_logFC > 1) {
  logFC_t <- 1
}

if (cut_logFC < 1) {
  logFC_t <- cut_logFC
}
logFC_t

pvalue_t <- 0.05


deg2_up <- DEG_limma_voom[with(DEG_limma_voom,logFC > logFC_t & P.Value < pvalue_t),]
deg2_up <- deg2_up[order(deg2_up$logFC,decreasing = T),]
deg2_up1 <- rownames(deg2_up)

deg2_down <- DEG_limma_voom[with(DEG_limma_voom,logFC < -logFC_t & P.Value < pvalue_t),]
deg2_down <- deg2_down[order(deg2_down$logFC,decreasing = F),]
deg2_down1 <- rownames(deg2_down)

# 1.3 load data AD1----
load(file = path_AD1)
head(DEG_DESeq2)

# 1. write deg result to csv----
DEG_limma_voom = DEG_DESeq2
names(DEG_limma_voom)
colnames(DEG_limma_voom)[2]='logFC'
colnames(DEG_limma_voom)[5]='P.Value'
DEG_limma_voom=DEG_limma_voom[order(DEG_limma_voom$logFC),]
head(DEG_limma_voom)

## set threshold to filter up and down genes
## in general, we choose 1 as the logFC threshold, 0.05 as the pvalue threshold
cut_logFC <- with(DEG_limma_voom,mean(abs(DEG_limma_voom$logFC)) + 2*sd(abs(DEG_limma_voom$logFC)) )
cut_logFC

if (cut_logFC > 1) {
  logFC_t <- 1
}

if (cut_logFC < 1) {
  logFC_t <- cut_logFC
}
logFC_t

pvalue_t <- 0.05


deg3_up <- DEG_limma_voom[with(DEG_limma_voom,logFC > logFC_t & P.Value < pvalue_t),]
deg3_up <- deg3_up[order(deg3_up$logFC,decreasing = T),]
deg3_up1 <- rownames(deg3_up)

deg3_down <- DEG_limma_voom[with(DEG_limma_voom,logFC < -logFC_t & P.Value < pvalue_t),]
deg3_down <- deg3_down[order(deg3_down$logFC,decreasing = F),]
deg3_down1 <- rownames(deg3_down)

# 1.4 load data AD2----
load(file = path_AD2)
head(DEG_DESeq2)

# 1. write deg result to csv----
DEG_limma_voom = DEG_DESeq2
names(DEG_limma_voom)
colnames(DEG_limma_voom)[2]='logFC'
colnames(DEG_limma_voom)[5]='P.Value'
DEG_limma_voom=DEG_limma_voom[order(DEG_limma_voom$logFC),]
head(DEG_limma_voom)

## set threshold to filter up and down genes
## in general, we choose 1 as the logFC threshold, 0.05 as the pvalue threshold
cut_logFC <- with(DEG_limma_voom,mean(abs(DEG_limma_voom$logFC)) + 2*sd(abs(DEG_limma_voom$logFC)) )
cut_logFC

if (cut_logFC > 1) {
  logFC_t <- 1
}

if (cut_logFC < 1) {
  logFC_t <- cut_logFC
}
logFC_t

pvalue_t <- 0.05


deg4_up <- DEG_limma_voom[with(DEG_limma_voom,logFC > logFC_t & P.Value < pvalue_t),]
deg4_up <- deg4_up[order(deg4_up$logFC,decreasing = T),]
deg4_up1 <- rownames(deg4_up)

deg4_down <- DEG_limma_voom[with(DEG_limma_voom,logFC < -logFC_t & P.Value < pvalue_t),]
deg4_down <- deg4_down[order(deg4_down$logFC,decreasing = F),]
deg4_down1 <- rownames(deg4_down)


# 2 RRA ---------------------------------------------------------------------


library(RobustRankAggreg)

library(clusterProfiler)
set.seed(123456789)

glist=list(deg1_up1,deg2_up1,deg3_up1,deg4_up1)

ups=aggregateRanks(glist = glist, N = length(unique(unlist(glist))))

tmp=as.data.frame(table(unlist(glist)))
ups$Freq=tmp[match(ups$Name,tmp[,1]),2]
head(ups)


gs=ups[ups$Score < 0.05 ,1]
gs=ups[ups$Score < 0.05 & ups$Freq>1,1]
gs
updat=data.frame(deg1=deg1_up[gs,'logFC'],
                 deg2=deg2_up[gs,'logFC'],
                 deg3=deg3_up[gs,'logFC'],
                 deg4=deg4_up[gs,'logFC'])
rownames(updat)=gs

updat <- updat[,-3]
updat <- na.omit(updat)

library(pheatmap)
head(updat)
colnames(updat)
colnames(updat) <- c("AR1","AA1","AD2")
pheatmap(updat,display_numbers=T)


glist=list(deg1_down1,deg2_down1,deg3_down1,deg4_down1)

downs=aggregateRanks(glist = glist, N = length(unique(unlist(glist))))

tmp=as.data.frame(table(unlist(glist)))
downs$Freq=tmp[match(downs$Name,tmp[,1]),2]
head(downs)


gs=downs[downs$Score < 0.05,1]
gs=downs[downs$Score < 0.05 & downs$Freq > 1,1]
gs=downs[downs$Score < 0.05 & downs$Freq > 2,1]

gs
downdat=data.frame(deg1=deg1_down[gs,'logFC'],
                 deg2=deg2_down[gs,'logFC'],
                 deg3=deg3_down[gs,'logFC'],
                 deg4=deg4_down[gs,'logFC'])
rownames(downdat)=gs

downdat[1,3]=downdat[1,4]
downdat <- downdat[,-4]

library(pheatmap)
head(downdat)
pheatmap(downdat,display_numbers=T)

colnames(downdat)
colnames(downdat) <- c("AR","AA","AD")
pheatmap(downdat,display_numbers=T)




