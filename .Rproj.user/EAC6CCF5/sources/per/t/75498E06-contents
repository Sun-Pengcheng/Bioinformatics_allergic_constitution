rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

all_path <- read.table("./output/meta_path_validation.txt",sep = '\t',header = T)

tizhi_gwas <- str_split(all_path$Input.a,"\\|") %>% unlist() %>% unique()
tizhi_tran <- str_split(all_path$Input.b,"\\|") %>% unlist() %>% unique()
tizhi_gmk <- str_split(all_path$Input.c,"\\|") %>% unlist() %>% unique()
tizhi_pro <- str_split(all_path$Input.d,"\\|") %>% unlist() %>% unique()

tizhi_gsea <- str_split(all_path$core_enrichment,"\\/") %>% unlist() %>% unique()

tizhi <- c(tizhi_gwas,tizhi_tran,tizhi_gmk,tizhi_pro) %>% unique()


tizhi_genelist <- bitr(tizhi,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")






#1 MALA CARDS --------------------------------------------------------------

ar_mala <- read.table("../../target/MALA/AR.txt",sep = '\t',header = F)
ar_mala <- ar_mala$V1 %>% unique()

ar <- bitr(ar_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

ra_mala <- read.table("../../target/MALA/RA.txt",sep = '\t',header = F)
ra_mala <- ra_mala$V1 %>% unique()

ra <- bitr(ra_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

aa_mala <- read.table("../../target/MALA/AA.txt",sep = '\t',header = F)
aa_mala <- aa_mala$V1 %>% unique()
aa <- bitr(aa_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

ad_mala <- read.table("../../target/MALA/AD.txt",sep = '\t',header = F)
ad_mala <- ad_mala$V1 %>% unique()
ad <- bitr(ad_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

ac_mala <- read.table("../../target/MALA/AC.txt",sep = '\t',header = F)
ac_mala <- ac_mala$V1 %>% unique()
ac <- bitr(ac_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

acd_mala <- read.table("../../target/MALA/ACD.txt",sep = '\t',header = F)
acd_mala <- acd_mala$V1 %>% unique()
acd <- bitr(acd_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

au_mala <- read.table("../../target/MALA/AU.txt",sep = '\t',header = F)
au_mala <- au_mala$V1 %>% unique()
au <- bitr(au_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

hv_mala <- read.table("../../target/MALA/HV.txt",sep = '\t',header = F)
hv_mala <- hv_mala$V1 %>% unique()
hv <- bitr(hv_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

eg_mala <- read.table("../../target/MALA/EG.txt",sep = '\t',header = F)
eg_mala <- eg_mala$V1 %>% unique()
eg <- bitr(eg_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

fa_mala <- read.table("../../target/MALA/FA.txt",sep = '\t',header = F)
fa_mala <- fa_mala$V1 %>% unique()
fa <- bitr(fa_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

da_mala <- read.table("../../target/MALA/DA.txt",sep = '\t',header = F)
da_mala <- da_mala$V1 %>% unique()
da <- bitr(da_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

pa_mala <- read.table("../../target/MALA/PA.txt",sep = '\t',header = F)
pa_mala <- pa_mala$V1 %>% unique()
pa <- bitr(pa_mala,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = "org.Hs.eg.db")

ige_mala <- read.table("../../target/MALA/IgE.txt",sep = '\t',header = F)
ige_mala <- ige_mala$V1 %>% unique()
ige <- bitr(ige_mala,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")




TIZHI <- data.frame(Term=rep("Allergic Constitution (Tizhi)",length(tizhi)),GeneSymbol=tizhi)

AR <- data.frame(Term=rep("Allergic Rhinitis",length(ar_mala)),GeneSymbol=ar_mala)
AA <- data.frame(Term=rep("Allergic Asthma",length(aa_mala)),GeneSymbol=aa_mala)
AD <- data.frame(Term=rep("Atopic Dermatitis",length(ad_mala)),GeneSymbol=ad_mala)
AC <- data.frame(Term=rep("Allergic Conjunctivitis",length(ac_mala)),GeneSymbol=ac_mala)
HV <- data.frame(Term=rep("Hypersensitivity Vasculitis",length(hv_mala)),GeneSymbol=hv_mala)
EG <- data.frame(Term=rep("Eosinophilic Gastroenteritis",length(eg_mala)),GeneSymbol=eg_mala)
FA <- data.frame(Term=rep("Food Allergy",length(fa_mala)),GeneSymbol=fa_mala)
DA <- data.frame(Term=rep("Drug Allergy",length(da_mala)),GeneSymbol=da_mala)
PA <- data.frame(Term=rep("Pollen Allergy",length(pa_mala)),GeneSymbol=pa_mala)
IRA <- data.frame(Term=rep("Ige Responsiveness, Atopic",length(ige_mala)),GeneSymbol=ige_mala)

AU <- data.frame(Term=rep("Allergic Urticaria",length(au_mala)),GeneSymbol=au_mala)
RA <- data.frame(Term=rep("Respiratory Allergy",length(ra_mala)),GeneSymbol=ra_mala)
ACD <- data.frame(Term=rep("Allergic Contact Dermatitis",length(acd_mala)),GeneSymbol=acd_mala)
MiA <- data.frame(Term=rep("Milk Allergy",length(ma_mala)),GeneSymbol=ma_mala)
MeA <- data.frame(Term=rep("Mental Allergy",length(metal_mala)),GeneSymbol=metal_mala)
EA <- data.frame(Term=rep("Egg Allergy",length(ea_mala)),GeneSymbol=ea_mala)
RA <- data.frame(Term=rep("Respiratory Allergy",length(ra_mala)),GeneSymbol=ra_mala)
PeA <- data.frame(Term=rep("Peanut Allergy",length(Peanut_mala)),GeneSymbol=Peanut_mala)
WA <- data.frame(Term=rep("Wheat Allergy",length(wh_mala)),GeneSymbol=wh_mala)
FishA <- data.frame(Term=rep("Fish Allergy",length(fish_mala)),GeneSymbol=fish_mala)
LA <- data.frame(Term=rep("Latex Allergy",length(la_mala)),GeneSymbol=la_mala)
PeA <- data.frame(Term=rep("Penicillin Allergy",length(pe_mala)),GeneSymbol=pe_mala)


ALL <- rbind(TIZHI,AR,AA,AD,AC,HV,EG,FA,DA,PA,IRA,AU,RA,ACD,MiA,MeA,EA,RA,PeA,WA,FishA,LA,PeA)

ALL <- rbind(TIZHI,AR,AA,AD,AC,HV,EG,FA,DA,PA,IRA,AU,RA,ACD,RA)


library(readxl)
CH <- read_xlsx("/mnt/d/Onedrive/02_数据/过敏/09_过敏康数据/HIT数据库/H0650-target.xlsx")
CH <- CH$`Gene Symbol`
CH <- CH %>% str_split(.,";") %>% unlist() %>% unique()
CH_id <- bitr(CH,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")

CH_PATH <- enricher(CH,
                    TERM2GENE=ALL,
                    pvalueCutoff = 1)

res_CH <- CH_PATH@result

YCH <- read_xlsx("/mnt/d/Onedrive/02_数据/过敏/09_过敏康数据/HIT数据库/H0744-target.xlsx")
YCH <- YCH$`Gene Symbol`
YCH <- YCH %>% str_split(.,";") %>% unlist() %>% unique()

YCH_PATH <- enricher(YCH,
                    TERM2GENE=ALL,
                    pvalueCutoff = 1)

res_YCH <- YCH_PATH@result

WM <- read_xlsx("/mnt/d/Onedrive/02_数据/过敏/09_过敏康数据/HIT数据库/H0601-target.xlsx")
WM <- WM$`Gene Symbol`
WM <- WM %>% str_split(.,";") %>% unlist() %>% unique()

WM_PATH <- enricher(WM,
                    TERM2GENE=ALL,
                    pvalueCutoff = 1)

res_WM <- WM_PATH@result

WM <- read_xlsx("/mnt/d/Onedrive/02_数据/过敏/09_过敏康数据/HIT数据库/H0601-target.xlsx")
WM <- WM$`Gene Symbol`
WM <- WM %>% str_split(.,";") %>% unlist() %>% unique()

WM_PATH <- enricher(WM,
                    TERM2GENE=ALL,
                    pvalueCutoff = 1)

res_WM <- WM_PATH@result


JXT <- read_xlsx("/mnt/d/Onedrive/02_数据/过敏/09_过敏康数据/HIT数据库/H0017-target.xlsx")
JXT <- JXT$`Gene Symbol`
JXT <- JXT %>% str_split(.,";") %>% unlist() %>% unique()

JXT_PATH <- enricher(JXT,
                    TERM2GENE=ALL,
                    pvalueCutoff = 1)

res_JXT <- JXT_PATH@result



# write.table(res,file = "./result/CEZS_final_enrich.txt",sep = '\t',col.names = T,row.names = F,quote = F)


GMK_PATH <- enricher(GMK,
                     TERM2GENE=ALL,
                     pvalueCutoff = 1)
res2 <- GMK_PATH@result
write.table(res2,file = "./result/GMK_final_enrich.txt",sep = '\t',col.names = T,row.names = F,quote = F)
