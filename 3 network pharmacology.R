rm(list = ls())
library(tidyverse)

load("./rdata/NP/GMK_4herb.Rda")

# GMK swiss target ------------------------------------------------------------

GMK_swiss_filter <- GMK_swiss[GMK_swiss$Probability.>0.1,]
length(GMK_swiss_filter$Uniprot.ID %>% unique())

GMK_swiss_sig <- GMK_swiss[GMK_swiss$Probability.>=0.9,]

# write.table(GMK_swiss$Uniprot.ID %>% unique(),
#             file = "/rdata/NP/Swiss_target_all.txt",
#             col.names = F,row.names = F,quote = F)

swiss_target <- read.table("./rdata/NP/Swiss_uniprot-2022.10.27.tsv",sep = '\t',
                           header = T, quote = "")
swiss_target <- swiss_target[swiss_target$From%in%GMK_swiss_filter$Uniprot.ID,]
setdiff(GMK_swiss_sig$Uniprot.ID,swiss_target$From)

swiss_gene <- swiss_target$To %>% unique()
swiss_gene_sig <- swiss_target[swiss_target$From%in%GMK_swiss_sig$Uniprot.ID,]$To %>% unique()

# GMK pharmMapper target ---------------------------------------------------------

GMK_pharm <- GMK_pharm[GMK_pharm$zscore>0,]
length(GMK_pharm$Uniplot %>% unique() %>% na.omit())

GMK_pharm_sig <- GMK_pharm[GMK_pharm$zscore>4 & GMK_pharm$Norm.Fit>0.4,]
length(GMK_pharm_sig$Uniplot %>% unique())

# write.table(GMK_pharm$Uniplot %>% unique() %>% na.omit(),
#             file = "./rdata/NP/GMK_pharm_zover0.txt",
#             col.names = F, row.names = F, quote = F)

pharm_target <- read.table("./rdata/NP/Pharm_uniprot-2022.10.27.tsv",sep = '\t',
                           header = T, quote = "")
pharm_target <- pharm_target[pharm_target$From%in%GMK_pharm$Uniplot,]
pharm_gene <- pharm_target$To %>% unique()
setdiff(GMK_pharm_sig$Uniplot,pharm_target$From)
pharm_gene_sig <- pharm_target[pharm_target$From%in%GMK_pharm_sig$Uniplot,]$To


# HERB target --------------------------------------------------------------------

herb_paper <- readxl::read_xlsx("./rdata/NP/HERB_paper_target_2022_10_18.xlsx")
herb_paper <- herb_paper$`Target name` %>% unique()

herb_base <- read.table("./rdata/NP/HERB_target_2022_10_18.txt", sep = '\t', header = T)
max(herb_base$FDR_BH)
min(herb_base$FDR_BH)

herb_gene <- c(herb_paper,herb_base$Target.name) %>% unique()


# SymMap ------------------------------------------------------------------

symap <- read.table("./rdata/NP/SymMap_target_2022_10_18.txt", sep = '\t', header = T, quote = "")
table(symap$herb)
length(unique(symap$Gene.symbol) %>% unique())
symap <- symap[symap$FDR.BH.<0.05,]
sym_gene <- symap$Gene.symbol %>% unique()

# Venn plot ---------------------------------------------------------------
library(VennDiagram)

A <- pharm_gene
B <- swiss_gene
C <- herb_gene
D <- sym_gene

venn.plot <- venn.diagram(
  x = list(B,C,A,D),
  category.names = c("Swiss", "HERB","PharmMapper",  "SymMap"),
  filename = "./figure/VennPlot_target_GMK_4herb.tiff",    
  euler.d = F,scaled = F,
  lty = "blank",
  fill = c("#FEB26D", "purple", "lightyellow", "greenyellow"),  
  alpha = 0.50,     
  cex = 1.2,   
  cat.cex = 1.2,   
  cat.fontfamily = "serif"
)

#两两
a1 <- intersect(A,B)
a2 <- intersect(A,C)
a3 <- intersect(A,D)
a4 <- intersect(B,C)
a5 <- intersect(B,D)
a6 <- intersect(C,D)

all <- c(a1,a2,a3,a4,a5,a6) %>% unique()
rm(a1,a2,a3,a4,a5,a6)

all_suppl <- c(pharm_gene_sig,swiss_gene_sig,herb_paper) %>% unique()
all_all <- c(all,all_suppl) %>% unique()
GMK_gene <- all_all


# validation --------------------------------------------------------------------

GMK_deg <- read.table("./rdata/NP/GMK_intervention_DEG.txt", sep = '\t', header = T)
im_deg <- GMK_deg[GMK_deg$type=="im_GMK",]$gene
m_deg <- GMK_deg[GMK_deg$type=="m_GMK",]$gene

GMK_gene_validation <- intersect(GMK_gene, union(im_deg,m_deg))

write.table(GMK_gene_validation, file = "./output/GMK_gene_4herb.txt",sep = '\t',col.names = F,row.names = F,quote = F)


