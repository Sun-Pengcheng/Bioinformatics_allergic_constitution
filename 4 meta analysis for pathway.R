library(tidyverse)

#----load data----
rm(list=ls())

gwas <- read.table("./rdata/Pathway/output_identify_20230903211440_GWAS.txt",
                   sep = '\t', skip = 4, comment.char = "", header = T, stringsAsFactors = F,
                   quote = "")
gwas <- gwas[gwas$Input.number>=2 & gwas$Background.number<200,]
gwas <- gwas[gwas$P.Value<0.05,]
gwas <- gwas[gwas$Database!="BioCyc" & gwas$Database!="PANTHER",]
table(gwas$Database)

rra <- read.table("./rdata/Pathway/output_identify_20230903211611_RRA.txt",
                  sep = '\t', skip = 4, comment.char = "", header = T, stringsAsFactors = F,
                  quote = "")
rra <- rra[rra$Input.number>=2 & rra$Background.number<200,]
rra <- rra[rra$P.Value<0.05,]
rra <- rra[rra$Database!="BioCyc" & rra$Database!="PANTHER",]
table(rra$Database)

blue <- read.table("./rdata/Pathway/output_identify_20230903212000_blue.txt",
                    sep = '\t', skip = 4, comment.char = "", header = T, stringsAsFactors = F,
                    quote = "")
blue <- blue[blue$Input.number>=2 & blue$Background.number<200,]
blue <- blue[blue$P.Value<0.05,]
blue <- blue[blue$Database!="BioCyc" & blue$Database!="PANTHER",]
table(blue$Database)

turquoise <- read.table("./rdata/Pathway/output_identify_20230903212205_turquoise.txt",
                         sep = '\t', skip = 4, comment.char = "", header = T, stringsAsFactors = F,
                         quote = "")
turquoise <- turquoise[turquoise$Input.number>=2 & turquoise$Background.number<200,]
turquoise <- turquoise[turquoise$P.Value<0.05,]
turquoise <- turquoise[turquoise$Database!="BioCyc" & turquoise$Database!="PANTHER",]
table(turquoise$Database)

dep <- read.table("./rdata/Pathway/output_identify_20230903212614_PROTEIN.txt",
                  sep = '\t', skip = 4, comment.char = "", header = T, stringsAsFactors = F,
                  quote = "")
dep <- dep[dep$Input.number>=2 & dep$Background.number<200,]
dep <- dep[dep$P.Value<0.05,]
dep <- dep[dep$Database!="BioCyc" & dep$Database!="PANTHER",]
table(dep$Database)

gmk <- read.table("./rdata/Pathway/output_identify_20230903131859_GMK_107.txt",
                  sep = '\t', skip = 4, comment.char = "", header = T, stringsAsFactors = F,
                  quote = "")
gmk <- gmk[gmk$Input.number>=2 & gmk$Background.number<200,]
gmk <- gmk[gmk$P.Value<0.05,]
gmk <- gmk[gmk$Database!="BioCyc" & gmk$Database!="PANTHER",]
table(gmk$Database)


# meta1 GWAS+RRA+DEP+GMK ---------------------------------------------------------------------

a <- gwas
b <- rra
c <- gmk
d <- dep

colnames(d)
d <- d[,c(3,4,6,7,8)]
colnames(d)
colnames(d) <- c("ID","Input.number.d","P.Value.d","Corrected.P.Value.d","Input.d")

tmp <- merge(a,b,by="ID")
tmp2 <- merge(tmp, c, by="ID")
tmp2 <- merge(tmp2,d,by="ID")

tmp2$Pfisher <- pchisq(log(tmp2$P.Value.x * tmp2$P.Value.y * tmp2$P.Value * tmp2$P.Value.d)*(-2),
                       df=8, lower.tail = F)
tmp2$Pfisher_bon <- tmp2$Pfisher * (18437+333+2226)
tmp2 <- tmp2[tmp2$Pfisher_bon<0.05,]

colnames(tmp2)
tmp2 <- tmp2[,c("ID","X.Term.x","Database.x","Input.x","Input.y","Input","Input.d","Pfisher","Pfisher_bon")]
colnames(tmp2) <- c("ID","Term","Database","Input.a","Input.b","Input.c","Input.d","Pfisher","Pfisher_bon")

table(tmp2$Database)
tmp2$Database <- factor(tmp2$Database,levels = c("KEGG PATHWAY","Reactome","BioCyc","PANTHER","Gene Ontology"))
tmp2 <- tmp2[order(tmp2$Database,tmp2$Pfisher),]

tmp2$Source <- rep("RRA",nrow(tmp2))

all_path <- data.frame()
all_path <- rbind(all_path,tmp2)

write.table(tmp2,file = './rdata/Pathway/meta/meta1_rra.txt',sep = '\t',col.names = T,row.names = F,quote = F)


# meta2 GWAS+blue+DEP+GMK ---------------------------------------------------------------------

rm(a,b,c,tmp,tmp2)
a <- gwas
b <- blue
c <- gmk

tmp <- merge(a,b,by="ID")
tmp2 <- merge(tmp, c, by="ID")

tmp2 <- merge(tmp2,d,by="ID")

tmp2$Pfisher <- pchisq(log(tmp2$P.Value.x * tmp2$P.Value.y * tmp2$P.Value * tmp2$P.Value.d)*(-2),
                       df=8, lower.tail = F)
tmp2$Pfisher_bon <- tmp2$Pfisher * (18437+333+2226)
length(which(tmp2$Pfisher_bon<0.01))
tmp2 <- tmp2[tmp2$Pfisher_bon<0.05,]

colnames(tmp2)
tmp2 <- tmp2[,c("ID","X.Term.x","Database.x","Input.x","Input.y","Input","Input.d","Pfisher","Pfisher_bon")]
colnames(tmp2) <- c("ID","Term","Database","Input.a","Input.b","Input.c","Input.d","Pfisher","Pfisher_bon")

delete <- c("hsa05160","hsa05167","hsa05162")
tmp2 <- tmp2[!tmp2$ID%in%delete,]
table(tmp2$Database)
tmp2$Database <- factor(tmp2$Database,levels = c("KEGG PATHWAY","Reactome","BioCyc","PANTHER","Gene Ontology"))
tmp2 <- tmp2[order(tmp2$Database,tmp2$Pfisher),]

tmp2$Source <- rep("blue",nrow(tmp2))
all_path <- rbind(all_path,tmp2)

write.table(tmp2,file = './rdata/Pathway/meta/meta2_blue.txt',sep = '\t',col.names = T,row.names = F,quote = F)

# meta3 GWAS+turquoise+DEP+GMK ---------------------------------------------------------------------

rm(a,b,c,tmp,tmp2)
a <- gwas
b <- turquoise
c <- gmk

tmp <- merge(a,b,by="ID")
tmp2 <- merge(tmp, c, by="ID")
tmp2 <- merge(tmp2,d,by="ID")
tmp2$Pfisher <- pchisq(log(tmp2$P.Value.x * tmp2$P.Value.y * tmp2$P.Value * tmp2$P.Value.d)*(-2),
                       df=8, lower.tail = F)
tmp2$Pfisher_bon <- tmp2$Pfisher * (18437+333+2226)
tmp2 <- tmp2[tmp2$Pfisher_bon<0.05,]

colnames(tmp2)
tmp2 <- tmp2[,c("ID","X.Term.x","Database.x","Input.x","Input.y","Input","Input.d","Pfisher","Pfisher_bon")]
colnames(tmp2) <- c("ID","Term","Database","Input.a","Input.b","Input.c","Input.d","Pfisher","Pfisher_bon")

delete <- c("hsa05142","hsa05152","hsa05160","hsa05162","hsa05164","hsa05167")
tmp2 <- tmp2[!tmp2$ID%in%delete,]
table(tmp2$Database)
tmp2$Database <- factor(tmp2$Database,levels = c("KEGG PATHWAY","Reactome","BioCyc","PANTHER","Gene Ontology"))
tmp2 <- tmp2[order(tmp2$Database,tmp2$Pfisher),]

tmp2$Source <- rep("turquoise",nrow(tmp2))
all_path <- rbind(all_path,tmp2)

write.table(tmp2,file = './rdata/Pathway/meta/meta3_turquoise.txt',sep = '\t',col.names = T,row.names = F,quote = F)

# all path ----------------------------------------------------------------

length(unique(all_path$ID))
colnames(all_path)
write.table(all_path,file = "./output/meta_all_path.txt",sep = '\t',col.names = T,row.names = F,quote = F)

# reload
rm(list = ls())
all_path <- read.table("./output/meta_all_path.txt",sep = '\t',header = T)
table(all_path$Source)

# Validation -----------------------------------------------------------------
rm(list = ls())

GMK_path <- read.table("./rdata/Pathway/imDC_mDC_GSEA.txt",header = T,sep = '\t')

all <- read.table("./output/meta_all_path.txt",sep = '\t',header = T,quote = "")

all$id <- c(1:nrow(all))

all_validation <- merge(all,GMK_path,by="ID")

write.table(all_validation,"./output/meta_path_validation.txt",sep = '\t',col.names = T,row.names = F,quote = F)

length(unique(all_validation$Term))
