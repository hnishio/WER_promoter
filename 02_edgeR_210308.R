
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)
library(VennDiagram)

setwd("data")


##### Extraction of genes with seasonal difference #####

# H3K4me3
K4_TSSup1000.rep1 <- read.csv("K4_TSSup1000_rep1.csv",header=T,sep=",")
K4_TSSup1000.rep2 <- read.csv("K4_TSSup1000_rep1.csv",header=T,sep=",")

colnames(K4_TSSup1000.rep1) <- c("K4_TSSup1000_chr","K4_TSSup1000_start","K4_TSSup1000_end", 
                       "K4_TSSup1000_rep1_11.6","K4_TSSup1000_rep1_12.4","K4_TSSup1000_rep1_1.8","K4_TSSup1000_rep1_2.5",
                       "K4_TSSup1000_rep1_3.5","K4_TSSup1000_rep1_4.2","K4_TSSup1000_rep1_4.30","K4_TSSup1000_rep1_5.28",
                       "K4_TSSup1000_rep1_7.2","K4_TSSup1000_rep1_7.30","K4_TSSup1000_rep1_8.27","K4_TSSup1000_rep1_9.24")
colnames(K4_TSSup1000.rep2) <- c("K4_TSSup1000_chr","K4_TSSup1000_start", "K4_TSSup1000_end", 
                       "K4_TSSup1000_rep2_11.6","K4_TSSup1000_rep2_12.4","K4_TSSup1000_rep2_1.8","K4_TSSup1000_rep2_2.5",
                       "K4_TSSup1000_rep2_3.5","K4_TSSup1000_rep2_4.2","K4_TSSup1000_rep2_4.30","K4_TSSup1000_rep2_5.28",
                       "K4_TSSup1000_rep2_7.2","K4_TSSup1000_rep2_7.30","K4_TSSup1000_rep2_8.27","K4_TSSup1000_rep2_9.24")
K4_TSSup1000.rep1.rep2 <- merge(K4_TSSup1000.rep1, K4_TSSup1000.rep2, 
                      by=c("K4_TSSup1000_chr","K4_TSSup1000_start","K4_TSSup1000_end"), 
                      all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup1000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K4_TSSup1000_chr","K4_TSSup1000_start","K4_TSSup1000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K4_TSSup1000.rep1.rep2.gene <- merge(K4_TSSup1000.rep1.rep2, gene, by=c("K4_TSSup1000_chr","K4_TSSup1000_start","K4_TSSup1000_end"), 
                           all=F, sort=F)
K4_TSSup1000.rep1.rep2.gene <- K4_TSSup1000.rep1.rep2.gene[!duplicated(paste(K4_TSSup1000.rep1.rep2.gene$K4_TSSup1000_chr,K4_TSSup1000.rep1.rep2.gene$K4_TSSup1000_start,
                                                         K4_TSSup1000.rep1.rep2.gene$K4_TSSup1000_end,K4_TSSup1000.rep1.rep2.gene$Ahal_ID),sep=","),]
K4_TSSup1000.rep1.rep2.gene <- K4_TSSup1000.rep1.rep2.gene[,-((ncol(K4_TSSup1000.rep1.rep2.gene)-5):ncol(K4_TSSup1000.rep1.rep2.gene))]

K4_TSSup1000.mat <- K4_TSSup1000.rep1.rep2.gene[,4:27]
row.names(K4_TSSup1000.mat) <- K4_TSSup1000.rep1.rep2.gene[,28]
head(K4_TSSup1000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K4_TSSup1000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK4_TSSup1000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K4_TSSup1000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK4_TSSup1000_CPMover1_edgeR_result.csv", row.names=F)



### H3K27me3
K27_TSSup1000.rep1 <- read.csv("K27_TSSup1000_rep1.csv",header=T,sep=",")
K27_TSSup1000.rep2 <- read.csv("K27_TSSup1000_rep2.csv",header=T,sep=",")

colnames(K27_TSSup1000.rep1) <- c("K27_TSSup1000_chr","K27_TSSup1000_start","K27_TSSup1000_end", 
                       "K27_TSSup1000_rep1_11.6","K27_TSSup1000_rep1_12.4","K27_TSSup1000_rep1_1.8","K27_TSSup1000_rep1_2.5",
                       "K27_TSSup1000_rep1_3.5","K27_TSSup1000_rep1_4.2","K27_TSSup1000_rep1_4.30","K27_TSSup1000_rep1_5.28",
                       "K27_TSSup1000_rep1_7.2","K27_TSSup1000_rep1_7.30","K27_TSSup1000_rep1_8.27","K27_TSSup1000_rep1_9.24")
colnames(K27_TSSup1000.rep2) <- c("K27_TSSup1000_chr","K27_TSSup1000_start", "K27_TSSup1000_end", 
                       "K27_TSSup1000_rep2_11.6","K27_TSSup1000_rep2_12.4","K27_TSSup1000_rep2_1.8","K27_TSSup1000_rep2_2.5",
                       "K27_TSSup1000_rep2_3.5","K27_TSSup1000_rep2_4.2","K27_TSSup1000_rep2_4.30","K27_TSSup1000_rep2_5.28",
                       "K27_TSSup1000_rep2_7.2","K27_TSSup1000_rep2_7.30","K27_TSSup1000_rep2_8.27","K27_TSSup1000_rep2_9.24")
K27_TSSup1000.rep1.rep2 <- merge(K27_TSSup1000.rep1, K27_TSSup1000.rep2, 
                      by=c("K27_TSSup1000_chr","K27_TSSup1000_start","K27_TSSup1000_end"), 
                      all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup1000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K27_TSSup1000_chr","K27_TSSup1000_start","K27_TSSup1000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K27_TSSup1000.rep1.rep2.gene <- merge(K27_TSSup1000.rep1.rep2, gene, by=c("K27_TSSup1000_chr","K27_TSSup1000_start","K27_TSSup1000_end"), 
                           all=F, sort=F)
K27_TSSup1000.rep1.rep2.gene <- K27_TSSup1000.rep1.rep2.gene[!duplicated(paste(K27_TSSup1000.rep1.rep2.gene$K27_TSSup1000_chr,K27_TSSup1000.rep1.rep2.gene$K27_TSSup1000_start,
                                                         K27_TSSup1000.rep1.rep2.gene$K27_TSSup1000_end,K27_TSSup1000.rep1.rep2.gene$Ahal_ID),sep=","),]
K27_TSSup1000.rep1.rep2.gene <- K27_TSSup1000.rep1.rep2.gene[,-((ncol(K27_TSSup1000.rep1.rep2.gene)-5):ncol(K27_TSSup1000.rep1.rep2.gene))]

K27_TSSup1000.mat <- K27_TSSup1000.rep1.rep2.gene[,4:27]
row.names(K27_TSSup1000.mat) <- K27_TSSup1000.rep1.rep2.gene[,28]
head(K27_TSSup1000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K27_TSSup1000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK27_TSSup1000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K27_TSSup1000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK27_TSSup1000_CPMover1_edgeR_result.csv", row.names=F)



### Venn diagram
dir.create("../figs/edgeR")
K4_TSSup1000.mat <- read.csv("seasonalK4_TSSup1000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27_TSSup1000.mat <- read.csv("seasonalK27_TSSup1000_CPMover1_edgeR_result.csv",header=T,sep=",")

K4_TSSup1000_season_osci <- K4_TSSup1000.mat[which(K4_TSSup1000.mat$FDR<0.05),]
K27_TSSup1000_season_osci <- K27_TSSup1000.mat[which(K27_TSSup1000.mat$FDR<0.05),]
K4TSSup1000_K27TSSup1000_season_osci <- merge(K4_TSSup1000_season_osci,K27_TSSup1000_season_osci,by=c("Ahal_ID"),all=F,sort=F)


pdf("../figs/edgeR/venn_K4TSSup1000_K27TSSup1000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4_TSSup1000_season_osci), area2=nrow(K27_TSSup1000_season_osci), cross.area=nrow(K4TSSup1000_K27TSSup1000_season_osci), 
  cex=1, category=c("K4_TSSup1000 seasonal","K27_TSSup1000 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans", 
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()





##### Extraction of genes with seasonal difference #####

# H3K4me3
K4_TSSup2000.rep1 <- read.csv("K4_TSSup2000_rep1.csv",header=T,sep=",")
K4_TSSup2000.rep2 <- read.csv("K4_TSSup2000_rep1.csv",header=T,sep=",")

colnames(K4_TSSup2000.rep1) <- c("K4_TSSup2000_chr","K4_TSSup2000_start","K4_TSSup2000_end", 
                                 "K4_TSSup2000_rep1_11.6","K4_TSSup2000_rep1_12.4","K4_TSSup2000_rep1_1.8","K4_TSSup2000_rep1_2.5",
                                 "K4_TSSup2000_rep1_3.5","K4_TSSup2000_rep1_4.2","K4_TSSup2000_rep1_4.30","K4_TSSup2000_rep1_5.28",
                                 "K4_TSSup2000_rep1_7.2","K4_TSSup2000_rep1_7.30","K4_TSSup2000_rep1_8.27","K4_TSSup2000_rep1_9.24")
colnames(K4_TSSup2000.rep2) <- c("K4_TSSup2000_chr","K4_TSSup2000_start", "K4_TSSup2000_end", 
                                 "K4_TSSup2000_rep2_11.6","K4_TSSup2000_rep2_12.4","K4_TSSup2000_rep2_1.8","K4_TSSup2000_rep2_2.5",
                                 "K4_TSSup2000_rep2_3.5","K4_TSSup2000_rep2_4.2","K4_TSSup2000_rep2_4.30","K4_TSSup2000_rep2_5.28",
                                 "K4_TSSup2000_rep2_7.2","K4_TSSup2000_rep2_7.30","K4_TSSup2000_rep2_8.27","K4_TSSup2000_rep2_9.24")
K4_TSSup2000.rep1.rep2 <- merge(K4_TSSup2000.rep1, K4_TSSup2000.rep2, 
                                by=c("K4_TSSup2000_chr","K4_TSSup2000_start","K4_TSSup2000_end"), 
                                all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup2000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K4_TSSup2000_chr","K4_TSSup2000_start","K4_TSSup2000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K4_TSSup2000.rep1.rep2.gene <- merge(K4_TSSup2000.rep1.rep2, gene, by=c("K4_TSSup2000_chr","K4_TSSup2000_start","K4_TSSup2000_end"), 
                                     all=F, sort=F)
K4_TSSup2000.rep1.rep2.gene <- K4_TSSup2000.rep1.rep2.gene[!duplicated(paste(K4_TSSup2000.rep1.rep2.gene$K4_TSSup2000_chr,K4_TSSup2000.rep1.rep2.gene$K4_TSSup2000_start,
                                                                             K4_TSSup2000.rep1.rep2.gene$K4_TSSup2000_end,K4_TSSup2000.rep1.rep2.gene$Ahal_ID),sep=","),]
K4_TSSup2000.rep1.rep2.gene <- K4_TSSup2000.rep1.rep2.gene[,-((ncol(K4_TSSup2000.rep1.rep2.gene)-5):ncol(K4_TSSup2000.rep1.rep2.gene))]

K4_TSSup2000.mat <- K4_TSSup2000.rep1.rep2.gene[,4:27]
row.names(K4_TSSup2000.mat) <- K4_TSSup2000.rep1.rep2.gene[,28]
head(K4_TSSup2000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K4_TSSup2000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK4_TSSup2000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K4_TSSup2000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK4_TSSup2000_CPMover1_edgeR_result.csv", row.names=F)



### H3K27me3
K27_TSSup2000.rep1 <- read.csv("K27_TSSup2000_rep1.csv",header=T,sep=",")
K27_TSSup2000.rep2 <- read.csv("K27_TSSup2000_rep2.csv",header=T,sep=",")

colnames(K27_TSSup2000.rep1) <- c("K27_TSSup2000_chr","K27_TSSup2000_start","K27_TSSup2000_end", 
                                  "K27_TSSup2000_rep1_11.6","K27_TSSup2000_rep1_12.4","K27_TSSup2000_rep1_1.8","K27_TSSup2000_rep1_2.5",
                                  "K27_TSSup2000_rep1_3.5","K27_TSSup2000_rep1_4.2","K27_TSSup2000_rep1_4.30","K27_TSSup2000_rep1_5.28",
                                  "K27_TSSup2000_rep1_7.2","K27_TSSup2000_rep1_7.30","K27_TSSup2000_rep1_8.27","K27_TSSup2000_rep1_9.24")
colnames(K27_TSSup2000.rep2) <- c("K27_TSSup2000_chr","K27_TSSup2000_start", "K27_TSSup2000_end", 
                                  "K27_TSSup2000_rep2_11.6","K27_TSSup2000_rep2_12.4","K27_TSSup2000_rep2_1.8","K27_TSSup2000_rep2_2.5",
                                  "K27_TSSup2000_rep2_3.5","K27_TSSup2000_rep2_4.2","K27_TSSup2000_rep2_4.30","K27_TSSup2000_rep2_5.28",
                                  "K27_TSSup2000_rep2_7.2","K27_TSSup2000_rep2_7.30","K27_TSSup2000_rep2_8.27","K27_TSSup2000_rep2_9.24")
K27_TSSup2000.rep1.rep2 <- merge(K27_TSSup2000.rep1, K27_TSSup2000.rep2, 
                                 by=c("K27_TSSup2000_chr","K27_TSSup2000_start","K27_TSSup2000_end"), 
                                 all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup2000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K27_TSSup2000_chr","K27_TSSup2000_start","K27_TSSup2000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K27_TSSup2000.rep1.rep2.gene <- merge(K27_TSSup2000.rep1.rep2, gene, by=c("K27_TSSup2000_chr","K27_TSSup2000_start","K27_TSSup2000_end"), 
                                      all=F, sort=F)
K27_TSSup2000.rep1.rep2.gene <- K27_TSSup2000.rep1.rep2.gene[!duplicated(paste(K27_TSSup2000.rep1.rep2.gene$K27_TSSup2000_chr,K27_TSSup2000.rep1.rep2.gene$K27_TSSup2000_start,
                                                                               K27_TSSup2000.rep1.rep2.gene$K27_TSSup2000_end,K27_TSSup2000.rep1.rep2.gene$Ahal_ID),sep=","),]
K27_TSSup2000.rep1.rep2.gene <- K27_TSSup2000.rep1.rep2.gene[,-((ncol(K27_TSSup2000.rep1.rep2.gene)-5):ncol(K27_TSSup2000.rep1.rep2.gene))]

K27_TSSup2000.mat <- K27_TSSup2000.rep1.rep2.gene[,4:27]
row.names(K27_TSSup2000.mat) <- K27_TSSup2000.rep1.rep2.gene[,28]
head(K27_TSSup2000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K27_TSSup2000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK27_TSSup2000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K27_TSSup2000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK27_TSSup2000_CPMover1_edgeR_result.csv", row.names=F)



### Venn diagram
dir.create("../figs/edgeR")
K4_TSSup2000.mat <- read.csv("seasonalK4_TSSup2000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27_TSSup2000.mat <- read.csv("seasonalK27_TSSup2000_CPMover1_edgeR_result.csv",header=T,sep=",")

K4_TSSup2000_season_osci <- K4_TSSup2000.mat[which(K4_TSSup2000.mat$FDR<0.05),]
K27_TSSup2000_season_osci <- K27_TSSup2000.mat[which(K27_TSSup2000.mat$FDR<0.05),]
K4TSSup2000_K27TSSup2000_season_osci <- merge(K4_TSSup2000_season_osci,K27_TSSup2000_season_osci,by=c("Ahal_ID"),all=F,sort=F)


pdf("../figs/edgeR/venn_K4TSSup2000_K27TSSup2000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4_TSSup2000_season_osci), area2=nrow(K27_TSSup2000_season_osci), cross.area=nrow(K4TSSup2000_K27TSSup2000_season_osci), 
  cex=1, category=c("K4_TSSup2000 seasonal","K27_TSSup2000 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans", 
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()





##### Extraction of genes with seasonal difference #####

# H3K4me3
K4_TSSup3000.rep1 <- read.csv("K4_TSSup3000_rep1.csv",header=T,sep=",")
K4_TSSup3000.rep2 <- read.csv("K4_TSSup3000_rep1.csv",header=T,sep=",")

colnames(K4_TSSup3000.rep1) <- c("K4_TSSup3000_chr","K4_TSSup3000_start","K4_TSSup3000_end", 
                                 "K4_TSSup3000_rep1_11.6","K4_TSSup3000_rep1_12.4","K4_TSSup3000_rep1_1.8","K4_TSSup3000_rep1_2.5",
                                 "K4_TSSup3000_rep1_3.5","K4_TSSup3000_rep1_4.2","K4_TSSup3000_rep1_4.30","K4_TSSup3000_rep1_5.28",
                                 "K4_TSSup3000_rep1_7.2","K4_TSSup3000_rep1_7.30","K4_TSSup3000_rep1_8.27","K4_TSSup3000_rep1_9.24")
colnames(K4_TSSup3000.rep2) <- c("K4_TSSup3000_chr","K4_TSSup3000_start", "K4_TSSup3000_end", 
                                 "K4_TSSup3000_rep2_11.6","K4_TSSup3000_rep2_12.4","K4_TSSup3000_rep2_1.8","K4_TSSup3000_rep2_2.5",
                                 "K4_TSSup3000_rep2_3.5","K4_TSSup3000_rep2_4.2","K4_TSSup3000_rep2_4.30","K4_TSSup3000_rep2_5.28",
                                 "K4_TSSup3000_rep2_7.2","K4_TSSup3000_rep2_7.30","K4_TSSup3000_rep2_8.27","K4_TSSup3000_rep2_9.24")
K4_TSSup3000.rep1.rep2 <- merge(K4_TSSup3000.rep1, K4_TSSup3000.rep2, 
                                by=c("K4_TSSup3000_chr","K4_TSSup3000_start","K4_TSSup3000_end"), 
                                all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup3000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K4_TSSup3000_chr","K4_TSSup3000_start","K4_TSSup3000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K4_TSSup3000.rep1.rep2.gene <- merge(K4_TSSup3000.rep1.rep2, gene, by=c("K4_TSSup3000_chr","K4_TSSup3000_start","K4_TSSup3000_end"), 
                                     all=F, sort=F)
K4_TSSup3000.rep1.rep2.gene <- K4_TSSup3000.rep1.rep2.gene[!duplicated(paste(K4_TSSup3000.rep1.rep2.gene$K4_TSSup3000_chr,K4_TSSup3000.rep1.rep2.gene$K4_TSSup3000_start,
                                                                             K4_TSSup3000.rep1.rep2.gene$K4_TSSup3000_end,K4_TSSup3000.rep1.rep2.gene$Ahal_ID),sep=","),]
K4_TSSup3000.rep1.rep2.gene <- K4_TSSup3000.rep1.rep2.gene[,-((ncol(K4_TSSup3000.rep1.rep2.gene)-5):ncol(K4_TSSup3000.rep1.rep2.gene))]

K4_TSSup3000.mat <- K4_TSSup3000.rep1.rep2.gene[,4:27]
row.names(K4_TSSup3000.mat) <- K4_TSSup3000.rep1.rep2.gene[,28]
head(K4_TSSup3000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K4_TSSup3000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK4_TSSup3000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K4_TSSup3000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK4_TSSup3000_CPMover1_edgeR_result.csv", row.names=F)



### H3K27me3
K27_TSSup3000.rep1 <- read.csv("K27_TSSup3000_rep1.csv",header=T,sep=",")
K27_TSSup3000.rep2 <- read.csv("K27_TSSup3000_rep2.csv",header=T,sep=",")

colnames(K27_TSSup3000.rep1) <- c("K27_TSSup3000_chr","K27_TSSup3000_start","K27_TSSup3000_end", 
                                  "K27_TSSup3000_rep1_11.6","K27_TSSup3000_rep1_12.4","K27_TSSup3000_rep1_1.8","K27_TSSup3000_rep1_2.5",
                                  "K27_TSSup3000_rep1_3.5","K27_TSSup3000_rep1_4.2","K27_TSSup3000_rep1_4.30","K27_TSSup3000_rep1_5.28",
                                  "K27_TSSup3000_rep1_7.2","K27_TSSup3000_rep1_7.30","K27_TSSup3000_rep1_8.27","K27_TSSup3000_rep1_9.24")
colnames(K27_TSSup3000.rep2) <- c("K27_TSSup3000_chr","K27_TSSup3000_start", "K27_TSSup3000_end", 
                                  "K27_TSSup3000_rep2_11.6","K27_TSSup3000_rep2_12.4","K27_TSSup3000_rep2_1.8","K27_TSSup3000_rep2_2.5",
                                  "K27_TSSup3000_rep2_3.5","K27_TSSup3000_rep2_4.2","K27_TSSup3000_rep2_4.30","K27_TSSup3000_rep2_5.28",
                                  "K27_TSSup3000_rep2_7.2","K27_TSSup3000_rep2_7.30","K27_TSSup3000_rep2_8.27","K27_TSSup3000_rep2_9.24")
K27_TSSup3000.rep1.rep2 <- merge(K27_TSSup3000.rep1, K27_TSSup3000.rep2, 
                                 by=c("K27_TSSup3000_chr","K27_TSSup3000_start","K27_TSSup3000_end"), 
                                 all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup3000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K27_TSSup3000_chr","K27_TSSup3000_start","K27_TSSup3000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K27_TSSup3000.rep1.rep2.gene <- merge(K27_TSSup3000.rep1.rep2, gene, by=c("K27_TSSup3000_chr","K27_TSSup3000_start","K27_TSSup3000_end"), 
                                      all=F, sort=F)
K27_TSSup3000.rep1.rep2.gene <- K27_TSSup3000.rep1.rep2.gene[!duplicated(paste(K27_TSSup3000.rep1.rep2.gene$K27_TSSup3000_chr,K27_TSSup3000.rep1.rep2.gene$K27_TSSup3000_start,
                                                                               K27_TSSup3000.rep1.rep2.gene$K27_TSSup3000_end,K27_TSSup3000.rep1.rep2.gene$Ahal_ID),sep=","),]
K27_TSSup3000.rep1.rep2.gene <- K27_TSSup3000.rep1.rep2.gene[,-((ncol(K27_TSSup3000.rep1.rep2.gene)-5):ncol(K27_TSSup3000.rep1.rep2.gene))]

K27_TSSup3000.mat <- K27_TSSup3000.rep1.rep2.gene[,4:27]
row.names(K27_TSSup3000.mat) <- K27_TSSup3000.rep1.rep2.gene[,28]
head(K27_TSSup3000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K27_TSSup3000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK27_TSSup3000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K27_TSSup3000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK27_TSSup3000_CPMover1_edgeR_result.csv", row.names=F)



### Venn diagram
dir.create("../figs/edgeR")
K4_TSSup3000.mat <- read.csv("seasonalK4_TSSup3000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27_TSSup3000.mat <- read.csv("seasonalK27_TSSup3000_CPMover1_edgeR_result.csv",header=T,sep=",")

K4_TSSup3000_season_osci <- K4_TSSup3000.mat[which(K4_TSSup3000.mat$FDR<0.05),]
K27_TSSup3000_season_osci <- K27_TSSup3000.mat[which(K27_TSSup3000.mat$FDR<0.05),]
K4TSSup3000_K27TSSup3000_season_osci <- merge(K4_TSSup3000_season_osci,K27_TSSup3000_season_osci,by=c("Ahal_ID"),all=F,sort=F)


pdf("../figs/edgeR/venn_K4TSSup3000_K27TSSup3000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4_TSSup3000_season_osci), area2=nrow(K27_TSSup3000_season_osci), cross.area=nrow(K4TSSup3000_K27TSSup3000_season_osci), 
  cex=1, category=c("K4_TSSup3000 seasonal","K27_TSSup3000 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans", 
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()





##### Extraction of genes with seasonal difference #####

# H3K4me3
K4_TSSup6000.rep1 <- read.csv("K4_TSSup6000_rep1.csv",header=T,sep=",")
K4_TSSup6000.rep2 <- read.csv("K4_TSSup6000_rep1.csv",header=T,sep=",")

colnames(K4_TSSup6000.rep1) <- c("K4_TSSup6000_chr","K4_TSSup6000_start","K4_TSSup6000_end", 
                                 "K4_TSSup6000_rep1_11.6","K4_TSSup6000_rep1_12.4","K4_TSSup6000_rep1_1.8","K4_TSSup6000_rep1_2.5",
                                 "K4_TSSup6000_rep1_3.5","K4_TSSup6000_rep1_4.2","K4_TSSup6000_rep1_4.30","K4_TSSup6000_rep1_5.28",
                                 "K4_TSSup6000_rep1_7.2","K4_TSSup6000_rep1_7.30","K4_TSSup6000_rep1_8.27","K4_TSSup6000_rep1_9.24")
colnames(K4_TSSup6000.rep2) <- c("K4_TSSup6000_chr","K4_TSSup6000_start", "K4_TSSup6000_end", 
                                 "K4_TSSup6000_rep2_11.6","K4_TSSup6000_rep2_12.4","K4_TSSup6000_rep2_1.8","K4_TSSup6000_rep2_2.5",
                                 "K4_TSSup6000_rep2_3.5","K4_TSSup6000_rep2_4.2","K4_TSSup6000_rep2_4.30","K4_TSSup6000_rep2_5.28",
                                 "K4_TSSup6000_rep2_7.2","K4_TSSup6000_rep2_7.30","K4_TSSup6000_rep2_8.27","K4_TSSup6000_rep2_9.24")
K4_TSSup6000.rep1.rep2 <- merge(K4_TSSup6000.rep1, K4_TSSup6000.rep2, 
                                by=c("K4_TSSup6000_chr","K4_TSSup6000_start","K4_TSSup6000_end"), 
                                all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup6000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K4_TSSup6000_chr","K4_TSSup6000_start","K4_TSSup6000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K4_TSSup6000.rep1.rep2.gene <- merge(K4_TSSup6000.rep1.rep2, gene, by=c("K4_TSSup6000_chr","K4_TSSup6000_start","K4_TSSup6000_end"), 
                                     all=F, sort=F)
K4_TSSup6000.rep1.rep2.gene <- K4_TSSup6000.rep1.rep2.gene[!duplicated(paste(K4_TSSup6000.rep1.rep2.gene$K4_TSSup6000_chr,K4_TSSup6000.rep1.rep2.gene$K4_TSSup6000_start,
                                                                             K4_TSSup6000.rep1.rep2.gene$K4_TSSup6000_end,K4_TSSup6000.rep1.rep2.gene$Ahal_ID),sep=","),]
K4_TSSup6000.rep1.rep2.gene <- K4_TSSup6000.rep1.rep2.gene[,-((ncol(K4_TSSup6000.rep1.rep2.gene)-5):ncol(K4_TSSup6000.rep1.rep2.gene))]

K4_TSSup6000.mat <- K4_TSSup6000.rep1.rep2.gene[,4:27]
row.names(K4_TSSup6000.mat) <- K4_TSSup6000.rep1.rep2.gene[,28]
head(K4_TSSup6000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K4_TSSup6000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK4_TSSup6000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K4_TSSup6000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK4_TSSup6000_CPMover1_edgeR_result.csv", row.names=F)



### H3K27me3
K27_TSSup6000.rep1 <- read.csv("K27_TSSup6000_rep1.csv",header=T,sep=",")
K27_TSSup6000.rep2 <- read.csv("K27_TSSup6000_rep2.csv",header=T,sep=",")

colnames(K27_TSSup6000.rep1) <- c("K27_TSSup6000_chr","K27_TSSup6000_start","K27_TSSup6000_end", 
                                  "K27_TSSup6000_rep1_11.6","K27_TSSup6000_rep1_12.4","K27_TSSup6000_rep1_1.8","K27_TSSup6000_rep1_2.5",
                                  "K27_TSSup6000_rep1_3.5","K27_TSSup6000_rep1_4.2","K27_TSSup6000_rep1_4.30","K27_TSSup6000_rep1_5.28",
                                  "K27_TSSup6000_rep1_7.2","K27_TSSup6000_rep1_7.30","K27_TSSup6000_rep1_8.27","K27_TSSup6000_rep1_9.24")
colnames(K27_TSSup6000.rep2) <- c("K27_TSSup6000_chr","K27_TSSup6000_start", "K27_TSSup6000_end", 
                                  "K27_TSSup6000_rep2_11.6","K27_TSSup6000_rep2_12.4","K27_TSSup6000_rep2_1.8","K27_TSSup6000_rep2_2.5",
                                  "K27_TSSup6000_rep2_3.5","K27_TSSup6000_rep2_4.2","K27_TSSup6000_rep2_4.30","K27_TSSup6000_rep2_5.28",
                                  "K27_TSSup6000_rep2_7.2","K27_TSSup6000_rep2_7.30","K27_TSSup6000_rep2_8.27","K27_TSSup6000_rep2_9.24")
K27_TSSup6000.rep1.rep2 <- merge(K27_TSSup6000.rep1, K27_TSSup6000.rep2, 
                                 by=c("K27_TSSup6000_chr","K27_TSSup6000_start","K27_TSSup6000_end"), 
                                 all=F, sort=F)

gene <- read.table("~/Ahal_v2_2.2/TSSup6000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K27_TSSup6000_chr","K27_TSSup6000_start","K27_TSSup6000_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K27_TSSup6000.rep1.rep2.gene <- merge(K27_TSSup6000.rep1.rep2, gene, by=c("K27_TSSup6000_chr","K27_TSSup6000_start","K27_TSSup6000_end"), 
                                      all=F, sort=F)
K27_TSSup6000.rep1.rep2.gene <- K27_TSSup6000.rep1.rep2.gene[!duplicated(paste(K27_TSSup6000.rep1.rep2.gene$K27_TSSup6000_chr,K27_TSSup6000.rep1.rep2.gene$K27_TSSup6000_start,
                                                                               K27_TSSup6000.rep1.rep2.gene$K27_TSSup6000_end,K27_TSSup6000.rep1.rep2.gene$Ahal_ID),sep=","),]
K27_TSSup6000.rep1.rep2.gene <- K27_TSSup6000.rep1.rep2.gene[,-((ncol(K27_TSSup6000.rep1.rep2.gene)-5):ncol(K27_TSSup6000.rep1.rep2.gene))]

K27_TSSup6000.mat <- K27_TSSup6000.rep1.rep2.gene[,4:27]
row.names(K27_TSSup6000.mat) <- K27_TSSup6000.rep1.rep2.gene[,28]
head(K27_TSSup6000.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K27_TSSup6000.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK27_TSSup6000_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K27_TSSup6000.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK27_TSSup6000_CPMover1_edgeR_result.csv", row.names=F)



### Venn diagram
dir.create("../figs/edgeR")
K4_TSSup6000.mat <- read.csv("seasonalK4_TSSup6000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27_TSSup6000.mat <- read.csv("seasonalK27_TSSup6000_CPMover1_edgeR_result.csv",header=T,sep=",")

K4_TSSup6000_season_osci <- K4_TSSup6000.mat[which(K4_TSSup6000.mat$FDR<0.05),]
K27_TSSup6000_season_osci <- K27_TSSup6000.mat[which(K27_TSSup6000.mat$FDR<0.05),]
K4TSSup6000_K27TSSup6000_season_osci <- merge(K4_TSSup6000_season_osci,K27_TSSup6000_season_osci,by=c("Ahal_ID"),all=F,sort=F)


pdf("../figs/edgeR/venn_K4TSSup6000_K27TSSup6000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4_TSSup6000_season_osci), area2=nrow(K27_TSSup6000_season_osci), cross.area=nrow(K4TSSup6000_K27TSSup6000_season_osci), 
  cex=1, category=c("K4_TSSup6000 seasonal","K27_TSSup6000 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans", 
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()


