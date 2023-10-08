
library(season)
library(VennDiagram)
library(fields)

source("functions/GOanalysis_functions.R")

setwd("data")


###### Integration of H3K4me3 and H3K27me3 data #####

# Preparation for cosinor
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup1000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup1000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c("2012-11-06", "2012-12-04", "2013-01-08", "2013-02-05", "2013-03-05", "2013-04-02",
          "2013-04-30", "2013-05-28", "2013-07-02", "2013-07-30", "2013-08-27", "2013-09-24")
date <- as.Date(c(date,date))
title <- K4.K27.duplicate[,1:7]


K4.mat <- matrix(nrow=nrow(K4),ncol=9)
K27.mat <- matrix(nrow=nrow(K4),ncol=9)


for(i in 1:nrow(K4)){
  his.data <- data.frame(date=date,K4=K4[i,],K27=K27[i,])
  
  res.K4 <- cosinor(K4~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K4)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
 
  lphase <- gsub("Month = ", "2012-", summary(res.K4)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K4.mat[i,] <- c(summary(res.K4)$amp*2,peak,trough,
                  summary(res.K4$glm)$coefficients[,"Estimate"],
                  summary(res.K4$glm)$coefficients[,"Pr(>|t|)"])
  
  res.K27 <- cosinor(K27~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K27)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K27)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K27.mat[i,] <- c(summary(res.K27)$amp*2,peak,trough,
                  summary(res.K27$glm)$coefficients[,"Estimate"],
                  summary(res.K27$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(K4.mat) <- c("amp.K4","peak.K4","trough.K4",
                      "coef.inter.K4","coef.cosw.K4","coef.sinw.K4",
                      "pval.inter.K4","pval.cosw.K4","pval.sinw.K4")
colnames(K27.mat) <- c("amp.K27","peak.K27","trough.K27",
                      "coef.inter.K27","coef.cosw.K27","coef.sinw.K27",
                      "pval.inter.K27","pval.cosw.K27","pval.sinw.K27")

K4.mat <- as.data.frame(K4.mat)
K27.mat <- as.data.frame(K27.mat)
K4.mat$pval.inter.K4 <- ng.BHFDR(K4.mat$pval.inter.K4)
K4.mat$pval.cosw.K4 <- ng.BHFDR(K4.mat$pval.cosw.K4)
K4.mat$pval.sinw.K4 <- ng.BHFDR(K4.mat$pval.sinw.K4)
K27.mat$pval.inter.K27 <- ng.BHFDR(K27.mat$pval.inter.K27)
K27.mat$pval.cosw.K27 <- ng.BHFDR(K27.mat$pval.cosw.K27)
K27.mat$pval.sinw.K27 <- ng.BHFDR(K27.mat$pval.sinw.K27)

K4.K27.mat <- cbind(title,K4.mat,K27.mat)
write.csv(K4.K27.mat, file="K4K27_TSSup1000_cosinor.csv", quote=F, row.names=F)


# Venn diagram
dir.create("../figs/cosinor_venn")

K4.K27.mat <- read.csv("K4K27_TSSup1000_cosinor.csv",header=T,sep=",")
K4.CPMover1 <- read.csv("seasonalK4_TSSup1000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.CPMover1 <- read.csv("seasonalK27_TSSup1000_CPMover1_edgeR_result.csv",header=T,sep=",")
K4 <- merge(K4.K27.mat,K4.CPMover1,by="Ahal_ID",all=F,sort=F)
K27 <- merge(K4.K27.mat,K27.CPMover1,by="Ahal_ID",all=F,sort=F)
K4.K27 <- merge(K4,K27,by="Ahal_ID",all=F,sort=F)
write.csv(K4[,1:ncol(K4.K27.mat)], file="K4K27_TSSup1000_seasonal_cosinor_K4_CPMover1.csv", quote=F, row.names=F)
write.csv(K27[,1:ncol(K4.K27.mat)], file="K4K27_TSSup1000_seasonal_cosinor_K27_CPMover1.csv", quote=F, row.names=F)

K4_season_osci <- length(which(K4$pval.cosw.K4<0.025 | K4$pval.sinw.K4<0.025))
K27_season_osci <- length(which(K27$pval.cosw.K27<0.025 | K27$pval.sinw.K27<0.025))
K4K27_season_osci <- length(which((K4.K27$pval.cosw.K4.x<0.025 | K4.K27$pval.sinw.K4.x<0.025) &
                                    (K4.K27$pval.cosw.K27.x<0.025 | K4.K27$pval.sinw.K27.x<0.025)))

pdf("../figs/cosinor_venn/venn_K4K27_TSSup1000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=K4_season_osci, area2=K27_season_osci, cross.area=K4K27_season_osci, 
  cex=1, category=c("K4 seasonal","K27 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()





###### Integration of H3K4me3 and H3K27me3 data #####

# Preparation for cosinor
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup2000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup2000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c("2012-11-06", "2012-12-04", "2013-01-08", "2013-02-05", "2013-03-05", "2013-04-02",
          "2013-04-30", "2013-05-28", "2013-07-02", "2013-07-30", "2013-08-27", "2013-09-24")
date <- as.Date(c(date,date))
title <- K4.K27.duplicate[,1:7]


K4.mat <- matrix(nrow=nrow(K4),ncol=9)
K27.mat <- matrix(nrow=nrow(K4),ncol=9)


for(i in 1:nrow(K4)){
  his.data <- data.frame(date=date,K4=K4[i,],K27=K27[i,])
  
  res.K4 <- cosinor(K4~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K4)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K4)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K4.mat[i,] <- c(summary(res.K4)$amp*2,peak,trough,
                  summary(res.K4$glm)$coefficients[,"Estimate"],
                  summary(res.K4$glm)$coefficients[,"Pr(>|t|)"])
  
  res.K27 <- cosinor(K27~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K27)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K27)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K27.mat[i,] <- c(summary(res.K27)$amp*2,peak,trough,
                   summary(res.K27$glm)$coefficients[,"Estimate"],
                   summary(res.K27$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(K4.mat) <- c("amp.K4","peak.K4","trough.K4",
                      "coef.inter.K4","coef.cosw.K4","coef.sinw.K4",
                      "pval.inter.K4","pval.cosw.K4","pval.sinw.K4")
colnames(K27.mat) <- c("amp.K27","peak.K27","trough.K27",
                       "coef.inter.K27","coef.cosw.K27","coef.sinw.K27",
                       "pval.inter.K27","pval.cosw.K27","pval.sinw.K27")

K4.mat <- as.data.frame(K4.mat)
K27.mat <- as.data.frame(K27.mat)
K4.mat$pval.inter.K4 <- ng.BHFDR(K4.mat$pval.inter.K4)
K4.mat$pval.cosw.K4 <- ng.BHFDR(K4.mat$pval.cosw.K4)
K4.mat$pval.sinw.K4 <- ng.BHFDR(K4.mat$pval.sinw.K4)
K27.mat$pval.inter.K27 <- ng.BHFDR(K27.mat$pval.inter.K27)
K27.mat$pval.cosw.K27 <- ng.BHFDR(K27.mat$pval.cosw.K27)
K27.mat$pval.sinw.K27 <- ng.BHFDR(K27.mat$pval.sinw.K27)

K4.K27.mat <- cbind(title,K4.mat,K27.mat)
write.csv(K4.K27.mat, file="K4K27_TSSup2000_cosinor.csv", quote=F, row.names=F)


# Venn diagram
dir.create("../figs/cosinor_venn")

K4.K27.mat <- read.csv("K4K27_TSSup2000_cosinor.csv",header=T,sep=",")
K4.CPMover1 <- read.csv("seasonalK4_TSSup2000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.CPMover1 <- read.csv("seasonalK27_TSSup2000_CPMover1_edgeR_result.csv",header=T,sep=",")
K4 <- merge(K4.K27.mat,K4.CPMover1,by="Ahal_ID",all=F,sort=F)
K27 <- merge(K4.K27.mat,K27.CPMover1,by="Ahal_ID",all=F,sort=F)
K4.K27 <- merge(K4,K27,by="Ahal_ID",all=F,sort=F)
write.csv(K4[,1:ncol(K4.K27.mat)], file="K4K27_TSSup2000_seasonal_cosinor_K4_CPMover1.csv", quote=F, row.names=F)
write.csv(K27[,1:ncol(K4.K27.mat)], file="K4K27_TSSup2000_seasonal_cosinor_K27_CPMover1.csv", quote=F, row.names=F)

K4_season_osci <- length(which(K4$pval.cosw.K4<0.025 | K4$pval.sinw.K4<0.025))
K27_season_osci <- length(which(K27$pval.cosw.K27<0.025 | K27$pval.sinw.K27<0.025))
K4K27_season_osci <- length(which((K4.K27$pval.cosw.K4.x<0.025 | K4.K27$pval.sinw.K4.x<0.025) &
                                    (K4.K27$pval.cosw.K27.x<0.025 | K4.K27$pval.sinw.K27.x<0.025)))

pdf("../figs/cosinor_venn/venn_K4K27_TSSup2000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=K4_season_osci, area2=K27_season_osci, cross.area=K4K27_season_osci, 
  cex=1, category=c("K4 seasonal","K27 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()





###### Integration of H3K4me3 and H3K27me3 data #####

# Preparation for cosinor
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup3000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup3000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c("2012-11-06", "2012-12-04", "2013-01-08", "2013-02-05", "2013-03-05", "2013-04-02",
          "2013-04-30", "2013-05-28", "2013-07-02", "2013-07-30", "2013-08-27", "2013-09-24")
date <- as.Date(c(date,date))
title <- K4.K27.duplicate[,1:7]


K4.mat <- matrix(nrow=nrow(K4),ncol=9)
K27.mat <- matrix(nrow=nrow(K4),ncol=9)


for(i in 1:nrow(K4)){
  his.data <- data.frame(date=date,K4=K4[i,],K27=K27[i,])
  
  res.K4 <- cosinor(K4~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K4)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K4)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K4.mat[i,] <- c(summary(res.K4)$amp*2,peak,trough,
                  summary(res.K4$glm)$coefficients[,"Estimate"],
                  summary(res.K4$glm)$coefficients[,"Pr(>|t|)"])
  
  res.K27 <- cosinor(K27~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K27)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K27)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K27.mat[i,] <- c(summary(res.K27)$amp*2,peak,trough,
                   summary(res.K27$glm)$coefficients[,"Estimate"],
                   summary(res.K27$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(K4.mat) <- c("amp.K4","peak.K4","trough.K4",
                      "coef.inter.K4","coef.cosw.K4","coef.sinw.K4",
                      "pval.inter.K4","pval.cosw.K4","pval.sinw.K4")
colnames(K27.mat) <- c("amp.K27","peak.K27","trough.K27",
                       "coef.inter.K27","coef.cosw.K27","coef.sinw.K27",
                       "pval.inter.K27","pval.cosw.K27","pval.sinw.K27")

K4.mat <- as.data.frame(K4.mat)
K27.mat <- as.data.frame(K27.mat)
K4.mat$pval.inter.K4 <- ng.BHFDR(K4.mat$pval.inter.K4)
K4.mat$pval.cosw.K4 <- ng.BHFDR(K4.mat$pval.cosw.K4)
K4.mat$pval.sinw.K4 <- ng.BHFDR(K4.mat$pval.sinw.K4)
K27.mat$pval.inter.K27 <- ng.BHFDR(K27.mat$pval.inter.K27)
K27.mat$pval.cosw.K27 <- ng.BHFDR(K27.mat$pval.cosw.K27)
K27.mat$pval.sinw.K27 <- ng.BHFDR(K27.mat$pval.sinw.K27)

K4.K27.mat <- cbind(title,K4.mat,K27.mat)
write.csv(K4.K27.mat, file="K4K27_TSSup3000_cosinor.csv", quote=F, row.names=F)


# Venn diagram
dir.create("../figs/cosinor_venn")

K4.K27.mat <- read.csv("K4K27_TSSup3000_cosinor.csv",header=T,sep=",")
K4.CPMover1 <- read.csv("seasonalK4_TSSup3000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.CPMover1 <- read.csv("seasonalK27_TSSup3000_CPMover1_edgeR_result.csv",header=T,sep=",")
K4 <- merge(K4.K27.mat,K4.CPMover1,by="Ahal_ID",all=F,sort=F)
K27 <- merge(K4.K27.mat,K27.CPMover1,by="Ahal_ID",all=F,sort=F)
K4.K27 <- merge(K4,K27,by="Ahal_ID",all=F,sort=F)
write.csv(K4[,1:ncol(K4.K27.mat)], file="K4K27_TSSup3000_seasonal_cosinor_K4_CPMover1.csv", quote=F, row.names=F)
write.csv(K27[,1:ncol(K4.K27.mat)], file="K4K27_TSSup3000_seasonal_cosinor_K27_CPMover1.csv", quote=F, row.names=F)

K4_season_osci <- length(which(K4$pval.cosw.K4<0.025 | K4$pval.sinw.K4<0.025))
K27_season_osci <- length(which(K27$pval.cosw.K27<0.025 | K27$pval.sinw.K27<0.025))
K4K27_season_osci <- length(which((K4.K27$pval.cosw.K4.x<0.025 | K4.K27$pval.sinw.K4.x<0.025) &
                                    (K4.K27$pval.cosw.K27.x<0.025 | K4.K27$pval.sinw.K27.x<0.025)))

pdf("../figs/cosinor_venn/venn_K4K27_TSSup3000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=K4_season_osci, area2=K27_season_osci, cross.area=K4K27_season_osci, 
  cex=1, category=c("K4 seasonal","K27 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()





###### Integration of H3K4me3 and H3K27me3 data #####

# Preparation for cosinor
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup6000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup6000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c("2012-11-06", "2012-12-04", "2013-01-08", "2013-02-05", "2013-03-05", "2013-04-02",
          "2013-04-30", "2013-05-28", "2013-07-02", "2013-07-30", "2013-08-27", "2013-09-24")
date <- as.Date(c(date,date))
title <- K4.K27.duplicate[,1:7]


K4.mat <- matrix(nrow=nrow(K4),ncol=9)
K27.mat <- matrix(nrow=nrow(K4),ncol=9)


for(i in 1:nrow(K4)){
  his.data <- data.frame(date=date,K4=K4[i,],K27=K27[i,])
  
  res.K4 <- cosinor(K4~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K4)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K4)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K4.mat[i,] <- c(summary(res.K4)$amp*2,peak,trough,
                  summary(res.K4$glm)$coefficients[,"Estimate"],
                  summary(res.K4$glm)$coefficients[,"Pr(>|t|)"])
  
  res.K27 <- cosinor(K27~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K27)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K27)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K27.mat[i,] <- c(summary(res.K27)$amp*2,peak,trough,
                   summary(res.K27$glm)$coefficients[,"Estimate"],
                   summary(res.K27$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(K4.mat) <- c("amp.K4","peak.K4","trough.K4",
                      "coef.inter.K4","coef.cosw.K4","coef.sinw.K4",
                      "pval.inter.K4","pval.cosw.K4","pval.sinw.K4")
colnames(K27.mat) <- c("amp.K27","peak.K27","trough.K27",
                       "coef.inter.K27","coef.cosw.K27","coef.sinw.K27",
                       "pval.inter.K27","pval.cosw.K27","pval.sinw.K27")

K4.mat <- as.data.frame(K4.mat)
K27.mat <- as.data.frame(K27.mat)
K4.mat$pval.inter.K4 <- ng.BHFDR(K4.mat$pval.inter.K4)
K4.mat$pval.cosw.K4 <- ng.BHFDR(K4.mat$pval.cosw.K4)
K4.mat$pval.sinw.K4 <- ng.BHFDR(K4.mat$pval.sinw.K4)
K27.mat$pval.inter.K27 <- ng.BHFDR(K27.mat$pval.inter.K27)
K27.mat$pval.cosw.K27 <- ng.BHFDR(K27.mat$pval.cosw.K27)
K27.mat$pval.sinw.K27 <- ng.BHFDR(K27.mat$pval.sinw.K27)

K4.K27.mat <- cbind(title,K4.mat,K27.mat)
write.csv(K4.K27.mat, file="K4K27_TSSup6000_cosinor.csv", quote=F, row.names=F)


# Venn diagram
dir.create("../figs/cosinor_venn")

K4.K27.mat <- read.csv("K4K27_TSSup6000_cosinor.csv",header=T,sep=",")
K4.CPMover1 <- read.csv("seasonalK4_TSSup6000_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.CPMover1 <- read.csv("seasonalK27_TSSup6000_CPMover1_edgeR_result.csv",header=T,sep=",")
K4 <- merge(K4.K27.mat,K4.CPMover1,by="Ahal_ID",all=F,sort=F)
K27 <- merge(K4.K27.mat,K27.CPMover1,by="Ahal_ID",all=F,sort=F)
K4.K27 <- merge(K4,K27,by="Ahal_ID",all=F,sort=F)
write.csv(K4[,1:ncol(K4.K27.mat)], file="K4K27_TSSup6000_seasonal_cosinor_K4_CPMover1.csv", quote=F, row.names=F)
write.csv(K27[,1:ncol(K4.K27.mat)], file="K4K27_TSSup6000_seasonal_cosinor_K27_CPMover1.csv", quote=F, row.names=F)

K4_season_osci <- length(which(K4$pval.cosw.K4<0.025 | K4$pval.sinw.K4<0.025))
K27_season_osci <- length(which(K27$pval.cosw.K27<0.025 | K27$pval.sinw.K27<0.025))
K4K27_season_osci <- length(which((K4.K27$pval.cosw.K4.x<0.025 | K4.K27$pval.sinw.K4.x<0.025) &
                                    (K4.K27$pval.cosw.K27.x<0.025 | K4.K27$pval.sinw.K27.x<0.025)))

pdf("../figs/cosinor_venn/venn_K4K27_TSSup6000_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=K4_season_osci, area2=K27_season_osci, cross.area=K4K27_season_osci, 
  cex=1, category=c("K4 seasonal","K27 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()


