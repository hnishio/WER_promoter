
library(fields)
library(LSD)
library(calibrate)
library(cocor)

source("functions/ng.Colors.R")

setwd("data")


###### Integration of H3K4me3 and H3K27me3 data #####

# rep1
K4 <- read.csv("K4_TSSup1000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup1000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup1000_rep1_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

#rep2
K4 <- read.csv("K4_TSSup1000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup1000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup1000_rep2_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup1000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup1000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22],K4.K27.rep1[,11:22],K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37],K4.K27.rep1[,26:37],K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22],K4.K27.rep2[,11:22],K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37],K4.K27.rep2[,26:37],K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27.duplicate[,1:7]



##### Calculation of max, min, mean, sd, amp #####
maxK4 <- rep(NA,length=nrow(K4))
minK4 <- rep(NA,length=nrow(K4))
meanK4 <- rep(NA,length=nrow(K4))
sdK4 <- rep(NA,length=nrow(K4))
ampK4 <- rep(NA,length=nrow(K4))
maxK27 <- rep(NA,length=nrow(K4))
minK27 <- rep(NA,length=nrow(K4))
meanK27 <- rep(NA,length=nrow(K4))
sdK27 <- rep(NA,length=nrow(K4))
ampK27 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(K4.K27.duplicate,maxK4,minK4,meanK4,sdK4,ampK4,
                maxK27,minK27,meanK27,sdK27,ampK27)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  K4.K27$maxK4[i] <- max(predK4$y)
  K4.K27$minK4[i] <- min(predK4$y)
  K4.K27$meanK4[i] <- mean(predK4$y)
  K4.K27$sdK4[i] <- sd(predK4$y)
  K4.K27$ampK4[i] <- max(predK4$y)-min(predK4$y)
  K4.K27$maxK27[i] <- max(predK27$y)
  K4.K27$minK27[i] <- min(predK27$y)
  K4.K27$meanK27[i] <- mean(predK27$y)
  K4.K27$sdK27[i] <- sd(predK27$y)
  K4.K27$ampK27[i] <- max(predK27$y)-min(predK27$y)
}
write.csv(K4.K27, file="K4K27_TSSup1000_data_stat.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
K4.K27 <- read.csv("K4K27_TSSup1000_data_stat.csv",header=T,sep=",")

K4.over2 <- subset(K4.K27,maxK4>2)
K27.over2 <- subset(K4.K27,maxK27>2)
K4.K27.over2 <- subset(K4.K27,maxK4>2&maxK27>2)

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.K27.over2, 
          file="K4K27_TSSup1000_data_stat_maxK4K27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.over2, 
          file="K4K27_TSSup1000_data_stat_maxK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup1000_data_stat_maxK27over2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup1000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup1000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:77]

colnames <- vector(length=ncol(K4.K27.over1))
list <- strsplit(colnames(K4.K27.over1), "\\.x")
for(i in 1:ncol(K4.K27.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over1) <- colnames

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)

write.csv(K4.over1, 
          file="K4K27_TSSup1000_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27.over1, 
          file="K4K27_TSSup1000_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over1, 
          file="K4K27_TSSup1000_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)


##### Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup1000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup1000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over2 <- subset(K4,ampK4>2)
K27.over2 <- subset(K27,ampK27>2)
K4.K27.over2 <- 
  merge(K4.over2, K27.over2, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over2 <- K4.K27.over2[,1:77]

colnames <- vector(length=ncol(K4.K27.over2))
list <- strsplit(colnames(K4.K27.over2), "\\.x")
for(i in 1:ncol(K4.K27.over2)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over2) <- colnames

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.over2, 
          file="K4K27_TSSup1000_data_stat_maxK4over2_ampK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup1000_data_stat_maxK27over2_ampK27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over2, 
          file="K4K27_TSSup1000_data_stat_maxK4K27over2_ampK4K27over2.csv", 
          quote=F, row.names=F)





###### Integration of H3K4me3 and H3K27me3 data #####

# rep1
K4 <- read.csv("K4_TSSup2000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup2000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup2000_rep1_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

#rep2
K4 <- read.csv("K4_TSSup2000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup2000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup2000_rep2_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup2000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup2000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22],K4.K27.rep1[,11:22],K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37],K4.K27.rep1[,26:37],K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22],K4.K27.rep2[,11:22],K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37],K4.K27.rep2[,26:37],K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27.duplicate[,1:7]



##### Calculation of max, min, mean, sd, amp #####
maxK4 <- rep(NA,length=nrow(K4))
minK4 <- rep(NA,length=nrow(K4))
meanK4 <- rep(NA,length=nrow(K4))
sdK4 <- rep(NA,length=nrow(K4))
ampK4 <- rep(NA,length=nrow(K4))
maxK27 <- rep(NA,length=nrow(K4))
minK27 <- rep(NA,length=nrow(K4))
meanK27 <- rep(NA,length=nrow(K4))
sdK27 <- rep(NA,length=nrow(K4))
ampK27 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(K4.K27.duplicate,maxK4,minK4,meanK4,sdK4,ampK4,
                maxK27,minK27,meanK27,sdK27,ampK27)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  K4.K27$maxK4[i] <- max(predK4$y)
  K4.K27$minK4[i] <- min(predK4$y)
  K4.K27$meanK4[i] <- mean(predK4$y)
  K4.K27$sdK4[i] <- sd(predK4$y)
  K4.K27$ampK4[i] <- max(predK4$y)-min(predK4$y)
  K4.K27$maxK27[i] <- max(predK27$y)
  K4.K27$minK27[i] <- min(predK27$y)
  K4.K27$meanK27[i] <- mean(predK27$y)
  K4.K27$sdK27[i] <- sd(predK27$y)
  K4.K27$ampK27[i] <- max(predK27$y)-min(predK27$y)
}
write.csv(K4.K27, file="K4K27_TSSup2000_data_stat.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
K4.K27 <- read.csv("K4K27_TSSup2000_data_stat.csv",header=T,sep=",")

K4.over2 <- subset(K4.K27,maxK4>2)
K27.over2 <- subset(K4.K27,maxK27>2)
K4.K27.over2 <- subset(K4.K27,maxK4>2&maxK27>2)

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.K27.over2, 
          file="K4K27_TSSup2000_data_stat_maxK4K27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.over2, 
          file="K4K27_TSSup2000_data_stat_maxK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup2000_data_stat_maxK27over2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup2000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup2000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:77]

colnames <- vector(length=ncol(K4.K27.over1))
list <- strsplit(colnames(K4.K27.over1), "\\.x")
for(i in 1:ncol(K4.K27.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over1) <- colnames

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)

write.csv(K4.over1, 
          file="K4K27_TSSup2000_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27.over1, 
          file="K4K27_TSSup2000_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over1, 
          file="K4K27_TSSup2000_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup2000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup2000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over2 <- subset(K4,ampK4>2)
K27.over2 <- subset(K27,ampK27>2)
K4.K27.over2 <- 
  merge(K4.over2, K27.over2, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over2 <- K4.K27.over2[,1:77]

colnames <- vector(length=ncol(K4.K27.over2))
list <- strsplit(colnames(K4.K27.over2), "\\.x")
for(i in 1:ncol(K4.K27.over2)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over2) <- colnames

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.over2, 
          file="K4K27_TSSup2000_data_stat_maxK4over2_ampK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup2000_data_stat_maxK27over2_ampK27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over2, 
          file="K4K27_TSSup2000_data_stat_maxK4K27over2_ampK4K27over2.csv", 
          quote=F, row.names=F)





###### Integration of H3K4me3 and H3K27me3 data #####

# rep1
K4 <- read.csv("K4_TSSup3000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup3000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup3000_rep1_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

#rep2
K4 <- read.csv("K4_TSSup3000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup3000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup3000_rep2_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup3000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup3000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22],K4.K27.rep1[,11:22],K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37],K4.K27.rep1[,26:37],K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22],K4.K27.rep2[,11:22],K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37],K4.K27.rep2[,26:37],K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27.duplicate[,1:7]



##### Calculation of max, min, mean, sd, amp #####
maxK4 <- rep(NA,length=nrow(K4))
minK4 <- rep(NA,length=nrow(K4))
meanK4 <- rep(NA,length=nrow(K4))
sdK4 <- rep(NA,length=nrow(K4))
ampK4 <- rep(NA,length=nrow(K4))
maxK27 <- rep(NA,length=nrow(K4))
minK27 <- rep(NA,length=nrow(K4))
meanK27 <- rep(NA,length=nrow(K4))
sdK27 <- rep(NA,length=nrow(K4))
ampK27 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(K4.K27.duplicate,maxK4,minK4,meanK4,sdK4,ampK4,
                maxK27,minK27,meanK27,sdK27,ampK27)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  K4.K27$maxK4[i] <- max(predK4$y)
  K4.K27$minK4[i] <- min(predK4$y)
  K4.K27$meanK4[i] <- mean(predK4$y)
  K4.K27$sdK4[i] <- sd(predK4$y)
  K4.K27$ampK4[i] <- max(predK4$y)-min(predK4$y)
  K4.K27$maxK27[i] <- max(predK27$y)
  K4.K27$minK27[i] <- min(predK27$y)
  K4.K27$meanK27[i] <- mean(predK27$y)
  K4.K27$sdK27[i] <- sd(predK27$y)
  K4.K27$ampK27[i] <- max(predK27$y)-min(predK27$y)
}
write.csv(K4.K27, file="K4K27_TSSup3000_data_stat.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
K4.K27 <- read.csv("K4K27_TSSup3000_data_stat.csv",header=T,sep=",")

K4.over2 <- subset(K4.K27,maxK4>2)
K27.over2 <- subset(K4.K27,maxK27>2)
K4.K27.over2 <- subset(K4.K27,maxK4>2&maxK27>2)

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.K27.over2, 
          file="K4K27_TSSup3000_data_stat_maxK4K27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.over2, 
          file="K4K27_TSSup3000_data_stat_maxK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup3000_data_stat_maxK27over2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup3000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup3000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:77]

colnames <- vector(length=ncol(K4.K27.over1))
list <- strsplit(colnames(K4.K27.over1), "\\.x")
for(i in 1:ncol(K4.K27.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over1) <- colnames

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)

write.csv(K4.over1, 
          file="K4K27_TSSup3000_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27.over1, 
          file="K4K27_TSSup3000_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over1, 
          file="K4K27_TSSup3000_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)


##### Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup3000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup3000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over2 <- subset(K4,ampK4>2)
K27.over2 <- subset(K27,ampK27>2)
K4.K27.over2 <- 
  merge(K4.over2, K27.over2, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over2 <- K4.K27.over2[,1:77]

colnames <- vector(length=ncol(K4.K27.over2))
list <- strsplit(colnames(K4.K27.over2), "\\.x")
for(i in 1:ncol(K4.K27.over2)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over2) <- colnames

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.over2, 
          file="K4K27_TSSup3000_data_stat_maxK4over2_ampK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup3000_data_stat_maxK27over2_ampK27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over2, 
          file="K4K27_TSSup3000_data_stat_maxK4K27over2_ampK4K27over2.csv", 
          quote=F, row.names=F)





###### Integration of H3K4me3 and H3K27me3 data #####

# rep1
K4 <- read.csv("K4_TSSup6000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup6000_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup6000_rep1_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

#rep2
K4 <- read.csv("K4_TSSup6000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_TSSup6000_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_TSSup6000_rep2_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_TSSup6000_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_TSSup6000_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22],K4.K27.rep1[,11:22],K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37],K4.K27.rep1[,26:37],K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22],K4.K27.rep2[,11:22],K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37],K4.K27.rep2[,26:37],K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27.duplicate[,1:7]



##### Calculation of max, min, mean, sd, amp #####
maxK4 <- rep(NA,length=nrow(K4))
minK4 <- rep(NA,length=nrow(K4))
meanK4 <- rep(NA,length=nrow(K4))
sdK4 <- rep(NA,length=nrow(K4))
ampK4 <- rep(NA,length=nrow(K4))
maxK27 <- rep(NA,length=nrow(K4))
minK27 <- rep(NA,length=nrow(K4))
meanK27 <- rep(NA,length=nrow(K4))
sdK27 <- rep(NA,length=nrow(K4))
ampK27 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(K4.K27.duplicate,maxK4,minK4,meanK4,sdK4,ampK4,
                maxK27,minK27,meanK27,sdK27,ampK27)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  K4.K27$maxK4[i] <- max(predK4$y)
  K4.K27$minK4[i] <- min(predK4$y)
  K4.K27$meanK4[i] <- mean(predK4$y)
  K4.K27$sdK4[i] <- sd(predK4$y)
  K4.K27$ampK4[i] <- max(predK4$y)-min(predK4$y)
  K4.K27$maxK27[i] <- max(predK27$y)
  K4.K27$minK27[i] <- min(predK27$y)
  K4.K27$meanK27[i] <- mean(predK27$y)
  K4.K27$sdK27[i] <- sd(predK27$y)
  K4.K27$ampK27[i] <- max(predK27$y)-min(predK27$y)
}
write.csv(K4.K27, file="K4K27_TSSup6000_data_stat.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
K4.K27 <- read.csv("K4K27_TSSup6000_data_stat.csv",header=T,sep=",")

K4.over2 <- subset(K4.K27,maxK4>2)
K27.over2 <- subset(K4.K27,maxK27>2)
K4.K27.over2 <- subset(K4.K27,maxK4>2&maxK27>2)

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.K27.over2, 
          file="K4K27_TSSup6000_data_stat_maxK4K27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.over2, 
          file="K4K27_TSSup6000_data_stat_maxK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup6000_data_stat_maxK27over2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup6000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup6000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:77]

colnames <- vector(length=ncol(K4.K27.over1))
list <- strsplit(colnames(K4.K27.over1), "\\.x")
for(i in 1:ncol(K4.K27.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over1) <- colnames

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)

write.csv(K4.over1, 
          file="K4K27_TSSup6000_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27.over1, 
          file="K4K27_TSSup6000_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over1, 
          file="K4K27_TSSup6000_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup6000_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_TSSup6000_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over2 <- subset(K4,ampK4>2)
K27.over2 <- subset(K27,ampK27>2)
K4.K27.over2 <- 
  merge(K4.over2, K27.over2, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over2 <- K4.K27.over2[,1:77]

colnames <- vector(length=ncol(K4.K27.over2))
list <- strsplit(colnames(K4.K27.over2), "\\.x")
for(i in 1:ncol(K4.K27.over2)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over2) <- colnames

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.over2, 
          file="K4K27_TSSup6000_data_stat_maxK4over2_ampK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_TSSup6000_data_stat_maxK27over2_ampK27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over2, 
          file="K4K27_TSSup6000_data_stat_maxK4K27over2_ampK4K27over2.csv", 
          quote=F, row.names=F)
