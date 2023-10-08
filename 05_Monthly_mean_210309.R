

setwd("data")

##### Ordering of genes according to amplitude of K4 #####
K4 <- read.csv("K4K27_TSSup1000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup1000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)

K4 <- read.csv("K4K27_TSSup2000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup2000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)

K4 <- read.csv("K4K27_TSSup3000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup3000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)

K4 <- read.csv("K4K27_TSSup6000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup6000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)


##### Ordering of genes according to ampK27 #####
K27 <- read.csv("K4K27_TSSup1000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup1000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)

K27 <- read.csv("K4K27_TSSup2000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup2000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)

K27 <- read.csv("K4K27_TSSup3000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup3000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)

K27 <- read.csv("K4K27_TSSup6000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup6000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)





##### Ordering of genes according to ampK4 #####
K4 <- read.csv("K4K27_TSSup1000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup1000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)

K4 <- read.csv("K4K27_TSSup2000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup2000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)

K4 <- read.csv("K4K27_TSSup3000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup3000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)

K4 <- read.csv("K4K27_TSSup6000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor.csv")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
nrow(K4.amp.order)
write.csv(K4.amp.order, 
          file="K4K27_TSSup6000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv",
          quote=F, row.names=F)


##### Ordering of genes according to ampK27 #####
K27 <- read.csv("K4K27_TSSup1000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup1000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)

K27 <- read.csv("K4K27_TSSup2000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup2000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)

K27 <- read.csv("K4K27_TSSup3000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup3000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)

K27 <- read.csv("K4K27_TSSup6000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor.csv")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
nrow(K27.amp.order)
write.csv(K27.amp.order, 
          file="K4K27_TSSup6000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv",
          quote=F, row.names=F)






##### Calculation of monthly mean #####
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


K4.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K4.month) <- 
  c("K4.Nov", "K4.Dec", "K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", 
    "K4.May", "K4.Jun", "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct")
K27.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K27.month) <- 
  c("K27.Nov", "K27.Dec", "K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", 
    "K27.May", "K27.Jun", "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct")

K4.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K4.min.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.min.date.fromJan1 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(title, K4.month, K27.month, 
                K4.max.date.fromJan1, K4.min.date.fromJan1, 
                K27.max.date.fromJan1, K27.min.date.fromJan1)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  # monthly mean
  K4.K27$K4.Nov[i] <- mean(predK4$y[365<predK4$x & predK4$x<=365+30])
  K4.K27$K4.Dec[i] <- mean(predK4$y[365+30<predK4$x & predK4$x<=365+61])
  K4.K27$K4.Jan[i] <- mean(predK4$y[365+61<predK4$x & predK4$x<=365+92])
  K4.K27$K4.Feb[i] <- mean(predK4$y[365+92<predK4$x & predK4$x<=365+120])
  K4.K27$K4.Mar[i] <- mean(predK4$y[365+120<predK4$x & predK4$x<=365+151])
  K4.K27$K4.Apr[i] <- mean(predK4$y[365+151<predK4$x & predK4$x<=365+181])
  K4.K27$K4.May[i] <- mean(predK4$y[365+181<predK4$x & predK4$x<=365+212])
  K4.K27$K4.Jun[i] <- mean(predK4$y[365+212<predK4$x & predK4$x<=365+242])
  K4.K27$K4.Jul[i] <- mean(predK4$y[365+242<predK4$x & predK4$x<=365+273])
  K4.K27$K4.Aug[i] <- mean(predK4$y[365+273<predK4$x & predK4$x<=365+304])
  K4.K27$K4.Sep[i] <- mean(predK4$y[365+304<predK4$x & predK4$x<=365+334])
  K4.K27$K4.Oct[i] <- mean(predK4$y[365+334<predK4$x & predK4$x<=365+365])
  
  K4.K27$K27.Nov[i] <- mean(predK27$y[365<predK27$x & predK27$x<=365+30])
  K4.K27$K27.Dec[i] <- mean(predK27$y[365+30<predK27$x & predK27$x<=365+61])
  K4.K27$K27.Jan[i] <- mean(predK27$y[365+61<predK27$x & predK27$x<=365+92])
  K4.K27$K27.Feb[i] <- mean(predK27$y[365+92<predK27$x & predK27$x<=365+120])
  K4.K27$K27.Mar[i] <- mean(predK27$y[365+120<predK27$x & predK27$x<=365+151])
  K4.K27$K27.Apr[i] <- mean(predK27$y[365+151<predK27$x & predK27$x<=365+181])
  K4.K27$K27.May[i] <- mean(predK27$y[365+181<predK27$x & predK27$x<=365+212])
  K4.K27$K27.Jun[i] <- mean(predK27$y[365+212<predK27$x & predK27$x<=365+242])
  K4.K27$K27.Jul[i] <- mean(predK27$y[365+242<predK27$x & predK27$x<=365+273])
  K4.K27$K27.Aug[i] <- mean(predK27$y[365+273<predK27$x & predK27$x<=365+304])
  K4.K27$K27.Sep[i] <- mean(predK27$y[365+304<predK27$x & predK27$x<=365+334])
  K4.K27$K27.Oct[i] <- mean(predK27$y[365+334<predK27$x & predK27$x<=365+365])
  
  # Days from 1 Jan.
  if(max(predK4$y)>0){
    if((365+61 < predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <=365+365)){
      K4.K27$K4.max.date.fromJan1[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 365 -61
    }else{K4.K27$K4.max.date.fromJan1[i] <- 
      predK4$x[predK4$y == max(predK4$y)] + 365 - 365 -61}
    if((365+61 < predK4$x[predK4$y == min(predK4$y)]) && 
       (predK4$x[predK4$y == min(predK4$y)] <=365+365)){
      K4.K27$K4.min.date.fromJan1[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 365 -61
    }else{K4.K27$K4.min.date.fromJan1[i] <- 
      predK4$x[predK4$y == min(predK4$y)] + 365 - 365 -61}	
  }else{}
  
  if(max(predK27$y)>0){
    if((365+61 < predK27$x[predK27$y == max(predK27$y)]) && 
       (predK27$x[predK27$y == max(predK27$y)] <=365+365)){
      K4.K27$K27.max.date.fromJan1[i] <- 
        predK27$x[predK27$y == max(predK27$y)] - 365 -61
    }else{K4.K27$K27.max.date.fromJan1[i] <- 
      predK27$x[predK27$y == max(predK27$y)] + 365 - 365 -61}
    if((365+61 < predK27$x[predK27$y == min(predK27$y)]) && 
       (predK27$x[predK27$y == min(predK27$y)] <=365+365)){
      K4.K27$K27.min.date.fromJan1[i] <- 
        predK27$x[predK27$y == min(predK27$y)] - 365 -61
    }else{K4.K27$K27.min.date.fromJan1[i] <- 
      predK27$x[predK27$y == min(predK27$y)] + 365 - 365 -61}
  }else{}
}

K4.K27.2 <- K4.K27[,c(1:7,10:19,8,9,22:31,20,21,32:35)]
write.csv(K4.K27.2, file="K4K27_TSSup1000_monthmean.csv", quote=F, row.names=F)



#####  Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup1000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup1000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup1000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup1000_monthmean_ampK4over1_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup1000_monthmean_ampK27over1_edgeR_cosinor.csv", quote=F, row.names=F)


#####  Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup1000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup1000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup1000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup1000_monthmean_ampK4over2_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup1000_monthmean_ampK27over2_edgeR_cosinor.csv", quote=F, row.names=F)





##### Calculation of monthly mean #####
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


K4.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K4.month) <- 
  c("K4.Nov", "K4.Dec", "K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", 
    "K4.May", "K4.Jun", "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct")
K27.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K27.month) <- 
  c("K27.Nov", "K27.Dec", "K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", 
    "K27.May", "K27.Jun", "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct")

K4.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K4.min.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.min.date.fromJan1 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(title, K4.month, K27.month, 
                K4.max.date.fromJan1, K4.min.date.fromJan1, 
                K27.max.date.fromJan1, K27.min.date.fromJan1)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  # monthly mean
  K4.K27$K4.Nov[i] <- mean(predK4$y[365<predK4$x & predK4$x<=365+30])
  K4.K27$K4.Dec[i] <- mean(predK4$y[365+30<predK4$x & predK4$x<=365+61])
  K4.K27$K4.Jan[i] <- mean(predK4$y[365+61<predK4$x & predK4$x<=365+92])
  K4.K27$K4.Feb[i] <- mean(predK4$y[365+92<predK4$x & predK4$x<=365+120])
  K4.K27$K4.Mar[i] <- mean(predK4$y[365+120<predK4$x & predK4$x<=365+151])
  K4.K27$K4.Apr[i] <- mean(predK4$y[365+151<predK4$x & predK4$x<=365+181])
  K4.K27$K4.May[i] <- mean(predK4$y[365+181<predK4$x & predK4$x<=365+212])
  K4.K27$K4.Jun[i] <- mean(predK4$y[365+212<predK4$x & predK4$x<=365+242])
  K4.K27$K4.Jul[i] <- mean(predK4$y[365+242<predK4$x & predK4$x<=365+273])
  K4.K27$K4.Aug[i] <- mean(predK4$y[365+273<predK4$x & predK4$x<=365+304])
  K4.K27$K4.Sep[i] <- mean(predK4$y[365+304<predK4$x & predK4$x<=365+334])
  K4.K27$K4.Oct[i] <- mean(predK4$y[365+334<predK4$x & predK4$x<=365+365])
  
  K4.K27$K27.Nov[i] <- mean(predK27$y[365<predK27$x & predK27$x<=365+30])
  K4.K27$K27.Dec[i] <- mean(predK27$y[365+30<predK27$x & predK27$x<=365+61])
  K4.K27$K27.Jan[i] <- mean(predK27$y[365+61<predK27$x & predK27$x<=365+92])
  K4.K27$K27.Feb[i] <- mean(predK27$y[365+92<predK27$x & predK27$x<=365+120])
  K4.K27$K27.Mar[i] <- mean(predK27$y[365+120<predK27$x & predK27$x<=365+151])
  K4.K27$K27.Apr[i] <- mean(predK27$y[365+151<predK27$x & predK27$x<=365+181])
  K4.K27$K27.May[i] <- mean(predK27$y[365+181<predK27$x & predK27$x<=365+212])
  K4.K27$K27.Jun[i] <- mean(predK27$y[365+212<predK27$x & predK27$x<=365+242])
  K4.K27$K27.Jul[i] <- mean(predK27$y[365+242<predK27$x & predK27$x<=365+273])
  K4.K27$K27.Aug[i] <- mean(predK27$y[365+273<predK27$x & predK27$x<=365+304])
  K4.K27$K27.Sep[i] <- mean(predK27$y[365+304<predK27$x & predK27$x<=365+334])
  K4.K27$K27.Oct[i] <- mean(predK27$y[365+334<predK27$x & predK27$x<=365+365])
  
  # Days from 1 Jan.
  if(max(predK4$y)>0){
    if((365+61 < predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <=365+365)){
      K4.K27$K4.max.date.fromJan1[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 365 -61
    }else{K4.K27$K4.max.date.fromJan1[i] <- 
      predK4$x[predK4$y == max(predK4$y)] + 365 - 365 -61}
    if((365+61 < predK4$x[predK4$y == min(predK4$y)]) && 
       (predK4$x[predK4$y == min(predK4$y)] <=365+365)){
      K4.K27$K4.min.date.fromJan1[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 365 -61
    }else{K4.K27$K4.min.date.fromJan1[i] <- 
      predK4$x[predK4$y == min(predK4$y)] + 365 - 365 -61}	
  }else{}
  
  if(max(predK27$y)>0){
    if((365+61 < predK27$x[predK27$y == max(predK27$y)]) && 
       (predK27$x[predK27$y == max(predK27$y)] <=365+365)){
      K4.K27$K27.max.date.fromJan1[i] <- 
        predK27$x[predK27$y == max(predK27$y)] - 365 -61
    }else{K4.K27$K27.max.date.fromJan1[i] <- 
      predK27$x[predK27$y == max(predK27$y)] + 365 - 365 -61}
    if((365+61 < predK27$x[predK27$y == min(predK27$y)]) && 
       (predK27$x[predK27$y == min(predK27$y)] <=365+365)){
      K4.K27$K27.min.date.fromJan1[i] <- 
        predK27$x[predK27$y == min(predK27$y)] - 365 -61
    }else{K4.K27$K27.min.date.fromJan1[i] <- 
      predK27$x[predK27$y == min(predK27$y)] + 365 - 365 -61}
  }else{}
}

K4.K27.2 <- K4.K27[,c(1:7,10:19,8,9,22:31,20,21,32:35)]
write.csv(K4.K27.2, file="K4K27_TSSup2000_monthmean.csv", quote=F, row.names=F)



#####  Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup2000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup2000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup2000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup2000_monthmean_ampK4over1_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup2000_monthmean_ampK27over1_edgeR_cosinor.csv", quote=F, row.names=F)


#####  Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup2000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup2000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup2000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup2000_monthmean_ampK4over2_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup2000_monthmean_ampK27over2_edgeR_cosinor.csv", quote=F, row.names=F)





##### Calculation of monthly mean #####
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


K4.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K4.month) <- 
  c("K4.Nov", "K4.Dec", "K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", 
    "K4.May", "K4.Jun", "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct")
K27.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K27.month) <- 
  c("K27.Nov", "K27.Dec", "K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", 
    "K27.May", "K27.Jun", "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct")

K4.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K4.min.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.min.date.fromJan1 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(title, K4.month, K27.month, 
                K4.max.date.fromJan1, K4.min.date.fromJan1, 
                K27.max.date.fromJan1, K27.min.date.fromJan1)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  # monthly mean
  K4.K27$K4.Nov[i] <- mean(predK4$y[365<predK4$x & predK4$x<=365+30])
  K4.K27$K4.Dec[i] <- mean(predK4$y[365+30<predK4$x & predK4$x<=365+61])
  K4.K27$K4.Jan[i] <- mean(predK4$y[365+61<predK4$x & predK4$x<=365+92])
  K4.K27$K4.Feb[i] <- mean(predK4$y[365+92<predK4$x & predK4$x<=365+120])
  K4.K27$K4.Mar[i] <- mean(predK4$y[365+120<predK4$x & predK4$x<=365+151])
  K4.K27$K4.Apr[i] <- mean(predK4$y[365+151<predK4$x & predK4$x<=365+181])
  K4.K27$K4.May[i] <- mean(predK4$y[365+181<predK4$x & predK4$x<=365+212])
  K4.K27$K4.Jun[i] <- mean(predK4$y[365+212<predK4$x & predK4$x<=365+242])
  K4.K27$K4.Jul[i] <- mean(predK4$y[365+242<predK4$x & predK4$x<=365+273])
  K4.K27$K4.Aug[i] <- mean(predK4$y[365+273<predK4$x & predK4$x<=365+304])
  K4.K27$K4.Sep[i] <- mean(predK4$y[365+304<predK4$x & predK4$x<=365+334])
  K4.K27$K4.Oct[i] <- mean(predK4$y[365+334<predK4$x & predK4$x<=365+365])
  
  K4.K27$K27.Nov[i] <- mean(predK27$y[365<predK27$x & predK27$x<=365+30])
  K4.K27$K27.Dec[i] <- mean(predK27$y[365+30<predK27$x & predK27$x<=365+61])
  K4.K27$K27.Jan[i] <- mean(predK27$y[365+61<predK27$x & predK27$x<=365+92])
  K4.K27$K27.Feb[i] <- mean(predK27$y[365+92<predK27$x & predK27$x<=365+120])
  K4.K27$K27.Mar[i] <- mean(predK27$y[365+120<predK27$x & predK27$x<=365+151])
  K4.K27$K27.Apr[i] <- mean(predK27$y[365+151<predK27$x & predK27$x<=365+181])
  K4.K27$K27.May[i] <- mean(predK27$y[365+181<predK27$x & predK27$x<=365+212])
  K4.K27$K27.Jun[i] <- mean(predK27$y[365+212<predK27$x & predK27$x<=365+242])
  K4.K27$K27.Jul[i] <- mean(predK27$y[365+242<predK27$x & predK27$x<=365+273])
  K4.K27$K27.Aug[i] <- mean(predK27$y[365+273<predK27$x & predK27$x<=365+304])
  K4.K27$K27.Sep[i] <- mean(predK27$y[365+304<predK27$x & predK27$x<=365+334])
  K4.K27$K27.Oct[i] <- mean(predK27$y[365+334<predK27$x & predK27$x<=365+365])
  
  # Days from 1 Jan.
  if(max(predK4$y)>0){
    if((365+61 < predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <=365+365)){
      K4.K27$K4.max.date.fromJan1[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 365 -61
    }else{K4.K27$K4.max.date.fromJan1[i] <- 
      predK4$x[predK4$y == max(predK4$y)] + 365 - 365 -61}
    if((365+61 < predK4$x[predK4$y == min(predK4$y)]) && 
       (predK4$x[predK4$y == min(predK4$y)] <=365+365)){
      K4.K27$K4.min.date.fromJan1[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 365 -61
    }else{K4.K27$K4.min.date.fromJan1[i] <- 
      predK4$x[predK4$y == min(predK4$y)] + 365 - 365 -61}	
  }else{}
  
  if(max(predK27$y)>0){
    if((365+61 < predK27$x[predK27$y == max(predK27$y)]) && 
       (predK27$x[predK27$y == max(predK27$y)] <=365+365)){
      K4.K27$K27.max.date.fromJan1[i] <- 
        predK27$x[predK27$y == max(predK27$y)] - 365 -61
    }else{K4.K27$K27.max.date.fromJan1[i] <- 
      predK27$x[predK27$y == max(predK27$y)] + 365 - 365 -61}
    if((365+61 < predK27$x[predK27$y == min(predK27$y)]) && 
       (predK27$x[predK27$y == min(predK27$y)] <=365+365)){
      K4.K27$K27.min.date.fromJan1[i] <- 
        predK27$x[predK27$y == min(predK27$y)] - 365 -61
    }else{K4.K27$K27.min.date.fromJan1[i] <- 
      predK27$x[predK27$y == min(predK27$y)] + 365 - 365 -61}
  }else{}
}

K4.K27.2 <- K4.K27[,c(1:7,10:19,8,9,22:31,20,21,32:35)]
write.csv(K4.K27.2, file="K4K27_TSSup3000_monthmean.csv", quote=F, row.names=F)



#####  Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup3000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup3000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup3000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup3000_monthmean_ampK4over1_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup3000_monthmean_ampK27over1_edgeR_cosinor.csv", quote=F, row.names=F)


#####  Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup3000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup3000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup3000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup3000_monthmean_ampK4over2_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup3000_monthmean_ampK27over2_edgeR_cosinor.csv", quote=F, row.names=F)





##### Calculation of monthly mean #####
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


K4.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K4.month) <- 
  c("K4.Nov", "K4.Dec", "K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", 
    "K4.May", "K4.Jun", "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct")
K27.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K27.month) <- 
  c("K27.Nov", "K27.Dec", "K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", 
    "K27.May", "K27.Jun", "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct")

K4.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K4.min.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.min.date.fromJan1 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(title, K4.month, K27.month, 
                K4.max.date.fromJan1, K4.min.date.fromJan1, 
                K27.max.date.fromJan1, K27.min.date.fromJan1)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  # monthly mean
  K4.K27$K4.Nov[i] <- mean(predK4$y[365<predK4$x & predK4$x<=365+30])
  K4.K27$K4.Dec[i] <- mean(predK4$y[365+30<predK4$x & predK4$x<=365+61])
  K4.K27$K4.Jan[i] <- mean(predK4$y[365+61<predK4$x & predK4$x<=365+92])
  K4.K27$K4.Feb[i] <- mean(predK4$y[365+92<predK4$x & predK4$x<=365+120])
  K4.K27$K4.Mar[i] <- mean(predK4$y[365+120<predK4$x & predK4$x<=365+151])
  K4.K27$K4.Apr[i] <- mean(predK4$y[365+151<predK4$x & predK4$x<=365+181])
  K4.K27$K4.May[i] <- mean(predK4$y[365+181<predK4$x & predK4$x<=365+212])
  K4.K27$K4.Jun[i] <- mean(predK4$y[365+212<predK4$x & predK4$x<=365+242])
  K4.K27$K4.Jul[i] <- mean(predK4$y[365+242<predK4$x & predK4$x<=365+273])
  K4.K27$K4.Aug[i] <- mean(predK4$y[365+273<predK4$x & predK4$x<=365+304])
  K4.K27$K4.Sep[i] <- mean(predK4$y[365+304<predK4$x & predK4$x<=365+334])
  K4.K27$K4.Oct[i] <- mean(predK4$y[365+334<predK4$x & predK4$x<=365+365])
  
  K4.K27$K27.Nov[i] <- mean(predK27$y[365<predK27$x & predK27$x<=365+30])
  K4.K27$K27.Dec[i] <- mean(predK27$y[365+30<predK27$x & predK27$x<=365+61])
  K4.K27$K27.Jan[i] <- mean(predK27$y[365+61<predK27$x & predK27$x<=365+92])
  K4.K27$K27.Feb[i] <- mean(predK27$y[365+92<predK27$x & predK27$x<=365+120])
  K4.K27$K27.Mar[i] <- mean(predK27$y[365+120<predK27$x & predK27$x<=365+151])
  K4.K27$K27.Apr[i] <- mean(predK27$y[365+151<predK27$x & predK27$x<=365+181])
  K4.K27$K27.May[i] <- mean(predK27$y[365+181<predK27$x & predK27$x<=365+212])
  K4.K27$K27.Jun[i] <- mean(predK27$y[365+212<predK27$x & predK27$x<=365+242])
  K4.K27$K27.Jul[i] <- mean(predK27$y[365+242<predK27$x & predK27$x<=365+273])
  K4.K27$K27.Aug[i] <- mean(predK27$y[365+273<predK27$x & predK27$x<=365+304])
  K4.K27$K27.Sep[i] <- mean(predK27$y[365+304<predK27$x & predK27$x<=365+334])
  K4.K27$K27.Oct[i] <- mean(predK27$y[365+334<predK27$x & predK27$x<=365+365])
  
  # Days from 1 Jan.
  if(max(predK4$y)>0){
    if((365+61 < predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <=365+365)){
      K4.K27$K4.max.date.fromJan1[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 365 -61
    }else{K4.K27$K4.max.date.fromJan1[i] <- 
      predK4$x[predK4$y == max(predK4$y)] + 365 - 365 -61}
    if((365+61 < predK4$x[predK4$y == min(predK4$y)]) && 
       (predK4$x[predK4$y == min(predK4$y)] <=365+365)){
      K4.K27$K4.min.date.fromJan1[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 365 -61
    }else{K4.K27$K4.min.date.fromJan1[i] <- 
      predK4$x[predK4$y == min(predK4$y)] + 365 - 365 -61}	
  }else{}
  
  if(max(predK27$y)>0){
    if((365+61 < predK27$x[predK27$y == max(predK27$y)]) && 
       (predK27$x[predK27$y == max(predK27$y)] <=365+365)){
      K4.K27$K27.max.date.fromJan1[i] <- 
        predK27$x[predK27$y == max(predK27$y)] - 365 -61
    }else{K4.K27$K27.max.date.fromJan1[i] <- 
      predK27$x[predK27$y == max(predK27$y)] + 365 - 365 -61}
    if((365+61 < predK27$x[predK27$y == min(predK27$y)]) && 
       (predK27$x[predK27$y == min(predK27$y)] <=365+365)){
      K4.K27$K27.min.date.fromJan1[i] <- 
        predK27$x[predK27$y == min(predK27$y)] - 365 -61
    }else{K4.K27$K27.min.date.fromJan1[i] <- 
      predK27$x[predK27$y == min(predK27$y)] + 365 - 365 -61}
  }else{}
}

K4.K27.2 <- K4.K27[,c(1:7,10:19,8,9,22:31,20,21,32:35)]
write.csv(K4.K27.2, file="K4K27_TSSup6000_monthmean.csv", quote=F, row.names=F)



#####  Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_TSSup6000_data_stat_maxK4over2_ampK4over1_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup6000_data_stat_maxK27over2_ampK27over1_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup6000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup6000_monthmean_ampK4over1_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup6000_monthmean_ampK27over1_edgeR_cosinor.csv", quote=F, row.names=F)


#####  Extraction of genes with amp > 2 #####
K4 <- read.csv("K4K27_TSSup6000_data_stat_maxK4over2_ampK4over2_edgeR_cosinor_ampK4order.csv")
K27 <- read.csv("K4K27_TSSup6000_data_stat_maxK27over2_ampK27over2_edgeR_cosinor_ampK27order.csv")
nrow(K4)
nrow(K27)

K4K27mean <- read.csv("K4K27_TSSup6000_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_TSSup6000_monthmean_ampK4over2_edgeR_cosinor.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_TSSup6000_monthmean_ampK27over2_edgeR_cosinor.csv", quote=F, row.names=F)


