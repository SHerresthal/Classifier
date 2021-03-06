# loading libraries
library("data.table")
library("class")
library("e1071")
library("ROCR")
library("MASS")
library("randomForest")
library("glmnet")
library("foreach")
library("doParallel")
library("multtest")
library("affy")
library("pamr")
library("dplyr")

# define directory with files 
dir <- "/data/Results/Test/"

# load data
load(paste(dir, "Datasets.RData", sep = ""))


# source functions 
source(paste(dir, "Functions_V2.R", sep = ""))

# define number of permutations
nperm_random <- 2 # for random sampling
nperm_crossstudy <- 2 # for cross study sampling
nperm_crossplatform <- 2 # for cross platform sampling
linux_run <- T # for linux-machine = T, for windows = F
ncores <- 2 # number of cores for parallel processing


# create output directories
dir.create(paste(dir,"Randomsampling/Dataset_1/all/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Randomsampling/Dataset_1/leukemia/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Randomsampling/Dataset_2/all/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Randomsampling/Dataset_2/leukemia/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Randomsampling/Dataset_3/all/", sep = ""), showWarnings = F, recursive = T)

dir.create(paste(dir,"Crossstudy/Dataset_1/all/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossstudy/Dataset_1/leukemia/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossstudy/Dataset_2/all/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossstudy/Dataset_2/leukemia/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossstudy/Dataset_3/all/", sep = ""), showWarnings = F, recursive = T)

dir.create(paste(dir,"Crossplatform/raw/D1_D2/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossplatform/raw/D1_D3/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossplatform/raw/D2_D3/", sep = ""), showWarnings = F, recursive = T)

dir.create(paste(dir,"Crossplatform/rank/D1_D2/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossplatform/rank/D1_D3/", sep = ""), showWarnings = F, recursive = T)
dir.create(paste(dir,"Crossplatform/rank/D2_D3/", sep = ""), showWarnings = F, recursive = T)


# Calculations
# 1.A Random Sampling Dataset 1, all samples
ntest <- 500
for(ntrain in c(100, 250, 500, 1000, 1500, 2000)){
  print(ntrain)
  randomsampling(nperm = nperm_random, # number of permutations
                       info = info.1, # metadata 
                       data = data.1, # data
                       server = linux_run, # for linux-machine = T, for windows = F
                       size.ts = ntrain, # training set size
                       size.vs = ntest, # test set size
                       cores = ncores, # number of cores for parallel processing
                       leukemia =  F,       
                       dir = paste(dir,"Randomsampling/Dataset_1/all/", sep = "")) # output directory
}

# 1.B Random Sampling Dataset 1, leukemia samples
nsamples <- sum(info.1$Disease %in% c("AML", "AMKL", "ALL", "CLL", "CML", "MDS", "DS transient myeloproliferative disorder", "T.ALL"))
ntest <- floor(0.2 * nsamples)

for(ntrain in c(100, 250, 500, 1000, 1500, 1748)){
  print(ntrain)
  randomsampling(nperm = nperm_random, 
                 info = info.1, 
                 data = data.1, 
                 size.ts = ntrain, 
                 size.vs = ntest, 
                 server = linux_run,
                 cores = ncores,
                 dir = paste(dir,"Randomsampling/Dataset_1/leukemia/", sep = ""),
                 leukemia = T) 
}

# 2.A Random Sampling Dataset 2, all samples

nsamples <- sum(info.2$Condition == "CASE") + sum(info.2$Condition == "CONTROL")
ntest <- floor(0.2 * nsamples)
nsamples-ntest
for(ntrain in c(100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000, nsamples-ntest)){
  print(ntrain)
  randomsampling(nperm = nperm_random, 
                 info = info.2, 
                 data = data.2, 
                 size.ts = ntrain, 
                 size.vs = ntest, 
                 leukemia = F,
                 server = linux_run,
                 cores = ncores,
                 dir =  paste(dir,"Randomsampling/Dataset_2/all/", sep = ""))
}

# 2.B Random Sampling Dataset 2, leukemia samples
nsamples <- sum(info.2$Disease %in% c("AML", "AMKL", "ALL", "CLL", "CML", "MDS", "DS transient myeloproliferative disorder", "T.ALL"))
ntest <- floor(0.2 * nsamples)

# blade 10
for(ntrain in c(100, 250, 500, 1000, 2000, 3000, 4000, nsamples-ntest)){
  print(ntrain)
  randomsampling(nperm = nperm_random, 
                 info = info.2, 
                 data = data.2, 
                 size.ts = ntrain, 
                 size.vs = ntest, 
                 leukemia = T,
                 server = linux_run,
                 cores = ncores,
                 dir = paste(dir,"Randomsampling/Dataset_2/leukemia/", sep = ""))
}

# 3.A. Random Sampling Dataset 3, all samples
nsamples <- sum(info.3$Condition == "CASE") + sum(info.3$Condition == "CONTROL")
ntest <- floor(0.2 * nsamples)
for(ntrain in c(100, 250, 500, nsamples-ntest)){
  print(ntrain)
  randomsampling(nperm = nperm_random, 
                 info = info.3, 
                 data = data.3, 
                 size.ts = ntrain, 
                 size.vs = ntest, 
                 classifiers. = c("SVM_linear", "SVM_radial", "SVM_polynomial", "SVM_sigmoid", "PAM", "KNN", "LASSO", "RTF", "LDA"),
                 leukemia = F,
                 server = linux_run,
                 cores = ncores,
                 dir =  paste(dir,"Randomsampling/Dataset_3/all/", sep = ""))
  }

# 2.1.A Cross-study sampling, Dataset 1, all samples

# prepare indices of test and training data
info.1.num <- info.1
info.1.num$Filename <- NULL
info.1.num[] <- lapply(info.1.num[,-8], as.numeric)
info.1.num <- t(info.1.num)
info.1.num <- as.data.frame(info.1.num)

info.1.num<-rbind(info.1.num[3,],info.1.num[-3,]) # change column order
info.1.num<-t(info.1.num) # transpose --> samples in rows
info.1.num[info.1.num[,1]==2,1]<-0 # control = 0, case = 1
n<-nrow(info.1.num) # n = number of samples
n.test<-round(n*0.2) # number of testing samples
study.id<-as.numeric(names(table(info.1.num[,3]))) # study IDs
study.size<-as.numeric(table(info.1.num[,3])) # study sizes
tables<-rbind(study.id,study.size)
colnames(tables)<-as.character(tables[1,])
times<-2000
counts<-matrix(NA,times,23)
shuffle<-matrix(NA,times,23)
for(i in 1:times){
  set.seed(i);shuffle[i,]<-sample(as.numeric(names(table(info.1.num[,3]))),23) # shuffle study IDs
  counts[i,]<-cumsum(tables[2,as.character(shuffle[i,])]) # cumulative sums for each shuffle
}
bingo<-apply(counts<n.test,1,sum)+1 # at what study is the n.test-threshold reached?
study.train.indices<-list() 
study.test.indices<-list()
indices.test<-list()
indices.train<-list()
for(i in 1:times){ # training and testing indices for studies and samples
  study.test.indices[[i]]<-shuffle[i,1:bingo[i]]
  study.train.indices[[i]]<-tables[1,!tables[1,]%in%as.character(study.test.indices[[i]])]
  indices.test[[i]]<-which(info.1.num[,3]%in%study.test.indices[[i]])
  indices.train[[i]]<-which(info.1.num[,3]%in%study.train.indices[[i]])
}
surplus<-lengths(indices.test)-n.test # testing samples bigger than ntest
for(i in which(surplus!=0)){
  indices.test[[i]]<-indices.test[[i]][-sample(1:length(indices.test[[i]]),surplus[i])] # downsample test indices to ntest
}
p<-c()
for(i in 1:times){
  p[i]<-sum(info.1.num[indices.test[[i]],1])/n.test # prevalence of testing set
}
Keep<-sample(which((p>0.3)&(p<0.7)),nperm_crossstudy)  # indices for samples with prevalence between 0.3 and 0.7
indices.test<-indices.test[Keep]
indices.train<-indices.train[Keep]
y.train<-info.1.num[,1]

for(ntrain in c(100, 250, 500, 1000, 1500, 2000)){ 
  print(ntrain)
  crossstudy(nperm=rep, 
             server = linux_run, 
             info = info.1, 
             data = data.1,
             size.ts = ntrain,
             cores = ncores, 
             indices.test = indices.test,
             indices.train = indices.train,
             dir = paste(dir,"Crosstudy/Dataset_1/all/", sep = ""))
}



# 2.1.A Cross-study Sampling Dataset 1, leukemia samples
info.1.num <- info.1[rownames(info.1[info.1$Disease %in% c("MDS", "AMKL", "CLL", "CML", "ALL", "AML", "DS transient myeloproliferative disorder", "T.ALL"),]),]
info.1.num$Filename <- NULL
info.1.num[] <- lapply(info.1.num[,-8], as.numeric)
info.1.num <- t(info.1.num)
info.1.num <- as.data.frame(info.1.num)

info.1.num<-rbind(info.1.num[3,],info.1.num[-3,])
info.1.num<-t(info.1.num)
info.1.num[info.1.num[,1]==2,1]<-0
table(info.1.num[,4])[2]==table(info.1.num[,1])[2]                      ## Check if TRUE
sum(info.1.num[info.1.num[,1]==1,4]==3)==table(info.1.num[,1])[2]            ## Check if TRUE
n<-nrow(info.1.num)
n.test<-round(n*0.2)
study.id<-as.numeric(names(table(info.1.num[,3])))
study.size<-as.numeric(table(info.1.num[,3]))
tables<-rbind(study.id,study.size)
colnames(tables)<-as.character(tables[1,])
times<-2000
counts<-matrix(NA,times,19)
shuffle<-matrix(NA,times,19)
for(i in 1:times){
  set.seed(i);shuffle[i,]<-sample(as.numeric(names(table(info.1.num[,3]))),19)
  counts[i,]<-cumsum(tables[2,as.character(shuffle[i,])])
}
bingo<-apply(counts<n.test,1,sum)+1
study.train.indices<-list()
study.test.indices<-list()
indices.test<-list()
indices.train<-list()
for(i in 1:times){
  study.test.indices[[i]]<-shuffle[i,1:bingo[i]]
  study.train.indices[[i]]<-tables[1,!tables[1,]%in%as.character(study.test.indices[[i]])]
  indices.test[[i]]<-which(info.1.num[,3]%in%study.test.indices[[i]])
  indices.train[[i]]<-which(info.1.num[,3]%in%study.train.indices[[i]])
}
surplus<-lengths(indices.test)-n.test
for(i in which(surplus!=0)){
  indices.test[[i]]<-indices.test[[i]][-sample(1:length(indices.test[[i]]),surplus[i])]
}
p<-c()
for(i in 1:times){
  p[i]<-sum(info.1.num[indices.test[[i]],1])/n.test
}
Keep<-sample(which((p>0.3)&(p<0.7)),nperm_crossstudy)  # indices for samples with prevalence between 0.3 and 0.7
indices.test<-indices.test[Keep]
indices.train<-indices.train[Keep]
y.train<-info.1.num[,1]

for(ntrain in c(50, 100, 250, 500, 1000, 1500)){ 
  print(ntrain)
  crossstudy(nperm=nperm_crossstudy, 
                           server = linux_run, 
                           info = info.1, 
                           data = data.1,
                           size.ts = ntrain,
                           print_values = T, 
                           cores = ncores, 
                           indices.test = indices.test,
                           indices.train = indices.train,
                           dir = paste(dir, "Crossstudy/Dataset_1/leukemia", sep = ""))
}


# 2.2.B Cross-study Sampling Dataset 2, all samples

info.2.num <- info.2
info.2.num$Filename <- NULL
info.2.num[] <- lapply(info.2.num[,-8], as.numeric)
info.2.num <- t(info.2.num)
info.2.num <- as.data.frame(info.2.num)


info.2.num<-rbind(info.2.num[3,],info.2.num[-3,]) # change order: 
info.2.num<-t(info.2.num) # transpose --> samples in rows
info.2.num[info.2.num[,1]==2,1]<-0 # control = 0, case = 1
n<-nrow(info.2.num) # n = number of samples
n.test<-round(n*0.2) # number of testing samples
study.id<-as.numeric(names(table(info.2.num[,3]))) # study IDs
study.size<-as.numeric(table(info.2.num[,3])) # study sizes
tables<-rbind(study.id,study.size)
colnames(tables)<-as.character(tables[1,])
times<-2000
counts<-matrix(NA,times,64)
shuffle<-matrix(NA,times,64)
for(i in 1:times){
  set.seed(i);shuffle[i,]<-sample(as.numeric(names(table(info.2.num[,3]))),64) # shuffle study IDs
  counts[i,]<-cumsum(tables[2,as.character(shuffle[i,])]) # cumulative sums for each shuffle
}
bingo<-apply(counts<n.test,1,sum)+1 # at what study is the n.test-threshold reached?
study.train.indices<-list() 
study.test.indices<-list()
indices.test<-list()
indices.train<-list()
for(i in 1:times){ # training and testing indices for studies and samples
  study.test.indices[[i]]<-shuffle[i,1:bingo[i]]
  study.train.indices[[i]]<-tables[1,!tables[1,]%in%as.character(study.test.indices[[i]])]
  indices.test[[i]]<-which(info.2.num[,3]%in%study.test.indices[[i]])
  indices.train[[i]]<-which(info.2.num[,3]%in%study.train.indices[[i]])
}
surplus<-lengths(indices.test)-n.test # testing samples bigger than ntest
for(i in which(surplus!=0)){
  indices.test[[i]]<-indices.test[[i]][-sample(1:length(indices.test[[i]]),surplus[i])] # downsample test indices to ntest
}
p<-c()
for(i in 1:times){
  p[i]<-sum(info.2.num[indices.test[[i]],1])/n.test # prevalence of testing set
}
Keep<-sample(which((p>0.3)&(p<0.7)),nperm_crossstudy)  # indices for samples with prevalence between 0.3 and 0.7
indices.test<-indices.test[Keep]
indices.train<-indices.train[Keep]
y.train<-info.2.num[,1]


for(ntrain in c(100, 250, 500, 1000, 1500, 2000, 4000, 6000)){ 
  print(ntrain)
  crossstudy(nperm = nperm_crossstudy, 
                           server = linux_run, 
                           classifiers. = c("SVM_linear", "SVM_radial", "SVM_polynomial", "SVM_sigmoid", "PAM", "LDA", "KNN", "LASSO", "RTF"),
                           info = info.2, 
                           data = data.2,
                           size.ts = ntrain,
                           cores = ncores, 
                           indices.test = indices.test,
                           indices.train = indices.train,
                           dir = paste(dir, "Crossstudy/Dataset_2/all", sep = ""))
}


# 2.2.B Cross-study Sampling Dataset 2, leukemia samples

info.2.num <- info.2[rownames(info.2[info.2$Disease %in% c("MDS", "AMKL", "CLL", "CML", "ALL", "AML", "DS transient myeloproliferative disorder", "T.ALL"),]),]
info.2.num$Filename <- NULL
info.2.num[] <- lapply(info.2.num[,-8], as.numeric)
info.2.num <- t(info.2.num)
info.2.num <- as.data.frame(info.2.num)

info.2.num<-rbind(info.2.num[3,],info.2.num[-3,]) # change order
info.2.num<-t(info.2.num) # transpose --> samples in rows
info.2.num[info.2.num[,1]==2,1]<-0 # control = 0, case = 1
n<-nrow(info.2.num) # n = number of samples
n.test<-round(n*0.2) # number of testing samples
study.id<-as.numeric(names(table(info.2.num[,3]))) # study IDs
study.size<-as.numeric(table(info.2.num[,3])) # study sizes
tables<-rbind(study.id,study.size)
colnames(tables)<-as.character(tables[1,])
times<-2000
counts<-matrix(NA,times,42)
shuffle<-matrix(NA,times,42)
for(i in 1:times){
  set.seed(i);shuffle[i,]<-sample(as.numeric(names(table(info.2.num[,3]))),42) # shuffle study IDs
  counts[i,]<-cumsum(tables[2,as.character(shuffle[i,])]) # cumulative sums for each shuffle
}
bingo<-apply(counts<n.test,1,sum)+1 # at what study is the n.test-threshold reached?
study.train.indices<-list() 
study.test.indices<-list()
indices.test<-list()
indices.train<-list()
for(i in 1:times){ # training and testing indices for studies and samples
  study.test.indices[[i]]<-shuffle[i,1:bingo[i]]
  study.train.indices[[i]]<-tables[1,!tables[1,]%in%as.character(study.test.indices[[i]])]
  indices.test[[i]]<-which(info.2.num[,3]%in%study.test.indices[[i]])
  indices.train[[i]]<-which(info.2.num[,3]%in%study.train.indices[[i]])
}
surplus<-lengths(indices.test)-n.test # testing samples bigger than ntest
for(i in which(surplus!=0)){
  indices.test[[i]]<-indices.test[[i]][-sample(1:length(indices.test[[i]]),surplus[i])] # downsample test indices to ntest
}
p<-c()
for(i in 1:times){
  p[i]<-sum(info.2.num[indices.test[[i]],1])/n.test # prevalence of testing set
}
Keep<-sample(which((p>0.3)&(p<0.7)),nperm_crossstudy)  # indices for samples with prevalence between 0.3 and 0.7
indices.test<-indices.test[Keep]
indices.train<-indices.train[Keep]
y.train<-info.2.num[,1]

for(ntrain in c(100, 250, 500, 1000, 1500, 2000, 4000)){ 
  print(ntrain)
  crossstudy(nperm=nperm_crossstudy, 
             server = linux_run, 
             info = info.2, 
             data = data.2,
             size.ts = ntrain,
             cores = ncores, 
             indices.test = indices.test,
             indices.train = indices.train,
             dir = paste(dir, "Crossstudy/Dataset_2/leukemia", sep = ""))
}

# 2.3. Cross-study Sampling Dataset 3, all samples

info.3.num <- info.3
info.3.num$Filename <- NULL
info.3.num[] <- lapply(info.3.num[,-8], as.numeric)
info.3.num <- t(info.3.num)
info.3.num <- as.data.frame(info.3.num)

info.3.num<-t(info.3.num)
info.3.num[info.3.num[,3]==2,3]<-0
n<-nrow(info.3.num)
n.test<-round(n*0.2)
study.id<-as.numeric(names(table(info.3.num[,2])))
study.size<-as.numeric(table(info.3.num[,2]))
tables<-rbind(study.id,study.size)
colnames(tables)<-as.character(tables[1,])
times<-20000
counts<-matrix(NA,times,23)
shuffle<-matrix(NA,times,23)
for(i in 1:times){
  set.seed(i);shuffle[i,]<-sample(as.numeric(names(table(info.3.num[,2]))),23)
  counts[i,]<-cumsum(tables[2,as.character(shuffle[i,])])
}
bingo<-apply(counts<n.test,1,sum)+1
study.train.indices<-list()
study.test.indices<-list()
indices.test<-list()
indices.train<-list()
for(i in 1:times){
  study.test.indices[[i]]<-shuffle[i,1:bingo[i]]
  study.train.indices[[i]]<-tables[1,!tables[1,]%in%as.character(study.test.indices[[i]])]
  indices.test[[i]]<-which(info.3.num[,2]%in%study.test.indices[[i]])
  indices.train[[i]]<-which(info.3.num[,2]%in%study.train.indices[[i]])
}
surplus<-lengths(indices.test)-n.test
for(i in which(surplus!=0)){
  indices.test[[i]]<-indices.test[[i]][-sample(1:length(indices.test[[i]]),surplus[i])]
}
times<-2000
Keep<-sample(which(lengths(indices.train)>700),times)
indices.test<-indices.test[Keep]
indices.train<-indices.train[Keep]
p<-c()
for(i in 1:times){
  p[i]<-sum(info.3.num[indices.test[[i]],3])/n.test
}
Keep<-sample(which((p>0.3)&(p<0.7)),nperm_crossstudy) # indices for samples with prevalence between 0.3 and 0.7
indices.test<-indices.test[Keep]
indices.train<-indices.train[Keep]
y.train<-info.3.num[,1]


for(ntrain in c(100, 250, 500, 1000)){ 
  print(ntrain)
  permute_ts_vs_crossstudy(nperm = nperm_crossstudy, 
                           server = linux_run, 
                           info = info.3, 
                           data = data.3,
                           size.ts = ntrain,
                           cores = ncores, 
                           indices.test = indices.test,
                           indices.train = indices.train,
                           dir = paste(dir, "Crossstudy/Dataset_3/all", sep = ""))
}

# 3.Cross-Platform Sampling
# 3A on raw (normalized) data
for(ntrain in c(100, 250, 500, 1000, 2000)){
  print(ntrain)
  crossplatform(nperm = nperm_crossplatform,  
                              classifiers = c("LDA"),
                              size.ts = ntrain, 
                              info.ts = info.1, 
                              data.ts = data.1,
                              info.vs = info.2, 
                              data.vs = data.2,
                              size.vs = ncol(data.2),
                              server = linux_run,
                              cores = ncores,
                              dir = paste(dir, "Crossplatform/raw/D1_D2", sep = ""))
}

crossplatform(nperm = 1,  
                            size.ts = 2500, 
                            info.ts = info.1, 
                            data.ts = data.1,
                            info.vs = info.2, 
                            data.vs = data.2,
                            size.vs = ncol(data.2),
                            server = linux_run,
                            cores = 1,
                            dir = paste(dir, "Crossplatform/raw/D1_D2", sep = ""))


### TS 3.2. Dataset 2, VS Dataset 3 

for(ntrain in c(100, 250, 500, 1000, 2000, 4000, 6000, 8000)){
  print(ntrain)
  crossplatform(nperm = nperm_crossplatform,  
                              size.ts = ntrain, 
                              info.ts = info.2, 
                              data.ts = data.2,
                              info.vs = info.3, 
                              data.vs = data.3,
                              size.vs = ncol(data.3),
                              server = linux_run,
                              cores = ncores,
                              dir = paste(dir, "Crossplatform/raw/D2_D3", sep = ""))
}

permute_ts_vs_crossplatform(nperm = 1,  
                            size.ts = ncol(data.2), 
                            info.ts = info.2, 
                            data.ts = data.2,
                            info.vs = info.3, 
                            data.vs = data.3,
                            size.vs = ncol(data.3),
                            server = linux_run,
                            cores = 1,
                            dir = paste(dir, "Crossplatform/raw/D2_D3", sep = ""))

### 3.3.  TS Dataset 1, VS Dataset 3 
for(ntrain in c(100, 250, 500, 1000, 2000)){
  print(ntrain)
  permute_ts_vs_crossplatform(nperm = nperm_crossplatform,  
                              size.ts = ntrain, 
                              info.ts = info.1, 
                              data.ts = data.1,
                              info.vs = info.3, 
                              data.vs = data.3,
                              size.vs = ncol(data.3),
                              server = linux_run,
                              cores = ncores,
                              dir = paste(dir, "Crossplatform/raw/D1_D3", sep = ""))
}

permute_ts_vs_crossplatform(nperm = 1,  
                            size.ts = 2500, 
                            info.ts = info.1, 
                            data.ts = data.1,
                            info.vs = info.3, 
                            data.vs = data.3,
                            size.vs = ncol(data.3),
                            server = linux_run,
                            cores = 1,
                            dir = paste(dir, "Crossplatform/raw/D1_D3", sep = ""))





## Prepare Data
data.1.rank <- t(data.1)
data.1.rank <- apply(data.1.rank,2,rank)/(nrow(data.1.rank)+1)
data.1.rank <- apply(data.1.rank,2,qnorm)
data.1.rank <- t(data.1.rank)

data.2.rank <- t(data.2)
data.2.rank <- apply(data.2.rank,2,rank)/(nrow(data.2.rank)+1)
data.2.rank <- apply(data.2.rank,2,qnorm)
data.2.rank <- t(data.2.rank)

data.3.rank <- t(data.3)
data.3.rank <- apply(data.3.rank,2,rank)/(nrow(data.3.rank)+1)
data.3.rank <- apply(data.3.rank,2,qnorm)
data.3.rank <- t(data.3.rank)

## run calculations
### 3.1. TS Dataset 1, VS Dataset 2 

for(ntrain in c(100, 250, 500, 1000, 2000)){
  print(ntrain)
  permute_ts_vs_crossplatform(nperm = nperm_crossplatform,  
                              size.ts = ntrain, 
                              info.ts = info.1, 
                              data.ts = data.1.rank,
                              info.vs = info.2, 
                              data.vs = data.2.rank,
                              size.vs = ncol(data.2.rank),
                              server = linux_run,
                              cores = ncores,
                              dir = paste(dir, "Crossplatform/rank/D1_D2", sep = ""))
}

permute_ts_vs_crossplatform(nperm = 1,  
                            size.ts = 2500, 
                            info.ts = info.1, 
                            data.ts = data.1.rank,
                            info.vs = info.2, 
                            data.vs = data.2.rank,
                            size.vs = ncol(data.2.rank),
                            server = linux_run,
                            cores = 1,
                            dir = paste(dir, "Crossplatform/rank/D1_D2", sep = ""))


### TS 3.2. Dataset 2, VS Dataset 3 

for(ntrain in c(100, 250, 500, 1000, 2000, 4000, 6000, 8000)){
  print(ntrain)
  permute_ts_vs_crossplatform(nperm = nperm_crossplatform,  
                              size.ts = ntrain, 
                              info.ts = info.2, 
                              data.ts = data.2.rank,
                              info.vs = info.3, 
                              data.vs = data.3.rank,
                              size.vs = ncol(data.3.rank),
                              server = linux_run,
                              cores = ncores,
                              dir = paste(dir, "Crossplatform/rank/D2_D3", sep = ""))
}

permute_ts_vs_crossplatform(nperm = 1,  
                            size.ts = ncol(data.2.rank), 
                            info.ts = info.2, 
                            data.ts = data.2.rank,
                            info.vs = info.3, 
                            data.vs = data.3.rank,
                            size.vs = ncol(data.3.rank),
                            server = linux_run,
                            cores = 1,
                            dir = paste(dir, "Crossplatform/rank/D2_D3", sep = ""))

### 3.3.  TS Dataset 1, VS Dataset 3 
for(ntrain in c(100, 250, 500, 1000, 2000)){
  print(ntrain)
  permute_ts_vs_crossplatform(nperm = nperm_crossplatform,  
                              size.ts = ntrain, 
                              info.ts = info.1, 
                              data.ts = data.1.rank,
                              info.vs = info.3, 
                              data.vs = data.3.rank,
                              size.vs = ncol(data.3.rank),
                              server = linux_run,
                              cores = ncores,
                              dir = paste(dir, "Crossplatform/rank/D1_D3", sep = ""))
}

permute_ts_vs_crossplatform(nperm = 1,  
                            size.ts = 2500, 
                            info.ts = info.1, 
                            data.ts = data.1.rank,
                            info.vs = info.3, 
                            data.vs = data.3.rank,
                            size.vs = ncol(data.3.rank),
                            server = linux_run,
                            cores = 1,
                            dir = paste(dir, "Crossplatform/rank/D1_D3", sep = ""))





