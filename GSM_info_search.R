library(GEOquery)
library(dplyr)
library(plyr)

load("E:/Stefanie/Classifier/Paper/Scripts/Datasets.RData")

# Dataset 1
info.1.char <- info.1 %>%
  mutate_all(as.character)

info.1.char$Filename <- unlist(lapply(strsplit(info.1.char$Filename, split = "\\."), `[[`, 1)) # remove .CEL


for(i in 1:length(levels(droplevels(info.1$GSE)))){
  gse <- getGEO(levels(info.1$GSE)[i],GSEMatrix=FALSE)
  
  for(j in 1:length(names(GSMList(gse)))){
    tmp <- as.matrix(lapply(GSMList(gse)[[j]]@header, paste, collapse = ";"))
    tmp <- as.data.frame(t(matrix(unlist(tmp), nrow=length(tmp), byrow=T)))
    colnames(tmp) <- attr(as.matrix(lapply(GSMList(gse)[[j]]@header, paste, collapse = ";")), "dimnames")[[1]]
    res <- rbind.fill(tmp, res)
  }

res <- merge(info.1.char[info.1.char$GSE %in% levels(droplevels(info.1$GSE))[i],], res, by.x = "Filename", by.y = "geo_accession")  
  
write.table(res, paste("E:/Stefanie/Classifier/AML_subtypes/Series_info/", levels(droplevels(info.1$GSE))[i], ".txt", sep = ""), sep = "\t")
  
}




