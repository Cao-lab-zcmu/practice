
args<-commandArgs(TRUE)

setwd(args[1])

datas=list.files(path = ".", pattern = "*pca.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)

for(data in datas){

savename=strsplit(data,split="_pca.tsv")

data <- read.table(file=data,header= T,row.names= 1,sep= "\t",check.names=F)

pca<- prcomp(data, scale. = TRUE)

pca_x=data.frame(pca$x)

rownames=rownames(pca$x)

pca_x=cbind(rownames,pca_x)

pca=summary(pca)

summary=c("summary",round(pca$importance[2,],4))

pca_x=rbind(pca_x,summary)

write.table(pca_x, file = paste0(savename,".prcomp"), quote = FALSE, append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

}
           
           
           
