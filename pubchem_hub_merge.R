library(tidyverse)
data_list <- list.files(path = ".", pattern = "*smiles.csv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
df <- data.frame()
for(i in data_list){
  data <- read.csv(file=i, header=T, check.names=F, sep=",")
  df <- rbind(df, data)
}
file <- read.csv(file="cid_metadata.tsv", check.names=F, sep="\t", header=T)
df <- merge(df, file, by.x="CID", by.y="cid", all.x=T, sort=T)
#file2 <- read.csv(file="../fingerid_first_score.tsv", check.names=F, sep="\t", header=T)
file2 <- read.csv(file="fingerid_first_score.tsv", check.names=F, sep="\t", header=T)
df <- merge(df, file2, all.x=T, by="id", sort=T)
df <- df[!duplicated(df$id), ]
write.table(df, file="merge.tsv", col.names=T, row.names=F, quote=F, sep="\t")

