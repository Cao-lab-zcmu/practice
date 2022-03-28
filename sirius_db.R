library(tidyverse)
path1="cnumber_cid.tsv"
path2="merge_smiles.tsv"
df1 <- read.csv(file=path1, sep="\t", header=T)
df2 <- read.csv(file=path2, sep="\t", header=T)
df <- merge(df1, df2, all.y=T, by.x="pubchem.id", by.y="CID", sort=T)
df <- df[, c(3,2,1)]
write.table(df, file="sirius_db.tsv", sep="\t", col.names=F, row.names=F, quote=F)

