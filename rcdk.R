###########  R rcdk structure 

####  wd ~/operation/back/0703_all/results

library(tidyverse)

library(rcdk)

data <- read.csv(file="fingerid_first_score.tsv", header=T,sep="\t")

data <- data[which(is.na(data$smiles)==F), colnames(data) %in% c("id", "smiles")]

for(i in 1:nrow(data)){

id <- data[i, 1]

stru <- data[i, 2]

mols <- parse.smiles(stru)

#mols <- generate.2d.coordinates(mols[[1]])

#mols <- get.smiles(mols,smiles.flavors(c('CxSmiles')))

write.molecules(mols, paste0("structure_2d/smiles_draw/", id, "_"), together = F, write.props = F)

}
