library(ggsci)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggforce)
library(ggpmisc)

setwd("linear_regression")

files=list.files(path = ".", pattern = "*.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)

for(file in files){

savename=strsplit(file,split=".tsv")

data <- read.csv(file=file,header=T,sep="\t")

p <- ggplot(data, aes(x=chem, y=id_content)) +

  geom_smooth(method = "lm", se=FALSE, color = "black", alpha=0.3, size=0.4, formula = y ~ x) +
  
  geom_point(shape=21, stroke=0, size=4, aes(fill=subgroup)) +
  
  labs(fill="Group") +
  
  coord_cartesian(clip = "off") +
  
  annotate("text", x = min(data$chem), y = max(data$id_content)*(16/20), 
   	    label = paste0("Linear correlation: ",as.character(savename)),
            color="black",size = 4, fontface="bold", family="Times", hjust = 0, alpha=0.8) +
  
  annotate("text", x = min(data$chem), y = max(data$id_content)*(15/20), 
   	    label = paste0("r(pearson): ",as.character(round(data[1,colnames(data) %in% c("r.pearson.")],4))),
            color="black",size = 4, fontface="bold", family="Times", hjust = 0, alpha=0.8) +
  
  annotate("text", x = min(data$chem), y = max(data$id_content)*(14/20), 
   	    label = paste0("p-value: ",as.character(round(data[1,colnames(data) %in% c("p.value")],4))),
            color="black",size = 4, fontface="bold", family="Times", hjust = 0, alpha=0.8) +
  
  scale_size_continuous( trans="exp", range=c(2, 5)) +
  
  scale_fill_manual(values = c("control"="#4A6990FF","drug"="#95CC5EFF","model"="black",
  
   	 	 	 	"pro_low"="#FDAE6BFF","pro_medium"="#FD8D3CFF","pro_high"="#E6550DFF",
   	 	 	 	
   	 	 	 	"raw_low"="#9ECAE1FF","raw_medium"="#6BAED6FF","raw_high"="#3182BDFF")) +
  
  theme(text=element_text(family="serif"))
  
ggsave(p,file=paste0(savename,".svg"),width=8,height=6.5)

print(file)

}
