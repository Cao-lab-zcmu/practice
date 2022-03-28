library(ggbiplot)
library(ggsci)
library(scales)
library(ggrepel)

args<-commandArgs(TRUE)

setwd(args[1])

datas=list.files(path = ".", pattern = "*pca.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)

metadatas=list.files(path = ".", pattern = "^metadata*", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
n=length(datas)

for(i in 1:n){
savename=strsplit(datas[i],split=".tsv")
source<- read.table(file= datas[i],header= T,row.names= 1,sep= "\t")
class<- read.table(file= metadatas[i],header= T,row.names= 1,sep= "\t")
pca<- prcomp(source, scale. = TRUE)
pca_anno <- as.data.frame(pca[5])

output <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = class$subgroup, ellipse = TRUE, circle = TRUE,
 	 	 	varname.size=0, var.axes = F) +
 	 	geom_label_repel(data=pca_anno, aes(x=x.PC1, y=x.PC2, label=class$name),
	  	 	          color="black", alpha=0.5, fontface="bold", size=2, angle= 0, 
	  	 	          direction="both", segment.size = 0.2, segment.alpha = 0.3,
	  	 	          inherit.aes = FALSE, hjust = 0,
	  	 	          family="Times") +
 	 	scale_color_npg() +
 	 	scale_fill_npg() +
 	 	theme(legend.position = "right",text=element_text(family="serif"))

ggsave(output,file=paste0(savename,".svg"))

}
