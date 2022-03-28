
library(ggplot2)
library(ggsci)

args<-commandArgs(TRUE)

setwd(args[1])

df <- read.csv(file= "pca_facet.tsv",header= T,sep= "\t")

PC1 <- df[which(df$PC1=="PC1"),]

PC1$figure <- paste0("PC1(",PC1$figure*100,"%)")

PC2 <- df[which(df$PC1=="PC2"),]

PC2$figure <- paste0("PC2(",PC2$figure*100,"%)")

df <- df[which(df$sample!="anno"),c(1:6)]

anno_x=min(as.numeric(df$PC1))*(20/20)

anno_y=max(as.numeric(df$PC2))*(22/20)

p <- ggplot(df, aes(x=as.numeric(PC1), y=as.numeric(PC2), fill=subgroup)) +
 	geom_point(alpha=0.8, size=3, shape=21, stroke=0.1) +
 	stat_ellipse(aes(color=subgroup), level = 0.7) +
 	#scale_color_npg() +
 	#scale_fill_npg() +
 	scale_color_manual(values = c("control"="grey","model"="#374E55FF",
  
   	 	 	 	"pro_low"="#FDAE6BFF","pro_medium"="#FD8D3CFF","pro_high"="#E6550DFF",
   	 	 	 	
   	 	 	 	"raw_low"="#9ECAE1FF","raw_medium"="#6BAED6FF","high_raw"="#3182BDFF")) +
   	
   	scale_fill_manual(values = c("control"="grey","model"="#374E55FF",
  
   	 	 	 	"pro_low"="#FDAE6BFF","pro_medium"="#FD8D3CFF","pro_high"="#E6550DFF",
   	 	 	 	
   	 	 	 	"raw_low"="#9ECAE1FF","raw_medium"="#6BAED6FF","raw_high"="#3182BDFF")) +
   	
   	guides(color= "none") +
   	
   	geom_text(data=PC1, aes(x=anno_x, y=anno_y, label=figure), 
   	 	 	hjust=0, color="black", fontface="bold",alpha=0.6, size=1.5, inherit.aes = FALSE, family="Times") +
   	
   	geom_text(data=PC2, aes(x=anno_x, y=anno_y*(18/22), label=figure), 
   	 	 	hjust=0, color="black", fontface="bold",alpha=0.6, size=1.5, inherit.aes = FALSE, family="Times") +
   	
   	labs(y="PC2", x="PC1", fill="Group") +
   	
   	#scale_x_continuous(limits=c(-60 ,60)) +
   	
   	#scale_y_continuous(limits=c(-60 ,60)) +
   	
 	facet_grid(g2_deep ~ g1_from) +
 	
 	theme(legend.position = "right",text=element_text(family="serif"))

ggsave(p,file="pca_facet.svg",width=8,height=6.5)



