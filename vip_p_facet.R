
library(ggplot2)
library(ggsci)

args<-commandArgs(TRUE)

setwd(args[1])

df <- read.csv(file="vip_p_facet.tsv",header= T,sep= "\t",check.names=FALSE)

colnames(df)[3]="p"

anno=unique(df[,colnames(df) %in% c("group","g1_from","g2_deep")])

p <- ggplot(df, aes(x=p, y=vip, color=FC)) +

 	geom_point(alpha=0.8, size=1.5, stroke=0) +
 	
 	xlim(0,0.5) +
 	
 	scale_color_gradientn(limits = c(-5, +5), 
                       breaks = c(-3, 0, +3), 
                       colours = c("#E6550DFF", "#79AF97FF", "#3182BDFF")) +
 	#scale_color_npg() +
 	#scale_fill_npg() +
   	
   	labs(y="VIP", x="p-value", color="FC") +
   	
   	geom_hline(yintercept = 1,linetype=4,size=0.8) +
  
   	geom_vline(xintercept = 0.05,linetype=4,size=0.8) + 
   	
   	geom_text(data=anno, aes(x=max(df$p)/4, y=max(df$vip), label=group), 
   	 	 	hjust=0.5, color="black", fontface="bold",alpha=0.6, size=2, inherit.aes = FALSE, family="Times") + 
   	
   	ggtitle("VIP-P") +
 	
 	facet_grid(g2_deep ~ g1_from) +
 	
 	theme(legend.position = "right",text=element_text(family="serif"), plot.title = element_text(hjust = 0.5))

ggsave(p,file="vip_p_facet.svg")

