
library(ggplot2)
library(ggsci)

args<-commandArgs(TRUE)

setwd(args[1])

datas=list.files(path = ".", pattern = "*.vip_p$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)

for(i in 1:length(datas)){

savename=strsplit(datas[i],split=".vip_p")

df <- read.csv(file=datas[i],header= T,sep= "\t",check.names=FALSE)

truenames=colnames(df)

set=unlist(strsplit(truenames[4],split=": "))

legend=paste0("VIP-P-FC: ",set[2])

colnames(df)=c("id","vip","p","fold")

p <- ggplot(df, aes(x=p, y=vip, color=fold)) +

 	geom_point(alpha=0.8, size=3, stroke=0) +
 	
 	scale_color_gradientn(limits = c(-5, +5), 
                       breaks = c(-3, 0, +3),
                       colours = c("#E6550DFF", "#79AF97FF", "#3182BDFF")) +
 	#scale_color_npg() +
 	#scale_fill_npg() +
   	
   	labs(y="VIP", x="p-value", color="FC") +
   	
   	geom_hline(yintercept = 1,linetype=4,size=0.8) +
  
   	geom_vline(xintercept = 0.05,linetype=4,size=0.8) + 
   	
   	ggtitle(legend) +
   	
   	#xlim(0,0.5) +
   	
   	#scale_x_continuous(limits=c(-60 ,60)) +
   	
   	#scale_y_continuous(limits=c(-60 ,60)) +
   	
 	#facet_grid(g2_deep ~ g1_from) +
 	
 	theme(legend.position = "right",text=element_text(family="serif"), plot.title = element_text(hjust = 0.5))

ggsave(p,file=paste0(savename,"_vip_p.svg"))
}

