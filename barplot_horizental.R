 library(ggplot2)
 library(ggrepel)
 library(ggsci)
 library(Hmisc)
 args<-commandArgs(TRUE)
 setwd(args[1])
 file=args[2]
 savename=strsplit(file,split=".tsv")
 df <- read.csv(file=file,header= T,sep= "\t")
 df[,1] <- capitalize(df[,1])
 annonation <- df[which(df[,1]=="Annonation"),]
 df <- df[which(df[,1]!="Annonation"),]
 df[,2] <- as.numeric(df[,2])
 #order=df[order(df[,2],decreasing=F),]
 p <- ggplot(data=df, aes(x=reorder(df[,1],df[,2]), y=df[,2], fill=df[,2])) +
   geom_col(width = 0.7) +
   scale_fill_gradientn(colours = c("#FFBB78FF","#FF7F0EFF")) +
   coord_flip() +
   guides(fill="none") +
   labs(y=paste0("Number(", annonation[,2]), x="FC group", ")", title="Stat features") +
   geom_text(data=df, aes(x=df[,1], y=df[,2]+max(df[,2])*(1/30), label=df[,2]), 
             hjust=0, color="black", family="Times", alpha=0.6, size=5, inherit.aes = FALSE ) +
theme(legend.position = "right",text=element_text(family="serif", size=20), plot.title = element_text(hjust = 0.3))
 ggsave(p,file=paste0(savename,".svg"), width=15, height=7)

