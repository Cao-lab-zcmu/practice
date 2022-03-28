library(ropls)
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

savename=strsplit(datas[i],split="_pca.tsv")
data <- read.table(file=datas[i],header= T,row.names= 1,sep= "\t",check.names=F)
class <- read.table(file=metadatas[i],header= T,row.names= 1,sep= "\t",check.names=F)
oplsda <- opls(x = data, y = class[, "group"], predI = 1, orthoI = NA)
df <- cbind(oplsda@scoreMN[, 1], oplsda@orthoScoreMN[, 1])
colnames(df) <- c("h1", paste0("o", 1))
x_lab <- paste0("T score[1](", oplsda@modelDF[1, "R2X"] * 100, "%)")
y_lab <- paste0("Orthogonal T score[1](", oplsda@modelDF[2, "R2X"] * 100, "%)")
df <- as.data.frame(df)
vip=data.frame(oplsda@vipVn)
vip=cbind(rownames(vip),vip)
colnames(vip)=c("id","vip")
p <- ggplot(df, aes(x=h1, y=o1)) +
 	geom_point(alpha=0.8, size=3, shape=21, stroke=0.1,aes(fill=class$subgroup)) +
 	stat_ellipse(aes(color=class$group), level = 0.95) +
 	geom_label_repel(data=df, aes(x=h1, y=o1, label=class$name),
	  	 	          color="black", alpha=0.5, fontface="bold", size=2, angle= 0, 
	  	 	          direction="both", segment.size = 0.2, segment.alpha = 0.3,
	  	 	          inherit.aes = FALSE, hjust = 0,
	  	 	          family="Times") +
 	scale_color_npg() +
 	scale_fill_npg() +
 	labs(x=x_lab,y=y_lab,title="OPLS-DA") +
 	theme(plot.title = element_text(hjust = 0.5),text=element_text(family="serif")) 
 	
ggsave(p,file=paste0(savename,"_oplsda.svg"), height=6, width=10)

df=cbind(rownames(df),df)

colnames(df)=c("sample","h1","o1")

df=rbind(df,c("summary",paste0("t1(", oplsda@modelDF[1, "R2X"] * 100, "%)"),paste0("ot1(", oplsda@modelDF[2, "R2X"] * 100, "%)")))

write.table(vip, file = paste0(savename,".vip"), quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)

write.table(df, file = paste0(savename,".ropls"), quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
}


