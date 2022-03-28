 library(ggplot2)
 library(ggrepel)
 library(ggsci)
 library(tidyr)
 library(ggforce)
 library(reshape2)
 library(ggthemes)
 library(stringr)
 file="for_violin_0.5.tsv"
 dfpath="fingerid_first_score.tsv"
 df <- read.csv(file=dfpath, header=T, sep="\t")
 df <- df[,colnames(df) %in% c("id", "similarity")]
 data <- read.csv(file=file,header=T,sep="\t")
 colnames(data)[3:4]=c("Raw-Eucommia", "Pro-Eucommia")
 data <- melt(data, id.vars=c("id","classification"), variable.name="condition", value.name="expr")
 data$expr[which(data$expr==-Inf)] <- 0
 data <- merge(data, df, all.x=T, by="id", sort=F)
 data <- data[which(data$similarity>0.5), ]
 p <- ggplot(data, aes(x=classification, y=expr, fill=ifelse(condition=="Raw-Eucommia", "1", "2"))) +
  	geom_violin(trim=F,color="transparent") +
  	geom_boxplot(width=0.2,position=position_dodge(0.9)) +
  	guides(fill="none") +
  	coord_flip() +
  	scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF")) +
  	scale_x_discrete(limits = rev(levels(as.factor(data$classification)))) +
  	labs(x="Classification", y="Log10(peak area)(Tanimoto similarity > 0.5; Posterior probability > 0.5)") +
  	theme(text=element_text(family="Times"), 
  	 	axis.title = element_text(face="bold"),
  	 	axis.title.x = element_text(hjust=1.7),
  	 	axis.text.x = element_text(angle=0),
  	 	axis.text = element_text(size=12)) +
  	facet_wrap(~condition, ncol=2)
  ggsave(p, file="violin.svg", width=8, height=18)
  #ggsave(p, file="violin.svg", width=10, height=15)
  ####################  geom_density
  list=grep("lignan|Iridoid", data$classification, ignore.case=T)
  df <- data[list,]
  df$classification <- str_wrap(df$classification, width=25)
  pp <- ggplot() +
   	geom_density(data=df, aes(x=expr, fill=condition), alpha=0.3) + 
   	facet_grid(classification~., 
   	 	   scales="free"
   	) +
   	scale_fill_manual(values = c("#4DBBD5FF", "#DC0000FF")) +
   	theme_classic() +
   	labs(x="Log10(peak area)", y="Density") +
   	guides(ncol=2) +
   	theme(text=element_text(family="Times"), 
  	 	axis.title = element_text(face="bold", size=15),
  	 	axis.text.x = element_text(angle=0),
  	 	axis.text = element_text(size=15),
  	 	legend.text = element_text(size=15), 
  	 	legend.title = element_blank(), 
  	 	legend.position = "bottom",
  	 	strip.text = element_text(size = 12, face = "bold"), 
  	 	strip.background=element_rect(color="white")
  	 	) 
  ggsave(pp, file="density_free.svg", width=8, height=14)
