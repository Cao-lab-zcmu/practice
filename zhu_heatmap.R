###  ggplot2_heatmap.R

library(ggplot2)
library(ggtree)
library(aplot)
library(tidyr)
library(rayshader)
library(ggsci)
library(stringr)

file="heatmap_1.tsv"

savename="zhu_heatmap.svg"

data <- read.csv(file=file,header=T,sep="\t",row.names=1, check.name=F)

colnames(data) <- str_wrap(colnames(data), width=30)

phr <- hclust(dist(data)) %>% ggtree(layout="rectangular", branch.length="none")

phc <- hclust(dist(t(data))) %>% ggtree(layout="rectangular", branch.length="none") + layout_dendrogram()

data$names <- rownames(data)

p1 <- gather(data, 1:(ncol(data)-1), key="condition", value='expr')  ## keep the last col

pp <- ggplot(p1,aes(x=condition,y=names,fill=expr)) + 
  geom_tile() +
  #theme_minimal()+
  scale_fill_gradient2() +
  #scale_fill_npg() +
  scale_y_discrete(position="right") +
  labs(x="Compounds", y="Sample", fill="Log10(Peak area)") +
  theme(
   	axis.text.x = element_text(angle=90, hjust=1),
  	axis.text = element_text(size=20),
  	axis.title = element_text(size=20, face="bold"),
  	#axis.text = element_blank(),
  	legend.title = element_text(size=20, face="bold"),
  	legend.text = element_text(size=20),
  	legend.key.width = unit(2, "cm"),
    	legend.key.height = unit(4, "cm"),
    	text=element_text(family="serif")
   	#axis.title.y = element_text(size = 14),
   	#plot.title = element_text(hjust = 1,vjust=-40,size=14)
   	)

pp_com <- pp %>% 
  insert_left(phr, width=.1) %>%
  insert_top(phc, height=.1) 

ggsave(pp_com,file=savename, width=20, height=20)

#plot_gg(pp_com, multicore = TRUE, width = 20 ,height=20, scale=250) #加载图形
#render_depth(focallength=100,focus=0.72)

library(ggupset)

file="instance_norm_random.tsv"

data <- read.csv(file=file,header=T,sep="\t", check.name=F)

data$upset_x=paste0(data$names, "_", data$condition)

p <- ggplot(data, aes(x=factor(upset_x), y=expr)) +
  geom_col() +
  #scale_x_mergelist(sep = "_") +
  axis_combmatrix(sep = "_") +
  #coord_flip() +
  labs(x="Feature link", y="Normalized ftalign similarity") +
  theme(
   	#axis.text.x = element_text(angle=90),
  	axis.text.y = element_text(size=20, angle=90),
  	axis.title.x = element_text(size=30, face="bold", angle=180),
  	axis.title.y = element_text(size=30, face="bold"),
  	#axis.text = element_blank(),
  	legend.title = element_text(size=20, face="bold"),
  	legend.text = element_text(size=20),
  	#legend.key.width = unit(2, "cm"),
    	#legend.key.height = unit(4, "cm"),
    	text=element_text(family="serif")
   	#axis.title.y = element_text(size = 14),
   	#plot.title = element_text(hjust = 1,vjust=-40,size=14)
   	) +
  theme_combmatrix(combmatrix.panel.point.size = 5, 
   	 	   #combmatrix.panel.margin = unit(c(0.5,0.5),"pt"),
   	 	   combmatrix.panel.line.size = 2,
   	 	   combmatrix.label.height = unit(500, "pt"), 
   	 	   combmatrix.label.text = element_text(family="serif", angle=145, size=12, face="bold", hjust=0, vjust=-0.2),
   	 	   combmatrix.label.make_space = F)

ggsave(p, file="test1.pdf", width=20, height=20)



