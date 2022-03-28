library(pheatmap)
map<- read.table(file= "hot_image-1.tsv",header= T,row.names= 1,sep= "\t")
test=pheatmap(map,
	fontsize=5,cellheight=10,cellwidth=12,border_color="#0a0404",
	treeheight_row=15,treeheight_col=15,
	clustering_method = "complete",
	color=colorRampPalette(c("#e60612","#e60612","white","#4a85c5","#4a85c5"))(100))
pdf("re_heatmap.pdf")
test
dev.off()

################echo /media/wizard/Seagate/MDMN/sep_neg0628/*_*_*/ | xargs -n 1 cp -v compound.config
################
################
#heatmap
## clustering_method参数设定不同聚类方法，默认为"complete",可以设定为'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
## clustering_distance_rows = "correlation"参数设定行聚类距离方法为Pearson corralation，默认为欧氏距离"euclidean"
#clustering_distance_rows = "correlation")
#clustering_method = "ward")
library(pheatmap)
map<- read.table(file= "hot_image.tsv",header= T,row.names= 1,sep= "\t")
group= read.table(file= "group.tsv",header= T,sep="\t")
test=pheatmap(map,annotation_row=group[1],
	fontsize=5,cellheight=10,cellwidth=12,border_color="#0a0404",
	treeheight_row=15,treeheight_col=15,
	clustering_method = "complete",
	color=colorRampPalette(c("#e60612","#e60612","white","#4a85c5","#4a85c5"))(100))
pdf("heatmap.pdf")
test
dev.off()

#################
#################
#################
# pca ggbiplot



################################

library(ggbiplot)
library(ggsci)
library(scales)
library(ggrepel)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggforce)
library(ggpmisc)

datas=list.files(path = ".", pattern = "*pca.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
savename=strsplit(datas[1],split=".tsv")
source<- read.table(file= datas[1],header= T,row.names= 1,sep= "\t")
class<- read.table(file= "metadata_level.tsv",header= T,row.names= 1,sep= "\t")
pca<- prcomp(source, scale. = TRUE)
pca_anno <- as.data.frame(pca[5])
pca<-cbind(pca_anno,class)

p <- ggplot(pca, aes(x=x.PC1, y=x.PC2, fill=level)) +
 	geom_point() +
 	#stat_ellipse(aes(color=class$subgroup), level = 0.95) +
 	scale_color_npg() +
 	scale_fill_npg() +
 	guides(fill= guide_colourbar(order = 1)) +
 	
 	theme(legend.position = "right",text=element_text(family="serif"))

ggsave(p,file=paste0(savename,"_level.svg"))

#################
#################
##pca  ------------ropls
library(ropls)
library(ggbiplot)
library(ggsci)
library(scales)
source<- read.table(file= "fecal_pos_YH_oplsda.tsv",header= T,row.names= 1,sep= "\t")
class<- read.table(file= "metadata_YH_oplsda.tsv",header= T,row.names= 1,sep= "\t")
pca <- opls(x = source)
df <- pca@scoreMN
df <- as.data.frame(df)
p <- ggplot(df, aes(x=p1, y=p2, fill=class$subgroup)) +
 	geom_point(alpha=0.8, size=3, shape=21, stroke=0.1) +
 	stat_ellipse(aes(color=class$subgroup), level = 0.95) +
 	scale_color_npg() +
 	scale_fill_npg() +
 	theme(legend.position = "right",text=element_text(family="serif"))

svg("raw_and_pro_YH_pca_Nversion.svg", height=6, width=10)
p
dev.off()	

#################
#################
##pls-da

library(ropls)
library(ggbiplot)
library(ggsci)
library(scales)
data <- read.table(file= "fecal_pos_pca.tsv",header= T,row.names= 1,sep= "\t",check.names=F)
class <- read.table(file= "metadata.tsv",header= T,row.names= 1,sep= "\t",check.names=F)
plsda <- opls(x = data, y = class[, "subgroup"], orthoI = 0)
df <- plsda@scoreMN
df <- as.data.frame(df)
p <- ggplot(df, aes(x=p1, y=p2, fill=class$subgroup)) +
 	geom_point(alpha=0.8, size=3, shape=21, stroke=0.1) +
 	stat_ellipse(aes(color=class$subgroup), level = 0.95) +
 	scale_color_npg() +
 	scale_fill_npg() +
 	theme(legend.position = "right",text=element_text(family="serif"))

svg("pls_da_ggplot.svg", height=6, width=10)
p
dev.off()	

#################
#################
## opls-da
##


#################
#################
#################




#################
#################
#################
#boxplot
library(ggplot2)
source<-read.table(file="chang",header=T,sep="\t")
boxplot<-ggplot(source,aes(x=days,y=length,fill=days))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(width=0.7)+
  geom_dotplot(binaxis ="y",
               stackdir="center",
               #position="jitter",
               dotsize = 0.4,)+
  stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="grey")
pdf("output_boxplot.pdf")
boxplot
dev.off()

#################
#################
#################
#bar
library(ggplot2)
source<-read.table(file="bar.tsv",header=T,sep="\t")
mzscale=1200
colwidth=mzscale/200
barplot<-ggplot(source, aes(x=mz, y=abun)) + 
geom_bar(stat="identity",width=colwidth,fill=ifelse(source$abun>0,'black','red')) + 
geom_point(size=0.95,color=ifelse(source$abun>0,'black','red')) + 
xlim(0,mzscale) +
theme(panel.background=element_rect(fill='transparent',color='white'),
 	panel.grid=element_line(color='grey85'), plot.margin = unit(c(3, 1, 3, 1), "cm"))
pdf("output_bar.pdf")
barplot
dev.off()

################
################
################
行列布局
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 3)))
print(p.scatter, vp=viewport(layout.pos.row=2:3, layout.pos.col=1:2))
print(p.hist.len, vp=viewport(layout.pos.row=1, layout.pos.col=1:2))
print(p.hist.wid, vp=viewport(layout.pos.row=2:3, layout.pos.col=3))

################
################
################

library(grid)
library(ggplot2)
# prepare ggplot charts
p.hist.len <- ggplot(iris) + geom_histogram(aes(x=Sepal.Length))
p.hist.wid <- ggplot(iris) + geom_histogram(aes(x=Sepal.Width)) + coord_flip()
p.scatter <- ggplot(iris) + geom_point(aes(x=Sepal.Length, y=Sepal.Width))
 
 
 
 

# create viewports
library(grid)
grid.newpage()
vp1 <- viewport(x=0, y=0, width=5, height=5)
vp2 <- viewport(x=0, y=0, width=2, height=2)
 
# direct the charts into the specified viewport
s1 <- print(p, vp=vp1)
s2 <- print(l, vp=vp2)

################
################

library("grImport2")
readPicture("file")
grid.picture(var)

################
################
################
library(visNetwork)
lists = read.table(file="net.csv",header=T,sep="\t")
nodes = read.table()
edges = read.table()

network = visNetwork(nodes, edges, height = "100%", width = "100%")
htmlwidgets::saveWidget(network, "network.html")

################
################
################
#Radial bar plot

library(tidyverse)
data <- read.table(file=args[1],header=T,sep="\t")
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)),
  value1=sample( seq(10,100), 60, replace=T),
  value2=sample( seq(10,100), 60, replace=T),
  value3=sample( seq(10,100), 60, replace=T)
)
data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
nObsType <- nlevels(as.factor(data$observation))
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType), each=nObsType)
################
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
p <- ggplot(data) +
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  scale_fill_brewer(palette = "Paired") +
  ylim(-150,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE )

png('tr_tst2.png',width=300,height=300,units="px",bg="transparent")
p
dev.off()

################
nodes(imagepath,shape"image")
nodes <- data.frame(id = 1:4, 
                    shape = c("image", "circularImage"),
                    label = "I'm an image")



visNetwork(nodes, edges) %>%
  visPhysics(stabilization = FALSE) %>%
  visNodes(shapeProperties = list(useBorderWithImage = FALSE)) %>%
  visEdges(smooth = FALSE) %>%
  #visIgraphLayout() %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group", manipulation = TRUE) %>%
  visLayout(randomSeed = 123)

########################
#bubble plot
#######################

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggforce)
library(ggpmisc)
data <- read.csv(file="fecal_pos_mzmine.tsv",header=T,sep="\t")
p <- ggplot(data, aes(x=rt, y=xlogp, size=similarity, fill=m.z, color=str_wrap(classification,30))) +
  geom_smooth(method = "lm", se=TRUE, color = "black", fill = "skyblue", alpha=0.3, size=0.4, formula = y ~ x) +
  geom_point(alpha=0.4, shape=21, stroke=0.1) +
  geom_mark_ellipse(aes(color=str_wrap(classification,30)), alpha=0.1, size=0.2, expand=0.02, fill=NA) +
  labs(size="Similarity") +
  labs(colour="Classification") +
  coord_cartesian(clip = "off") +
  scale_size_continuous( trans="exp", range=c(2, 5)) +
  theme_minimal() +
  theme(text=element_text(family="serif"))
pdf("bubble.pdf",width=7,height=5)
p
dev.off()

####################

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggforce)
library(ggpmisc)
data <- read.csv(file="0.5_com_compound.tsv",header=T,sep="\t")
p <- ggplot(data, aes(x=rt, y=xlogp, size=similarity, fill=m.z, color=classification)) +
  geom_smooth(method = "lm", se=TRUE, color = "black", fill = "skyblue", alpha=0.3, size=0.4, formula = y ~ x) +
  geom_point(alpha=0.4, shape=21, stroke=0.1) +
  scale_size_continuous(trans="exp", range=c(2, 5)) +
   guides(
   colour = "none", 
   fill = guide_colourbar(order = 1),
   size = guide_legend(order = 2)
 ) +
  theme_minimal() +
  theme(text=element_text(family="serif"))
pdf("sep_bubble.pdf",width=7,height=5)
p
dev.off()

###################

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggforce)
library(ggpmisc)
data <- read.csv(file="0.5_com_compound.tsv",header=T,sep="\t")
data_sep <- data[which(data$m.z<=800 & data$m.z>=400),]
p <- ggplot(data_sep, aes(x=rt, y=xlogp, size=similarity, fill=m.z, color=classification)) +
  geom_smooth(method = "lm", se=TRUE, color = "black", fill = "skyblue", alpha=0.3, size=0.4, formula = y ~ x) +
  
  stat_cor(method = "pearson",label.x = 3, label.y = 30) +
  
  stat_poly_eq(
    aes(label = ..eq.label..),
    formula = formula,parse = TRUE, geom = "text",label.x = 3,label.y = 28, hjust = 0) +
     
  geom_point(alpha=0.4, shape=21, stroke=0.1) +
  scale_size_continuous(trans="exp", range=c(2, 5)) +
   guides(
   colour = "none", 
   fill = guide_colourbar(order = 1),
   size = guide_legend(order = 2)
 ) +
  theme_minimal() +
  theme(text=element_text(family="serif"))
pdf("sep_bubble_400_800.pdf",width=7,height=5)
p
dev.off()

###################
#rank
library(tidyverse)
data <- read.csv(file="rank.tsv",header=T,sep="\t")
p <- ggplot(data) +
  geom_point(aes(x=rank, y=log10_raw, size=similarity),fill="#99CCFF",alpha=0.4, shape=21, color="black",stroke=0.001) +
  geom_point(aes(x=rank, y=log10_pro, size=similarity),fill="#FF6666",alpha=0.4, shape=21, color="black",stroke=0.001) +
  scale_size_continuous(trans="exp", range=c(1, 3)) +
   guides(
   colour = "none",
   size = guide_legend(order = 2)
 ) +
  theme_minimal() +
  labs(y = "log10(area)")
  theme(text=element_text(family="serif"))
pdf("rank.pdf",width=7,height=5)
p
dev.off()

##################

library(ggplot2)

library(ggforce)

library(ggrepel)

n=list.files(path = ".", pattern = "*.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)

for (x in n) {

file <- strsplit(x, split=".tsv")

savename=paste0(file,".svg")

source<-read.csv(file=paste0(file,".tsv"),header=T,sep="\t")

source_anno=source[which(source$rel.intensity<0),]

barplot<-ggplot(source, aes(x=mz, y=rel.intensity)) + 

  geom_bar(stat="identity",width=max(source$mz)/150,fill=ifelse(source$rel.intensity>0,"black","red")) + 

  geom_point(size=1.3,color=ifelse(source$rel.intensity>0,"black","red"),alpha=ifelse(source$match==0,0,1)) + 
  
  xlim(0,max(source$mz, na.rm=T)+50) +
  
  annotate("text", x = 10, y = 90, label = paste0("Precursor m/z: ",source[1,colnames(source) %in% c("precursor.m.z")]),
            color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
            
  annotate("text", x = 10, y = 78, label = paste0("RT (min): ", source[1,colnames(source) %in% c("RT..min.")]),
            color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
            
  annotate("text", x = 10, y = 66, label = paste0("Tanimoto similarity: ", source[1,colnames(source) %in% c("similarity")]),
            color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
  
  geom_label_repel(data=source_anno, aes(x=mz, y=rel.intensity, label=round(mz,2)),
	  	 	          color="black", alpha=0.5, fontface="bold", size=2, angle= 0, 
	  	 	          direction="both", ylim=c(0, -100), segment.size = 0.2, segment.alpha = 0.3,
	  	 	          inherit.aes = FALSE, hjust = 0,
	  	 	          family="Times") +
  
  # annotate("curve", x = source_anno$mz+1, xend = source_anno$mz+4, y = source_anno$rel.intensity-3, yend = source_anno$y_label,
  #          colour = "black", curvature =0, size=0.1, alpha=0.5, arrow = arrow(length = unit(0.5, "mm"))) +
  
  theme(text=element_text(family="serif"),

  panel.background=element_rect(fill="transparent",color="white"),

  panel.grid=element_line(color="grey85"), 

  plot.margin = unit(c(3, 1, 3, 1), "cm"))
  
ggsave(barplot,file=savename,width=8,height=6.5)

print(file)}

##################


 
 
   
 #########################
 ### line
 
 library(ggplot2)
 
 library(ggrepel)
 
 library(ggsci)
 
 list <- list.files(path = ".", pattern = "*.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 
 for(file in list){
 
 savename <- strsplit(file, split=".tsv")
 
 data <- read.csv(file=file,header=T,sep="\t")
 
 data <- data[which(data$group!=""),]
 
 data_anno_mz <- data[1,colnames(data) %in% c("mz")]
 
 data_anno_rt <- data[1,colnames(data) %in% c("center_rt")]
 
 tolerance = 0.005
 
 data_label <- data[which(data$label==1),]
 
 anno_x_range <- c(0, max(data$rt))
 
 anno_y_range <- c(0, max(data$intensity))
 
 delta <- max(data$rt)-min(data$rt)
 
 p <- ggplot(data,aes(x=rt, y=intensity, group=sample, colour=group)) +
 
  geom_line(alpha=0.8) +
  
  geom_label_repel(data=data_label, aes(x=rt, y=intensity, label=sample),
	  	 	          color="black", alpha=0.5, fontface="bold", size=2, angle= 0, xlim=anno_x_range, ylim=anno_y_range,
	  	 	          direction="both", segment.size = 0.2, segment.alpha = 0.3, force = 1,
	  	 	          nudge_x = runif(1, min=delta*(1/18), max=delta*(1/10)), nudge_y = max(data$intensity)*(1/10),
	  	 	          inherit.aes = FALSE, hjust = 0,
	  	 	          family="Times") +
	  	 	          
  scale_color_npg() +
  
  scale_fill_npg() +
  
  labs(color="Peak attribution", x="RT (min)", y="Intensity") +
  
  annotate("text", x = min(data$rt), y = max(data$intensity)*(16/20), 
   	    label = paste0("Precursor m/z: ",data_anno_mz-tolerance," ~ ",data_anno_mz+tolerance),
            color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
  
  annotate("text", x = min(data$rt), y = max(data$intensity)*(15/20), 
   	    label = paste0("RT (min): ",data_anno_rt),
            color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
  
  #theme_minimal() +
  
  theme(text=element_text(family="serif"),
  
    legend.position = c(0.85,0.75),
    
    #axis.line = element_line(colour = "black", size=0.2),
    
    legend.background = element_rect(fill = "transparent", color = "transparent"),

    plot.margin = unit(c(3, 1, 3, 1), "cm"))

 ggsave(p,file=paste0(savename,".svg"),width=8,height=6.5)
 
 print(paste0("Finish >>> ",file))}

 #######################################
 # python Draw smiles
 
 ######################
 ######################
 #change svg to pdf
 
 library(rsvg)
 
 files=list.files(path = ".", pattern = "*.svg$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
           
 for(file in files){
 
   name_set=strsplit(file, split=".svg")
   
   id=name_set[1]
   
   savename=paste0(id,".pdf")
   
   rsvg_pdf(file, savename)
   }

 ########################
 
 library(staplr)
 
 ms1=list.files(path = "EIC_rt_during", pattern = "*.pdf$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 
 for(file in ms1){
 
 id=strsplit(file, split=".pdf")
 
 ms1=paste0("EIC_rt_during/",id,".pdf")
 
 ms2=paste0("ms2_figures_label/",id,".pdf")
 
 structure=paste0("structure_2d/",id,".pdf")
 
 gather=c(ms1, ms2, structure)
 
 savename=paste0("report_pdf/",file)
 
 staple_pdf(input_files=gather, output_filepath= savename)
 }


########################
###t.test and volcano


########################
### Pearson correlation coefficient



#########################
### linear regression



#########################
### volcano plot



######## network

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)



















