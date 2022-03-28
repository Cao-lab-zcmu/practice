 ### line
 library(xcms)
 path="/media/wizard/back/thermo_mzML_0518"
 data <- read.csv(file="com_lignans_and_iridoids.tsv",header=T,sep="\t")
 data <- data[!duplicated(data$id),]
 data <- data[which(data$similarity >= 0.5), ]
 # data <- data[which(data$id==1746|data$id==2081),]
 metadata <- data[which(log2(data$pro.raw)>1 | log2(data$pro.raw)<(-1)), colnames(data) %in% c("id", "m.z")]
 # metadata <- data[, colnames(data) %in% c("id", "m.z")]
 dda_file=list.files(path = path, pattern = "*.mzML$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 dir.create("/media/wizard/back/thermo_mzML_0518/EIC")
 tolerance=0.005
 for(filename in dda_file){
   dda_data <- readMSData(paste0(path,"/",filename), mode = "onDisk")
   dir.create(paste0(path,"/EIC/EIC_",filename))
   for(number in 1:nrow(metadata)){
     id <- metadata[number,colnames(metadata) %in% c("id")]
     mz <- metadata[number,colnames(metadata) %in% c("m.z")]
     mzrange <- c(as.numeric(mz)-tolerance,as.numeric(mz)+tolerance)
     ex_data <- chromatogram(dda_data, msLevel = 1L, mz = mzrange, aggregationFun = "max")
     ex_data_1 <- ex_data[1,1]
     if(number==1){
       write.table(rtime(ex_data_1),paste0(path, "/EIC/EIC_",filename,"/rt",".tsv"),col.names = FALSE,sep="\t")}
     write.table(intensity(ex_data_1),paste0(path, "/EIC/EIC_",filename,"/",id,"_intensity",".tsv"),col.names = FALSE,sep="\t")
     print(paste0(filename," >>> ",number,"/",nrow(metadata)))
   }
 }
 write.table(metadata, paste0(path, "/metadata.tsv"), col.names = T, row.names=F, sep="\t")
 ### run bash script colorful_line_eucommia_bash.sh
 library(ggplot2)
 library(ggrepel)
 library(ggsci)
 library(ggalt)
 path="results/EIC_rt_during"
 list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                    full.names = FALSE, recursive = FALSE,
                    ignore.case = FALSE, include.dirs = FALSE)
 data <- read.csv(file="com_lignans_and_iridoids.tsv",header=T,sep="\t")
 data <- data[which(data$similarity >= 0.5), ]
 draw_list <- data[which(log2(data$pro.raw)>1 | log2(data$pro.raw)<(-1)), colnames(data) %in% c("id")]
 draw_list <- paste0(as.character(draw_list), ".tsv")
 list <- list[list %in% draw_list]
 for(file in list){
   savename <- strsplit(file, split=".tsv")
   data <- read.csv(file=paste0(path, "/", file),header=T,sep="\t")
   data <- data[which(data$group!=""),]
   data_anno_mz <- data[1,colnames(data) %in% c("mz")]
   data_anno_rt <- data[1,colnames(data) %in% c("center_rt")]
   tolerance = 0.005
   data_label <- data[which(data$label==1),]
   data_label$sample <- strsplit(data_label$sample, split=".mzML")
   anno_x_range <- c(0, max(data$rt))
   anno_y_range <- c(0, max(data$intensity))
   delta <- max(data$rt)-min(data$rt)
   if(file=="3938.tsv"){data <- data[which(data$rt<=12.08),] }
   p <- ggplot(data,aes(x=rt, y=intensity, group=sample, colour=color)) +
     geom_line() +
     geom_point(size=0.5, stroke=0) +
     #  geom_label_repel(data=data_label, aes(x=rt, y=intensity, label=sample),
     #	  	 	          color="black", alpha=0.5, fontface="bold", size=2, angle= 0, xlim=anno_x_range, ylim=anno_y_range,
     #	  	 	          direction="both", segment.size = 0.2, segment.alpha = 0.3, force = 1,
     #	  	 	          #nudge_x = runif(1, min=delta*(1/18), max=delta*(1/10))*sample(c(1,-1),1), 
     #	  	 	          nudge_y = max(data$intensity)*(1/10),
     #	  	 	          inherit.aes = FALSE, hjust = 0,
     #	  	 	          family="Times") +
     scale_color_manual(values = c("Blank"="#4A6990FF",
                                   "Non feature"="#B8B8B8FF",
                                   "Pro_Eucommia"="#E6550DFF",
                                   "Raw_Eucommia"="#3182BDFF")) +
#drug green
labs(color="Peak attribution", x="RT (min)", y="Intensity") +
annotate("text", x = min(data$rt), y = max(data$intensity)*(17/20), 
         label = paste0("ID: ",strsplit(file, split=".tsv")),
         color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
annotate("text", x = min(data$rt), y = max(data$intensity)*(16/20), 
         label = paste0("Precursor m/z: ",data_anno_mz-tolerance," ~ ",data_anno_mz+tolerance),
         color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
annotate("text", x = min(data$rt), y = max(data$intensity)*(15/20), 
         label = paste0("RT (min): ",data_anno_rt),
         color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
#theme_minimal() +
theme(text=element_text(family="Times"),
      legend.position = c(0.85,0.75),
      #axis.line = element_line(colour = "black", size=0.2),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      plot.margin = unit(c(3, 1, 3, 1), "cm"))
ggsave(p,file=paste0(path, "/", savename,".svg"),width=8,height=6.5)
print(paste0("Finish >>> ",file))
 }
 #########################
 #########################
 #########################
 library(ggplot2)
 library(ggrepel)
 library(ggsci)
 library(ggalt)
 library(hrbrthemes)
 path="results/EIC_rt_during"
 list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 data <- read.csv(file="com_lignans_and_iridoids.tsv",header=T,sep="\t")
 data <- data[which(data$similarity >= 0.5), ]
 draw_list <- data[which(log2(data$pro.raw)>1 | log2(data$pro.raw)<(-1)), colnames(data) %in% c("id")]
 draw_list <- paste0(as.character(draw_list), ".tsv")
 list <- list[list %in% draw_list]
 n=0
 for(file in list){
 id <- strsplit(file, split=".tsv")
 n=n+1
 savename <- strsplit(file, split=".tsv")
 data <- read.csv(file=paste0(path, "/", file),header=T,sep="\t")
 data <- data[which(data$group!=""),]
 data_anno_mz <- data[1,colnames(data) %in% c("mz")]
 data_anno_rt <- data[1,colnames(data) %in% c("center_rt")]
 tolerance = 0.005
 data_label <- data[which(data$label==1),]
 data_label$sample <- strsplit(data_label$sample, split=".mzML")
 if(file=="3938.tsv"){data <- data[which(data$rt<=12.08),] }
 if(file=="3380.tsv"){data <- data[which(data$rt>=12.55),] }
 data$id <- id
 if(n==1){sum_data=data[,c(1:6,9)]}else{sum_data<-rbind(sum_data, data[,c(1:6,9)])}
 #text1=c("x"=min(data$rt), "y"=max(data$intensity)*(17/20), "label"=paste0("ID: ", id))
 #
 text2=c("x"=min(data$rt), "y"=max(data$intensity)*(17/20), 
  	 	"label"=paste0("Precursor m/z: ",data_anno_mz-tolerance," ~ ",data_anno_mz+tolerance))
 text3=c("x"=min(data$rt), "y"=max(data$intensity)*(15/20), "label"=paste0("RT (min): ",data_anno_rt) )
 text=data.frame(rbind(text2, text3))
 text$id <- id
 if(n==1){sum_text=text}else{sum_text<-rbind(sum_text, text)}
 print(paste0("Finish >>> ",file))}
 #######################################  plot
 #######################################
 select <- c(279, 458, 574, 1107, 1445, 2227, 2529, 2664, 2824, 3380, 3918, 3938)
 sum_data <- sum_data[sum_data$id %in% select,]
 sum_data$id <- paste0("ID:", sum_data$id)
 sum_text <- sum_text[sum_text$id %in% select,]
 sum_text$id <- paste0("ID:", sum_text$id)
 p <- ggplot(sum_data, aes(x=rt, y=intensity, group=sample, colour=color)) +
  geom_line() +
  geom_point(size=0.5, stroke=0) +
  scale_color_manual(values = c("Blank"="#4A6990FF",
   	 	 	 	"Non feature"="#B8B8B8FF",
   	 	 	 	"Pro_Eucommia"="#E6550DFF",
   	 	 	 	"Raw_Eucommia"="#3182BDFF")) +
  #drug green
  labs(color="Peak attribution", x="RT (min)", y="Intensity") +
  geom_text(data=sum_text, aes(x=as.numeric(x), y=as.numeric(y), label=label), 
   	 	hjust=0, fontface="bold", alpha=1, size=3.5, inherit.aes = FALSE, family="Times") +
  #theme_minimal() +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(~as.character(id), scales="free", ncol=3) +
  guides(color=guide_legend(nrow=1)) +
  theme_ipsum() +
  theme(text=element_text(family="Times"),
    plot.background = element_rect(fill = "white", size=0),
    legend.position = "bottom",
    axis.title.y = element_text(face="bold", size=20, hjust=0.5, family="Times"),
    axis.title.x = element_text(face="bold", size=20, hjust=0.5, family="Times"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size=15),
    legend.title = element_text(size=15, face="bold"),
    #axis.line = element_line(colour = "black", size=0.2),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size=15, face="bold", family="Times"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    #panel.spacing = unit(c(0,0,0,0),"cm")
    )
 ggsave(p,file=paste0("results", "/", "ms1_line.svg"), width=12, height=15)
 #######################################  plot
 #######################################
library(ggplot2)
library(ggforce)
library(ggrepel)
path="results/ms2_figures_label"
n=list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE)
select <- c(279, 458, 574, 1107, 1445, 2227, 2529, 2664, 2824, 3380, 3918, 3938)
# select <- c(1746, 2081)
n <- n[n %in% paste0(select, ".tsv")]
m=0
for(x in n){
  m=m+1
  file <- strsplit(x, split=".tsv")
  id <- as.character(file)
  savename=paste0(file,".svg")
  data <- read.csv(file=paste0(path, "/", file,".tsv"),header=T,sep="\t")
  data$id <- id
  if(m==1){
    sum_data=data[,c(1:3,10)]
  }else{
    sum_data<-rbind(sum_data, data[,c(1:3,10)])
  }
  text1=c("x"=max(data$mz), "y"=90, "label"=paste0("Precursor m/z: ",data[1,colnames(data) %in% c("precursor.m.z")]) )
  text2=c("x"=max(data$mz), "y"=72, 
          "label"=paste0("RT (min): ", data[1,colnames(data) %in% c("RT..min.")]) )
  text3=c("x"=max(data$mz), "y"=54, "label"=paste0("Tanimoto similarity: ", data[1,colnames(data) %in% c("similarity")]) )
  text=data.frame(rbind(text1, text2, text3))
  text$id <- id
  if(m==1){
    sum_text=text
  }else{
    sum_text<-rbind(sum_text, text)
  }
  print(file) 
}
 #####################################
 #####################################
 data=sum_data
 data$id <- paste0("ID:", data$id)
 sum_text$id <- paste0("ID:", sum_text$id)
 p <- ggplot(data) + 
  geom_segment(aes(x=mz, xend=mz, y=0, yend=rel.intensity), color=ifelse(data$rel.intensity>0,"black","#E6550DFF"), size=0.8) +
  geom_point(data=data[which(data$match>=1),], size=0.9, aes(x=mz, y=rel.intensity),
   	 	 	color=ifelse(data[which(data$match>=1),]$rel.intensity>0,"black","red")) + 
  geom_text(data=sum_text, aes(x=as.numeric(x)*3/5, y=as.numeric(y), label=label), 
   	 	hjust=0, fontface="bold", alpha=1, size=3.5, inherit.aes = FALSE, family="Times") +
  facet_wrap(~as.character(id), scales="free_x", ncol=3) +
  labs(x="m/z", y="Relative intensity") +
  theme(text=element_text(family="Times"),
  panel.background=element_rect(fill="white", size=0),
  panel.grid=element_line(color="grey85"), 
  axis.title.y = element_text(face="bold", size=20, hjust=0.5, family="Times"),
  axis.title.x = element_text(face="bold", size=20, hjust=0.5, family="Times"),
  strip.text = element_text(size=15, face="bold", family="Times", hjust=0),
  strip.background  = element_rect(fill = "white", size = 0), 
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )
ggsave(p,file=paste0("results", "/", "ms2_segment.svg"),width=12,height=15)
###################################
###################################
library(ggplot2)
library(ggtree)
library(aplot)
library(tidyr)
library(rayshader)
library(ggsci)
library(reshape2)
 datapath="results/re_neg_RT.tsv"
 select <- c(279, 458, 574, 1107, 1445, 2227, 2529, 2664, 2824, 3380, 3918, 3938)
 data <- read.csv(file=datapath, header=T,sep="\t")
 grouppath="results/stat_classification.tsv"
 group <- read.csv(file=grouppath, header=T,sep="\t")
 col <- grep("row.ID|area", colnames(data), ignore.case=T)
 data <- data[data$row.ID %in% select, c(col)]
 colnames(data) <- c("id", "blank", "Pro-Eu1", "Pro-Eu2", "Pro-Eu3", "Raw-Eu1", "Raw-Eu2", "Raw-Eu3")
 rownames(data) <- paste0("ID: ",data$id)
 group_anno <- merge(data, group[, colnames(group) %in% c("id", "specific", "specific_pp")], by="id", all.x=T, sort=T)
 group_anno$id <- paste0("ID: ",group_anno$id)
 data <- data[, c(-1,-2)]
 data <- log2(data)
 ############### heatmap
 savename="results/features_heatmap.svg"
 phr <- hclust(dist(data)) %>% 
  	ggtree(layout="rectangular", branch.length="none", size=1.2) +
  	theme(
    	plot.margin = unit(c(0, 0, 0, 0), "cm")
   	)
phc <- hclust(dist(t(data))) %>% 
  	ggtree(layout="rectangular", branch.length="none", size=1.2) + layout_dendrogram() +
  	theme(
    	plot.margin = unit(c(0, 0, 0, 0), "cm")
   	)
data$names <- rownames(data)
 p1 <- gather(data, 1:(ncol(data)-1), key="condition", value='expr')  ## keep the last col
 p1 <- merge(p1, group_anno[, colnames(group_anno) %in% c("id", "specific")], by.x="names", by.y="id", all.x=T, sort=T)
 anno <- ggplot(p1, aes(x="group", y=names)) +
  geom_tile(size=2, aes(fill=specific), color="black") +
  scale_fill_npg() +
  labs(x="", y="", fill="Classes") +
  #guides(fill="none") +
  theme_minimal() +
  theme(
   	axis.text.x = element_text(angle=90, face="bold"),
   	axis.text.y = element_blank(),
  	legend.title = element_text(face="bold"),
  	#legend.text = element_text(size=20),
  	legend.key.width = unit(0.5, "cm"),
    	legend.key.height = unit(1, "cm"),
    	text=element_text(family="Times", size=15, face="bold"),
    	plot.margin = unit(c(0, 0, 0, 0), "cm")
   	)
 pp <- ggplot(p1,aes(x=condition,y=names,fill=expr)) + 
  geom_tile(size=1.5, color="black") +
  #theme_minimal()+
  scale_fill_gradientn(colors=c("#3182BDFF", "white", "#E6550DFF")) +
  #scale_color_gradientn(colors=c( "black", "#8F7700FF", "black")) +
  #scale_fill_gradient2(low="#3182BDFF", mid="black", high="#E6550DFF", midpoint=21) +
  scale_y_discrete(position="right") +
  guides(color="none") +
  labs(x="Sample name", y="", fill="Log2(peak area)") +
  theme_minimal()+
  theme(
   	axis.text.x = element_text(angle=90),
   	#axis.text.y = element_blank(),
  	legend.title = element_text(face="bold"),
  	#legend.text = element_text(size=20),
  	legend.key.width = unit(0.5, "cm"),
    	legend.key.height = unit(1, "cm"),
    	text=element_text(family="Times", size=15, face="bold"),
    	plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")
   	)
 pp_com <- pp %>% 
  insert_right(anno, width=.15) %>%
  insert_left(phr, width=.1) %>%
  insert_top(phc, height=.1) 
 ggsave(pp_com,file=savename, width=10, height=12)
 ###############

# extended data
 library(ggplot2)
 library(ggrepel)
 library(ggsci)
 library(ggalt)
 library(hrbrthemes)
 path="EIC_rt_during"
 list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 data <- read.csv(file="com_lignans_and_iridoids.tsv",header=T,sep="\t")
 data <- data[which(data$similarity >= 0.5), ]
 n=0
 for(file in list){
 id <- strsplit(file, split=".tsv")
 n=n+1
 savename <- strsplit(file, split=".tsv")
 data <- read.csv(file=paste0(path, "/", file),header=T,sep="\t")
 data <- data[which(data$group!=""),]
 data_anno_mz <- data[1,colnames(data) %in% c("mz")]
 data_anno_rt <- data[1,colnames(data) %in% c("center_rt")]
 tolerance = 0.005
 data_label <- data[which(data$label==1),]
 data_label$sample <- strsplit(data_label$sample, split=".mzML")
 if(file=="3938.tsv"){data <- data[which(data$rt<=12.08),] }
 if(file=="3380.tsv"){data <- data[which(data$rt>=12.55),] }
 data$id <- id
 if(n==1){sum_data=data[,c(1:6,9)]}else{sum_data<-rbind(sum_data, data[,c(1:6,9)])}
 #text1=c("x"=min(data$rt), "y"=max(data$intensity)*(17/20), "label"=paste0("ID: ", id))
 #
 text2=c("x"=min(data$rt), "y"=max(data$intensity)*(17/20), 
  	 	"label"=paste0("Precursor m/z: ",data_anno_mz-tolerance," ~ ",data_anno_mz+tolerance))
 text3=c("x"=min(data$rt), "y"=max(data$intensity)*(15/20), "label"=paste0("RT (min): ",data_anno_rt) )
 text=data.frame(rbind(text2, text3))
 text$id <- id
 if(n==1){sum_text=text}else{sum_text<-rbind(sum_text, text)}
 print(paste0("Finish >>> ",file))}
 #######################################  plot
 #######################################
 select <- c(1746,2081)
 sum_data <- sum_data[sum_data$id %in% select,]
 sum_data$id <- paste0("ID:", sum_data$id)
 sum_text <- sum_text[sum_text$id %in% select,]
 sum_text$id <- paste0("ID:", sum_text$id)
 p <- ggplot(sum_data, aes(x=rt, y=intensity, group=sample, colour=color)) +
  geom_line() +
  geom_point(size=0.5, stroke=0) +
  scale_color_manual(values = c("Blank"="#4A6990FF",
   	 	 	 	"Non feature"="#B8B8B8FF",
   	 	 	 	"Pro_Eucommia"="#E6550DFF",
   	 	 	 	"Raw_Eucommia"="#3182BDFF")) +
  #drug green
  labs(color="Peak attribution", x="RT (min)", y="Intensity") +
  geom_text(data=sum_text, aes(x=as.numeric(x), y=as.numeric(y), label=label), 
   	 	hjust=0, fontface="bold", alpha=1, size=5, inherit.aes = FALSE, family="Times") +
  #theme_minimal() +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(~as.character(id), scales="free", ncol=3) +
  guides(color=guide_legend(nrow=1)) +
  theme_ipsum() +
  theme(text=element_text(family="Times"),
    plot.background = element_rect(fill = "white", size=0),
    legend.position = "bottom",
    axis.title.y = element_text(face="bold", size=10, hjust=0.5, family="Times"),
    axis.title.x = element_text(face="bold", size=10, hjust=0.5, family="Times"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size=15),
    legend.title = element_text(size=15, face="bold"),
    #axis.line = element_line(colour = "black", size=0.2),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size=15, face="bold", family="Times"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    #panel.spacing = unit(c(0,0,0,0),"cm")
    )
 ggsave(p,file=paste0("extended_ms1_line.svg"), width=12, height=5)

