 ### line
 
 library(ggplot2)
 
 library(ggrepel)
 
 library(ggsci)
 
 library(ggalt)
 
 path="results/EIC_rt_during"
 
 list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 
 for(file in list){
 
 savename <- strsplit(file, split=".tsv")
 
 data <- read.csv(file=paste0(path, "/", file),header=T,sep="\t")
 
 data <- data[which(data$group!=""),]
 
 data_anno_mz <- data[1,colnames(data) %in% c("mz")]
 
 data_anno_rt <- data[1,colnames(data) %in% c("center_rt")]
 
 tolerance = 0.005
 
 data_label <- data[which(data$label==1),]
 
 anno_x_range <- c(0, max(data$rt))
 
 anno_y_range <- c(0, max(data$intensity))
 
 delta <- max(data$rt)-min(data$rt)
 
 
 
 p <- ggplot(data,aes(x=rt, y=intensity, group=sample, colour=color)) +
 
  geom_line() +
  
  geom_point(size=0.5, stroke=0) +
  
  geom_label_repel(data=data_label, aes(x=rt, y=intensity, label=sample),
	  	 	          color="black", alpha=0.5, fontface="bold", size=2, angle= 0, xlim=anno_x_range, ylim=anno_y_range,
	  	 	          direction="both", segment.size = 0.2, segment.alpha = 0.3, force = 1,
	  	 	          nudge_x = runif(1, min=delta*(1/18), max=delta*(1/10)), nudge_y = max(data$intensity)*(1/10),
	  	 	          inherit.aes = FALSE, hjust = 0,
	  	 	          family="Times") +
	  	 	          
  scale_color_manual(values = c("control"="#4A6990FF","drug"="#95CC5EFF","model"="#374E55FF",
  
   	 	 	 	"Non feature"="#B8B8B8FF",
  
   	 	 	 	"pro_low"="#FDAE6BFF","pro_medium"="#FD8D3CFF","pro_high"="#E6550DFF",
   	 	 	 	
   	 	 	 	"raw_low"="#9ECAE1FF","raw_medium"="#6BAED6FF","high_raw"="#3182BDFF")) +
  #drug green
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

 ggsave(p,file=paste0(path, "/", savename,".svg"),width=8,height=6.5)
 
 print(paste0("Finish >>> ",file))}

