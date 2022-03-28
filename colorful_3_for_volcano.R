 	 	 	library(tidyverse)
 	 	 	
 	 	 	library(RColorBrewer)
 	 	 	
 	 	 	library(ggsci)
 	 	 	
 	 	 	data <- read.csv(file="pearson_filter_volcano_sun.tsv",header=T,sep="\t")
	 	 	
	 	 	data$group_comple<-ifelse(data$FC>0,"up","down")
	 	 	
	 	 	label_data <- data %>% group_by(colnames(data), group) %>% summarize(tot=sum(FC))
	 	 	
	 	 	label_data=label_data[order(label_data$group),]
	 	 	
	 	 	label_data$num=rownames(label_data)
	 	 	
	 	 	number_of_bar <- nrow(label_data)
	 	 	
	 	 	angle <- 90 - 360 * (as.numeric(label_data$num)-0.5) /number_of_bar
	 	 	
	 	 	label_data$hjust <- ifelse( angle < -90, 1, 0)
	 	 	
	 	 	label_data$angle <- ifelse(angle < -90, angle+180, angle)
	 	 	
	 	 	colourCount = length(unique(data$id))
	 	 	
 	 	 	getPalette = colorRampPalette(brewer.pal(12, "Paired"))
	 	 	
	 	 	test1=label_data[nrow(label_data),]
	 	 	
	 	 	p <- 
	 	 	 ggplot(data) +
	  	 	 geom_bar(aes(x=as.factor(name), y=fold_change, fill=as.factor(id)),
	  	 	          stat="identity", alpha=0.5
	  	 	          ) +
	  	 	 geom_hline(yintercept = 0, linetype = "twodash", size=0.5, alpha=0.5) +
	  	 	 geom_text(data=label_data, aes(x=name, y=max(data$fold_change), label=name, hjust=hjust),
	  	 	          color="black", alpha=0.6, fontface="bold", size=3, angle= label_data$angle,
	  	 	          inherit.aes = FALSE,
	  	 	          family="Times") +
	  	 	 geom_text(data=test1, aes(x=name, y=-40, label="p-value < 0.05"),
	  	 	          color="black", alpha=0.6, fontface="bold", size=3, hjust=0.5,
	  	 	          inherit.aes = FALSE,
	  	 	          family="Times") +
	  	 	 geom_text(data=test1, aes(x=name, y=-50, label="r(pearson) >= 0.4"),
	  	 	          color="black", alpha=0.6, fontface="bold", size=3, hjust=0.5,
	  	 	          inherit.aes = FALSE,
	  	 	          family="Times") +
	  	 	 scale_fill_manual(values = getPalette(colourCount)) +
	  	 	 theme_minimal() + 
	  	 	 ylim(-50,50) +
	  	 	 theme(
	  	 	       text=element_text(family="serif"),
	    	 	       axis.text = element_blank(),
	    	 	       axis.title.y = element_blank(),
	    	 	       panel.grid = element_blank(),
	    	 	       legend.position = c(0.5,0.1),
	    	 	       legend.title= element_blank(),
	    	 	       legend.key.width = unit(0.2, "cm"),
   	 	 	       legend.key.height = unit(0.2, "cm"),
   	 	 	       legend.text=element_text(size=6)
	  	 	      ) +
	  	 	 labs(x = "", y = "")+
	  	 	 coord_polar() +
	 	 	 guides(fill=guide_legend(nrow=11))
	  	 	ggsave(p,file="test.pdf")
	  	 	
	  	 	
	  	 	
	  	 	          
