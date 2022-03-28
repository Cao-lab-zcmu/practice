 	 	 	library(tidyverse)
 	 	 	
 	 	 	library(RColorBrewer)
 	 	 	
 	 	 	data <- read.csv(file="Sun_0827-1.tsv",header=T,sep="\t")
	 	 	
	 	 	data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
	 	 	
	 	 	nObsType <- nlevels(as.factor(data$observation))
	 	 	
	 	 	data <- data %>% arrange(group, name)
	 	 	
	 	 	data$id <- rep( seq(1, nrow(data)/nObsType), each=nObsType)
	 	 	
	 	 	label_data <- data %>% group_by(id, name) %>% summarize(tot=sum(value))
	 	 	
	 	 	number_of_bar <- nrow(label_data)
	 	 	
	 	 	angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
	 	 	
	 	 	label_data$hjust <- ifelse( angle < -90, 1, 0)
	 	 	
	 	 	label_data$angle <- ifelse(angle < -90, angle+180, angle)
	 	 	
	 	 	colourCount = length(unique(data$observation))
	 	 	
 	 	 	getPalette = colorRampPalette(brewer.pal(12, "Paired"))
	 	 	
	 	 	p <- 
	 	 	 ggplot(data) +
	  	 	 geom_bar(aes(x=as.factor(id), y=value, fill=observation), 
	  	 	          stat="identity", alpha=0.5
	  	 	          ) +
	  	 	 scale_fill_manual(values = getPalette(colourCount)) +
	  	 	 ylim(-10,13) +
	  	 	 theme_minimal() +
	  	 	 theme(
	  	 	       text=element_text(family="serif"),
	    	 	       axis.text = element_blank(),
	    	 	       axis.title.y = element_blank(),
	    	 	       panel.grid = element_blank(),
	    	 	       legend.position = c(0.5,0.03),
	    	 	       legend.title= element_blank(),
	    	 	       legend.key.width = unit(0.2, "cm"),
   	 	 	       legend.key.height = unit(0.2, "cm"),
   	 	 	       legend.text=element_text(size=6)
	  	 	      ) +
	  	 	labs(x = "", y = "")+
	  	 	guides(fill=guide_legend(nrow=10))+
	  	 	coord_polar() +
	  	 	geom_text(data=label_data, aes(x=id, y=0, label=name, hjust=hjust),
	  	 	          color="black", alpha=1, fontface="bold", size=2, angle= label_data$angle, 
	  	 	          inherit.aes = FALSE,
	  	 	          family="Times")
	  	 	          
	  	 	pdf("test.pdf")
	  	 	p
	  	 	dev.off()
	  	 	
	  	 	
	  	 	
	  	 	
