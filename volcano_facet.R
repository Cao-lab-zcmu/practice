library(ggplot2)
library(ggrepel)

setwd("volcano")

data <- read.csv(file="volcano_facet",header=T,sep="\t")

#data <- na.omit(data)

savename="volcano_facet"

title="Volcano grid"

data$change <- factor(ifelse(data$p.value < 0.05 & abs(data$fold) >= 1,
                      ifelse(data$fold >= 1,"up","down"),"stable"),
                      levels = c("up","down","stable"))

p <- ggplot(data,aes(x=fold, y=-log10(p.value),color = change)) + 

  geom_point(alpha=0.8, stroke=0, size=1.5) + 
  
  scale_color_manual(values = c("down"="#4DBBD5FF","stable"="#8491B4FF","up"="#DC0000FF")) +
  
  #xlim(-10,10) + 
  
  ylim(1,max(-log10(data$p.value))) +
   
  geom_hline(yintercept = -log10(0.05),linetype=4,size=0.8) +
  
  geom_vline(xintercept = c(-1,1),linetype=4,size=0.8) + 
  
  labs(x = "log2(foldchange)",y="-log10(p-value)",title=title) + 
  
  geom_text_repel(data = data[data$p.value<0.01 & abs(data$fold) >= 2,],
                  aes(label = id),size = 3,family="Times") +
                  
  #annotate("text", x = -10, y = max(-log10(data$p.value))*0.8, label = paste0("Based on ",nrow(data)," features"),
  #          color="black",size = 3, fontface="bold", family="Times", hjust = 0 ) +
  
  theme(text=element_text(family="serif"),
    
    #axis.line = element_line(colour = "black", size=0.2),

    #plot.margin = unit(c(3, 1, 3, 1), "cm")
    
    plot.title = element_text(hjust = 0.5)) +
    
  facet_grid(g2_deep ~ g1_from)
    
 #svg(paste0(savename,".svg"),width=8,height=6.5)
 #p
 #dev.off()
 
 ggsave(p,file=paste0(savename,".svg"),width=8,height=6.5)
