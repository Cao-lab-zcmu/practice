library(tidyverse)
library(ggsci)
library(scales)
library(ggrepel)
pal= pal_d3("category20c")(20)
pal= data.frame(matrix(pal,5:4))
pal= c(pal[1,2:3],pal[2,3:1])
data <- read.csv(file="com_lignans_and_iridoids.tsv",header=T,sep="\t")
data <- data[!duplicated(data$id),]
data <- data[which(data$similarity >= 0.5), ]
data$id <- factor(data$id, levels=data[order(data$pro.raw, decreasing = T), colnames(data) %in% c("id")])
label_data <- data
label_data <- label_data[order(label_data$id), ]
label_data$sequ <- seq(nrow(label_data))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$sequ-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
p <- ggplot(data, aes(x=as.factor(id), y=log2(pro.raw), fill=log2(pro.raw))) + 
  geom_bar(stat="identity", alpha=1) +
  ylim(-7,7) +
  coord_polar(start = 0) +
  labs(fill="Log2(FC)") +
  scale_fill_gradientn(colours = pal, breaks=c(-4,-2,0,2,4)) +
  geom_text(data=label_data, aes(x=id, y=ifelse(log2(pro.raw)>0, log2(pro.raw)+0.3, log2(pro.raw)-0.3), angle=angle, 
                                 label=ifelse(log2(pro.raw)>1 | log2(pro.raw)<(-1), as.character(id), "")), hjust=0,
            size=ifelse(log2(label_data$pro.raw)>0, 1.5, 1), 
            color="black", fontface="bold",alpha=0.6, inherit.aes = FALSE, family="Times" ) +
annotate("text", x=data[which(data$pro.raw==max(data$pro.raw)), colnames(data) %in% c("id")], y=-5,
         label = "Lignans and iridoids\n(PPCP > 0.5;\nTanimoto similarity > 0.5)\nFC(peak area:\nEucommia-pro/raw)", 
         family="Times", fontface="bold") +
#ggtitle("Lignans and iridoids(PP >= 0.9)\nFC(peak area:\nEucommia-pro/raw)") +
theme_minimal() +
theme(
      #plot.title = element_text(hjust = 0.5, vjust=-190, face="bold", size=12),
      text=element_text(family="Times"),
      axis.ticks = element_blank(),
      #panel.background = element_rect(fill="white"), 
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = c(0.5,0.41),
      legend.title=element_text(face="bold", hjust= -0.5),
      #plot.margin = unit(rep(-1,4), "cm"), 
      plot.margin = unit(c(-2, -2, -4, -2), "cm")    # Adjust the margin to make in sort labels are not truncated!
) 
ggsave(p,file="fc.svg", width=8, height=8)
