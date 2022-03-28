library(tidyverse)
data <- read.csv(file="sun.tsv",header=T,sep="\t")
data <- data %>%
  gather(key = "observation", value="value", -c(1,2)) 
nObsType <- nlevels(as.factor(data$observation))
data <- data %>%
  arrange(group, name)
data$id <- rep( seq(1, nrow(data)/nObsType), each=nObsType)
data=data[which(data$value!=0),]
data$group_comple<-ifelse(data$value>0,"pro","raw")
label_data <- data %>%
  group_by(id, name) %>%
  summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
p <- 
  ggplot(data) +
  geom_bar(aes(x=as.factor(id), y=value, fill=group_comple),
           colour="#000033",
           stat="identity", alpha=0.5
           ) +
theme_minimal() +
theme(
      text=element_text(family="serif"),
      legend.key.size = unit(35, "pt")
      ) +
labs(x = "", y = "")+
guides(
       fill = guide_legend(order = 1)
       ) +
geom_text(data=label_data, aes(x=id, y=0, label=name, hjust=hjust),
          color="black", alpha=1, fontface="bold", size=15, angle= 90, 
          inherit.aes = FALSE,
          family="Times")
pdf("test.pdf",width=20,height=50)
p
dev.off()
