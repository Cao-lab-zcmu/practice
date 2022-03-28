library(tidyverse)
library(ggsci)
library(scales)
library(grid)
library(stringr)
library(reshape2)
library(hrbrthemes)
plimit=0.5
datapath1="stat_classification.tsv"
data <- read.csv(file=datapath1, header=T, sep="\t", quote = "")
df <- data[which(data$definition!="null"),]
datapath2="fingerid_first_score.tsv"
data2 <- read.csv(file=datapath2, header=T, sep="\t", quote = "")
df <- merge(df, data2[,colnames(data2) %in% c("id", "similarity")], by="id", all.x=T, sort=T)
df$similarity <- as.numeric(df$similarity)
instance <- c(2028, 274, 347, 495, 2268)
test <- df[df$id %in% instance, ]
## another: "specific_pp", "definition_pp"
facet <- c("specific_pp", "level_5_pp", "subclass_pp", "class_pp", "superclass_pp")
re_facet <- c("Most specific class", "Level 5", "Subclass", "Class", "Superclass")
colnames(test)[colnames(test) %in% facet]=re_facet
test <- melt(test, measure.vars=re_facet, variable.name="condition", value.name="expr")
p <- ggplot(test) + 
  geom_col(aes(x=2, y=expr, 
               fill=ifelse(expr>plimit, 
                           ifelse(condition==re_facet[1], "PPCP (filter)", "PPCP"), "filter")), 
           color="black") + # class
geom_col(aes(x=3, y=plimit, fill="PPCP threshold"), color="black") + # class expect
  geom_col(aes(x=1, y=0, 
               fill=ifelse(similarity>0.4, 
                           ifelse(condition==re_facet[1], "filter", "similarity"), "filter"))) + # similarity
geom_col(aes(x=4, y=0, fill="similarity threshold")) + # similarity expect
  theme_ipsum() +
    scale_fill_manual(values=c("PPCP (filter)"="#ADB6B6FF", "PPCP"="#4DBBD5FF",
                               "PPCP threshold"="#0073C2FF")) +
labs(fill="Status") +
guides(fill=guide_legend(nrow=1)) +
theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 6),
      panel.spacing = unit(0.1, "lines"),
      plot.background = element_rect(fill = "white", size=0),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom", 
      legend.key.height = unit(0.3, "cm"), 
      legend.key.width = unit(0.3, "cm"),
      legend.text = element_text(size=8), 
      legend.title = element_text(size=8),
      strip.text.x = element_text(size = 7, family="Times", face="bold", hjust=0.5),
      strip.text.y = element_blank(),
      text=element_text(family="Times"), 
      panel.border = element_blank(),
      plot.margin =unit(c(0,0.2,0.2,0.2),"cm")
      ) +
xlim(0, 5) +
ylim(0, 1.5) +
facet_grid(id ~ condition) 
 #ggsave(p, file="test.svg")
 ##########################
 ##########################
 eg <- cbind(data.frame(seq(6)), data.frame(c(0, 4, 3, 2, 1, 0)))
colnames(eg) <- c("x", "y")
p2 <- ggplot(eg, aes(x=x, y=y)) +
  geom_step(size=1) +
  theme_classic() +
  labs(y="Class priority") +
  theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(), 
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_blank(),
        plot.margin =unit(c(0.2,0.2,0,0.2),"cm"),
        #plot.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5, size=7, family="Times", face="bold"),
        legend.position = "bottom", 
        text=element_text(family="Times")
  )
  text <- df[df$id %in% instance, colnames(df) %in% c("id", "definition")]
  text$access <- "Access class"
  p3 <- ggplot() +
    geom_text(data=text, aes(x=5, y=5, label=str_wrap(definition, width=25)), 
              size=2, fontface="bold", family="Times") +
xlim(0, 10) +
ylim(0, 10) +
facet_grid(id~access) +
theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      panel.spacing = unit(0.1, "lines"),
      panel.border = element_blank(),
      plot.margin =unit(c(0,0.2,0,0),"cm"),
      #plot.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 7, family="Times", face="bold", hjust=0.5),
      axis.title.y = element_blank(),
      legend.position = "bottom", 
      text=element_text(family="Times")
)
svg("class_priority.svg", height=7, width=7)
#pdf("class_priority.svg", height=7, width=7)
grid.newpage()
pushViewport( viewport(layout = grid.layout(50, 60)) )
print( p2, vp=viewport(layout.pos.row=1:7, layout.pos.col=1:50 ))
print( p, vp=viewport(layout.pos.row=8:50, layout.pos.col=3:50 ))
print( p3, vp=viewport(layout.pos.row=8:46, layout.pos.col=51:60 ))
dev.off()
