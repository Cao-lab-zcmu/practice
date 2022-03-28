###  ggplot2_heatmap.R
library(ggplot2)
library(ggtree)
library(aplot)
library(tidyr)
library(rayshader)
library(ggsci)
# instance: 2655 3058 2631 4719 3591 2629 2835 3378 3188 2528 1121 525 196 1683 3983 2428 34 4322 2943 2053 4291 3038 403 1637 4684 2028 274 347 495 2268
###########################################################################
###########################################################################
file="0703_all/ftalign.tsv"
data <- read.csv(file=file,header=T,sep="\t",row.names=1, check.name=F)
name=data.frame(strsplit(colnames(data), split="initial_8_neg_"))[2,]
colnames(data)=name
rownames(data)=name
instance <- c(sample(colnames(data),27), c("347", "495", "2268"))
data <- data[colnames(data) %in% instance, colnames(data) %in% instance]
###########################################################################
###########################################################################
savename="instance.svg"
phr <- hclust(dist(data)) %>% ggtree(layout="rectangular", branch.length="none")
phc <- hclust(dist(t(data))) %>% ggtree(layout="rectangular", branch.length="none") + layout_dendrogram()
data$names <- rownames(data)
p1 <- gather(data, 1:(ncol(data)-1), key="condition", value='expr')  ## keep the last col
pp <- ggplot(p1,aes(x=condition,y=names,fill=expr)) + 
  geom_tile(size=0.5, color="black") +
  #theme_minimal()+
  scale_fill_viridis_c() +
  #scale_fill_npg() +
  scale_y_discrete(position="right") +
  labs(x="Feature ID", y="Feature ID", fill="FTAS") +
  theme(
        axis.text.x = element_text(angle=90),
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
  ggsave(pp_com,file=savename, width=16, height=15)
  data <- data[, c(ncol(data), 1:(ncol(data)-1))]
  write.table(data, file="instance_data.tsv", col.names=T, row.names=F, sep="\t")
  ##############################################################################
  ##############################################################################
  ###  normalized data
  file="norm_instance_data.tsv"
  data <- read.csv(file=file,header=T,sep="\t",row.names=1, check.name=F)
  ###########################################################################
  ###########################################################################
  savename="norm_instance.svg"
  phr <- hclust(dist(data)) %>% ggtree(layout="rectangular", branch.length="none")
  phc <- hclust(dist(t(data))) %>% ggtree(layout="rectangular", branch.length="none") + layout_dendrogram()
  data$names <- rownames(data)
  p1 <- gather(data, 1:(ncol(data)-1), key="condition", value='expr')  ## keep the last col
  pp <- ggplot(p1,aes(x=condition,y=names,fill=expr)) + 
    geom_tile(size=0.5, color="black") +
    #theme_minimal()+
    scale_fill_viridis_c() +
    #scale_fill_npg() +
    scale_y_discrete(position="right") +
    labs(x="Feature ID", y="Feature ID", fill="NFTAS") +
    theme(
          axis.text.x = element_text(angle=90),
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
    ggsave(pp_com,file=savename, width=16, height=15)
    #plot_gg(pp_com, multicore = TRUE, width = 20 ,height=20, scale=250) #加载图形
    #render_depth(focallength=100,focus=0.72)
    library(ggupset)
    library(reshape2)
    ### instance: 627 880 362 835 26 289 39 482 871 213 433 636 609 295 132 245 15 162 740 599 585 412 865 446 116 580 766 221 469 663 244 894 774 827 782 408 841 830 488 832 397 114 848 722 776 485 645 235 186 668
    file="norm_instance_data.tsv"
    data <- read.csv(file=file,header=T,sep="\t",row.names=1, check.name=F)
    data$names <- rownames(data)
    data <- melt(data, measure.vars=rownames(data), variable.name="condition", value.name="expr")
    instance <- sample(1:nrow(data), 50)
    data <- data[instance, ]
    data$upset_x=paste0(data$names, "_", data$condition)
    p <- ggplot(data, aes(x=factor(upset_x), y=expr, fill=expr)) +
      geom_col() +
      #scale_x_mergelist(sep = "_") +
      axis_combmatrix(sep = "_") +
      scale_fill_viridis_c() +
      #coord_flip() +
      labs(x="Feature link", y="NFTAS") +
      theme(
            #axis.text.x = element_text(angle=90),
            axis.text.y = element_text(size=20, angle=90),
            axis.title.x = element_text(size=30, face="bold", angle=180),
            axis.title.y = element_text(size=30, face="bold"),
            #axis.text = element_blank(),
            legend.title = element_text(size=20, face="bold"),
            legend.text = element_text(size=20),
            legend.position = "none",
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
######################
######################
library(ggplot2)
library(ggtree)
library(aplot)
library(tidyr)
library(rayshader)
library(ggsci)
file="0703_all/ftalign.tsv"
savename="small_heatmap.svg"
data <- read.csv(file=file,header=T,sep="\t",row.names=1, check.name=F)
name=data.frame(strsplit(colnames(data), split="initial_8_neg_"))[2,]
colnames(data)=name
rownames(data)=name
data <- data[colnames(data) %in% c("347", "495", "2268"), colnames(data) %in% c("347", "495", "2268")]
phr <- hclust(dist(data)) %>% ggtree(layout="rectangular", branch.length="none")
phc <- hclust(dist(t(data))) %>% ggtree(layout="rectangular", branch.length="none") + layout_dendrogram()
data$names <- rownames(data)
p1 <- gather(data, 1:(ncol(data)-1), key="condition", value='expr')  ## keep the last col
pp <- ggplot(p1,aes(x=condition,y=names,fill=expr)) + 
  geom_tile(color="black", size=2) +
  #theme_minimal()+
  scale_fill_viridis_c() +
  #scale_fill_npg() +
  scale_y_discrete(position="right") +
  labs(x="Feature ID", y="Feature ID", fill="Normalized\nftalign\nsimilarity\n") +
  theme_minimal() +
  theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        #axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        panel.grid = element_blank()
        #legend.text = element_text(size=20),
        #legend.key.width = unit(2, "cm"),
        #legend.key.height = unit(4, "cm"),
        #text=element_text(family="serif")
        #axis.title.y = element_text(size = 14),
        #plot.title = element_text(hjust = 1,vjust=-40,size=14)
  )
  pp_com <- pp %>% 
    insert_left(phr, width=.1) %>%
    insert_top(phc, height=.1) 
  ggsave(pp_com,file=savename, width=5, height=5)

