library(tidyverse)
library(gridtext)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggsci)
library(scales)
library(grid)
library(stringr)
library(ggimage)
library(gridExtra)
#### Eucommia analyses
path="network_facet_0.50"
edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
# facet_group <- read.csv(file="../for_violin.tsv", header=T,sep="\t")
# facet_group <- facet_group[order(facet_group$classification),] ### lead
nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t")
nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
nodes2$classification <- nodes2$definition
nodes2$classification[grep("null", nodes2$classification)] <- "Undifined"
nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=T)
nodes <- nodes[which(nodes$superclass!=""), ]
# edge <- list()
##### color palette
pal1= pal_simpsons()(16)
pal1= pal1[c(14,2,4,6,7,8,13,1,15,16)]
pal2= pal_d3("category20")(20)
pal3= pal_futurama()(12)
pal4= pal_uchicago("light")(9)
pal5= pal_jco()(10)
pal6= pal_startrek()(7)
pal7= pal_ucscgb()(26)
pal8= pal_locuszoom()(7)
pal9= pal_rickandmorty()(12)
palette= unique(c(pal1, pal2, pal3, pal4, pal5, pal6, pal7, pal8, pal9))
plot <- list()
for(i in 1:length(edge_list)){
  edges <- read.csv(file=paste0(path, "/", edge_list[i]), header=T,sep="\t")
  edges_id <- unique(c(edges$source, edges$target))
  cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
  ###########################################
  edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
  edges <- edges[!duplicated(edges[,1:2]), ]
  network_nodes <- as_tbl_graph(edges) %>%
    activate(nodes) %>% #as_tibble()
    mutate(deg = centrality_degree(mode='in')) %>%
      merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=F) %>%
      as_tibble()
    network_edges <- as_tbl_graph(edges) %>%
      activate(edges) %>%
      as_tibble()
    network <- tbl_graph(nodes = network_nodes, edges = network_edges)
    #########
    ######### candidate: graphopt kk fr mds
    layout=ifelse(nrow(network_nodes)>1000, "mds", "fr")
    plot[[i]] <- 
      ggraph(network, layout = layout) + 
      geom_edge_fan(aes(edge_width=ftalign_similarity), color="black", show.legend=F) + 
      geom_node_point(aes(size = as.numeric(similarity), fill=classification), shape=21) + 
      # geom_node_text(aes(label=name),size=3) +
      scale_color_manual(values=palette) +
      scale_fill_manual(values=palette) +
      scale_edge_width(range=c(0.1,0.7)) + 
      facet_edges(~str_wrap(facet, width=30)) + 
      guides(size="none", fill="none") +
      theme_grey() +
      theme(
            text=element_text(family="serif"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.background = element_rect(fill="white"), 
            #axis.line = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            strip.text = element_text(size=15, face="bold")
      )
      # ggsave(plot[[i]], file="test.svg", width=5, height=5)
      cat(i, edge_list[i], "\n") 
}
##
svg("test_child_nebula.svg", width=18*1.6, height=22*1.5)
n=length(edge_list)
s=n^(1/2); if(round(s)!=s){s=round(s); ss=s+1}else{ss=s}
grid.newpage()
pushViewport(viewport(layout = grid.layout(ss, s)))
r_ss=1
for(i in 1:n){
  c_s=i%%s
  if(c_s==0){c_s=s}
  print( plot[[i]], vp=viewport(layout.pos.row=r_ss, layout.pos.col=c_s ))
  cat("push view port of ",i,"\n")
  if(c_s==s){ r_ss=r_ss+1}
}
dev.off()
#############################
#############################
#############################
#############################
#############################
#############################
#############################
#############################
## network zoomed
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggsci)
library(scales)
library(grid)
library(stringr)
library(ggimage)
library(grImport2)
library(gridSVG)
library(rsvg)
library(gridExtra)
path="network_facet_0.50"
edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
# facet_group <- read.csv(file="../for_violin.tsv", header=T,sep="\t")
# facet_group <- facet_group[order(facet_group$classification),] ### lead
# edge <- list()
nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t")
nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
nodes2$classification <- nodes2$definition
nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=T)
nodes <- nodes[which(nodes$superclass!=""), ]
nodes$similarity <- as.numeric(nodes$similarity)
nodes$classification <- gsub("null", "Undefined", nodes$classification)
##### color palette
pal1= pal_simpsons()(16)
pal1= pal1[c(14,2,4,6,7,8,13,1,15,16)]
pal2= pal_d3("category20")(20)
pal3= pal_futurama()(12)
pal4= pal_uchicago("light")(9)
pal5= pal_jco()(10)
pal6= pal_startrek()(7)
pal7= pal_ucscgb()(26)
pal8= pal_locuszoom()(7)
pal9= pal_rickandmorty()(12)
palette= unique(c(pal1, pal2, pal3, pal4, pal5, pal6, pal7, pal8, pal9))
select1="Iridoids and derivatives.tsv"
select2="Lignans, neolignans and related compounds.tsv"
select_list=c(34,39)
plot <- list()
stru_path="structure_2d/smiles_draw"
structure_list <- list.files(path = stru_path, pattern = "*.mol.svg.cairo.svg$", all.files = FALSE,
                             full.names = FALSE, recursive = FALSE,
                             ignore.case = FALSE, include.dirs = FALSE)
matrix <- data.frame(structure_list) %>% mutate(catch="T")
canopus_set <- read.csv(file="canopus_pp_filter.tsv",header=T,sep="\t")
canopus_set <- t(canopus_set)
colnames(canopus_set)=as.character(canopus_set[1,])
canopus_set <- canopus_set[-1,]
######################################################################
metadata_path="../canopus_neg.tsv"
metadata <- read.csv(file=metadata_path, header=T, sep="\t", quote = "")
metadata <- metadata[,c(2,3,4)]
metadata$class <- paste0("C",metadata$absoluteIndex)
######################################################################
######################### start plot
#########################
for(i in select_list){
  edges <- read.csv(file=paste0(path, "/", edge_list[i]), header=T,sep="\t")
  edges_id <- unique(c(edges$source, edges$target))
  cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
  ###########################################
  edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
  edges <- edges[!duplicated(edges[,1:2]), ]
  network_nodes <- as_tbl_graph(edges) %>%
    activate(nodes) %>% #as_tibble()
    mutate(deg = centrality_degree(mode='in')) %>%
      merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=F) %>%
      as_tibble()
    network_edges <- as_tbl_graph(edges) %>%
      activate(edges) %>%
      as_tibble()
    network <- tbl_graph(nodes = network_nodes, edges = network_edges)
    #########
    ######### candidate: graphopt kk fr mds
    layout=ifelse(nrow(network_nodes)>1000, "mds", "fr")
    layout_n <- create_layout(network, layout = layout)
    ##########################
    ##########################
    ### figure elements
    ### structure
    elements <- layout_n[,colnames(layout_n) %in% c("name", "x", "y", "similarity", "classification")]
    elements$link_structure <- paste0(elements$name,".mol.svg.cairo.svg")
    grid_elements <- merge(elements, matrix, all.x=T, by.x="link_structure", by.y="structure_list", sort=T)
    ### grobify
    structure_p<-list()
    list_stru <- grid_elements[which(grid_elements$catch=="T"), colnames(grid_elements) %in% c("link_structure")]
    prefix=c()
    for(j in list_stru){
      id<-strsplit(j, split=".mol.svg.cairo.svg")[[1]]
      #structure_p[[as.numeric(id)]] = grobify(readPicture(paste0(path,"/",j)))
      assign(paste0("grob_",id), grobify(readPicture(paste0(stru_path,"/",j))))
      cat(id," >>> ", "structure_p\n")
    }
    aes_stru <- mutate(grid_elements[which(grid_elements$catch=="T"),], grob=paste0("grob_",name), id=name)
    ##########################
    ##########################
    ##########################
    ### ring_bar plot
    ring_data <- canopus_set[,colnames(canopus_set) %in% c(elements$name)]
    list_palette = data.frame(cbind(sort(unique(elements$classification)), palette[1:length(unique(elements$classification))]))
    colnames(list_palette) <- c("classification", "color")
    elements_palette <- merge(elements, list_palette, all.x=T, by="classification", sort=T)
    for(j in 1:ncol(ring_data)){
      id <- colnames(ring_data)[j]
      fill_in=elements_palette[which(elements_palette$name==id),]$color
      fill_border="black"
      df <- data.frame(ring_data[,colnames(ring_data) %in% c(id)]) 
      df <- mutate(df, class=rownames(df))
      colnames(df)=c("value","class")
      df$num <- seq(nrow(df))
      df <- merge(df, metadata, all.x=T, by="class", sort=T)
      df <- df[order(df$num),]
      df$fill <- paste0(df$class, ": ", df$name)
      df$fill <- factor(df$fill, levels=df[order(df$absoluteIndex), colnames(df) %in% c("fill")])
      p <- ggplot(df, aes(x=num, y=value)) +
        geom_ribbon(aes(x=ifelse(num==1, 0, ifelse(num==nrow(df), num+1, num)), ymin = -5, ymax = 0), fill =fill_in) +
        geom_ribbon(aes(x=ifelse(num==1, 0, ifelse(num==nrow(df), num+1, num)), ymin = 0, ymax = 1.1), fill=fill_border) +
        geom_col(alpha=1, aes(fill=fill), color="white", size=0.02) +
        ylim(-5,1.3) +
        coord_polar() +
        labs(fill="") +
        scale_fill_manual(values=palette) +
        #annotate("text", x=df[nrow(df), colnames(df) %in% "num"], y=-1, 
        # 	 	label=paste0("ID: ", id), hjust=0, family="Times", fontface="bold", size=0.7, alpha=0.5) +
        #scale_fill_gsea() +
        #geom_text(data=df, aes(x=class, y=value+0.1, label=ifelse(value>0.7, class, "")), 
        #	    color="black", fontface="bold",alpha=0.6, size=0.3, inherit.aes = FALSE ) +
        theme_minimal() +
        theme(
              text=element_text(family="Times"),
              axis.ticks = element_blank(),
              #plot.background = element_rect(fill = "transparent"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              panel.grid.major =element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              legend.title=element_text(face="bold", hjust= -0.5),
              panel.border = element_blank(),
              plot.margin =unit(c(0,0,0,0),"cm"),
              panel.spacing =unit(c(0,0,0,0),"cm") # Adjust the margin to make in sort labels are not truncated!
        )
        assign(paste0("ring_", id), p) 
    }
    ##########################
    ##########################
    plot[[i]] <- 
      ggraph(layout_n) + 
      geom_edge_fan(aes(edge_width=ftalign_similarity),
                    color="lightblue", show.legend=F) + 
geom_node_point(aes(fill=str_wrap(classification, width=25)), size=1, shape=21) + 
scale_color_manual(values=palette) +
scale_fill_manual(values=palette) +
scale_edge_width(range=c(0.1,0.7)) + 
facet_edges(~facet) + 
labs(fill="Access classes", size="Tanimoto\nsimilarity") +
guides(fill=guide_legend(override.aes=list(size=5))) +
theme_grey() +
theme(
      text=element_text(family="Times"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.title= element_text(face="bold"),
      #panel.background = element_rect(fill="white"), 
      legend.key.height = unit(0.6, "cm"),
      axis.line = element_blank(),
      #panel.grid = element_blank(),
      #legend.position = "none",
      strip.text = element_text(size=15, face="bold"),
      #plot.margin =unit(c(0,0,0,0),"cm"),
      #panel.spacing =unit(c(0,0,0,0),"cm")
)
assign("find_id", plot[[i]])
## size and posion adjust
min=min(elements$similarity)
delta=max(elements$similarity)-min
step=1/delta
if(i==39 | i==34){p_size=0.85; p_dist=0.105}else{p_size=0.65; p_dist=0.09}
## draw subview into the network
for(k in 1:ncol(ring_data)){
  id=colnames(ring_data)[k]
  aes_ring=elements[which(elements$name==id),]
  size=(aes_ring$similarity-min)*step+p_size
  ## ppcp datacet
  plot[[i]] <- plot[[i]] + geom_subview(x=aes_ring$x-p_dist, y=aes_ring$y-p_dist,
                                        subview=get(paste0("ring_", id)),
                                        width=size, height=size )  
  ## structure
  kk <- which(aes_stru$id==id)
  size=aes_stru[kk,]$similarity
  plot[[i]] <- plot[[i]] + geom_subview(x=aes_stru[kk,]$x, y=aes_stru[kk,]$y,
                                        subview=arrangeGrob( get(aes_stru[kk,]$grob) ),
                                        width=size*6/5, height=size*6/5)
}
#############################################
ggsave(plot[[i]], file=paste0("tt_zoom_",i,".svg"), width=7, height=5.5)
cat(i, edge_list[i], "\n") 
}

p <- find_id + geom_node_text(aes(label=name),size=3) 
#############################################
#############################################
#############################################
### ladder network
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggsci)
library(scales)
library(grid)
library(stringr)
library(ggimage)
library(grImport2)
library(gridSVG)
library(ggpubr)
path="network_facet_ladder1"
edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
# facet_group <- read.csv(file="../for_violin.tsv", header=T,sep="\t")
# facet_group <- facet_group[order(facet_group$classification),] ### lead
nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t")
nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
nodes2$classification <- nodes2$definition
nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=T)
nodes <- nodes[which(nodes$superclass!=""), ]
# edge <- list()
##### color palette
pal0= pal_npg()(10)
pal1= pal_simpsons()(16)
pal2= pal_d3("category20")(20)
pal3= pal_futurama()(12)
pal4= pal_uchicago("light")(9)
pal5= pal_jco()(10)
pal6= pal_startrek()(7)
pal7= pal_ucscgb()(26)
pal8= pal_locuszoom()(7)
pal9= pal_rickandmorty()(12)
palette= unique(c(pal0, pal1, pal2, pal3, pal4, pal5, pal6, pal7, pal8, pal9))
select1="Iridoids and derivatives.tsv"
select2="Lignans, neolignans and related compounds.tsv"
select_list=c(34,39)
plot <- list()
######################### start plot
#########################
#for(i in select_list){
i=39
edges <- read.csv(file=paste0(path, "/", edge_list[i]), header=T,sep="\t")
edges_id <- unique(c(edges$source, edges$target))
cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
###########################################
edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
edges <- edges[!duplicated(edges[,1:2]), ]
network_nodes <- as_tbl_graph(edges) %>%
  activate(nodes) %>% #as_tibble()
  mutate(deg = centrality_degree(mode='in')) %>%
    merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=F) %>%
    as_tibble()
  network_edges <- as_tbl_graph(edges) %>%
    activate(edges) %>%
    as_tibble()
  #########################################
  num_edges <- data.frame(matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7))) 
  colnames(num_edges) <- "tlimit"
  num_edges$num=NA
  for(tlimit in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)){
    network_edges_filter <- network_edges[which(network_edges$ftalign_similarity > tlimit),]
    list <- unique(c(network_edges_filter$from, network_edges_filter$to))
    network_nodes$rownames=rownames(network_nodes)
    num_edges[which(num_edges$tlimit==tlimit),]$num <- nrow(network_edges_filter)
    network <- tbl_graph(nodes = network_nodes, edges = network_edges_filter) %>%
      activate(nodes) %>% #as_tibble()
      mutate(deg = centrality_degree(mode='in'))
      #########
      ######### candidate: graphopt kk fr mds
      layout=ifelse(tlimit==0.2, "fr", "fr")
      layout_n <- create_layout(network, layout = layout)
      #########################################
      plot[[i]] <- 
        ggraph(layout_n) + 
        geom_edge_fan(aes(edge_width=ftalign_similarity), label_alpha=0.3, color="#4DBBD5FF", show.legend=F) + 
        #geom_node_point(aes(fill=str_wrap(classification, width=25), size=similarity), shape=21) + 
        geom_node_point(aes(fill=deg, size=as.numeric(similarity)), alpha=ifelse(layout_n$rownames %in% list, 1, 0.1), shape=21) + 
        #scale_fill_gradientn(labels=c("low", "", "", "", "high"), colors=c("#1B1919FF", "#4DBBD5FF")) +
        scale_fill_gradientn(colors=c("#1B1919FF", "#4DBBD5FF")) +
        scale_edge_width(range=c(0.1,0.7)) + 
        #facet_edges(~facet) + 
        guides(alpha="none") +
        labs(fill="Centrality\ndegree", size="Tanimoto\nsimilarity") +
        #guides(size="none", fill="none") +
        theme_grey() +
        theme(
              text=element_text(family="serif", size=20),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              legend.title= element_text(face="bold", size=15),
              #panel.background = element_rect(fill="white"), 
              legend.key.height = unit(0.3, "cm"),
              legend.key.width = unit(0.5, "cm"),
              axis.line = element_blank(),
              #panel.grid = element_blank(),
              #legend.position = "none",
              #strip.text = element_text(size=15, face="bold")
              plot.margin =unit(c(0,0,0,0),"cm"),
              panel.spacing =unit(c(0,0,0,0),"cm")
        )
        if(tlimit==0.7 & i==39){
          plot[[i]] <- plot[[i]] + scale_fill_gradientn(labels = c("low", "", "", "high"), 
                                                        colors=c("#1B1919FF", "#4DBBD5FF")) 
        }
        assign(paste0("plot_", tlimit), plot[[i]])
        #ggsave(plot[[i]], file=paste0("tlimit_", tlimit, ".svg"), width=7, height=5)
        cat(i, edge_list[i], "\n")  
  }
  #######################################################################
  p_step1 <- ggplot(num_edges, aes(x=as.factor(tlimit), y=num)) + 
    geom_col() +
    labs(x="NFTAS (PPCP ≥ 0.5)", y="Edges")+
    theme_classic() +
    theme(
          text=element_text(family="Times", size=20),
          axis.title = element_text(face="bold"),
          axis.title.y = element_text(hjust=1),
          plot.margin =unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.spacing =unit(c(0,0,0,0),"cm")
          ) +
      coord_cartesian(ylim = c(0,1500))
    p_step2 <- ggplot(num_edges, aes(x=as.factor(tlimit), y=num)) + 
      geom_col() +
      labs(x=NULL, y="number")+
      theme_classic() +
      ggtitle("Lignans, neolignans and related compounds") +
      theme(
            text=element_text(family="Times", size=20),
            axis.title.y = element_text(hjust=0, face="bold"), 
            plot.margin =unit(c(0.1,0.1,0.1,0.1),"cm"),
            plot.title = element_text(hjust=0.5, size=25, face="bold"),
            panel.spacing =unit(c(0,0,0,0),"cm"),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank() 
            ) +
      coord_cartesian(ylim = c(6500,7000))
    p_step <- ggarrange(p_step2, p_step1, heights=c(2/5, 3/5),ncol = 1, nrow = 2, common.legend = TRUE, legend="right", align = "v")
    # p2 <- gg.gap(plot = p_step,
    #           segments = c(1000, 3400),
    #           fontfamily = "serif",
    #           #tick_width = 10,
    #rel_heights = c(0.25, 0, 0.1),# 设置分隔为的三个部分的宽度
    #           ylim = c(0, 3800)
    #           ) +
    #       theme(
    #  	 axis.title = element_text(face="bold"),
    #  	 plot.margin =unit(c(0,0,0,0),"cm"),
    #         panel.spacing =unit(c(0,0,0,0),"cm"),
    #         text=element_text(family="serif")
    #        )
    #######################################################################
    #png("test.png", width=1000, height=1000)
    #tiff("step_network.tiff", width=1200*7/5, height=400*7/5, compression="lzw")
    svg("step_network.svg", width=24, height=8)
    #pdf("step_network.pdf", width=24, height=8)
    grid.newpage()
    n=0
    adjust=3.8
    posi_just=3
    pushViewport( viewport(layout = grid.layout(40, 60+adjust+posi_just)) )
    print( p_step, vp=viewport(layout.pos.row=1:20, layout.pos.col=1:(60+adjust) ))
    for(tlimit in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)){
      n=n+1
      print( get( paste0("plot_", tlimit) ), 
            vp=viewport(layout.pos.row=21:39, 
                        layout.pos.col=( ((n-1)*10+1+posi_just):(n*10+adjust+posi_just) ) 
            )
            ) }
    grid.text("Morphology", x=0.02, y=0.25, rot=90,
              gp = gpar( fontface = "bold", fontsize = 20, fontfamily = "Times", fontangle=90) )
    dev.off()
    #}
    ############################
    ############################
    ############################
    ############################
    ############################
    ############################
    library(tidyverse)
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(ggsci)
    library(scales)
    library(grid)
    library(stringr)
    library(ggimage)
    library(grImport2)
    library(gridSVG)
    library(ggpubr)
    library(reshape2)
    # facet_group <- read.csv(file="../for_violin.tsv", header=T,sep="\t")
    # facet_group <- facet_group[order(facet_group$classification),] ### lead
    nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t")
    nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
    nodes2$classification <- nodes2$definition
    nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=T)
    nodes <- nodes[which(nodes$superclass!=""), ]
    ##############################
    ##############################
    nums <- data.frame(matrix(c(0.3, 0.5, 0.7, 0.9, 0.95, 0.99))) 
    colnames(nums) <- "pp"
    nums$e_num=NA
    nums$n_num=NA
    ##############################
    ##############################
    for(pp in c(0.3,0.5,0.7,0.9,0.95,0.99)){
      path=paste0("network_facet_", pp)
      edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE)
      name="Lignans, neolignans and related compounds.tsv"
      edges <- read.csv(file=paste0(path, "/", edge_list[which(edge_list==name)]), header=T,sep="\t")
      edges_id <- unique(c(edges$source, edges$target))
      cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
      ###########################################
      edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
      edges <- edges[!duplicated(edges[,1:2]), ]
      network_nodes <- as_tbl_graph(edges) %>%
        activate(nodes) %>% #as_tibble()
        mutate(deg = centrality_degree(mode='in')) %>%
          merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=F) %>%
          as_tibble()
        network_edges <- as_tbl_graph(edges) %>%
          activate(edges) %>%
          as_tibble()
        list <- unique(c(network_edges[which(network_edges$ftalign_similarity!=1),]$from, 
                         network_edges[which(network_edges$ftalign_similarity!=1),]$to))
        network_nodes$rownames=rownames(network_nodes)
        ##############################
        nums[which(nums$pp==pp), ]$e_num= nrow(network_edges)
        nums[which(nums$pp==pp), ]$n_num= nrow(network_nodes)
        ##############################
        ##############################
        ##############################
        network <- tbl_graph( nodes = network_nodes, edges = network_edges ) %>%
          activate(nodes) %>% #as_tibble()
          mutate(deg = centrality_degree(mode='in'))
          #########
          ######### candidate: graphopt kk fr mds
          layout=ifelse(nrow(network_nodes)>1000, "mds", "fr")
          layout_n <- create_layout(network, layout = layout)
          p <- 
            ggraph(layout_n) + 
            geom_edge_fan(aes(edge_width=ftalign_similarity), color="#4DBBD5FF", show.legend=F) + 
            #geom_node_point(aes(fill=str_wrap(classification, width=25), size=similarity), shape=21) + 
            geom_node_point(aes(fill=deg, size=as.numeric(similarity)), 
                            alpha=ifelse(layout_n$rownames %in% list, 1, 0.1), shape=21) + 
      scale_edge_width(range=c(0.1,0.7)) + 
      scale_fill_gradient(low="#1B1919FF", high="#DC0000FF") +
      guides(alpha="none") +
      labs(fill="Centrality\ndegree", size="Tanimoto\nsimilarity") +
      theme_grey() +
      theme(
            text=element_text(family="Times", size=20),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.title= element_text(face="bold", size=15),
            #panel.background = element_rect(fill="white"), 
            legend.key.height = unit(0.3, "cm"),
            legend.key.width = unit(0.5, "cm"),
            axis.line = element_blank(),
            #panel.grid = element_blank(),
            #legend.position = "none",
            #strip.text = element_text(size=15, face="bold")
            plot.margin =unit(c(0,0,0,0),"cm"),
            panel.spacing =unit(c(0,0,0,0),"cm")
      )
      if(pp==0.99){
        p <- p + scale_fill_gradient(low="#1B1919FF", high="#DC0000FF", labels = c("low", "", "", "high", "")) 
      }
      assign(paste0("plot_", pp), p) 
    }
    nums1 <- melt(nums, measure.vars=c("e_num", "n_num"), variable.name="var", value.name="expr")
    p_step <- ggplot(nums1, aes(x=as.factor(pp), y=expr, fill=var)) + 
      geom_col(position=position_dodge(width = 0.9), width = 0.8) +
      labs(x="", y="")+
      theme_classic() +
      scale_fill_uchicago("dark", labels=c("nodes", "edges")) +
      labs(x="PPCP (NFTAS ≥ 0.4)", y="Number", fill="Type") +
      ggtitle("Lignans, neolignans and related compounds") +
      theme(
            text=element_text(family="Times", size=20),
            axis.title = element_text(face="bold"),
            axis.title.y = element_text(hjust=0.5),
            legend.title = element_text(size=15, face="bold"),
            plot.title = element_text(hjust=0.5, size=25, face="bold"), 
            legend.key.height = unit(0.3, "cm"),
            legend.key.width = unit(0.5, "cm"),
            legend.position = c(1.02, 0.5),
            plot.margin =unit(c(0.1,0.1,0.1,0.1),"cm"),
            panel.spacing =unit(c(0,0,0,0),"cm")
            ) +
      coord_cartesian(ylim = c(0,600))
    ##############################
    ##############################
    svg("step2_network.svg", width=24, height=8)
    grid.newpage()
    n=0
    adjust=3.8
    posi_just=3
    pushViewport( viewport(layout = grid.layout(40, 60+adjust+posi_just)) )
    print( p_step, vp=viewport(layout.pos.row=1:20, layout.pos.col=1:(60+adjust) ))
    for(pp in c(0.3, 0.5, 0.7, 0.9, 0.95, 0.99)){
      n=n+1
      print( get( paste0("plot_", pp) ), 
            vp=viewport(layout.pos.row=21:39, 
                        layout.pos.col=( ((n-1)*10+1+posi_just):(n*10+adjust+posi_just) ) 
            )
            ) }
    grid.text("Morphology", x=0.02, y=0.25, rot=90,
              gp = gpar( fontface = "bold", fontsize = 20, fontfamily = "Times", fontangle=90) )
    dev.off()
    #############################
    #############################
    # EX
    metadata_path="../canopus_neg.tsv"
    metadata <- read.csv(file=metadata_path, header=T, sep="\t", quote = "")
    metadata <- metadata[,c(2,3,4)]
    metadata$class <- paste0("C",metadata$absoluteIndex)
    id <- 527
    fill_in=elements_palette[which(elements_palette$name==id),]$color
    fill_border="black"
    df <- data.frame(ring_data[,colnames(ring_data) %in% c(id)]) 
    df <- mutate(df, class=rownames(df))
    colnames(df)=c("value","class")
    df$num <- seq(nrow(df))
    df <- merge(df, metadata, all.x=T, by="class", sort=T)
    df <- df[order(df$num),]
    df$fill <- paste0(df$class, ": ", df$name)
    df$fill <- factor(df$fill, levels=df[order(df$absoluteIndex), colnames(df) %in% c("fill")])
    label_data <- df
    number_of_bar <- nrow(label_data)
    label_data$id <- seq(1,nrow(label_data))
    angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar    
    label_data$hjust<-ifelse( angle < -90, 1, 0)
    label_data$angle<-ifelse(as.numeric(angle) < (-90), angle+180, angle)
    rm(angle)
    p <- ggplot(df, aes(x=num, y=value)) +
      geom_ribbon(aes(x=ifelse(num==1, 0, ifelse(num==nrow(df), num+1, num)), ymin = -5, ymax = 0), fill =fill_in) +
      geom_ribbon(aes(x=ifelse(num==1, 0, ifelse(num==nrow(df), num+1, num)), ymin = 0, ymax = 1.1), fill=fill_border) +
      geom_col(alpha=1, aes(fill=fill), color="white", size=0.1) +
      ylim(-5,1.3) +
      coord_polar() +
      labs(fill="") +
      scale_fill_manual(values=palette) +
      #scale_fill_gsea() +
      geom_text(data=label_data, aes(x=num, y=0.5, label=class, angle=angle), 
                color="white", fontface="bold",alpha=0.8, size=1.5,inherit.aes = FALSE ) +
      theme_minimal() +
      theme(
            text=element_text(family="Times", size=8),
            axis.ticks = element_blank(),
            #plot.background = element_rect(fill = "transparent"),
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major =element_blank(),
            panel.grid.minor = element_blank(),
            legend.key.width = unit(0.5,"cm"),
            legend.key.height = unit(0.5,"cm"),
            #legend.position = "none",
            legend.title=element_text(face="bold", hjust= -0.5),
            panel.border = element_blank(),
            plot.margin =unit(c(0,0,0,0),"cm"),
            panel.spacing =unit(c(0,0,0,0),"cm") # Adjust the margin to make in sort labels are not truncated!
      )
      complement <- readPicture(paste0(stru_path,"/", id,".mol.svg.cairo.svg"))
      svg("annotate.svg", width=14, height=7)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(50, 50)))
      #print( pal_grey, vp=viewport(layout.pos.row=1:50, layout.pos.col=1:20 ))
      print( p, vp=viewport(layout.pos.row=1:50, layout.pos.col=1:50 ))
      grid.picture(complement, x=0.18, y=0.5, width=0.5, height=0.5)
      dev.off()
      # ggsave(pal_grey, file="bg.svg")

      ## compare with gnps, so add noise
      library(tidyverse)
      library(igraph)
      library(ggraph)
      library(tidygraph)
      library(ggsci)
      library(scales)
      library(stringr)
      #### Eucommia analyses
      pal1= pal_simpsons()(16)
      pal2= pal_d3("category20")(20)
      pal3= pal_futurama()(12)
      pal4= pal_uchicago("light")(9)
      palette= unique(c(pal1[c(4, 2:3, 1, 5:8,10,11, 13:16)], pal2, pal3, pal4))
      nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t") 
      nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
      nodes2$classification <- nodes2$definition
      nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=F)
      nodes <- nodes[which(nodes$superclass!=""), ]
      edges <- read.csv(file="source_target_tree_0.4.tsv", header=T,sep="\t")
      edges_id <- unique(c(edges$source, edges$target))
      cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
      ###########################################
      edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
      edges <- edges[!duplicated(edges[,1:2]), ]
      ###########  only show pp > 0.9
      network_nodes <- as_tbl_graph(edges) %>%
        activate(nodes) %>% #as_tibble()
        mutate(deg = centrality_degree(mode='in')) %>%
          merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=F) %>%
          as_tibble()
        network_edges <- as_tbl_graph(edges) %>%
          activate(edges) %>%
          as_tibble()
        network <- tbl_graph(nodes = network_nodes, edges = network_edges)
        layout_n <- create_layout(network, layout = "mds")
        #### molecular network
        p <- ggraph(layout_n) + 
          geom_edge_fan(aes(edge_width=ftalign_similarity), color="lightblue", show.legend=F) + 
          geom_node_point(aes(size = as.numeric(similarity), fill=str_wrap(superclass, width=25)), shape=21) + 
          #geom_node_text(aes(filter= deg>12,label=name),size=1) +
          scale_color_manual(values=palette) +
          scale_fill_manual(values=palette) +
          scale_edge_width(range=c(0.1, 0.7)) + 
          #facet_nodes(~classification) + 
          #guides(size="none") +
          labs(fill="Superclass", size="Tanimoto similarity") +
          theme_grey() +
          theme(
                text=element_text(family="Times"),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                panel.background = element_rect(fill="white"), 
                #axis.line = element_blank(),
                legend.key.width = unit(1, "cm"),
                legend.key.height = unit(1.8, "cm"),
                legend.title = element_text(size=20, face="bold", hjust=0.2),
                legend.text = element_text(size=20),
                legend.background = element_rect(fill="transparent"),
                #legend.position = c(0.6,0.25),
                panel.grid = element_blank(),
                strip.text = element_text(size=20, face="bold")
          )
          # ggsave(p, file="parent_network.tiff", width=20, height=22)
          ggsave(p, file="parent_network.svg", width=20, height=16)



