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
## Many steps of this script may not be easy to implement, such as having to draw the compound structure image beforehand, and fine-tuning the position of the subgraphs.
path="network_facet_0.50"
edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
nodes1 <- read.csv(file="fingerid_first_score.tsv", quote="", header=T,sep="\t")
nodes2 <- read.csv(file="stat_classification.tsv", quote="", header=T,sep="\t")
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
## an instance to draw lignans and iridoids of E. ulmoides 
## select1="Iridoids and derivatives.tsv"
## select2="Lignans, neolignans and related compounds.tsv"
select_list=c(34,39)
plot <- list()
## Unfortunately, the structure of the compound must be drawn previously and converted into a Cairo SVG.
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
######################### start plot
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
        theme_minimal() +
        theme(
              text=element_text(family="Times"),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              panel.grid.major =element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              legend.title=element_text(face="bold", hjust= -0.5),
              panel.border = element_blank(),
              plot.margin =unit(c(0,0,0,0),"cm"),
              panel.spacing =unit(c(0,0,0,0),"cm")
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
      legend.key.height = unit(0.6, "cm"),
      axis.line = element_blank(),
      strip.text = element_text(size=15, face="bold"),
)
assign("find_id", plot[[i]])
## As we see it, there must be some bug in the ggimage package that causes the position to be misaligned, so it needs to be fine-tuned manually.
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

