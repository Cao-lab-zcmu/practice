library(tidyverse)
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
nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t")
nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
nodes2$classification <- nodes2$definition
nodes2$classification[grep("null", nodes2$classification)] <- "Undifined"
nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=T)
nodes <- nodes[which(nodes$superclass!=""), ]
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
## grid 
svg("child_nebula.svg", width=18*1.6, height=22*1.5)
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
 
