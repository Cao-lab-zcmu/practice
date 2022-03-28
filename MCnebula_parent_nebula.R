 library(tidyverse)
 library(igraph)
 library(ggraph)
 library(tidygraph)
 library(ggsci)
 library(scales)
 library(stringr)
 #### Eucommia analysis
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
 guides(fill = guide_legend(override.aes = list(size=5))) +
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

