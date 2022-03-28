 library(tidyverse)
 
 library(igraph)
 
 library(ggraph)
 
 library(tidygraph)
 
 library(ggsci)
 
 library(scales)
 
 library(grid)
 
 pal1= pal_jco()(10)
 
 pal2= pal_jama()(7)
 
 pal3= pal_uchicago("dark")(9)
 
 pal4= pal_igv("default")(51)
 
 palette= unique(c(pal1, pal2, pal3, pal4))
 
 #########################
 
 setwd("json_tree")
 
 for(i in c(495, 347, 2268)){
 
 id=i
 
 nodes <- read.csv(file=paste0("nodes_",id,".tsv"),header=T,sep="\t")
 
 edges <- read.csv(file=paste0("edges_",id,".tsv"),header=T,sep="\t")
 
 network_nodes <- as_tbl_graph(edges) %>%
 
  	activate(nodes) %>% #as_tibble()
  	
  	mutate(deg = centrality_degree(mode='in'), shape = ifelse(name==0, 16, 3)) %>%
  	
  	merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=TRUE) %>%
  	
  	as_tibble()
 
 network_edges <- as_tbl_graph(edges) %>%
  	
  	activate(edges) %>%
  	
  	as_tibble()
 
 network <- tbl_graph(nodes = network_nodes, edges = network_edges)
 
 ###### network
 
 ## layout="dendrogram"
 
 layout="tree"
 
 p <- ggraph(network, layout = layout) + 
 
 geom_edge_elbow(edge_width=1, color="black", show.legend=F, strength = 1) + 
 
 geom_node_point(size = 4, shape=3, color="red") + 
 
 #geom_node_label(aes(label=label, fill=label), color="white", size=7,
  	 	 #fill="transparent",
  	 	# label.padding = unit(0.5, "lines"), 
  	 	# label.r = unit(0.25, "lines"), 
  	 	# label.size = 2,
  	 	# fontface="bold",
  	 	# nudge_y = 0.1,
  	 	# family="Times",
  	 	# repel = T) +
 
 scale_color_manual(values=palette) +
 
 scale_fill_manual(values=palette) +
 
 # scale_edge_width(range=c(0.1,0.7)) + 
 
 #facet_edges(~facet) + 
 
 #guides(size="none", fill="none") +
 
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
       strip.text = element_text(size=20, face="bold")
      )
 
 ggsave(p, file=paste0("tree_",id,".svg"), width=4, height=8)  }
 
 
 
 
 
 
 
 
