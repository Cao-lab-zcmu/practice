 library(tidyverse)
 library(igraph)
 library(ggraph)
 library(tidygraph)
 library(ggsci)
 library(scales)
 #### Eucommia analyses
 pal1= pal_simpsons()(16)
 pal2= pal_d3("category20")(20)
 pal3= pal_futurama()(12)
 pal4= pal_uchicago("light")(9)
 palette= unique(c(pal1[c(4, 2:3, 1, 5:8,10,11, 13:16)], pal2, pal3, pal4))
 #lai; palette= unique(c(pal1[c(6:8,9:10,1:4,11:16)], pal2, pal3, pal4))
 nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t") 
 nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
 nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=F)
 nodes <- nodes[which(nodes$superclass!=""),]
 # from mcnebula
 edges <- read.csv(file="source_target_tree_0.4.tsv", header=T,sep="\t")
 edges_id <- unique(c(edges$source, edges$target))
 cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
 ###########################################
 edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
 edges <- edges[!duplicated(edges[,1:2]), ]
 allid <- unique(c(edges$source, edges$target))
 nrow(edges[edges$source!=edges$target,])
 # from gnps
 edges <- read.csv(file="gnps07_merge.tsv", header=T,sep="\t")
 colnames(edges)[1:5] <- c("source", "target",	"delta_m.z", "MEH", "ftalign_similarity")
 edges <- edges[which(edges$ftalign_similarity > 0.6808 & edges$source!=edges$target),]
 edges_id <- unique(c(edges$source, edges$target))
 # edges id which not find in nodes
 cont_edges_id <- edges_id[!(edges_id %in% allid)]
 edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
 edges <- edges[!duplicated(edges[,1:2]), ]
 nrow(edges)
 # insert the connectionless nodes
 edges_id <- unique(c(edges$source, edges$target))
 cont_edges_id <- allid[!(allid %in% edges_id)]
 cont_list <- data.frame(cont_edges_id)
 colnames(cont_list) <- c("source")
 cont_list <- mutate(cont_list, target=source, delta_m.z=0, MEH="N", ftalign_similarity=1)
 edges <- rbind(edges[, 1:5], cont_list)
 # write.table(edges[, c(1,2,5)], file="filter_net_0.6808", sep="\t", quote=F, row.names=F, col.names=F)
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
 layout_n <- create_layout(network, layout = 'stress')
 #### molecular network
 p <- ggraph(layout_n) + 
   geom_edge_fan(aes(edge_width=ftalign_similarity), color="lightblue", show.legend=F) + 
   geom_node_point(aes(size = as.numeric(similarity), fill=str_wrap(superclass, width=25)), shape=21) + 
   #geom_node_text(aes(filter= deg>12,label=name),size=1) +
   scale_color_manual(values=palette) +
   scale_fill_manual(values=palette) +
   scale_edge_width(range=c(0.1, 0.7)) + 
   #facet_nodes(~classification) + 
   labs(fill="Superclass", size="Tanimoto similarity") +
   guides(fill = guide_legend(override.aes = list(size=5))) +
   theme_grey() +
   theme(
         text=element_text(family="Times"),
         axis.ticks = element_blank(),
         axis.text = element_blank(),
         axis.title = element_blank(),
         panel.background = element_rect(fill="white"), 
         #axis.line = element_blank(),
         legend.key.width = unit(1.5, "cm"),
         legend.key.height = unit(1.8, "cm"),
         legend.title = element_text(size=20, face="bold", hjust=0.2),
         legend.text = element_text(size=20),
         legend.background = element_rect(fill="transparent"),
         #legend.position = c(0.6,0.25),
         panel.grid = element_blank(),
         strip.text = element_text(size=20, face="bold") 
   )
   # ggsave(p, file="parent_network.tiff", width=20, height=22)
   ggsave(p, file="gnps_network_merge07.svg", width=19, height=16)
## facet plot
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
 library(stringr)
 #### Eucommia analyses
 path="gnps_network_facet_0.5"
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
 cat(i, edge_list[i], "\n") }
 # theme_graph(base_family = "Times", 
 # 	     foreground = "#8491B4FF",
 # 	     strip_text_size = 20,
 # 	     fg_text_colour = "white",
 # 	     plot_margin = margin(10, 10, 10, 10))
 #############################
 #############################
 #png("test.png", width=1000, height=1000)
 #tiff("merge_grey.tiff", width=2000, height=2200, compression="lzw")
 svg("gnps_child_nebula.svg", width=18*1.6, height=22*1.5)
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
 ## zoom the specific class 
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
 path="gnps_network_facet_0.5"
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
 #grid.picture(readPicture(paste0(stru_path,"/",j)))
 for(j in list_stru){
   id<-strsplit(j, split=".mol.svg.cairo.svg")[[1]]
   #structure_p[[as.numeric(id)]] = grobify(readPicture(paste0(path,"/",j)))
   assign(paste0("grob_",id), grobify(readPicture(paste0(stru_path,"/",j))))
   cat(id," >>> ", "structure_p\n")
 }
 aes_stru <- mutate(grid_elements[which(grid_elements$catch=="T"),], grob=paste0("grob_",name))
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
 assign(paste0("ring_", id), p) }
 ##########################
 ##########################
 plot[[i]] <-
 ggraph(layout_n) +
 geom_edge_fan(aes(edge_width=ftalign_similarity, label=delta_m.z), label_size=0.7, label_alpha=0.3, color="black", show.legend=F) +
 geom_node_point(aes(fill=str_wrap(classification, width=25)), size=1, shape=21) +
 scale_color_manual(values=palette) +
 scale_fill_manual(values=palette) +
 scale_edge_width(range=c(0.1,0.7)) +
 facet_edges(~facet) +
 labs(fill="Access classes", size="Tanimoto\nsimilarity") +
 #guides(size="none", fill="none") +
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
 #plot[[i]] <- plot[[i]] + geom_subview(x=aes_ring$x-(size/4.7)*1.05, y=aes_ring$y-(size/4.7)*1.05,
  	 	# 	 	 	subview=get(paste0("ring_", id)), width=size*4.5*1.05, height=size*4.5*1.05 )
 #ggsave(plot[[i]], file="test.svg", width=7, height=5)
 min=min(elements$similarity)
 delta=max(elements$similarity)-min
 step=1/delta
 if(i==39 | i==34){p_size=0.85; p_dist=0.105}else{p_size=0.65; p_dist=0.09}
 for(k in 1:ncol(ring_data)){
 id=colnames(ring_data)[k]
 aes_ring=elements[which(elements$name==id),]
 size=(aes_ring$similarity-min)*step+p_size
 plot[[i]] <- plot[[i]] + geom_subview(x=aes_ring$x-p_dist, y=aes_ring$y-p_dist,
  	 	 	 	 	subview=get(paste0("ring_", id)), width=size, height=size )  }
 #############################################
 #############################################
 for(k in 1:nrow(aes_stru)){
 size=aes_stru[k,]$similarity
 plot[[i]] <- plot[[i]] + geom_subview(x=aes_stru[k,]$x, y=aes_stru[k,]$y,
  	 	 	 	 	subview=arrangeGrob( get(aes_stru[k,]$grob) ), width=size*6/5, height=size*6/5) }
 # plot[[i]] <- plot[[i]] + facet_zoom(xlim = c(2, 4), ylim=c(5,7))
 ggsave(plot[[i]], file=paste0("gnps_zoom_",i,".svg"), width=7, height=5)
 cat(i, edge_list[i], "\n")
 }
 p <- find_id + geom_node_text(aes(label=name),size=3)
 #############################################################################################
 ## add noise (other class)
 library(tidyverse)
 library(igraph)
 library(ggraph)
 library(tidygraph)
 library(ggsci)
 library(scales)
 #### Eucommia analyses
 pal1= pal_simpsons()(16)
 pal2= pal_d3("category20")(20)
 pal3= pal_futurama()(12)
 pal4= pal_uchicago("light")(9)
 palette= unique(c(pal1[c(4, 2:3, 1, 5:8,10,11, 13:16)], pal2, pal3, pal4))
 #lai; palette= unique(c(pal1[c(6:8,9:10,1:4,11:16)], pal2, pal3, pal4))
 nodes1 <- read.csv(file="fingerid_first_score.tsv",header=T,sep="\t")
 nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
 nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=F)
 nodes <- nodes[which(nodes$superclass!=""),]
 # from mcnebula
 edges <- read.csv(file="source_target_tree_0.4.tsv", header=T,sep="\t")
 edges_id <- unique(c(edges$source, edges$target))
 cont_edges_id <- edges_id[!(edges_id %in% nodes$id)]
 ###########################################
 edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
 edges <- edges[!duplicated(edges[,1:2]), ]
 allid <- unique(c(edges$source, edges$target))
 nrow(edges[edges$source!=edges$target,])
 # from gnps
 edges <- read.csv(file="gnps04.tsv", header=T,sep="\t")
 colnames(edges)[1:5] <- c("source", "target",	"delta_m.z", "MEH", "ftalign_similarity")
 edges <- edges[which(edges$ftalign_similarity > 0.6808 & edges$source!=edges$target),]
 edges_id <- unique(c(edges$source, edges$target))
 # edges id which not find in nodes
 cont_edges_id <- edges_id[!(edges_id %in% allid)]
 edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
 edges <- edges[!duplicated(edges[,1:2]), ]
 nrow(edges)
 # insert the connectionless nodes
 edges_id <- unique(c(edges$source, edges$target))
 cont_edges_id <- allid[!(allid %in% edges_id)]
 cont_list <- data.frame(cont_edges_id)
 colnames(cont_list) <- c("source")
 cont_list <- mutate(cont_list, target=source, delta_m.z=0, MEH="N", ftalign_similarity=1)
 edges <- rbind(edges[, 1:5], cont_list)
 # add lignans and iridoids
 path="gnps_network_facet_0.5"
 edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE)
 allid <- c()
 for(i in c(34,39,41)){
   ed <- read.csv(file=paste0(path, "/", edge_list[i]), header=T,sep="\t")
   ed_id <- unique(c(ed$source, ed$target))
   allid <- c(allid, ed_id)
   assign(paste0("from_", i), ed_id)
 }
 allid <- unique(allid)
 edges_id <- unique(c(edges$source, edges$target))
 # edges id which not find in nodes
 cont_edges_id <- edges_id[!(edges_id %in% allid)]
 edges <- edges[!(edges$source %in% cont_edges_id | edges$target %in% cont_edges_id), ]
 edges <- edges[!duplicated(edges[,1:2]), ]
 nodes <- mutate(nodes, from_class=
                 ifelse(nodes$id %in% from_34, "Iridoids",
                        ifelse(nodes$id %in% from_39, "Lignans", 
                               ifelse(nodes$id %in% from_41, "Long-chain fatty acids", NA)
                        )
                        ))
 
 # write.table(edges[, c(1,2,5)], file="filter_net_0.6808", sep="\t", quote=F, row.names=F, col.names=F)
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
 layout_n <- create_layout(network, layout = 'fr')
 #### molecular network
 p <- ggraph(layout_n) +
   geom_edge_fan(aes(edge_width=ftalign_similarity), color="lightblue", show.legend=F) +
   geom_node_point(aes(size = as.numeric(similarity), fill=str_wrap(from_class, width=25)), shape=21) +
   #geom_node_text(aes(filter= deg>12,label=name),size=1) +
   scale_color_manual(values=palette) +
   scale_fill_manual(values=palette) +
   scale_edge_width(range=c(0.1, 0.7)) +
   #facet_nodes(~classification) +
   labs(fill="From class", size="Tanimoto similarity") +
   guides(fill=guide_legend(override.aes=list(size=5))) +
   theme_grey() +
   theme(
         text=element_text(family="Times"),
         axis.ticks = element_blank(),
         axis.text = element_blank(),
         axis.title = element_blank(),
         panel.background = element_rect(fill="white"),
         #axis.line = element_blank(),
         legend.key.width = unit(1, "cm"),
         legend.key.height = unit(1, "cm"),
         legend.title = element_text(size=15, face="bold", hjust=0.2),
         legend.text = element_text(size=15),
         legend.background = element_rect(fill="transparent"),
         #legend.position = c(0.6,0.25),
         panel.grid = element_blank(),
         strip.text = element_text(size=20, face="bold")
   )
   # ggsave(p, file="parent_network.tiff", width=20, height=22)
   ggsave(p, file="gnps_add_noise_34_39_41.svg", width=12, height=10)
   #ggsave(p, file="test_gnps_add_noise_34_39.svg", width=12, height=10)

