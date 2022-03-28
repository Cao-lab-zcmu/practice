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
library(gridExtra)
path="network_facet_0.3"
edge_list <- list.files(path = path, pattern = "*.tsv$", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
# facet_group <- read.csv(file="../for_violin.tsv", header=T,sep="\t")
# facet_group <- facet_group[order(facet_group$classification),] ### lead
nodes1 <- read.csv(file="../com_compound.tsv",header=T,sep="\t")
nodes2 <- read.csv(file="stat_classification.tsv",header=T,sep="\t")
 nodes <- merge(nodes1, nodes2, by="id", all.x=T, sort=F)
 nodes <- nodes[which(nodes$superclass!=""), ]

################################
################################
stru_path="structure_2d/smiles_draw"
structure_list <- list.files(path = stru_path, pattern = "*.mol.svg.cairo.svg$", all.files = FALSE,
                             full.names = FALSE, recursive = FALSE,
                             ignore.case = FALSE, include.dirs = FALSE)
matrix <- data.frame(structure_list) %>%
  mutate(catch="T")
################################
################################
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
pal10= pal_igv("default")(51)
palette= unique(c(pal0, pal1, pal2, pal3, pal4, pal5, pal6, pal7, pal8, pal9))
select1="Iridoids and derivatives.tsv"
select2="Lignans, neolignans and related compounds.tsv"
select_list=c(34,39)
plot <- list()
i=34
############################################
edges <- read.csv(file=paste0(path, "/", edge_list[i]), header=T,sep="\t")
group_select2 <- c(3938, 2529, 1445)
group_select1 <- c(279, 458, 574, 1107, 1445, 2227, 2529, 2664, 2824)
if(i==34){
  select=group_select1
}else{
  select=group_select2
} ############### filter
colnames(edges)[colnames(edges) %in% "facet"] <- "facet_c"
# select <- c(3938, 2529, 1445) ##### lignans
edges <- edges[(edges$source %in% select | edges$target %in% select), ]
#edges2 <- edges[edges$target %in% select, ]
#edges <- rbind(edges1, edges2)
edges <- edges[!duplicated(edges[,1:2]),]
edges_set <- edges
########################################
tlimit <- ifelse(i==34, 0.2, 0.1)
#for(tlimit in c(0.2, 0.25, 0.3, 0.4)){
edges <- edges_set[which(edges_set$ftalign_similarity > tlimit),] %>% 
  mutate(group=ifelse(source %in% select, source, target), group_sub=ifelse(source %in% select, target, source))
edges <- merge(edges, nodes[,colnames(nodes) %in% c("id", "norm_delta")], by.x="group", by.y="id", all.x=T, sort=F)
colnames(edges)[which(colnames(edges)=="norm_delta")] <- "norm_delta_group"
edges <- merge(edges, nodes[,colnames(nodes) %in% c("id", "norm_delta", "similarity")], by.x="group_sub", by.y="id", all.x=T, sort=F)
colnames(edges)[colnames(edges) %in% c("norm_delta", "similarity")] <- c("norm_delta_group_sub", "sub_similarity")
edges <- edges[which(edges$norm_delta_group * edges$norm_delta_group_sub < 0),]
edges <- edges[, c(3:4, 1:2, 5:ncol(edges))]
s_limit <- ifelse(i==34, 0.5, 0.5)
edges <- edges[which(edges$sub_similarity > s_limit), ]
edges <- edges[order(edges$group, edges$norm_delta_group_sub), ]
edges <- edges[which(edges$ftalign_similarity > tlimit), ]
write.table(edges, file=paste0("structure_2d/candidate/", i, "_class.tsv"), sep="\t", col.names=T, row.names=F)
########################################
network_nodes <- as_tbl_graph(edges) %>%
  activate(nodes) %>% #as_tibble()
  mutate(deg = centrality_degree(mode='in')) %>%
    merge(nodes, by.x="name", by.y="id", all.x=TRUE, sort=F) %>%
    as_tibble()
  network_edges <- as_tbl_graph(edges) %>%
    activate(edges) %>%
    as_tibble()
  #########################################
  ###################
  network_edges_filter <- network_edges[which(network_edges$ftalign_similarity > tlimit),]
  list <- unique(c(network_edges_filter$from, network_edges_filter$to))
  network_nodes$rownames=rownames(network_nodes)
  network <- tbl_graph(nodes = network_nodes, edges = network_edges_filter) %>%
    activate(nodes) %>% #as_tibble()
    mutate(deg = centrality_degree(mode='in'))
    #########################################
    layout_n <- create_layout(network, layout = "circle")
    select <- unique(edges$group)
    outer_circle <- layout_n %>%
      filter(!(name %in% select)) %>%
      mutate(
             x = cos((row_number() - 1) / ( nrow(network_nodes) - length(select) ) * 2 * pi),
             y = sin((row_number() - 1) / ( nrow(network_nodes) - length(select) ) * 2 * pi)
      )
      #############
      angles <- seq(360, 0, -360/length(select) )
      angles <- angles[-1]
      radius <- rep(0.5, length(select))
      centers <- tibble(
                        x = radius * cos(angles / 180 * pi),
                        y = radius * sin(angles / 180 * pi)
      )
      inner_circle <- bind_cols(centers, select(filter(layout_n, name %in% select), -x, -y))
      #############
      layout_n[] <- bind_rows(outer_circle, inner_circle) %>%
        arrange(.ggraph.index)
      #########################################
      elements <- layout_n[,colnames(layout_n) %in% c("name", "x", "y", "similarity")]
      elements$link_structure <- paste0(elements$name,".mol.svg.cairo.svg")
      grid_elements <- merge(elements, matrix, all.x=T, by.x="link_structure", by.y="structure_list", sort=T)
      ### grobify
      structure_p<-list()
      list_stru <- grid_elements[which(grid_elements$catch=="T"), colnames(grid_elements) %in% c("link_structure")]
      # grid.picture(readPicture(paste0(stru_path,"/",j)))
      for(j in list_stru){
        id<-strsplit(j, split=".mol.svg.cairo.svg")[[1]]
        #structure_p[[as.numeric(id)]] = grobify(readPicture(paste0(path,"/",j)))
        assign(paste0("grob_",id), grobify(readPicture(paste0(stru_path,"/",j))))
        cat(id," >>> ", "structure_p\n")
      }
      aes_stru <- mutate(grid_elements[which(grid_elements$catch=="T"),], grob=paste0("grob_",name))
      #########################################
      ratio <- (nrow(network_nodes) - length(select))/length(select)
      plot[[i]] <- 
        ggraph(layout_n) + 
        #geom_edge_fan(aes(edge_width=ftalign_similarity), label_alpha=0.3, color="#4DBBD5FF", show.legend=F) + 
        geom_edge_diagonal(aes(edge_width=ftalign_similarity, edge_color = as.factor(group) ), 
                           label_alpha=0.3, show.legend=T, alpha=0.6) + 
#color="lightblue",
#geom_node_point(aes(fill=str_wrap(classification, width=25), size=similarity), shape=21) + 
geom_node_point(aes( fill=norm_delta*log2(10) ), size=ifelse( layout_n$name %in% select, 40, 40*4.5/ratio), 
                stroke=0,
                alpha=0.6, shape=21) + 
geom_node_text(aes(label=paste0("ID: ", name, "\n", "simi: ", round(similarity,2) )), 
               size=ifelse(layout_n$name %in% select, 5, 5*4.5/ratio), color="white", alpha=0.4) +
#scale_color_manual(values=palette) +
# blue green grey yellow red
scale_fill_gradientn(colors=c("#197EC0FF", "#B8B8B8FF", "#EFC000FF")) +
scale_color_npg() +
scale_edge_colour_manual(values=pal1) +
scale_edge_width(range=c(0.1,4)) + 
facet_edges(~facet_c) +
guides(alpha="none", guide_legend(override.aes = list(alpha = 0.6))) +
labs(fill="Log2(delta area)", edge_width="Normalized\nftalign similarity", edge_colour="From ID") +
theme_grey() +
theme(
      text=element_text(family="Times", size=ifelse(i==34, 10, 15) ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.title= element_text(face="bold"),
      #panel.background = element_rect(fill="white"), 
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.2, "cm"),
      axis.line = element_blank(),
      #panel.grid = element_blank(),
      #legend.position = "none",
      strip.text = element_text(size=ifelse(i==34, 15, 24), face="bold"),
      plot.margin =unit(c(0,0,0,0),"cm"),
      panel.spacing =unit(c(0,0,0,0),"cm")
) 
#assign(paste0("plot_", tlimit), plot[[i]])
#ggsave(plot_0.25, file=paste0("tlimit_", 0.25, ".svg"), width=7, height=5)
#if(tlimit==0.7 & i==28){plot[[i]] <- plot[[i]] + scale_fill_gradient(labels = c("low", "", "", "high")) }
#if(tlimit==0.7 & i==26){plot[[i]] <- plot[[i]] + scale_fill_gradient(labels = c("low", "", "", "","high")) }
n = ifelse(i==34, 3, 1.7)
for(k in 1:nrow(aes_stru)){
  size=aes_stru[k,]$similarity
  id <- aes_stru[k,]$name
  plot[[i]] <- plot[[i]] + geom_subview( x=aes_stru[k,]$x, y=aes_stru[k,]$y,
                                        subview=arrangeGrob( get(aes_stru[k,]$grob) ), 
                                        width=ifelse(id %in% select, size*n/5, size*n/ratio), 
                                        height=ifelse(id %in% select, size*n/5, size*n/ratio) ) }
assign(paste0("plot_", tlimit), plot[[i]])
#cat(i, edge_list[i], "\n")  }
ggsave(plot_0.2, file=paste0("iridoids_", 0.2, ".svg"), width=14, height=12)
#ggsave(plot_0.25, file=paste0("tlimit_", 0.25, ".svg"), width=6.5, height=5)
ggsave(plot_0.1, file=paste0("lignans_", 0.1, ".svg"), width=9, height=7.5)
################################
################################
######## import data
path <- "structure_2d/candidate"
structure_list <- list.files(path = path, pattern = "*.mol.svg.cairo.svg$", all.files = FALSE,
                             full.names = FALSE, recursive = FALSE,
                             ignore.case = FALSE, include.dirs = FALSE)
matrix <- data.frame(structure_list) %>% mutate(catch="T")
matrix <- mutate(matrix, idd=strsplit(structure_list, split=".mol.svg.cairo.svg")) %>%
  separate(c("idd"), c("id", "candidate"), sep="_can_", remove=F)
for(id in unique(matrix$id)){
  df <- matrix[which(matrix$id==id), ]
  df <- df[order(as.numeric(df$candidate)),]
  for(i in 1:nrow(df)){
    # for(i in 20:29){
    assign(paste0("stru_", df[i, colnames(df) %in% c("idd")]), readPicture(paste0(path,"/",df[i, colnames(df) %in% c("structure_list")])) )
    cat(i, "\n")  
  }
}
################################
################################
dir.create(paste0(path,"/merge"))
for(id in unique(matrix$id)){
  df <- matrix[which(matrix$id==id), ]
  df <- df[order(as.numeric(df$candidate)),]
  svg(paste0(path, "/merge/", id, "_merge.svg"), height=1.2, width=12)
  n=0
  for(i in 1:nrow(df)){
    n=n+1
    grid.picture( get(paste0("stru_", df[i, colnames(df) %in% c("idd")])) , x=n/11, y=0.5, height=1, width=1)
    cat("Info: the grid picture number is ", i, " of ID: ", id, "\n")
  }
  dev.off()
}
