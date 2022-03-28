library(tidyverse)

path="~/operation/re_fecal_neg"
setwd(path)
mzmine <- read.csv(file="fecal_neg_mzmine.csv", header=T, sep=",", check.name=F)
#filter the blank
mzmine <- mzmine[,!(1:ncol(mzmine) %in% grep("Blank",colnames(mzmine)))]
metadata <- read.csv(file="metadata.tsv", header=T, sep="\t")
data_area <- mzmine[, grep("ID|m/z|retention|area", colnames(mzmine))]
rownames(data_area) <- data_area$"row ID"
df_area <- data.frame(
		      t(
			data_area[,grep("area", colnames(data_area))]
		      )
)
df_area$file <- rownames(df_area)
df_area <- separate(df_area, col="file", into=c("file", "peak", "area"), remove=T, sep=" ")
# df_area[1:10,(ncol(df_area)-5):ncol(df_area)]
df_area <- df_area[, !(colnames(df_area) %in% c("peak", "area"))]
df_area <- merge(df_area, metadata, all.x=T, by.x="file", by.y="file")
# finish reformat
# compare group as follow
element1 <- c("pro", "raw", "model", "control")
element2 <- c("high", "medium", "low", "multi")
c_group <- data.frame(c(rep("model", 3), "pro", "multi"), c(element1[c(1,2,4)], "raw", "multi"))
colnames(c_group) <- c("n1", "n2")
# compare between multi-subgroup
for(i in 1:nrow(c_group)){
	for(j in element2){
		# grep group
		if(c_group[i,1]=="multi"){
			data <- df_area[grep("pro|raw|model|control|drug", df_area$subgroup), ]
		}else{
		 	data <- df_area[grep(paste0(c_group[i,1],"|",c_group[i,2]), df_area$subgroup), ]
	 	}
		# grep subgroup
		if(j!="multi"){
			data <- data[grep(paste0("model|control|", j), data$subgroup), ]
		}
		assign(paste0("list_", i, "_", j), data)
		# break the single comparation
		if(c_group[i,2]=="control"){
			break
		}
	}
}
# list all the compare group
team <- ls()[c(grep("list_", ls()))]
# single plot of pca
library(ggbiplot)
library(ggsci)
library(scales)
library(ggrepel)
dir.create("pca_plot")
for(i in team){
	df <- get(i)
	dfd <- df[, !(colnames(df) %in% c("file", "name", "subgroup", "group"))]	
	dfm <- df[, (colnames(df) %in% c("file", "name", "subgroup", "group"))]
	#compute pca matrix
  dfd <- dfd[, apply(dfd, 2, var) != 0]
	pca <- prcomp(dfd, scale.=T)
	pca_anno <- as.data.frame(pca[5])
	p <- ggbiplot(pca, obs.scale = 1, 
			   var.scale = 1, 
			   groups = dfm$subgroup, 
			   ellipse = TRUE, 
			   circle = TRUE,
 	 	 	   varname.size=0, 
			   var.axes = F) +
 	 	geom_label_repel(data=pca_anno, 
				   aes(x=x.PC1, y=x.PC2, label=dfm$name),
	  	 	          color="black", alpha=0.5, 
				  fontface="bold", size=2, angle= 0, 
	  	 	          direction="both", segment.size = 0.2, 
				  segment.alpha = 0.3,
	  	 	          inherit.aes = FALSE, hjust = 0,
	  	 	          family="Times") +
 	 	scale_color_npg() +
 	 	scale_fill_npg() +
 	 	theme(legend.position = "right",text=element_text(family="Times"))
	ggsave(p, file=paste0("pca_plot/", i, ".svg"))
}
# pca_facet
origin_team <- team
team <- team[!team %in% "list_3_high"]
for(i in team){
	df <- get(i)
	dfd <- df[, !(colnames(df) %in% c("file", "name", "subgroup", "group"))]
	dfm <- df[, (colnames(df) %in% c("file", "name", "subgroup", "group"))]
  dfd <- dfd[, apply(dfd, 2, var) != 0]
	pca<- prcomp(dfd, scale.=T)
	pca_x=data.frame(pca$x)
	pca_x <- pca_x[, 1:3]
	# join the facet group
	n <- strsplit(i, split="_")[[1]][2]
	deep <- strsplit(i, split="_")[[1]][3]
	if(c_group[n,1]!="multi"){
	 	pca_x$facet_col <- paste0(c_group[n, 1], "_", c_group[n, 2])
	}else{
		pca_x$facet_col <- "multi"
	}
	pca_x$facet_row <- deep 
	pca_x <- cbind(pca_x,dfm)
	# annotate PC value
	pca=summary(pca)
	summary=c(round(pca$importance[2,],4))[1:3]
	summary <- data.frame(summary)
	summary$annotate <- rownames(summary)
	if(c_group[n,1]!="multi"){
	 	summary$facet_col <- paste0(c_group[n,1], "_", c_group[n,2])
	}else{
		summary$facet_col <- "multi"
	}
	summary$facet_row <- deep
	if(i==team[1]){
		data_facet <- pca_x
		data_facet_anno <- summary
	}else{
		data_facet <- rbind(data_facet, pca_x)
		data_facet_anno <- rbind(data_facet_anno, summary)
	}
}

PC1 <- data_facet_anno[which(data_facet_anno$annotate=="PC1"),]
PC2 <- data_facet_anno[which(data_facet_anno$annotate=="PC2"),]
PC1$figure <- paste0("PC1(",PC1$summary*100,"%)")
PC2$figure <- paste0("PC2(",PC2$summary*100,"%)")

# calculate the annotate coord
anno_x=min(as.numeric(data_facet$PC1))*(20/20)
anno_y=max(as.numeric(data_facet$PC2))*(22/20)
# set the color
colors <- c("control"="grey",
  "model"="#374E55FF",
  "drug"="#00A087FF",
  "pro_low"="#FDAE6BFF",
  "pro_medium"="#FD8D3CFF",
  "pro_high"="#E6550DFF",
  "raw_low"="#9ECAE1FF",
  "raw_medium"="#6BAED6FF",
  "raw_high"="#3182BDFF")

p <- ggplot(data_facet, aes(x=as.numeric(PC1), y=as.numeric(PC2), fill=subgroup)) +
 	geom_point(alpha=0.8, size=3, shape=21, stroke=0.1) +
 	stat_ellipse(aes(color=subgroup), level = 0.95) +
 	#scale_color_npg() +
 	#scale_fill_npg() +
 	scale_color_manual(values = colors) +
   	scale_fill_manual(values = colors) +
   	guides(color= "none") +
   	geom_text(data=PC1, aes(x=anno_x*1.4, y=anno_y*1.4, label=figure), 
   	 	 	hjust=0, color="black", 
			fontface="bold",alpha=0.6, 
			size=1.5, inherit.aes = FALSE, 
			family="Times") +
   	geom_text(data=PC2, aes(x=anno_x*1.4, y=anno_y*(18/22)*1.4, label=figure), 
   	 	 	hjust=0, color="black", fontface="bold",
			alpha=0.6, size=1.5, inherit.aes = FALSE, 
			family="Times") +
   	labs(y="PC2", x="PC1", fill="Group") +
   	#scale_x_continuous(limits=c(-60 ,60)) +
   	#scale_y_continuous(limits=c(-60 ,60)) +
 	facet_grid(facet_row ~ facet_col) +
 	theme(legend.position = "right",text=element_text(family="Times"))
ggsave(p, file=paste0("pca_plot/pca_facet.svg"),width=8,height=6.5)

# opls_da plot
team <- origin_team 
escape <- grep(paste0("multi|",which(c_group$n1=="multi" | c_group$n2=="multi")), team)
team <- team[!team %in% team[c(escape)]]
# single plot
dir.create("opls_plot")
library(ropls)
library(scales)
for(i in team){
	df <- get(i)
	dfd <- df[, !(colnames(df) %in% c("file", "name", "subgroup", "group"))]
	dfm <- df[, (colnames(df) %in% c("file", "name", "subgroup", "group"))]
	oplsda <- opls(x = dfd, y = dfm[, "subgroup"], predI = 1, orthoI = NA)
	df <- cbind(oplsda@scoreMN[, 1], oplsda@orthoScoreMN[, 1])
	colnames(df) <- c("h1", paste0("o", 1))
	df <- as.data.frame(df)
	df <- cbind(df, dfm)
	# join the facet group
	n <- strsplit(i, split="_")[[1]][2]
	deep <- strsplit(i, split="_")[[1]][3]
	df$facet_col <- paste0(c_group[n, 1], "_", c_group[n, 2])
	df$facet_row <- deep
	# x lab and y lab text
	x_lab <- paste0("T score[1](", oplsda@modelDF[1, "R2X"] * 100, "%)")
	y_lab <- paste0("Orthogonal T score[1](", oplsda@modelDF[2, "R2X"] * 100, "%)")
	summary <- rbind(x_lab, y_lab)	
	summary <- data.frame(summary)
	# join the facet group
	summary$facet_col <- paste0(c_group[n, 1], "_", c_group[n, 2])
	summary$facet_row <- deep
	summary$annotate <- rownames(summary)
	# vip 
	vip=data.frame(oplsda@vipVn)
	vip=cbind(rownames(vip),vip)
	colnames(vip)=c("id","vip")
	# join the facet group
	vip$facet_col <- paste0(c_group[n, 1], "_", c_group[n, 2])
	vip$facet_row <- deep
	vip$team <- paste0(i, "_", vip$id)
	# gather the data into facet data.frame
	if(i==team[1]){
		data_facet <- df
		data_facet_anno <- summary
		data_vip <- vip
	}else{
		data_facet <- rbind(data_facet, df)
		data_facet_anno <- rbind(data_facet_anno, summary)
		data_vip <- rbind(data_vip, vip)
	}
	# ggplot plot the single plot
	p <- ggplot(df, aes(x=h1, y=o1)) +
		geom_point(alpha=0.8, size=3, shape=21, stroke=0.1,
			   aes(fill=subgroup)) +
		stat_ellipse(aes(color=subgroup), level = 0.95) +
		geom_label_repel(data=df, 
				 aes(x=h1, y=o1, label=name),
				 color="black", alpha=0.5, fontface="bold", 
				 size=2, angle= 0,
				 direction="both", segment.size = 0.2,
				 segment.alpha = 0.3,
				 inherit.aes = FALSE, hjust = 0,
				 family="Times") +
		scale_color_npg() +
		scale_fill_npg() +
		labs(x=x_lab,y=y_lab,title="OPLS-DA") +
		theme(plot.title = element_text(hjust = 0.5),
		      text=element_text(family="serif"))

	ggsave(p, file=paste0("opls_plot/", i,".svg"), height=6, width=8)
}
# opls facet plot 
df <- data_facet[which(data_facet$facet_col!="model_control"),]
df_anno <- data_facet_anno[which(data_facet_anno$facet_col!="model_control"),] 
df_xlab <- df_anno[which(df_anno$annotate=="x_lab"),]
df_ylab <- df_anno[which(df_anno$annotate=="y_lab"),]
anno_x=min(as.numeric(df$h1))*(20/20)
anno_y=max(as.numeric(df$o1))*(22/20)
colors <- c("control"="grey",
  "model"="#374E55FF",
  "drug"="#00A087FF",
  "pro_low"="#FDAE6BFF",
  "pro_medium"="#FD8D3CFF",
  "pro_high"="#E6550DFF",
  "raw_low"="#9ECAE1FF",
  "raw_medium"="#6BAED6FF",
  "raw_high"="#3182BDFF")
p <- ggplot(df, aes(x=as.numeric(h1), y=as.numeric(o1), fill=subgroup)) +
 	geom_point(alpha=0.8, size=3, shape=21, stroke=0.1) +
 	stat_ellipse(aes(color=subgroup), level = 0.95) +
 	#scale_color_npg() +
 	#scale_fill_npg() +
 	scale_color_manual(values = colors) +
   	scale_fill_manual(values = colors) + 
   	guides(color= "none") +
   	geom_text(data=df_xlab, 
		     aes(x=anno_x, y=anno_y*1.4, label=summary),
   	 	 hjust=0, color="black", fontface="bold",alpha=0.6, 
		 size=2, inherit.aes = FALSE, family="Times") +
   	geom_text(data=df_ylab, 
		     aes(x=anno_x, y=anno_y*(18/22)*1.4, label=summary),
   	 	 hjust=0, color="black", fontface="bold",alpha=0.6,
		 size=2, inherit.aes = FALSE, family="Times") +
   	labs(y="Orthogonal T score[1]", x="T score[1]", fill="Group") +
   	#scale_x_continuous(limits=c(-60 ,60)) +
   	#scale_y_continuous(limits=c(-60 ,60)) +
 	facet_grid(facet_row ~ facet_col) +
 	theme(legend.position = "right",text=element_text(family="Times"))

ggsave(p,file="opls_plot/opls_facet.svg",width=8,height=6.5)

# volcano_plot
team <- origin_team
escape <- grep(paste0("multi|",which(c_group$n1=="multi" | c_group$n2=="multi")), team)
team <- team[!team %in% team[c(escape)]]
# single plot
dir.create("volcano")
for(i in team){
	df <- get(i)
	dfd <- df[, !(colnames(df) %in% c("file", "name", "subgroup", "group"))]
	dfm <- df[, (colnames(df) %in% c("file", "name", "subgroup", "group"))]
	compare <- unique(dfm$group)
	xn <- which(dfm$group==compare[1])
	yn <- which(dfm$group==compare[2])
	fc_name <- paste0(compare[1], "_d_", compare[2])
	for(j in colnames(dfd)){
		# t.test calculate p.value
		x <- dfd[c(xn), colnames(dfd) %in% j]
		y <- dfd[c(yn), colnames(dfd) %in% j]
		stat=t.test(x, y, var.equal = T, paired = F)$p.value
		assign(paste0("p.value_", j), stat)
		# mean and fc
		fc=mean(x)/mean(y)
		assign(paste0("fc_", j), fc)
	}
	# gather the fc list
	fc_list <- ls()[c(grep("fc_X", ls()))]
	fc <- data.frame(fc_list)
	fc$fc <- NA
	for(k in 1:nrow(fc)){
		fc$fc[k] <- get(fc$fc_list[k])
	}
	fc <- separate(fc, col="fc_list", sep="fc_", 
		       into=c("m", "id"), remove=T)[, 2:3]
	# gather the p.value list
	p.value_list <- ls()[c(grep("p.value_X", ls()))] 
	p.value <- data.frame(p.value_list)
	p.value$p.value <- NA
	for(k in 1:nrow(p.value)){
		p.value$p.value[k] <- get(p.value$p.value_list[k])
	}
	p.value <- separate(p.value, col="p.value_list", sep="p.value_", 
			    into=c("m", "id"), remove=T)[, 2:3]
	# merge
	fc_p <- merge(fc, p.value, by="id", all.x=T, sort=T)
	fc_p$fc <- log2(fc_p$fc)
	fc_p$change <- factor(ifelse(fc_p$p.value < 0.05 & abs(fc_p$fc) >= 1,
                      ifelse(fc_p$fc >= 1,"up","down"),"stable"),
                      levels = c("up","down","stable"))
	#plot volcano
	data <- fc_p
	title <-  paste0(compare[1], "/", compare[2])
	p <- ggplot(data,aes(x=fc, y=-log10(p.value),color = change)) + 
		geom_point(alpha=0.8, stroke=0, size=3) + 
		scale_color_manual(values = c("down"="#4DBBD5FF",
					      "stable"="#8491B4FF",
					      "up"="#DC0000FF")) +
		ylim(1,max(-log10(data$p.value))) +
		geom_hline(yintercept = -log10(0.05), linetype=4, size=0.8) +
		geom_vline(xintercept = c(-1,1), linetype=4, size=0.8) + 
		labs(x = "log2(FC)", y="-log10(p-value)", title=title) + 
		geom_text_repel(data = data[data$p.value<0.01 & abs(data$fc) >= 2,],
				aes(label = substring(id, 2)),
				size = 3,family="Times") +
		theme(text=element_text(family="Times"),
			    #axis.line = element_line(colour = "black", size=0.2),
			    #plot.margin = unit(c(3, 1, 3, 1), "cm")
			    plot.title = element_text(hjust = 0.5))
	ggsave(p, file=paste0("volcano/", i, ".svg"),width=8,height=6.5)

	# for facet plot
	deep <- strsplit(i, split="_")[[1]][3]
	fc_p$facet_col <- paste0(compare[1], "/", compare[2])
	fc_p$facet_row <- deep
	fc_p$team <- paste0(i, "_", fc_p$id)
	if(i==team[1]){
		data_facet <- fc_p
	}else{
		data_facet <- rbind(data_facet, fc_p)
	}
}
data <- data_facet[!(1:nrow(data_facet) %in% grep("control", data_facet$facet_col)),]
p <- ggplot(data,aes(x=fc, y=-log10(p.value),color = change)) +
	geom_point(alpha=0.8, stroke=0, size=1.5) +
	scale_color_manual(values = c("down"="#4DBBD5FF",
				      "stable"="#8491B4FF",
				      "up"="#DC0000FF")) +
	ylim(1,max(-log10(data$p.value))) +
	geom_hline(yintercept = -log10(0.05),linetype=4,size=0.8) +
	geom_vline(xintercept = c(-1,1),linetype=4,size=0.8) +
	labs(x = "log2(FC)", y = "-log10(p-value)") +
	geom_text_repel(data = data[data$p.value<0.01 & abs(data$fc) >= 2,],
			aes(label = substring(id, 2)),
			size = 3,family="Times") +
	theme(text=element_text(family="Times"),
	      plot.title = element_text(hjust = 0.5)) +
	facet_grid(facet_row ~ facet_col)
ggsave(p, file=paste0("volcano/volcano_facet.svg"),width=8,height=6.5)

# plot vip-p plot
# merge the vip dataset and p.value dataset
dir.create("vip")
df <- merge(data, data_vip, all.x=T, by="team", sort=T)
p <- ggplot(df, aes(x=p.value, y=vip, color=fc)) +
	geom_point(alpha=0.8, size=1.5, stroke=0) +
	xlim(0,0.5) +
	#scale_color_gradientn(limits = c(-5, +5), 
	#		      breaks = c(-3, 0, +3), 
	#		      colours = c("#E6550DFF", "#79AF97FF", "#3182BDFF")) +
	scale_color_viridis_c() +
	labs(y="VIP", x="p-value", color="log2(FC)") +
	geom_hline(yintercept = 1,linetype=4,size=0.8) +
	geom_vline(xintercept = 0.05,linetype=4,size=0.8) +
	facet_grid(facet_row.x ~ facet_col.x) +
	theme(legend.position = "right",
	      text=element_text(family="Times"), 
	      plot.title = element_text(hjust = 0.5))

ggsave(p,file="vip/vip_p_facet.svg")
# filter the data according to vip, p.value, fc
# FDR revise
library(fdrtool)
df <- df[which(is.na(df$p.value)==F),] %>%
  mutate(id=substring(id.x, 2))
df$q.value <- fdrtool(df$p.value, statistic='pvalue', plot=F)$qval
write.table(df, file="discrepancy_raw.tsv", col.names=T, row.names=F, sep="\t")
# p.value filter
data <- df[which(df$fc>1 & df$vip>1 & df$p.value<0.05),]
data_u <- data[!duplicated(data$id.x),] 
write.table(data_u, file="discrepancy_fc1_vip1_p005.tsv", col.names=T, row.names=F, sep="\t")
# get sirius idenfication results
stru_path="0070_results/"
structure <- read.csv(file=paste0(stru_path, "fingerid_first_score.tsv"), sep="\t", header=T)
# merge the structure according to qdata
data1 <- merge(data_u, structure, by="id", all.x=T, sort=T)
write.table(data1, file="discrepancy_pdata_structure.tsv", sep="\t", col.names=T, row.names=F)

# q.value filter
qdata <- df[which(df$fc>1 & df$vip>1 & df$q.value<0.05),]
qdata_u <- qdata[!duplicated(qdata$id.x),]
# merge the structure according to qdata
data0 <- merge(qdata_u, structure, by="id", all.x=T, sort=T)
write.table(data0, file="discrepancy_qdata_structure.tsv", sep="\t", col.names=T, row.names=F)
# pathway enrichment analysis
library(FELLA)
# data0 <- read.csv(file="discrepancy_qdata_structure.tsv", sep="\t", header=T)
# data1 <- read.csv(file="discrepancy_pdata_structure.tsv", sep="\t", header=T)
# according to q.value
qdata_h <- data0[grep("high",data0$team), ]
qdata_h_simi <- qdata_h[which(qdata_h$similarity>=0.5), ]
qdata_h_simi <- qdata_h_simi[!duplicated(qdata_h_simi$smiles),]
# according to p.value
pdata_h <- data1[grep("high",data1$team), ]
pdata_h_simi <- pdata_h[which(pdata_h$similarity>=0.5), ]
pdata_h_simi <- pdata_h_simi[!duplicated(pdata_h_simi$smiles),]
# this table is summarized by fc, vip, p, dosage of high, idenfication of high tanimoto similarity, and is unique.
write.table(pdata_h_simi, file="p_cluster_pathway.tsv", sep="\t", col.names=T, row.names=F)
