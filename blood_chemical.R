 library(ggplot2)
 library(ggrepel)
 library(ggsci)
 library(Hmisc)
 library(outliers)
 library(ggthemes)
 library(stringr)
grubbs<-function(x){
  x<-round(x,4)
  grubbs_outliers<-c()
  grubbs_p.value<-c()
  grubbs_g.value<-c()
  grubbs_g<-c()
  grubbs_minormax<-c()
  grubbs_pvalue<-c()
  grubbs_p<-0
  while(grubbs_p<0.05){
    grubbs_outliers<-c(grubbs_outliers,grubbs_minormax)
    grubbs_p.value<-c(grubbs_p.value,grubbs_pvalue)
    grubbs_g<-c(grubbs_g,grubbs_g.value)
    if(sum(x==grubbs_minormax)!=0)x<-x[-which(x==grubbs_minormax)]
    if(sd(x)==0) break
    grubbs_test<-grubbs.test(x,type=10,opposite=F,two.sided=F)
    grubbs_p<-grubbs_test$p.value
    grubbs_pvalue<-grubbs_test$p.value
    grubbs_g.value<-grubbs_test$statistic[1]
    grubbs_a<-strsplit(grubbs_test$alternative," ",fixed=T)
    grubbs_minormax<-as.numeric(unlist(grubbs_a)[3])
  }
  outliner_res<-data.frame(outliers=grubbs_outliers,gvalue=grubbs_g,pvalue=grubbs_p.value)
  return(outliner_res)
}
dixon<-function(x){
    dixon_p.value<-c()
    dixon_q.value<-c()
    dixon_q<-c()
    dixon_pvalue<-c()
    dixon_outliers<-c()
    dixon_minormax<-c()
    dixon_p<-0
    while(dixon_p<0.05){
      dixon_outliers<-c(dixon_outliers,dixon_minormax)
      dixon_p.value<-c(dixon_p.value,dixon_pvalue)
      dixon_q<-c(dixon_q,dixon_q.value)
      if(sum(x==dixon_minormax)!=0)x<-x[-which(x==dixon_minormax)]
      if(sd(x)==0) break
      dixon_test<-dixon.test(x,type=0,opposite=F,two.sided=F)
      dixon_p<-dixon_test$p.value
      dixon_pvalue<-dixon_test$p.value
      dixon_q.value<-dixon_test$statistic[1]
      dixon_a<-strsplit(dixon_test$alternative," ",fixed=T)
      dixon_minormax<-as.numeric(unlist(dixon_a)[3])
    }
    outliner_res<-data.frame(outliers=dixon_outliers,qvalue=dixon_q,pvalue=dixon_p.value)
    return(outliner_res)
  }
 args<-commandArgs(TRUE)
 setwd(args[1])
 file=args[2]
 savename=strsplit(file,split=".tsv")
 df <- read.csv(file=file,header= T,sep= "\t")
 class=unique(df$index)
 if(exists("re_stat")==TRUE){rm(re_stat)}
 if(exists("stat_exclude")==TRUE){rm(stat_exclude)}
 for(i in class){
 data <- df[which(df$index==i),]
 subgroup=unique(data$subgroup)
  	for(k in subgroup){
  	group=data[which(data$subgroup==k),]
  	X=group$level
	#if(exists("x_filter")==TRUE){rm(x_filter)}
 	out1=grubbs(round(X,2))
	out2=dixon(round(X,2))
	out3=NULL
	for(t in X){
	 	if( t>(mean(X)+2*sd(X)) || t<(mean(X)-2*sd(X))){out3=c(out3,t)}}
	out=c(out1$outliers, out2$outliers, out3)
	#num_x=length(x_filter)
	if(is.null(out)){re_group=group}else{re_group=group[which(group$level!=out),]}
	if(is.null(out)){exclude=NULL}else{exclude=group[which(group$level==out),]}
  	if(exists("re_stat")==FALSE){re_stat=re_group}else{re_stat=rbind(re_stat,re_group)}
  	if(exists("stat_exclude")==FALSE){stat_exclude=exclude}else{stat_exclude=rbind(stat_exclude,exclude)}
  	}
 }
list=rbind(c("model","raw"),c("model","pro"),c("raw","pro"),c("model","control"))
t.test=array(c("compare1","compare2","p","index"),c(1,4))
for(chem in class){
	for(row in 1:nrow(list)){
	double=list[row,]
	escape=0
		for(k in c("low","medium","high")){
		if(double[1]!="model"){compare2=paste0(double[1],"_",k)}else{compare2=double[1]}
	 	if(double[2]!="model"){compare1=paste0(double[2],"_",k)}else{compare1=double[2]}
	 	if(double[1]=="model" & double[2]=="control"){compare2="model"; compare1="control"; escape=1}
	  	X=re_stat[which(re_stat$subgroup==compare2 & re_stat$index==chem),] #x~compare2
	  	Y=re_stat[which(re_stat$subgroup==compare1 & re_stat$index==chem),] #y~compare1
	  	stat_p=t.test(X$level, Y$level, var.equal = T, paired = F)
	  	stat_p$p.value=round(stat_p$p.value,5)
	  	#hx=fivenum(X$level)[4]+(fivenum(X$level)[4]-fivenum(X$level)[2])*(1.5)
	  	#hy=fivenum(Y$level)[4]+(fivenum(Y$level)[4]-fivenum(Y$level)[2])*(1.5)
	  	t.test=rbind(t.test,c(compare1, compare2, stat_p$p.value, chem),
	  	 	 	 	c(compare2, compare1, stat_p$p.value, chem))
	  	if(escape==1){break}
		}
	}
}
stat=data.frame(t.test[-1,])
colnames(stat)=t.test[1,]
write.table(stat, file = paste0("ttest",".tsv"), quote = FALSE, append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(re_stat, file = paste0("filter_chemical_index",".tsv"), quote = FALSE, 
 	 	append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(stat_exclude, file = paste0("exclude",".tsv"), quote = FALSE, 
 	 	append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
### loop
df=read.csv(file="filter_chemical_index.tsv",header= T,sep= "\t")
write.table(df, file = paste0("filter_chemical_index_",Sys.Date(),".tsv"), quote = FALSE, 
 	 	append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
class=unique(df$index)
 if(exists("re_stat")==TRUE){rm(re_stat)}
 if(exists("stat_exclude")==TRUE){rm(stat_exclude)}
 for(i in class){
 data <- df[which(df$index==i),]
 subgroup=unique(data$subgroup)
  	for(k in subgroup){
  	group=data[which(data$subgroup==k),]
  	X=group$level
	#if(exists("x_filter")==TRUE){rm(x_filter)}
 	out1=grubbs(round(X,2))
	out2=dixon(round(X,2))
	out3=NULL
	for(t in X){
	 	if( t>(mean(X)+2*sd(X)) || t<(mean(X)-2*sd(X))){out3=c(out3,t)}}
	out=c(out1$outliers, out2$outliers, out3)
	#num_x=length(x_filter)
	if(is.null(out)){re_group=group}else{re_group=group[which(group$level!=out),]}
	if(is.null(out)){exclude=NULL}else{exclude=group[which(group$level==out),]}
  	if(exists("re_stat")==FALSE){re_stat=re_group}else{re_stat=rbind(re_stat,re_group)}
  	if(exists("stat_exclude")==FALSE){stat_exclude=exclude}else{stat_exclude=rbind(stat_exclude,exclude)}
  	}
 }
list=rbind(c("model","raw"),c("model","pro"),c("raw","pro"),c("model","control"))
t.test=array(c("compare1","compare2","p","index"),c(1,4))
for(chem in class){
	for(row in 1:nrow(list)){
	double=list[row,]
	escape=0
		for(k in c("low","medium","high")){
		if(double[1]!="model"){compare2=paste0(double[1],"_",k)}else{compare2=double[1]}
	 	if(double[2]!="model"){compare1=paste0(double[2],"_",k)}else{compare1=double[2]}
	 	if(double[1]=="model" & double[2]=="control"){compare2="model"; compare1="control"; escape=1}
	  	X=re_stat[which(re_stat$subgroup==compare2 & re_stat$index==chem),] #x~compare2
	  	Y=re_stat[which(re_stat$subgroup==compare1 & re_stat$index==chem),] #y~compare1
	  	stat_p=t.test(X$level, Y$level, var.equal = T, paired = F)
	  	stat_p$p.value=round(stat_p$p.value,5)
	  	#hx=fivenum(X$level)[4]+(fivenum(X$level)[4]-fivenum(X$level)[2])*(1.5)
	  	#hy=fivenum(Y$level)[4]+(fivenum(Y$level)[4]-fivenum(Y$level)[2])*(1.5)
	  	t.test=rbind(t.test,c(compare1, compare2, stat_p$p.value, chem),
	  	 	 	 	c(compare2, compare1, stat_p$p.value, chem))
	  	if(escape==1){break}
		}
	}
}
stat=data.frame(t.test[-1,])
colnames(stat)=t.test[1,]
write.table(stat, file = paste0("ttest",".tsv"), quote = FALSE, append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(re_stat, file = paste0("re_filter_chemical_index",".tsv"), quote = FALSE, 
 	 	append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(stat_exclude, file = paste0("re_exclude",".tsv"), quote = FALSE, 
 	 	append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
### loop end here
re_stat$subgroup=capitalize(re_stat$subgroup)
stat$compare1=capitalize(stat$compare1)
t1=stat[which(stat$compare1!="Model" & stat$compare1!="Control"),]
	t1$anno=ifelse(t1$compare2=="model", 
	 	ifelse(t1$p<0.05, ifelse(t1$p<0.01, "**", "*"), ""), 
	 	ifelse(t1$p<0.05, ifelse(t1$p<0.01, "##", "#"), ""))
t2=stat[which(stat$compare1=="Model" & stat$compare2=="control"),]
	t2$anno=ifelse(t2$compare2=="control", 
	 	ifelse(t2$p<0.05, ifelse(t2$p<0.01, "**", "*"), ""), 
	 	ifelse(t2$p<0.05, ifelse(t2$p<0.01, "##", "#"), ""))
tt=rbind(t1[which(t1$anno!=""),], t2[which(t2$anno!=""),])
tt$h=NA
tt$h1=NA
tt$h2=NA
for(i in unique(re_stat$index)){
 	for(j in unique(re_stat$subgroup)){
 	calculate=re_stat[which(re_stat$index==i & re_stat$subgroup==j),]
 	rule=re_stat[which(re_stat$index==i),]
 	h=max(calculate$level) #fivenum(calculate$level)[4]+(fivenum(calculate$level)[4]-fivenum(calculate$level)[2])*(3/2)
 	h1=h+(max(rule$level)-min(rule$level))*(1/20)
 	h2=h+(max(rule$level)-min(rule$level))*(3/20)
 	#if(h > max(calculate$level)+(fivenum(calculate$level)[4]-fivenum(calculate$level)[2])){h=max(calculate$level)}
 	test=tt[which(tt$index==i & tt$compare1==j),]
 	if(nrow(test)!=0){tt[which(tt$index==i & tt$compare1==j),]$h=h; 
 	 	 	  tt[which(tt$index==i & tt$compare1==j),]$h1=h1
 	 	 	  tt[which(tt$index==i & tt$compare1==j),]$h2=h2}
 	}
}
for(i in class){
data=re_stat[which(re_stat$index==i),]
#delta=max(data$level)-min(data$level)
anno=tt[which(tt$index==i),]
savename=i
complement=anno[1,]
anno1=anno[c(which(anno$anno=="*"),which(anno$anno=="**")),]
anno2=anno[c(which(anno$anno=="#"),which(anno$anno=="##")),]
if(nrow(anno1)==0){anno1=complement; anno1$anno=""}
if(nrow(anno2)==0){anno2=complement; anno2$anno=""}
boxplot<-ggplot(data,aes(x=subgroup,y=level,fill=subgroup)) +
  stat_boxplot(geom="errorbar", width=0.4) +
  geom_boxplot(width=0.4) +
  #geom_point(aes(fill=subgroup), shape=21, color="black") +
  geom_jitter(aes(fill=subgroup, x=subgroup), shape=21, color="black", width=0.01, height=0, size=2) +
  stat_summary(fun="mean",geom="point",shape=23,size=3,fill="grey") +
  geom_label(data=anno1, aes(x=compare1, label=anno, y=h1),
  	 	hjust=0.5, color="black", family="Times", 
  	 	alpha=ifelse(anno1$anno[1]=="",0,0.7), label.size=ifelse(anno1$anno[1]=="",0,0.3),
  	 	size=5, inherit.aes = FALSE ) +
  geom_label(data=anno2, aes(x=compare1, label=anno, y=h2),
  	 	hjust=0.5, color="black", family="Times", 
  	 	alpha=ifelse(anno2$anno[1]=="",0,0.7), label.size=ifelse(anno2$anno[1]=="",0,0.3),
  	 	size=4, inherit.aes = FALSE ) +
  scale_fill_manual(values = c("Control"="#ADB6B6FF","Drug"="#95CC5EFF","Model"="#7E6148FF",
   	 	 	 	"Pro_low"="#FDAE6BFF","Pro_medium"="#FD8D3CFF","Pro_high"="#E6550DFF",
   	 	 	 	"Raw_low"="#9ECAE1FF","Raw_medium"="#6BAED6FF","Raw_high"="#3182BDFF")) +
  scale_color_manual(values = c("Control"="#ADB6B6FF","Drug"="#95CC5EFF","Model"="#7E6148FF",
   	 	 	 	"Pro_low"="#FDAE6BFF","Pro_medium"="#FD8D3CFF","Pro_high"="#E6550DFF",
   	 	 	 	"Raw_low"="#9ECAE1FF","Raw_medium"="#6BAED6FF","Raw_high"="#3182BDFF")) +
  labs(y=paste0(unique(data$index),"(", unique(data$unit), ")"), 
   	x="Classification(Control vs model or Eucommia vs model: * ~ p < 0.05; ** ~ p < 0.01. \
   	Pro- or raw- Eucommia in identical dosage compare with each other: # ~ p < 0.05; ## ~ p < 0.01)", 
   	title="Hematochemistry") +
  theme_minimal() +
  theme(legend.position = "right",text=element_text(family="serif", size=10), plot.title = element_text(hjust = 0.5),
   	axis.text.x = element_blank(), plot.background = element_rect(fill ="white", color="white"))
ggsave(boxplot, file=paste0(i,".svg"), width=10, height=6)
}
 #### facet_grid all plot
 data=re_stat
 data$index_unit=paste0(data$index, "(", data$unit, ")")
 anno=tt
 ct=unique(data[,colnames(data) %in% c("index","index_unit")])
 anno=merge(anno, ct, by="index", all.x=TRUE, sort=TRUE)
 anno1=anno[c(which(anno$anno=="*"),which(anno$anno=="**")),]
 anno2=anno[c(which(anno$anno=="#"),which(anno$anno=="##")),]
 complement=anno[1,]
 if(nrow(anno1)==0){anno1=complement; anno1$anno=""}
 if(nrow(anno2)==0){anno2=complement; anno2$anno=""}
 boxplot<-ggplot(data,aes(x=subgroup,y=level,fill=subgroup)) +
  stat_boxplot(geom="errorbar", width=0.4) +
  geom_boxplot(width=0.4) +
  #geom_point(aes(fill=subgroup), shape=21, color="black") +
  geom_jitter(aes(fill=subgroup, x=subgroup), shape=21, color="black", width=0.01, height=0, size=2) +
  stat_summary(fun="mean",geom="point",shape=23,size=3,fill="grey") +
  geom_label(data=anno1, aes(x=compare1, label=anno, y=h1),
  	 	hjust=0.5, color="black", family="Times", 
  	 	alpha=ifelse(anno1$anno[1]=="",0,0.7), label.size=ifelse(anno1$anno[1]=="",0,0.3),
  	 	size=5, inherit.aes = FALSE ) +
  geom_label(data=anno2, aes(x=compare1, label=anno, y=h2),
  	 	hjust=0.5, color="black", family="Times", 
  	 	alpha=ifelse(anno2$anno[1]=="",0,0.7), label.size=ifelse(anno2$anno[1]=="",0,0.3),
  	 	size=4, inherit.aes = FALSE ) +
  scale_fill_manual(values = c("Control"="#ADB6B6FF","Drug"="#95CC5EFF","Model"="#7E6148FF",
   	 	 	 	"Pro_low"="#FDAE6BFF","Pro_medium"="#FD8D3CFF","Pro_high"="#E6550DFF",
   	 	 	 	"Raw_low"="#9ECAE1FF","Raw_medium"="#6BAED6FF","Raw_high"="#3182BDFF")) +
  scale_color_manual(values = c("Control"="#ADB6B6FF","Drug"="#95CC5EFF","Model"="#7E6148FF",
   	 	 	 	"Pro_low"="#FDAE6BFF","Pro_medium"="#FD8D3CFF","Pro_high"="#E6550DFF",
   	 	 	 	"Raw_low"="#9ECAE1FF","Raw_medium"="#6BAED6FF","Raw_high"="#3182BDFF")) +
  labs(y="", 
   	x="Classification(Control vs model or Eucommia vs model: * ~ p < 0.05; ** ~ p < 0.01. \
   	Pro- or raw- Eucommia in identical dosage compare with each other: # ~ p < 0.05; ## ~ p < 0.01)", 
   	title="Hematochemistry") +
  theme_minimal() +
  facet_grid(index_unit~.,scales="free_y") +
  theme(legend.position = "right",text=element_text(family="serif", size=10), plot.title = element_text(hjust = 0.5),
   	axis.text.x = element_blank(), plot.background = element_rect(fill ="white", color="white"))
ggsave(boxplot, file=paste0("chemical_facet",".svg"), width=8, height=20)
 ##################### BUN and CR
 data=re_stat[c(which(re_stat$index=="Urea"), which(re_stat$index=="CR")),]
 data$index_unit=paste0(data$index, "(", data$unit, ")")
 anno=tt[c(which(tt$index=="Urea"), which(tt$index=="CR")),]
 ct=unique(data[,colnames(data) %in% c("index","index_unit")])
 anno=merge(anno, ct, by="index", all.x=TRUE, sort=TRUE)
 anno1=anno[c(which(anno$anno=="*"),which(anno$anno=="**")),]
 anno2=anno[c(which(anno$anno=="#"),which(anno$anno=="##")),]
 complement=anno[1,]
 if(nrow(anno1)==0){anno1=complement; anno1$anno=""}
 if(nrow(anno2)==0){anno2=complement; anno2$anno=""}
 boxplot<-ggplot(data,aes(x=subgroup,y=level,fill=subgroup)) +
  stat_boxplot(geom="errorbar", width=0.4) +
  geom_boxplot(width=0.4) +
  #geom_point(aes(fill=subgroup), shape=21, color="black") +
  geom_jitter(aes(fill=subgroup, x=subgroup), shape=21, color="black", width=0.01, height=0, size=2) +
  stat_summary(fun="mean",geom="point",shape=23,size=3,fill="grey") +
  geom_label(data=anno1, aes(x=compare1, label=anno, y=h1),
  	 	hjust=0.5, color="black", family="Times", 
  	 	alpha=ifelse(anno1$anno[1]=="",0,0.7), label.size=ifelse(anno1$anno[1]=="",0,0.3),
  	 	size=5, inherit.aes = FALSE ) +
  geom_label(data=anno2, aes(x=compare1, label=anno, y=h2),
  	 	hjust=0.5, color="black", family="Times", 
  	 	alpha=ifelse(anno2$anno[1]=="",0,0.7), label.size=ifelse(anno2$anno[1]=="",0,0.3),
  	 	size=4, inherit.aes = FALSE ) +
  scale_fill_manual(values = c("Control"="#ADB6B6FF","Drug"="#95CC5EFF","Model"="#7E6148FF",
   	 	 	 	"Pro_low"="#FDAE6BFF","Pro_medium"="#FD8D3CFF","Pro_high"="#E6550DFF",
   	 	 	 	"Raw_low"="#9ECAE1FF","Raw_medium"="#6BAED6FF","Raw_high"="#3182BDFF")) +
  scale_color_manual(values = c("Control"="#ADB6B6FF","Drug"="#95CC5EFF","Model"="#7E6148FF",
   	 	 	 	"Pro_low"="#FDAE6BFF","Pro_medium"="#FD8D3CFF","Pro_high"="#E6550DFF",
   	 	 	 	"Raw_low"="#9ECAE1FF","Raw_medium"="#6BAED6FF","Raw_high"="#3182BDFF")) +
  labs(y="", 
   	x="Classification(Control vs model or Eucommia vs model: * ~ p < 0.05; ** ~ p < 0.01. \
   	Pro- or raw- Eucommia in identical dosage compare with each other: # ~ p < 0.05; ## ~ p < 0.01)", 
   	title="Hematochemistry") +
  theme_minimal() +
  facet_grid(index_unit~.,scales="free_y") +
  theme(legend.position = "right",text=element_text(family="serif", size=10), plot.title = element_text(hjust = 0.5),
   	axis.text.x = element_blank(), plot.background = element_rect(fill ="white", color="white"))
ggsave(boxplot, file=paste0("BUN_CR_facet",".svg"), width=10, height=8)

