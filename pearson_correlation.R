 library(outliers)
 
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

savepath="correlation"

dir.create(savepath)

data <- read.csv(file="fecal_pos_correlation.tsv",header=T,sep="\t",check.names=F)

group=unique(data$subgroup)

data_set=array(NA,dim=c(nrow(data)*3,ncol(data)*3))

for(i in colnames(data)){
 	test<-try(as.numeric(i)); if(is.na(test)==F){begin=i;break}}

for(i in 1:ncol(data)){if(colnames(data)[i]==begin){begin=i}}

#for(i in 3:(begin-1)){data=data[which(is.na(data[,i])==F),]}

for(i in 3:ncol(data)){

 	for(j in group){
 
 	X=data[which(data$subgroup==j),i]
 	
 	n=try(out<-dixon(round(X,3)))
	
	if(class(n)=="try-error"){n=try(out<-grubbs(round(X,3)))}
	
	if(class(n)=="try-error"){out=NULL}
	
	if(is.null(out$outliers)){x_filter=X}else{x_filter=X[which(X!=c(out$outliers))]}
   	 	
 	if(exists("team_filter")==TRUE){
 	 	 	team_filter=c(team_filter,x_filter)
 	 	}else{
 	 	 	team_filter=c(x_filter)}
 	}
 	
 	data_set[1:length(team_filter),i]=team_filter
 	
 	rm(team_filter)}	

sink(paste0(savepath,"/","pearson_p_value"),append = FALSE, split = FALSE)

print("|factor|id|r|p_value|")

name_set=colnames(data)

for(i in 3:(begin-1)){

 	for(j in begin:ncol(data)){

 	stat=cor.test(data_set[,i],data_set[,j])
 	
 	print(paste0("|",name_set[i],"|",name_set[j],"|",as.vector(stat$estimate),"|",stat$p.value,"|"))}}

sink()
