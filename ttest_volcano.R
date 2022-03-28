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
savepath="volcano_ttest"
if(file.exists(savepath)==FALSE){dir.create(savepath)}
data <- read.csv(file="fecal_pos_volcano.tsv",header=T,sep="\t",check.names=F)
group <- unique(data$subgroup)
id_set=colnames(data)
list=rbind(c("model","raw"),c("model","pro"),c("raw","pro"),c("model","control"))
the_filter=0
for(row in 1:nrow(list)){
double=list[row,]
escape=0
	for(k in c("low","medium","high")){
	if(double[1]!="model"){compare2=paste0(double[1],"_",k)}else{compare2=double[1]}
 	if(double[2]!="model"){compare1=paste0(double[2],"_",k)}else{compare1=double[2]}
 	if(double[1]=="model" & double[2]=="control"){compare2="model"; compare1="control"; escape=1}
	sink(paste0(savepath,"/",compare2,"@",compare1),append = FALSE, split = FALSE)
	print(paste0("|id|p-value|number|fold:",compare2,"/",compare1,"|"))
	for(i in 2:ncol(data)){
	X=data[which(data$subgroup==compare1),i]
	n=try(out<-dixon(round(X,3)))
	if(class(n)=="try-error"){n=try(out<-grubbs(round(X,3)))}
	if(class(n)=="try-error"){out=NULL}
	if(is.null(out$outliers)){x_filter=X}else{x_filter=X[which(X!=c(out$outliers))];the_filter=the_filter+1}
	  num_x=length(x_filter)
	Y=data[which(data$subgroup==compare2),i]
	n=try(out<-dixon(round(Y,3)))
	if(class(n)=="try-error"){n=try(out<-grubbs(round(Y,3)))}
	if(class(n)=="try-error"){out=NULL}
	if(is.null(out$outliers)){y_filter=Y}else{y_filter=Y[which(Y!=c(out$outliers))];the_filter=the_filter+1}
	  num_y=length(y_filter)
	stat=t.test(X, Y, var.equal = T, paired = F)
	fold=mean(y_filter)/mean(x_filter)
	print(paste0("|",id_set[i],"|",stat$p.value,"|",num_x,"_",num_y,"|",fold,"|"))}
	sink()
 	print(paste0("fold change: ",compare2," divided by ",compare1))
 	if(escape==1){break}
	}
}
print(the_filter)
