 library(xcms)
 library(MSnbase)
 setwd("/media/wizard/back/thermo_mzML_0518")
 savepath="/home/wizard/operation/back/"
 dataname="EU-Pro2.mzML"
 com_data <- readMSData(file = dataname, mode = "onDisk")
   ####### ms1 in peaks
   data <- filterRt(com_data, c(600, 780))
   mzrange <- c(380, 390)  # build the eic model
   ex_data <- chromatogram(data, msLevel = 1L, mz = mzrange, aggregationFun = "max")
   ex_data_1 <- ex_data[1,1]
   rt=data.frame(rtime(ex_data_1))
   int=data.frame(intensity(ex_data_1))
   rt_int=cbind(rt,int)
   colnames(rt_int)=c("rt", "int")
   rt_int$rt=rt_int$rt/60
   ########
   n=3
   time=vector(mode="list", length=n)
   time[[1]]=c(10,11)
   time[[2]]=c(11,12)
   time[[3]]=c(12,13)
   peak=vector(mode="list", length=n)
   scan=vector(mode="list", length=n)
   for(i in 1:n){
     peak[[i]]=rt_int[which(rt_int$rt > time[[i]][1] & rt_int$rt < time[[i]][2]),] 
     scan[[i]]=rt_int[which(rt_int$int==max(peak[[i]]$int)),] }
   # rownames(scan1)
   ####### ms1 in scans
   mz=mz(data)
   inn=intensity(data)
   mz_get=vector(mode="list", length=n)
   int_get=vector(mode="list", length=n)
   scan_set=vector(mode="list", length=n)
   a=pi/3 ##### 45
   filter=50000
   xlim=c(10.6,12.7)
   mzfilter=c(360, 420)
   ratio = ((xlim[2])/(mzfilter[2]))  #((xlim[2]-xlim[1])/(mzfilter[2]-mzfilter[1]))
   ratio_x_y = max(rt_int$int)/(xlim[2]-xlim[1])
   lift=(max(rt_int$int))*(1/20)
   height_lift=lift*22
   base=(mzfilter[2]-mzfilter[1])*ratio *ratio_x_y*sin(a)
   #zoom=0.95
   for(i in 1:n){
   mz_get[[i]]=mz[which(names(mz)==rownames(scan[[i]]))] 
   int_get[[i]]=inn[which(names(inn)==rownames(scan[[i]]))] 
   scan_set[[i]]=data.frame(mz_get[[i]], int_get[[i]]);   colnames(scan_set[[i]])=c("mz", "int") 
   scan_set[[i]]=scan_set[[i]][which(scan_set[[i]]$int >= filter &
    	 	 	 	 	scan_set[[i]]$mz >=mzfilter[1] & 
    	 	 	 	 	scan_set[[i]]$mz <= mzfilter[2] ),] 
   scan_set[[i]]$x <- xlim[1] + ((scan_set[[i]]$mz-mzfilter[1])*ratio*cos(a) +(scan[[i]]$rt-xlim[1]))
   scan_set[[i]]$xend <- scan_set[[i]]$x
   scan_set[[i]]$s_x <- xlim[1] + (scan_set[[i]]$mz-mzfilter[1])*ratio*cos(a)
   scan_set[[i]]$s_xend <- scan_set[[i]]$s_x 
   scan_set[[i]]$y <- ((scan_set[[i]]$mz-mzfilter[1])*ratio *ratio_x_y*sin(a) + scan_set[[i]]$int) + lift
   scan_set[[i]]$yend <- ((scan_set[[i]]$mz-mzfilter[1])*ratio *ratio_x_y*sin(a)) + lift
   scan_set[[i]]$group=i
   scan_set[[i]][nrow(scan_set[[i]])+1,] <- c(0, 0, 
    	 	 	 	 	  	scan[[i]]$rt,  ## x
    	 	 	 	 	  	scan[[i]]$rt+(mzfilter[2]-mzfilter[1])*ratio*cos(a), ## xend
    	 	 	 	 	  	NA,
    	 	 	 	 	  	NA,
    	 	 	 	 	  	lift, ## y
    	 	 	 	 	  	base+lift, ## yend
    	 	 	 	 	  	0)  ## group
   scan_set[[i]] <- scan_set[[i]][order(scan_set[[i]]$mz),] }
  segment=scan_set[[1]]
  for(i in 2:n){segment <- rbind(segment, scan_set[[i]])}
  segment_bottom <- segment[which(segment$mz==0), ]
  segment <- segment[which(segment$mz!=0), ]
  segment_outlier=data.frame(rbind(
   	 	 	c(xlim[1], xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), 0, base ), # 1
   	 	 	c(xlim[2], xlim[2]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), 0, base ), # 2
   	 	 	c(xlim[1], xlim[2], 0, 0), # 3
   	 	 	#c(xlim[1], xlim[2], lift*21, lift*21), #top
   	 	 	c(xlim[1], xlim[1], 0, height_lift), # 4
   	 	 	c(xlim[2], xlim[2], 0, height_lift), # 5
   	 	 	c(xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), 
   	 	 	 	 	base, base+height_lift), # 6
   	 	 	# c(xlim[2]+(xlim[2]-xlim[1])*cos(a), xlim[2]+(xlim[2]-xlim[1])*cos(a), lift*21, lift*21*2),
   	 	 	c(xlim[1], xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), height_lift, base+height_lift), # 7
   	 	 	c(xlim[1], xlim[2], height_lift, height_lift)  # 8
   	 	 	))
  segment_outlier$group="0"
  colnames(segment_outlier) <- c("x", "xend", "y", "yend", "group")
 # segment_top=data.frame(rbind(c(xlim[1], xlim[2], base, base)))
 # segment_top$group="0"
 # colnames(segment_top) <- c("x", "xend", "y", "yend", "group")
  delta=(mzfilter[2]-mzfilter[1])*ratio*cos(a)
  point <- data.frame(rbind(c(xlim[1],0),
   	 	 	  c(xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), base),
   	 	 	  c(xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), base + height_lift),
   	 	 	  c(xlim[1], height_lift),
   	 	 	  c(xlim[2],0),
   	 	 	  c(xlim[2], height_lift),
   	 	 	  c(xlim[2]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), base)
   	 	 	  ))
  point$group="0"
  colnames(point)=c("x","y", "group")
  grey <- point[1:4,]
  white <- point[c(1,5,6,4),]
  grey_bottom <- point[c(1,5,7,2),]
  grey_bottom$y1 <- grey_bottom$y -lift*7
  grey_bottom$y2 <- grey_bottom$y -lift*6
  grey_bottom$y3 <- grey_bottom$y -lift*5
  grey_bottom$y4 <- grey_bottom$y -lift*4
  grey_bottom$y5 <- grey_bottom$y -lift*3
 library(ggalt)
 library(tidyverse)
 library(ggsci)
   data=raw_data=rt_int
   data=data[which(data$rt>=xlim[1] & data$rt<=xlim[2]),]
   data=data.frame(spline(data$rt, data$int, n=10000))
  ### beauty
  peak_to=0.15
  scan_tot=scan[[1]]
  for(i in 2:n){scan_tot <- rbind(scan_tot, scan[[i]])}
  for(i in 1:n){
   	scan_tot <- rbind(scan_tot, 
   	 	 	   c(scan_tot[i,1]-peak_to*4/5, scan_tot[i,2]*1/5),
   	 	 	   c(scan_tot[i,1]-peak_to, scan_tot[i,2]*1/10),
   	 	 	   c(scan_tot[i,1]+peak_to*4/5, scan_tot[i,2]*1/5),
   	 	 	   c(scan_tot[i,1]+peak_to, scan_tot[i,2]*1/10) ) 
   	if(i!=n){scan_tot <- rbind(scan_tot, 
   	 	 	   c(scan_tot[i,1]+(scan_tot[i+1,1]-scan_tot[i,1])*1/3, 0), 
   	 	 	   c(scan_tot[i,1]+(scan_tot[i+1,1]-scan_tot[i,1])*1/2, 0),
   	 	 	   c(scan_tot[i,1]+(scan_tot[i+1,1]-scan_tot[i,1])*2/3, 0)
   	 	 	    ) 
   	 	}else{scan_tot <- rbind(scan_tot,
   	 	 	   c(scan_tot[i,1]+peak_to*6/4, 0))} 
   	 	}
  data=rbind(scan_tot, c(xlim[1],0), c(xlim[2],0))
  data=data.frame(spline(data$rt, data$int, n=10000))
  peak_to=peak_to *1.2
   colnames(data)=c("rt", "int")
   group_set=vector(mode="list", length=n)
   for(i in 1:n){
   group_set[[i]] <- data[which(data$rt > as.numeric(scan[[i]][1]-peak_to) & data$rt < as.numeric(scan[[i]][1]+peak_to)),] %>% 
    	 	mutate(group=i) }
   set=group_set[[1]]
   for(i in 2:n){set <- rbind(set, group_set[[i]])}
   data=merge(data, set[,c(1,3)], by="rt", all.x=TRUE, sort=TRUE)
   data[which(is.na(data$group)==T),]$group=0
   color_set <- c("0"="black", "3"="#E64B35FF", "2"="#4DBBD5FF", "1"="#00A087FF")
   ###
   p <- ggplot(data, aes(x=rt, y=int, fill=as.character(group), color=as.character(group))) +
     geom_polygon(data=grey, aes(x=x, y=y), alpha=0.25, fill="grey") +
     geom_segment(data=segment, aes(x=s_x, y=y-lift, xend=s_xend, yend=yend-lift), color="#BDBDBDFF", size=4, alpha=0.8) +
    # geom_polygon(data=grey_bottom, aes(x=x, y=y1), alpha=0.15, fill="skyblue", color="#BDBDBDFF") +
     geom_polygon(data=grey_bottom, aes(x=x, y=y2), alpha=0.15, fill="skyblue", color="#BDBDBDFF") +
     geom_polygon(data=grey_bottom, aes(x=x, y=y3), alpha=0.15, fill="skyblue", color="#BDBDBDFF") +
     geom_polygon(data=grey_bottom, aes(x=x, y=y4), alpha=0.15, fill="skyblue", color="#BDBDBDFF") +
     geom_polygon(data=grey_bottom, aes(x=x, y=y5), alpha=0.15, fill="skyblue", color="#BDBDBDFF") +
     geom_polygon(data=grey_bottom, aes(x=x, y=y), alpha=0.15, fill="skyblue", color="#BDBDBDFF") +
     #geom_area(data=point[!duplicated(point$x),], aes(x=x, y=y), fill=grey, alpha=0.5) +
     geom_segment(data=segment_bottom, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.4) +
     annotate("rect", xmin=xlim[1]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), xmax=xlim[2]+(mzfilter[2]-mzfilter[1])*ratio*cos(a), 
      	 	      ymin=base, ymax=base+height_lift,
            alpha=1, color="#3C5488FF", fill="white", size=2) +
     geom_point(data=point, aes(x=x, y=y), size=1.35) +
     geom_segment(data=segment_outlier, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.95) +
     geom_segment(data=segment, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=4, alpha=0.8) +
     geom_line() +
     geom_polygon(data=white, aes(x=x, y=y), alpha=0.4, fill="white") +
     geom_area(alpha=1) +
     annotate("text", x=scan[[1]]$rt+delta, y=height_lift*8/10+base, label="LC-MS data set", family="Times", size=8, hjust=0) +
     annotate("text", x=scan[[1]]$rt, y=height_lift*8/10, label="The Sample", family="Times", size=8, hjust=0) +
     annotate("text", x=scan[[median(1:n)]]$rt, y=height_lift*1/10*(-1), label="Retention time", family="Times", size=8, hjust=0) +
     annotate("text", x=xlim[1]-peak_to/2, y=lift*7*(-1), label="Samples", family="Times", size=8, hjust=0, angle=90) +
     annotate("text", x=xlim[2]+delta/2+peak_to/2, y=base/2, label="m/z", family="Times", size=8, hjust=0, angle=0, vjust=1) +
   #  geom_segment(data=segment_top, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.95) +
     #geom_xspline(spline_shape = -0.5, size=5) +
     #geom_smooth(method = "loess", formula=y~poly(x), se=F) +
     #stat_smooth(method = "loess", geom="area", se=F, formula=y~poly(x)) +
     #xlim(xlim[1],xlim[2]*22/20) +
     scale_fill_manual(values=color_set) +
     scale_color_manual(values=color_set) +
     theme_classic() +
     labs(x="", y="") +
     theme(
	       text=element_text(family="serif"),
	       axis.ticks = element_blank(),
	       axis.text = element_blank(),
	       axis.title.y = element_blank(),
	       axis.line = element_blank(),
	       panel.grid = element_blank(),
	       legend.position = "none"
	      )
  ggsave(p, file=paste0(savepath, "test.svg"), width=12, height=7)
  #################################################################################################################################3
  #################################################################################################################################3
  #################################################################################################################################3
  ## card 
 library(ggalt)
 library(tidyverse)
 library(ggsci)
 library(grid)
 setwd("/media/wizard/back/thermo_mzML_0518")
 savepath="/home/wizard/operation/back/"
 dataname="EU-Pro2.mzML"
 com_data <- readMSData(file = dataname, mode = "onDisk")
   ####### ms1 in peaks
   data <- filterRt(com_data, c(600, 780))
   mzrange <- c(380, 390)  # build the eic model
   ex_data <- chromatogram(data, msLevel = 1L, mz = mzrange, aggregationFun = "max")
   ex_data_1 <- ex_data[1,1]
   rt=data.frame(rtime(ex_data_1))
   int=data.frame(intensity(ex_data_1))
   rt_int=cbind(rt,int)
   colnames(rt_int)=c("rt", "int")
   rt_int$rt=rt_int$rt/60
   ########
   n=3
   time=vector(mode="list", length=n)
   time[[1]]=c(10,11)
   time[[2]]=c(11,12)
   time[[3]]=c(12,13)
   peak=vector(mode="list", length=n)
   scan=vector(mode="list", length=n)
   for(i in 1:n){
     peak[[i]]=rt_int[which(rt_int$rt > time[[i]][1] & rt_int$rt < time[[i]][2]),] 
     scan[[i]]=rt_int[which(rt_int$int==max(peak[[i]]$int)),] }
   # rownames(scan1)
   ####### ms1 in scans
   xlim=c(10.6,12.7)
   data=raw_data=rt_int
   data=data[which(data$rt>=xlim[1] & data$rt<=xlim[2]),]
   data=data.frame(spline(data$rt, data$int, n=10000))
  ### beauty
  peak_to=0.15
  scan_tot=scan[[1]]
  for(i in 2:n){scan_tot <- rbind(scan_tot, scan[[i]])}
  for(i in 1:n){
   	scan_tot <- rbind(scan_tot, 
   	 	 	   c(scan_tot[i,1]-peak_to*4/5, scan_tot[i,2]*1/5),
   	 	 	   c(scan_tot[i,1]-peak_to, scan_tot[i,2]*1/10),
   	 	 	   c(scan_tot[i,1]+peak_to*4/5, scan_tot[i,2]*1/5),
   	 	 	   c(scan_tot[i,1]+peak_to, scan_tot[i,2]*1/10) ) 
   	if(i!=n){scan_tot <- rbind(scan_tot, 
   	 	 	   c(scan_tot[i,1]+(scan_tot[i+1,1]-scan_tot[i,1])*1/3, 0), 
   	 	 	   c(scan_tot[i,1]+(scan_tot[i+1,1]-scan_tot[i,1])*1/2, 0),
   	 	 	   c(scan_tot[i,1]+(scan_tot[i+1,1]-scan_tot[i,1])*2/3, 0)
   	 	 	    ) 
   	 	}else{scan_tot <- rbind(scan_tot,
   	 	 	   c(scan_tot[i,1]+peak_to*6/4, 0))} 
   	 	}
  data=rbind(scan_tot, c(xlim[1],0), c(xlim[2],0))
  data=data.frame(spline(data$rt, data$int, n=10000))
  peak_to=peak_to *1.2
   colnames(data)=c("rt", "int")
   group_set=vector(mode="list", length=n)
   for(i in 1:n){
   group_set[[i]] <- data[which(data$rt > as.numeric(scan[[i]][1]-peak_to) & data$rt < as.numeric(scan[[i]][1]+peak_to)),] %>% 
    	 	mutate(group=i) }
   set=group_set[[1]]
   for(i in 2:n){set <- rbind(set, group_set[[i]])}
   data=merge(data, set[,c(1,3)], by="rt", all.x=TRUE, sort=TRUE)
   data[which(is.na(data$group)==T),]$group=0
   color_set <- c("0"="black", "3"="#E64B35FF", "2"="#4DBBD5FF", "1"="#00A087FF")
  ####################################################
   segment <- data.frame(c(0, 10, 0, 0),
    	 	 	 c(0, 10, 10, 10),
    	 	 	 c(0, 0, 0, 10),
    	 	 	 c(10, 10, 0, 10)
    	 	 	)
   colnames(segment) <- c("x", "xend", "y", "yend")
  ####################################################
 setwd("/home/wizard/operation/back/0703_all/results")
  p1_1 <- ggplot(data[which(data$rt > 10.68 & data$rt < 11.12), ], 
   	 	 	aes(x=rt, y=int, fill=as.character(group), color=as.character(group))) +
    geom_area() +
    scale_fill_manual(values=color_set) +
    scale_color_manual(values=color_set) +
    #geom_segment(data=segment, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.4) +
    labs(x="") +
    xlim(10.6,12.6) +
    theme_minimal() +
    theme(
       text=element_text(family="serif"),
       axis.ticks = element_blank(),
       axis.text = element_blank(),
       axis.title.y = element_blank(),
       axis.line = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none"
      ) 
  p1_2 <- ggplot(data[which(data$rt > 11.32 & data$rt < 11.8), ], 
   	 	 	aes(x=rt, y=int, fill=as.character(group), color=as.character(group))) +
    geom_area() +
    scale_fill_manual(values=color_set) +
    scale_color_manual(values=color_set) +
    #geom_segment(data=segment, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.4) +
    labs(x="") +
    xlim(10.6,12.6) +
    theme_minimal() +
    theme(
       text=element_text(family="serif"),
       axis.ticks = element_blank(),
       axis.text = element_blank(),
       axis.title.y = element_blank(),
       axis.line = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none"
      ) 
  p1_3 <- ggplot(data[which(data$rt > 12.05 & data$rt < 12.55), ], 
   	 	 	aes(x=rt, y=int, fill=as.character(group), color=as.character(group))) +
    geom_area() +
    scale_fill_manual(values=color_set) +
    scale_color_manual(values=color_set) +
    #geom_segment(data=segment, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.4) +
    labs(x="") +
    xlim(10.6,12.6) +
    theme_minimal() +
    theme(
       text=element_text(family="serif"),
       axis.ticks = element_blank(),
       axis.text = element_blank(),
       axis.title.y = element_blank(),
       axis.line = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none"
      ) 
  p2 <- ggplot(segment) +
    #geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="#BDBDBDFF", size=4, alpha=0.8) +
    annotate("rect", xmin=0, xmax=10, ymin=0 , ymax=10, color="black", fill="white", size=3, alpha=0.5) +
    theme_minimal() +
    theme(
       text=element_text(family="serif"), 
       axis.ticks = element_blank(), 
       axis.text = element_blank(), 
       axis.title.y = element_blank(), 
       axis.line = element_blank(), 
       panel.grid = element_blank(), 
       legend.position = "none"
      )
  p3 <- ggplot(segment) +
    #geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="#BDBDBDFF", size=4, alpha=0.8) +
    annotate("rect", xmin=0, xmax=10, ymin=0 , ymax=10, color="black", fill="white", size=3, alpha=0) +
    theme_minimal() +
    theme(
       text=element_text(family="serif"), 
       axis.ticks = element_blank(), 
       axis.text = element_blank(), 
       axis.title.y = element_blank(), 
       axis.line = element_blank(), 
       panel.grid = element_blank(), 
       legend.position = "none"
      )
 svg("2d_ms.svg", width=12, height=8)
 grid.newpage()
 pushViewport( viewport(layout = grid.layout(100, 100) ))
 print( p2, vp=viewport(layout.pos.row=12:92, layout.pos.col=16:98 ))
 print( p1_1, vp=viewport(layout.pos.row=48:93, layout.pos.col=16:98 ))
 print( p3, vp=viewport(layout.pos.row=12:92, layout.pos.col=16:98 ))
 ############
 print( p2, vp=viewport(layout.pos.row=8:88, layout.pos.col=14:96 ))
 print( p1_2, vp=viewport(layout.pos.row=44:89, layout.pos.col=16:98 ))
 print( p3, vp=viewport(layout.pos.row=8:88, layout.pos.col=14:96 ))
 ############
 print( p2, vp=viewport(layout.pos.row=4:84, layout.pos.col=12:94 ))
 print( p1_3, vp=viewport(layout.pos.row=24:85, layout.pos.col=12:94 ))
 print( p3, vp=viewport(layout.pos.row=4:84, layout.pos.col=12:94 ))
 ############
 print( p3, vp=viewport(layout.pos.row=4:84, layout.pos.col=12:94 ))
 # grid.text("Morphology", x=0.02, y=0.25, rot=90, gp = gpar( fontface = "bold", fontsize = 20, fontfamily = "Times", fontangle=90) )
 dev.off()
 ######################################################################################################3
 #######################################################################################################
 path="ms2_figures"
 id_1=495
 id_2=347
 id_3=2268
 for(i in 1:3){
 id=get(paste0("id_", i))
 msms <- read.csv(file=paste0(path, "/", id, ".tsv"), header=T, sep="\t")
 set <- msms[which(msms$rel.intensity > 0),]
 pal=pal_npg()(10)
 width=max(set$mz)/150
 p <- ggplot(set, aes(x=mz, y=rel.intensity)) +
   geom_col(width=width, fill=pal[i]) +
    #geom_segment(data=segment, aes(x=x, y=y, xend=xend, yend=yend, color=as.character(group)), size=2, alpha=0.4) +
    labs(x="") +
    theme_minimal() +
    theme(
       text=element_text(family="serif"),
       axis.ticks = element_blank(),
       axis.text = element_blank(),
       axis.title.y = element_blank(),
       axis.line = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none"
      )
  assign(paste0("p1_", i), p)    } 
 p2 <- ggplot(segment) +
    #geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="#BDBDBDFF", size=4, alpha=0.8) +
    annotate("rect", xmin=0, xmax=10, ymin=0 , ymax=10, color="black", fill="white", size=3, alpha=0.5) +
    theme_minimal() +
    theme(
       text=element_text(family="serif"), 
       axis.ticks = element_blank(), 
       axis.text = element_blank(), 
       axis.title.y = element_blank(), 
       axis.line = element_blank(), 
       panel.grid = element_blank(), 
       legend.position = "none"
      )
 svg("msms.svg", width=12, height=8)
 grid.newpage()
 pushViewport( viewport(layout = grid.layout(100, 100) ))
 print( p2, vp=viewport(layout.pos.row=12:92, layout.pos.col=16:98 ))
 print( p1_3, vp=viewport(layout.pos.row=48:93, layout.pos.col=16:98 ))
 print( p3, vp=viewport(layout.pos.row=12:92, layout.pos.col=16:98 ))
 ############
 print( p2, vp=viewport(layout.pos.row=8:88, layout.pos.col=14:96 ))
 print( p1_2, vp=viewport(layout.pos.row=44:89, layout.pos.col=16:98 ))
 print( p3, vp=viewport(layout.pos.row=8:88, layout.pos.col=14:96 ))
 ############
 print( p2, vp=viewport(layout.pos.row=4:84, layout.pos.col=12:94 ))
 print( p1_1, vp=viewport(layout.pos.row=24:85, layout.pos.col=12:94 ))
 print( p3, vp=viewport(layout.pos.row=4:84, layout.pos.col=12:94 ))
 ############
 print( p3, vp=viewport(layout.pos.row=4:84, layout.pos.col=12:94 ))
 # grid.text("Morphology", x=0.02, y=0.25, rot=90, gp = gpar( fontface = "bold", fontsize = 20, fontfamily = "Times", fontangle=90) )
 dev.off()
 ######################################################################################################3
 #######################################################################################################
 path="ms2_figures_label"
 id_1=495
 id_2=347
 id_3=2268
 for(i in 1:3){
 id=get(paste0("id_", i))
 source <- read.csv(file=paste0(path, "/", id, ".tsv"), header=T, sep="\t")
 p <- ggplot(source, aes(x=mz, y=rel.intensity)) + 
  geom_bar(stat="identity",width=max(source$mz)/150,fill=ifelse(source$rel.intensity>0,"black","red")) + 
  geom_point(size=1.3,color=ifelse(source$rel.intensity>0,"black","red"),alpha=ifelse(source$match==0,0,1)) + 
  #xlim(0,max(source$mz, na.rm=T)) +
  theme(text=element_text(family="serif"),
    panel.background=element_rect(fill="transparent",color="white"),
    panel.grid=element_line(color="grey85"),
    panel.border = element_rect(fill=NA, color="black", size=3, linetype="solid"), 
       axis.ticks = element_blank(), 
       axis.text = element_blank(), 
       axis.title = element_blank(), 
       axis.line = element_blank() )
    #plot.margin = unit(c(3, 1, 3, 1), "cm"))
  assign(paste0("p1_", i), p)  
  ggsave(p, file=paste0("ms2_", id, ".svg"), width=6, height=4) }
 ######################################################################################################3
 #######################################################################################################
 library(ggforce)
 library(ggalt)
 library(tidyverse)
 library(ggsci)
 library(grid)
 p <- ggplot() +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = 7, b = 3, angle = 0), size=5, color="white", fill="#8491B4FF") +
  coord_fixed() +
  theme_minimal() +
  theme(text=element_text(family="serif"),
    #panel.background=element_rect(fill="transparent",color="white"),
    panel.grid=element_blank(),
    plot.margin =unit(c(0,0,0,0),"cm"),
    panel.spacing =unit(c(0,0,0,0),"cm"),
       axis.ticks = element_blank(), 
       axis.text = element_blank(), 
       axis.title = element_blank(), 
       axis.line = element_blank() )
 svg("database.svg", width=5, height=7)
 grid.newpage()
 pushViewport( viewport(layout = grid.layout(100, 100) ))
 print( p, vp=viewport(layout.pos.row=60:90, layout.pos.col=5:95 ) )
 print( p, vp=viewport(layout.pos.row=50:80, layout.pos.col=5:95 ) )
 print( p, vp=viewport(layout.pos.row=40:70, layout.pos.col=5:95 ) )
 print( p, vp=viewport(layout.pos.row=30:60, layout.pos.col=5:95 ) )
 print( p, vp=viewport(layout.pos.row=20:50, layout.pos.col=5:95 ) )
 dev.off()
 ######################################################################################################3
 #######################################################################################################
 library(ggforce)
 library(ggalt)
 library(tidyverse)
 library(ggsci)
 library(grid)
 data="/media/wizard/back/0703_all/490_initial_8_neg_495/fingerprints/C17H24O10_[M-H]-.fpt"
 fp <- read.csv(file=data, header=T, sep="\t")
 fp <- cbind(seq(nrow(fp)),fp)
 colnames(fp) <- c("num", "pp")
 fp$group <- "g1"
 eg_fp <- fp[1:15,]
 p <- ggplot(eg_fp, aes(x=as.factor(num), y=group, fill=pp)) +
   geom_tile(color="black", size=3) +
   scale_fill_gradient(low="white", high="#3C5488FF") +
   coord_fixed() +
    theme_minimal() +
    theme(text=element_text(family="serif"),
    panel.grid=element_blank(),
    plot.margin =unit(c(0,0,0,0),"cm"),
    panel.spacing =unit(c(0,0,0,0),"cm"),
    legend.position = "none", 
       axis.ticks = element_blank(), 
       axis.text = element_blank(), 
       axis.title = element_blank(), 
       axis.line = element_blank() )
 ggsave(p, file="495.fp.svg", width=8, height=1)
 ######################################################################################################3
 #######################################################################################################
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
 library(sunburstR)
 library(RColorBrewer)
 library(shiny)
 metadata_path="../canopus_neg.tsv"
 metadata <- read.csv(file=metadata_path, header=T, sep="\t", quote = "")
 ##### see in network
 edges <- metadata[, colnames(metadata) %in% c("id", "parentId")]
 network <- as_tbl_graph(edges) %>%
 	 activate(nodes) %>% #as_tibble()
  	 mutate(deg = centrality_degree(mode='in')) 
 layout=ifelse(nrow(metadata)>1000, "unrooted", "fr")
 layout_n <- create_layout(network, layout = layout)
 #############################################################
 #############################################################
 p <- 
 ggraph(layout_n) + 
 geom_edge_fan(edge_width=0.5, color="lightblue", show.legend=F) + 
 #geom_node_point(aes(fill=str_wrap(classification, width=25), size=similarity), shape=21) + 
 geom_node_point(shape=21, aes(fill=deg, size=deg)) + 
 scale_edge_width(range=c(0.1,0.7)) + 
 scale_fill_gradient(low="#1B1919FF", high="#DC0000FF") +
 guides(alpha="none") +
 #labs(fill="Centrality\ndegree", size="Tanimoto\nsimilarity") +
 theme_grey() +
 theme(
       text=element_text(family="serif", size=20),
       axis.ticks = element_blank(),
       axis.text = element_blank(),
       axis.title = element_blank(),
       legend.title= element_text(face="bold", size=15),
       #panel.background = element_rect(fill="white"), 
       legend.key.height = unit(0.3, "cm"),
       legend.key.width = unit(0.5, "cm"),
       axis.line = element_blank(),
       #panel.grid = element_blank(),
       #legend.position = "none",
       #strip.text = element_text(size=15, face="bold")
       plot.margin =unit(c(0,0,0,0),"cm"),
       panel.spacing =unit(c(0,0,0,0),"cm")
      )
 ggsave(p, file="class_network.svg", width=7, height=6)
 #############################################################
 #############################################################
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
 palette= unique(c(pal1, pal2, pal3, pal4, pal5, pal6, pal7, pal8, pal9))
 getPalette = colorRampPalette(pal0)
 #####
 # data=data[which(data$value > 0.01),]
 # data$color <- seq(nrow(data))
 # plot <- plot_ly( data=data, 
 # labels = ~id,
 # parents = ~parentid,
 # values = ~value,
 # color = ~color,
 # colors = c(palette),
  #colors = colorRamp(c("red", "blue")), 
  #colors = c(`4` = "red", `5` = "black", `6` = "blue", `8` = "green"), 
 # type = "sunburst", 
 # insidetextorientation="radial"
 # ) %>%
 # layout( margin = list(l = 0, r = 0, b = 0, t = 0))
 # https://stackoverflow.com/questions/12926779/how-to-make-a-sunburst-plot-in-r-or-pythonlibrary(sunburstR)sequences sunburst(sequences)
 metadata_parent_path="../canopus_parent_index.tsv"
 metadata_parent <- read.csv(file=metadata_parent_path, header=T, sep="\t", quote = "")
 data_path="/media/wizard/back/0703_all/490_initial_8_neg_495/canopus/C17H24O10_[M-H]-.fpt"
 data <- read.csv(file=data_path, header=F, sep="\t")
 guide_data <- cbind(metadata_parent[,c(1)], data)
 dataset <- cbind(metadata_parent[,c(2)], data)
 colnames(dataset) <- c("id", "value")
 dataset_s <- dataset[which(dataset$value>0.01), ]
 dataset_s <- dataset_s[order(dataset_s$id),]
 #pal<-palette[1:nrow(dataset_s)]
 pal_do <- getPalette(nrow(dataset_s))
 #sunburst(dataset_s, colors=list(range = RColorBrewer::brewer.pal(9, "Set3")))
 sund2b(dataset_s, colors=list(range = RColorBrewer::brewer.pal(9, "Set3")))
 ######################################################################################################3
 #######################################################################################################
 metadata_parent_path="../canopus_parent_index.tsv"
 metadata_parent <- read.csv(file=metadata_parent_path, header=T, sep="\t", quote = "")
 data <- metadata_parent[, 2:3]
 #data$num <- 1
 sund2b(data, colors=palette)
# sunburst(data, colors=palette)

