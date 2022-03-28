 library(xcms)
 setwd("/media/wizard/back/nanjing_sample/fecal/pos_mzml")
 metadata <- read.csv(file="EIC_metadata.tsv",header=T,sep="\t")
 dda_file=list.files(path = ".", pattern = "*.mzML$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE)
 dir.create("EIC")
 tolerance=0.005
 for(filename in dda_file){
 dda_data <- readMSData(filename, mode = "onDisk")
 dir.create(paste0("EIC/EIC_",filename))
   for(number in 1:nrow(metadata)){
   id <- metadata[number,colnames(metadata) %in% c("id")]
   mz <- metadata[number,colnames(metadata) %in% c("m.z")]
   mzrange <- c(as.numeric(mz)-tolerance,as.numeric(mz)+tolerance)
   ex_data <- chromatogram(dda_data, msLevel = 1L, mz = mzrange, aggregationFun = "max")
   ex_data_1 <- ex_data[1,1]
   if(number==1){
   	write.table(rtime(ex_data_1),paste0("EIC/EIC_",filename,"/rt",".tsv"),col.names = FALSE,sep="\t")}
   write.table(intensity(ex_data_1),paste0("EIC/EIC_",filename,"/",id,"_intensity",".tsv"),col.names = FALSE,sep="\t")
   print(paste0(filename," >>> ",number,"/",nrow(metadata)))}}
