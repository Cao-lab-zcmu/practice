library(xcms)
library(RColorBrewer)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
# data name list
path <- "."
dda_file=list.files(path = path, pattern = "*.mzML$", all.files = FALSE, 
                    full.names = FALSE, recursive = FALSE,
                    ignore.case = FALSE, include.dirs = FALSE)
# metadata
pd <- data.frame(sample_name = dda_file, 
                 sample_group = rep("test", length(dda_file)),
                 stringsAsFactors = FALSE)
# xcms read the data
raw_data <- readMSData(files = dda_file, pdata = new("NAnnotatedDataFrame", pd), 
                       mode = "onDisk")

