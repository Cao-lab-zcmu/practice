## pretty table
library(tidyverse)
library(reshape2)
library(gt)
library(Hmisc)
library(stringr)
data_save <- function(data_pip, file_name="df"){
  time_temp <- date()
  write.table(data_pip, file=paste0(file_name, "_", time_temp, ".tsv"),
              sep="\t", col.names=T, row.names=F, quote=F)
}

# datapath0="com_compound.tsv"
datapath0="merge.tsv"
#datapath1="indraw_lignans_and_iridoids.tsv"
datapath1="merge_lignans_and_iridoids.tsv"
#datapath11="iridoids_Mon Jan  3 19:56:31 2022.tsv"
#datapath1="com_carboxylic_acids.tsv"
datapath2="stat_classification.tsv"
# dataset <- read.csv(file=datapath0, header=T, sep="\t")
# dataset <- dataset[which(dataset$m.z > 193 & dataset$m.z <194),]
# write.table(dataset[,colnames(dataset) %in% c("id", "smiles")], file="test.tsv", sep="\t", col.names=T, row.names=F)
data0 <- read.csv(file=datapath0, header=T, sep="\t")
data0 <- data0[order(data0$id),]
data1 <- read.csv(file=datapath1, header=T, sep="\t")
data1 <- data1[order(data1$id),]
origin_order=colnames(data1)
data1 <- merge(data1, data0[,colnames(data0) %in% c("id", "IUPACName")], by="id", all.x=T, sort=T)
data1 <- mutate(data1, extra_name=ifelse(
                                         name!="null", name,
                                         ifelse(
                                                is.na(IUPACName)==T, "null", IUPACName
                                         )
)
)
data1 <- data1[, !(colnames(data1) %in% c("name", "IUPACName"))]
colnames(data1)[which(colnames(data1)=="extra_name")]="name"
data1 <- data1[, order(factor(colnames(data1), levels=origin_order))]
data_save(data1, "chem_supple")
#################################################################
#data11 <- read.csv(file=datapath11, header=T, sep="\t")
#data1 <- rbind(data1, data11)
data2 <- read.csv(file=datapath2, header=T, sep="\t")
data <- merge(data1, data2, by="id", all.x=T, sort=T)
data_branch <- data[which(data$similarity>0.45), ]
data_branch <- data_branch[, colnames(data_branch) %in% 
 	c("id", "rt", "m.z", "variety", "pro.raw", "similarity", "name", "formula", "inchikey2D", "Index")]
data_branch <- data_branch[,c(which(colnames(data_branch)!="Index"), which(colnames(data_branch)=="Index"))]
data_branch <- data_branch[order(data_branch$inchikey2D, -data_branch$similarity), ]
## data_branch <- data_branch[!duplicated(data_branch$inchikey2D), ] 
data_branch <- data_branch[order(data_branch$id, -data_branch$similarity), ]
data_branch <- data_branch[!duplicated(data_branch$id), ] 
colnames(data_branch) <- capitalize(colnames(data_branch))
colnames(data_branch)[c(1:6, 10)] <- c("ID", "RT", "m/z", "Processing Variations", "Pro/Raw(peak area)", "Tanimoto similarity", "Classes")
data_branch$Number <- seq(nrow(data_branch))
data_branch <- data_branch[c(ncol(data_branch), 1:(ncol(data_branch)-1))]
#data_branch <- data_branch[order(data_branch[,colnames(data_branch) %in% c("m/z")]),]
data_branch$Name <- str_wrap(data_branch$Name, width=40)
###################################################################
t <- gt(data_branch) %>% 
  opt_table_font(
    font=list(google_font(name="Times"))
    ) %>%
  tab_header(
    title = md("**Features of lignans and iridoids of _Eucommia ulmoides_ identified in LC-MS negative ion mode**")
    #subtitle = md("`gtcars` is an R dataset")
  ) %>%
  opt_align_table_header(align = "left") %>%
  tab_footnote(
               footnote = "These features (compounds) were mainly obtained by phase I clustering of MCnebula. The structures of these compounds were manually picked from the top 10 ranked CSI:FingerID structure prediction candidates based on phase I clustering similarity propagation. The compounds with Tanimoto similarity (≤ 0.45) of picked structure were filtered. Note that MS/MS spectrum is fail in unequivocally ascertaining the molecular scaffold (i.e. geometrical isomers, position isomers) and 3D structure. Therefore, the structures of these compounds may be the respective isomers.",
               locations = cells_title(
                                       groups = c("title")
               )
  ) %>%
  opt_table_lines(extent = c("none")) %>%
    cols_align(
               align = "left",
    columns = everything()
  ) %>%
  cols_width(
    Name ~ px(300)
    # ends_with("r") ~ px(100),
    # starts_with("date") ~ px(200),
    # everything() ~ px(60)
  ) %>%
  tab_style(
      style = cell_borders(
        sides = c("top", "bottom"),
        color = "black",
        weight = px(1.5),
        style = "solid"
      ),
      locations = cells_column_labels()
  ) %>%
  tab_style(
      style = cell_text(v_align="top"),
      locations = cells_column_labels(
        columns = everything()
        #rows = 1:nrow(data_branch)
    )
  ) %>%
  tab_style(
      style = cell_borders(
        sides = c("bottom"),
        color = "black",
        weight = px(1.5),
        style = "solid"
      ),
      locations = cells_body(
        columns=everything(),
        rows=nrow(data_branch)
      )
   ) %>%
  tab_style(
      style = cell_text(v_align="top"),
      locations = cells_body(
        columns = everything()
        #rows = 1:nrow(data_branch)
      )
   )
