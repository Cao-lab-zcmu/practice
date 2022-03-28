## library(igraph)
## library(ggraph)
## library(tidygraph)
## library(ggsci)
## library(scales)
library(tidyverse)
library(grid)
library(stringr)
library(ggimage)
library(grImport2)
library(gridSVG)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(gt)
library(Hmisc)
library(stringr)
## obabel, cairosvg, molconvert, these should be installed previously.
parse_structure <- function(smiles, id_n){
  system(paste0("molconvert mol \"", smiles, "\" -o Rtemp.mol"))
  system("obabel Rtemp.mol -imol -osvg -O Rtemp.svg")
  system("cairosvg Rtemp.svg -o Rtemp.cairo.svg")
  assign(id_n, readPicture("Rtemp.cairo.svg"))
  return(get(id_n))
}
## grid 10(n) structure candidate
top <- function(smiles_list, save_id_number=0){
  if(file.exists("structure_check")){
    cat("Dir exist\n")
  }else{
    dir.create("structure_check")
  }
  if(file.exists(paste0("structure_check/", save_id_number, ".svg"))){
    system(paste0("xdg-open ", "structure_check/", save_id_number, ".svg"))
    return(cat("Image exist >>>", paste0("structure_check/", save_id_number, ".svg") ,"\n"))
  }
  ## draw structure in grid
  structure_num=0
  for(i in smiles_list){
    structure_num=structure_num+1
    assign(paste0("project_structure_", structure_num), parse_structure(i, paste0("medium_p", structure_num)))
  }
  svg(paste0("structure_check/", save_id_number, ".svg"), height=2, width=2*structure_num)
  for(i in 1:structure_num){
    grid.picture( get(paste0("project_structure_", i)) , x=i/(structure_num+1), y=0.5, height=1, width=1)
  }
  dev.off()
  system(paste0("xdg-open ", "structure_check/", save_id_number, ".svg"))
  cat("Save Image >>>", paste0("structure_check/", save_id_number, ".svg"), "\n")
}
## html table get
tab <- function(data_branch, list_col=c("id", "similarity", "name", "formula")){
  t <- gt(data_branch[, colnames(data_branch) %in% list_col]) %>%
  opt_table_font(
    font=list(google_font(name="Times"))
    ) %>%
  tab_header(
    title = md("**Features structure list**")
    #subtitle = md("`gtcars` is an R dataset")
  ) %>%
  opt_align_table_header(align = "left") %>%
  opt_table_lines(extent = c("none")) %>%
  cols_align(
    align = "left",
    columns = everything()
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
  return(t)
}
#####################
#####################
database <- read.csv(file="fingerid_candidate_top10.tsv", header=T, sep="\t")
## catch and draw structure in db
smi <- function(id, n=id){
  smiles <- database[which(database$"id"==id), colnames(database) %in% "smiles"] 
  smiles <- top(smiles,n)
  return(smiles)
}
## smi_batch_save
smi_batch <- function(id, n=id){
  smiles <- database[which(database$"id"==id), colnames(database) %in% "smiles"]
  smiles <- top(smiles,n)
}
## catch and pretty table in db
sht <- function(id){
  table <- database[which(database$"id"==id),] %>%
    tab()
  return(table)
}
## replace the candidate structure annotation
rp <- function(data_pip, rp_id, rank_n,
               col=c("similarity",
                     "name",
                     "formula",
                     "xlogp",
                     "smiles",
                     "inchi",
                     "inchikey2D",
                     "links")){
  origin_order=colnames(data_pip)
  data_pip<-data_pip[,order(colnames(data_pip))]
  database<-database[,order(colnames(database))]
  data_pip[which(data_pip$"id"==rp_id),
       colnames(data_pip) %in% col] <- database[which(database$"id"==rp_id)[rank_n], colnames(database) %in% col]
  data_pip<-data_pip[, order(factor(colnames(data_pip), levels=origin_order))]
  return(data_pip)
}
cut_line <- function(data_pip, cut_id){
  data_pip <- data_pip[which(data_pip$"id"!=cut_id), ]
  return(data_pip)
}
data_save <- function(data_pip, file_name="df"){
  time_temp <- date()
  write.table(data_pip, file=paste0(file_name, "_", time_temp, ".tsv"),
              sep="\t", col.names=T, row.names=F, quote=F)
}
workflow <- function(id_list, sleep_time=5, data_pip=data){
  cut_id_list=c()
  system("echo > cut.list")
  for(i in id_list){
    smi_batch(i)
    cut_else=scan("", nlines=1, what="c")
    if(class(try(cut_else))=="try_error"){
      cat("Try again\n")
      cut_else=scan("", nlines=1, what="c")
    }
    if(cut_else=="s"){
      cat("Skip\n")
      next
    }
    if(cut_else=="c"){
      cut_id_list=c(cut_id_list, i)
      data_pip=cut_line(data_pip, i)
      system(paste0("echo ", i, " >> cut.list"))
      cat("Cut the line of ", i, "\n")
    }else if(is.na(as.numeric(cut_else))==F){
      data_pip=rp(data_pip, i, as.numeric(cut_else))
      cat("Replace id", i, "candidate\n")
    }
  }
  return(data_pip)
}
####################
####################
datapath="com_1.tsv"
data=read.csv(file=datapath, header=T, sep="\t")

