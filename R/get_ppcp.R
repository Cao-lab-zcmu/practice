get_ppcp <- 
  function(
           key_id = NULL,
           dir = NULL,
           precursor_formula = "method_pick_formula_excellent",
           adduct = NULL,
           reformat = T,
           filter = T,
           filter_threshold = 0.1,
           class_index = "canopus.tsv",
           ...
           ){
    ## get dir path
    if(is.null(dir) == T & is.null(key_id) == T){
      return()
    }else if(is.null(dir) == T){
      dir <- get_dir(key_id)
    }
    ## ---------------------------------------------------
    ## aquire formula via the method
    if( precursor_formula == "method_pick_formula_excellent" ){
      meta <- method_pick_formula_excellent(dir = dir)
      precursor_formula <- meta$precursorFormula
      adduct <- meta$adduct
    }
    ## ---------------------------------------------------
    ## read ppcp data
    file <- list.files(path = paste0(.MCn.sirius, "/", dir, "/", "canopus"),
                       pattern = paste0("^", precursor_formula, "(.*)", escape_ch(adduct), "(.*)", ".fpt$"),
                       full.names = T)
    ppcp <- read_fpt(file)
    ## ---------------------------------------------------
    ## reformat section
    if(reformat == F){
      return(ppcp)
    }
    ## check meta list
    if(exists(".MCn.class_tree_list") == F){
      build_classes_tree_list(class_index = class_index)
    }
    ## get the environment name for lapply function to invoke data
    assign(paste0("envir_", key_id), environment(), envir = parent.env(environment()))
    ## merge with meta table, and filter
    ppcp <- lapply(.MCn.class_tree_list, merge_class_ppcp,
                   ## parameter
                   key_id = key_id, filter = filter, filter_threshold = filter_threshold)
    return(ppcp)
  }
## a small function to get data of ppcp
read_fpt <- function(file){
  fpt = data.table::fread(input = file, header = F, quote = "")
  fpt$relativeIndex = seq(0, nrow(fpt) - 1)
  return(fpt)
}
## specific character in adduct description need to be revise, for pattern matching
escape_ch <- function(x){
  x <- gsub("\\[", "\\\\\\[", x)
  x <- gsub("\\]", "\\\\\\]", x)
  x <- gsub("\\+", "\\\\\\+", x)
  x <- gsub("\\-", "\\\\\\-", x)
  x <- gsub(" ", "", x)
  return(x)
}
## the function to merge raw ppcp with meta list
merge_class_ppcp <-
  function(
           class,
           filter = T,
           filter_threshold = 0.1,
           key_id = NULL,
           values = get("ppcp", envir = get(paste0("envir_", key_id))),
           filter_col = "V1"
           ){
    df <- merge(class, values, all.x = T, by = "relativeIndex", sort = F)
    df <- df[which(df[[filter_col]] > ifelse(filter == T, filter_threshold, 0)),] %>%
      dplyr::as_tibble()
    return(df)
  }
