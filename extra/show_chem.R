## ---------------------------------------------------------------------- 
## a function to fast show df
show_meta <- 
  function(
           x,
           df = meta
           ){
    prefix <- stringr::str_extract(df[1,]$".id", "^[a-z]{1,100}(?=[0-9])")
    df <- dplyr::filter(df, .id == paste0(prefix, x)) %>%
      data.table::data.table()
    return(df)
  } 
## ---------------------------------------------------------------------- 
getk <-
  function(
           x,
           col = "INCHIKEY"
           ){
    df <- show_meta(x)
    return(df[[col]])
  }
## ---------------------------------------------------------------------- 
show_stru <- 
  function(
           df,
           key = c("smiles", "SMILES")
           ){
    key <- key[key %in% colnames(df)][1]
    smiles <- df[[key]]
    molconvert_structure(smiles)
  }
## ---------------------------------------------------------------------- 
auto <- 
  function(
           id
           ){
    id <- as.character(substitute(id))
    show_meta(id) %>%
      show_stru()
  }
