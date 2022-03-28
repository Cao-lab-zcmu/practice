msp_to_mgf <-
  function(
           name,
           id_prefix,
           path = "~/Downloads/msp",
           write_meta_data = paste0(path, "/", name, ".meta.tsv")
           ){
    msp <- read_msp(paste0(path, "/", name))
    cache <- new.env()
    store <- new.env()
    assign("id", 0, envir = cache)
    mgf <- paste0(path, "/", name, ".mgf")
    assign("envir_meta", environment(), envir = parent.env(environment()))
    cat("", file = mgf)
    pbapply::pblapply(msp[[1]], deal_with_msp_record, 
                      id_prefix = id_prefix,
                      cache = cache,
                      store = store)
    set <- ls(envir = store)
    meta_data <- lapply(set, get_envir_df,
                                   envir = store)
    meta_data <- data.table::rbindlist(meta_data)
    if(is.null(write_meta_data) == F){
      write_tsv(meta_data, write_meta_data)
    }
    return(meta_data)
  }
read_msp <-
  function(
           filepath
           ){
    msp <- data.table::fread(filepath, sep = NULL, header = F)
  }
get_envir_df <-
  function(
           var,
           envir
           ){
    df <- get(var, envir = envir)
    return(df)
  }
