list_merge_df <- 
  function(
           list,
           df,
           ...
           ){
    assign("envir_meta", environment(), parent.env(environment()))
    list <- lapply(list, merge,
                   y = get("df", envir = get("envir_meta")),
                   ...)
    return(list)
  }
