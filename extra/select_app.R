select_app <- 
  function(
           df,
           col
           ){
    df <- dplyr::select(df, all_of(col))
    return(df)
  }
