auto_classy <- 
  function(
           df,
           ...
           ){
    ## classyfireR
    list <- by_group_as_list(df, ".id")
    pbapply::pblapply(list, base_auto_classy,
                        ...)
  }
base_auto_classy <- 
  function(
           df
           ){
    .id <- df[1,][[".id"]]
    lapply(df[["InChIKey"]], base2_classy,
                      .id = .id)
  }
base2_classy <- 
  function(
           inchi,
           .id
           ){
    ch <- try(read_tsv(paste0(.id)), silent = T)
    if(class(ch) == "try-error"){
      ch <- classyfireR::get_classification(inchi)
    }else{
      return()
    }
    if(is.null(ch)){
      return()
    }else{
      ch <- classyfireR::classification(ch)
      write_tsv(ch, paste0(.id))
    }
  }
