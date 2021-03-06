## main function
collate_ppcp <- 
  function(
           formula_adduct = .MCn.formula_set,
           path = .MCn.sirius,
           dirs = "all",
           output = paste0(.MCn.output, "/", .MCn.results),
           write_output = T,
           nebula_class = T,
           nebula_index = T,
           ...
           ){
    cat("[INFO] MCnebula run: collate_ppcp\n")
    ## ---------------------------------------------------------------------- 
    ## check dirs ---- canopus
    cat("## collate_ppcp: check_dir\n")
    if(dirs == "all"){
      dirs <- list.files(path = path, pattern="^[0-9](.*)_(.*)_(.*)$", full.names = F)
      check <- pbapply::pbsapply(dirs, check_dir, file = "canopus") %>% unname
    }else{
      check <- pbapply::pbsapply(dirs, check_dir, file = "canopus") %>% unname
    }
    ## ---------------------------------------------------------------------- 
    ## lock on file location
    meta_dir <- dirs[which(check == T)] %>%
      data.frame() %>%
      dplyr::rename(dir = ".") %>%
      dplyr::mutate(.id = sapply(dir, grep_id)) %>%
      merge(formula_adduct, by = ".id", all.x = T, sort = F) %>%
      dplyr::mutate(adduct_trans = gsub(" ", "", adduct),
             target = paste0(precursorFormula, "_", adduct_trans, ".fpt"), 
             full.name = paste0(path, "/", dir, "/", "canopus", "/", target), 
             ## these files need to be check and filter (whether exist)
             ## note that some formula is no fingerprint computed
             ppcp = file.exists(full.name))
      meta_dir_filter <- dplyr::filter(meta_dir, ppcp == T)
      cat("## STAT of PPCP dataset:",
          paste0(nrow(meta_dir_filter), "(formula with PPCP)", "/", nrow(meta_dir), "(all formula)"), 
          "\n")
    ## ---------------------------------------------------------------------- 
    ## load all ppcp dataset
    ppcp_dataset <- pbapply::pblapply(meta_dir_filter$full.name, read_fpt)
    names(ppcp_dataset) <- meta_dir_filter$".id"
    .MCn.ppcp_dataset <<- ppcp_dataset
    ## ---------------------------------------------------------------------- 
    ## summarize nebula_class
    if(nebula_class == T){
      cat("# Collate_ppcp |", date(), "| USE Method: method_summarize_nebula_class.\n")
      metadata <- data.table::rbindlist(.MCn.class_tree_list, idcol = T) %>%
        dplyr::rename(hierarchy = .id)
      .MCn.class_tree_data <<- dplyr::as_tibble(metadata)
      ## transmit environment
      assign("envir_meta", environment(), envir = parent.env(environment()))
      ## get nebula classes
      nebula_class <- pbapply::pblapply(ppcp_dataset, method_summarize_nebula_class, 
                             class_data_type = "classes_tree_data",
                             ...)
      .MCn.nebula_class <<- nebula_class
    }
    ## ---------------------------------------------------------------------- 
    if(nebula_index == T){
      cat("# Collate_ppcp |", date(), "| USE Method: method_summarize_nebula_index.\n")
      ## gather all nebula classes
      nebula_index <- method_summarize_nebula_index(ppcp_dataset,
                                                    ...)
      .MCn.nebula_index <<- nebula_index
    ## ---------------------------------------------------------------------- 
      if(write_output == T){
        write_tsv(nebula_index, file = paste0(output, "/", "nebula_index.tsv"))
      }
    }
    cat("[INFO] MCnebula Job Done: collate_ppcp.\n")
    return(nebula_index)
  }
