nebula_re_rank <-
  function(
           nebula_name,
           top_n = 10,
           match_pattern = c("precursorFormula"), ## or c("precursorFormula", "adduct")
           collate_factor = 0.85,
           cluster_method = "method_rerank_binning_cluster",
           revise_MCn_formula_set = T,
           revise_MCn_structure_set = T,
           ...
           ){
    cat("[INFO] MCnebula run: nebula_re_rank\n")
    ## get formula
    id_set <- dplyr::filter(.MCn.nebula_index, name == nebula_name)
    formula_adduct <- dplyr::filter(.MCn.formula_set, .id %in% id_set$".id")
    ## ---------------------------------------------------------------------- 
    ## match patern
    if("precursorFormula" %in% match_pattern == F){
      formula_adduct$precursorFormula = NULL
    }
    if("adduct" %in% match_pattern == F){
      formula_adduct$adduct = NULL
    }
    ## ---------------------------------------------------------------------- 
    ## catch file
    ## due to the return data type of mapply, lapply is used instead
    ## first, as list
    assign("envir_nebula", environment(), envir = parent.env(environment()))
    formula_adduct <- lapply(formula_adduct$".id", by_group_for_list,
                             col = ".id",
                             df = get("formula_adduct", envir = get("envir_nebula")))
    ## then, use lapply match file
    cat("## netbula_re_rank: get_structure\n")
    structure_set <- pbapply::pblapply(formula_adduct, df_get_structure,
                                       top_n = top_n,
                                       collate_factor = collate_factor,
                                       ...)
    structure_set <- data.table::rbindlist(structure_set, fill = T)
    cat("## STAT of structure_set:",
        paste0(nrow(structure_set), " (structure sum)/", length(unique(structure_set$".id")), "(.id sum)"), "\n")
    ## ---------------------------------------------------------------------- 
    ## get sdfset, and further get apset via ChemmineR
    cat("## convert data: SMILES_set -> SDF_set -> AP_set\n")
    sdfset <- smiles_to_sdfset(structure_set)
    apset <- sdf2ap(sdfset)
    ## ---------------------------------------------------------------------- 
    ## cluster method
    method_fun <- match.fun(cluster_method)
    meta_rank <- method_fun(apset, ...)
    print(meta_rank)
    ## ---------------------------------------------------------------------- 
    structure_set <- merge(structure_set,
                           meta_rank[, c(".id", "structure_rank", "size")],
                           by = c(".id", "structure_rank"), all.x = T) %>%
      dplyr::filter(is.na(size) == F) %>%
      dplyr::as_tibble()
    ## ---------------------------------------------------------------------- 
    ## revise .GlobalVar .MCn.formula_set
    if(revise_MCn_formula_set == T){
      ## prepare replace data
      rp <- dplyr::arrange(structure_set, .id) %>%
        tidyr::separate(col = "file_name", sep = "_", into = c("precursorFormula", "adduct")) %>%
        dplyr::mutate(adduct = gsub("\\+(?!$)", " \\+ ", adduct, perl = T),
                      adduct = gsub("\\-(?!$)", " \\- ", adduct, perl = T)) %>%
        dplyr::select(.id, precursorFormula, adduct, molecularFormula)
      ## replace
      fset <- dplyr::arrange(.MCn.formula_set, .id)
      fset[fset$".id" %in% rp$".id", c(".id", "precursorFormula", "adduct", "molecularFormula")] <- rp
      .MCn.formula_set <<- fset
    }
    ## ---------------------------------------------------------------------- 
    ## revise .GlobalVar .MCn.structure_set -------
    if(revise_MCn_structure_set == T){
      sset <- dplyr::arrange(.MCn.structure_set, .id)
      ## prepare replace data
      rp <- dplyr::arrange(structure_set, .id) %>%
        dplyr::select(colnames(sset))
      ## replace
      sset <- dplyr::distinct(rbind(rp, sset), .id, .keep_all = T)
      .MCn.structure_set <<- sset
      ## rename exist structure picture -------
      tmp_stru <- paste0(.MCn.output, "/", .MCn.results, "/tmp/structure")
      if(file.exists(tmp_stru) == T){
        lapply(paste0(tmp_stru, "/", rp$".id", ".svg"), rename_file)
      }
    }
    ## ---------------------------------------------------------------------- 
    cat("[INFO] MCnebula Job Done: nebula_re_rank\n")
    return(structure_set)
  }
smiles_to_sdfset <-
  function(
           structure_set
           ){
    ##
    smiles_set <- structure_set$smiles
    names(smiles_set) <- paste0(structure_set$".id", "_", structure_set$structure_rank)
    ## this function automaticly set the vector name as name of each subset
    sdf_set <- ChemmineR::smiles2sdf(smiles_set)
    return(sdf_set)
  }
df_get_structure <-
  function(
           x,
           top_n = 10,
           collate_factor = 0.85,
           ...
           ){
    df <- get_structure(
                        x[[".id"]],
                        x[["precursorFormula"]],
                        x[["adduct"]],
                        return_row = 1:top_n,
                        ...)
    if(nrow(df) == 0){
      return(df)
    }
    df <- dplyr::mutate(df, .id = x[[".id"]]) ## add key_id
    top_simi <- df[1, "tanimotoSimilarity"]
    df <- dplyr::filter(df, tanimotoSimilarity >= top_simi * collate_factor)
    return(df)
  }
rename_file <-
  function(
           file,
           suffix = "prefix"
           ){
    if(file.exists(file) == T){
      file.rename(file, paste0(file, ".", suffix))
    }
  }
