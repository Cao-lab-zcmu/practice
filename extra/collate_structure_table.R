collate_structure_table <- 
  function(
           file = "method_pick_formula_excellent.structure.tsv",
           class = "../canopus_summary.tsv",
           cut_tanimoto = 0.4,
           delete_null = T,
           ...
           ){
    order = c(".id", "name", "molecularFormula", "adduct",
              "most specific class", "inchikey2D", "tanimotoSimilarity")
    df <- read_tsv(file)
    df <- dplyr::filter(df, tanimotoSimilarity >= cut_tanimoto) %>%
      dplyr::select(.id, name, tanimotoSimilarity, molecularFormula, inchikey2D) %>%
      dplyr::mutate(.id = as.character(.id))
    class <- read_tsv(class) %>%
      dplyr::select(name, `most specific class`, adduct) %>%
      dplyr::mutate(.id = stringr::str_extract(name, "(?<=_)[0-9]{1,4}$")) %>%
      dplyr::select(.id, `most specific class`, adduct)
    df <- merge(df, class, by = ".id", all.x = T, sort = T) %>%
      dplyr::select(all_of(order)) %>%
      dplyr::distinct(name, .keep_all = T)
    if(delete_null == T)
      df <- dplyr::filter(df, name != "null")
    write_tsv(df, "table_of_pdf.tsv")
    return(df)
  }
