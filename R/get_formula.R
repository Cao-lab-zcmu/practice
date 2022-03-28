get_formula <-
  function(
           key_id,
           exclude_element = NULL, ## e.g., c("S", "B", "P", "Si")
           formula_method = "top_zodiac",
           rank = 1:5, # or "all"
           ppm_error = 20,
           return_col = c("rank", "precursorFormula", "molecularFormula",
                        "adduct", "ZodiacScore", "massErrorPrecursor(ppm)"),
           ...
           ){
    path <- list.files(path = .MCn.sirius, pattern=paste0("*_", key_id, "$"), full.names=T)
    file <- read_tsv(paste0(path, "/", "formula_candidates.tsv"))
    file$rank <- as.numeric(file$rank)
    ## ---------------------------------------------------------------------- 
    if("ZodiacScore" %in% colnames(file) == F){
      file$ZodiacScore = 0
    }
    ## ---------------------------------------------------------------------- 
    if(is.null(exclude_element) == F){
      file <- file[!unname(sapply(file$precursorFormula, grep_element,
                                  exclude_element = exclude_element)), ]
    }
    ## ---------------------------------------------------------------------- 
    if(formula_method == "top_zodiac"){
      if(rank[1] == "all"){
        rank <- unique(file$rank)
      }
      df <- file[which(file$rank %in% rank & abs(file$"massErrorPrecursor(ppm)") <= ppm_error), c(return_col)]
    }
    return(df)
}
## ---------------------------------------------------------------------- 
grep_element <-
  function(
           formula,
           exclude_element = c("S", "P", "B")
           ){
    if(length(grep(paste(exclude_element, collapse = "|"), formula)) == 1){
      return(T)
    }else{
      return(F)
    }
  }
