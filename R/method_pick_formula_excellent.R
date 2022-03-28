method_pick_formula_excellent <- 
  function(
           dir = NULL,
           key_id = NULL, 
           return_formula = T,
           return_structure = F,
           formula_cache = NULL,
           structure_cache = NULL,
           exclude_element = NULL,
           ppm_error = 20,
           fc = 1.5,
           formula_info = c("precursorFormula", "molecularFormula", "adduct", "ZodiacScore")
           ){
    ## check parameter
    if(is.null(key_id) == T & is.null(dir) == T){
      return(0)
    }
    if(is.null(key_id) == F){
      dir = get_dir(key_id)
    }else{
      key_id = grep_id(dir)
    }
    ## -------------------------------------------------------------------
    ## first, test formula
    ## check whether.MCn.sirius compute formula for this id
    if(file.exists(paste0(.MCn.sirius, "/", dir, "/", "formula_candidates.tsv")) == F){
      return(0)
    }
    ## -------------------------------------------------------------------
    formula_df <- get_formula(key_id,
                              rank = "all",
                              ppm_error = ppm_error,
                              exclude_element = exclude_element)
    ## sometimes the top formula is not rank 1, but must be top
    n = formula_df[1, "rank"] 
    formula_zodiac_rank1 <- formula_df[ which(formula_df$rank == n), c(formula_info) ]
    ## check whether there are fingerid structure files exist.
    structure_df = ""
    if(file.exists(paste0(.MCn.sirius, "/", dir, "/", "fingerid")) == T){
      check = T
      structure_df <- try(get_structure(key_id, return_row= "all"), silent = T)
    }else{
      check = F
    }
    if(class(structure_df)[1] != "try-error" & check == T){
      ## In sirius workflow, if some enforce setting (e.g., adduct enfoce) have done, max mass error (ppm) may too large. here to filter.
      structure_df <- structure_df[structure_df$molecularFormula %in% formula_df$molecularFormula, ]
      ## get formula and adduct type (select from zodiac top score formula or formula of top score structure)
      if(nrow(structure_df) > 0){
        ## get top score (of structure) formula
        top_struc_formula <- structure_df[1,]$molecularFormula
        ## ------------------- get score
        formula_structure_rank1 <-
          formula_df[ which(formula_df$molecularFormula == top_struc_formula), c(formula_info)]
        score_rank1_zodiac <- formula_zodiac_rank1[1,]$ZodiacScore %>% # get score
          as.numeric()
        score_rank1_structure <- formula_structure_rank1[1,]$ZodiacScore %>% # get score
          as.numeric()
        ## ------------------- comparation
        if( score_rank1_zodiac >= (score_rank1_structure * fc) ){
          use_zodiac = T
          ## sometimes rank 1 zodiac formulae are plural, due to complex adduct type
          ## hence based on structure score to filter them
          structure_df <- structure_df[structure_df$molecularFormula %in% formula_zodiac_rank1$molecularFormula, ]
        }else{
          use_zodiac = F
        }
        ## ------------------- aquisition
        if( nrow(structure_df) > 0 ){
          ## in rank 1 zodiac formula, top score (structure) formula were picked
          structure_pick <- structure_df[1,]
          formula_adduct <- 
            formula_df[which(formula_df$molecularFormula == structure_pick$molecularFormula), c(formula_info)]
        }else{
          structure_pick <- data.frame()
          formula_adduct <- formula_zodiac_rank1[1,]
        }
      }else{ # nrow(structure_df) == 0, no results
        use_zodiac = T
        structure_pick <- data.frame()
        formula_adduct <- formula_zodiac_rank1[1,]
      }
    }else{ # try-error or check == F, no such files
      use_zodiac = T
      structure_pick <- data.frame()
      formula_adduct <- formula_zodiac_rank1[1,]
    }
    ## -------------------------------------------------------------------
    ## add annotation col
    formula_adduct <- dplyr::mutate(formula_adduct, use_zodiac = use_zodiac)
    ## -------------------------------------------------------------------
    if(is.null(formula_cache) == F & is.null(structure_cache) == F){
      assign(paste0(key_id), formula_adduct, envir = formula_cache)
      assign(paste0(key_id), structure_pick, envir = structure_cache)
    }
    ## return data
    if(return_formula == T){
      return(formula_adduct)
    }else if(return_structure == T){
      return(structure_pick)
    }
  }
## get dir from key_id
get_dir <- function(key_id, path = .MCn.sirius){
  dir <- list.files(path = path, pattern=paste0("^[0-9](.*)_(.*)_", key_id, "$"), full.names = F)
  check <- check_dir(dir)
  if(check == T){
    return(dir)
  }
}

