#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dir PARAM_DESCRIPTION, Default: NULL
#' @param key_id PARAM_DESCRIPTION, Default: NULL
#' @param return_formula PARAM_DESCRIPTION, Default: T
#' @param return_structure PARAM_DESCRIPTION, Default: F
#' @param formula_cache PARAM_DESCRIPTION, Default: NULL
#' @param structure_cache PARAM_DESCRIPTION, Default: NULL
#' @param exclude_element PARAM_DESCRIPTION, Default: NULL
#' @param ppm_error PARAM_DESCRIPTION, Default: 20
#' @param fc PARAM_DESCRIPTION, Default: 1.5
#' @param formula_info PARAM_DESCRIPTION, Default: c("precursorFormula", "molecularFormula", "adduct", "ZodiacScore")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}
#' @rdname method_pick_formula_excellent
#' @export 
#' @importFrom dplyr mutate
method_pick_formula_excellent <- 
  function(
           key_id = NULL, 
           dir = NULL,
           return_formula = T,
           return_structure = F,
           formula_cache = NULL,
           structure_cache = NULL,
           exclude_element = NULL,
           ppm_error = 20,
           fc = 1.5,
           formula_info = c("precursorFormula", "molecularFormula", "adduct", "ZodiacScore")
           ){
    ## ----------------------------------------------------------------------
    if(is.null(dir) == F){
      if(length(dir) >= 1){
        meta <- data.table::data.table(dir = dir)
        meta <- dplyr::mutate(meta, key_id = grep_id(dir))
      }
    }else{
      meta <- data.table::data.table(key_id = key_id, dir = get_dir(key_id))
    }
    ## ---------------------------------------------------------------------- 
    meta <- dplyr::mutate(meta, formula_dir = paste0(.MCn.sirius, "/", dir, "/", "formula_candidates.tsv"),
                          fingerid = paste0(.MCn.sirius, "/", dir, "/", "fingerid"))
    ## ------------------------------------- 
    cat("## method part: check_dir: formula set\n")
    check <- pbapply::pblapply(meta$formula_dir, file.exists)
    meta$check_fm <- check
    ## ------------------------------------- 
    cat("## method part: check_dir: fingerid\n")
    check <- pbapply::pblapply(meta$fingerid, file.exists)
    meta$check_fd <- check
    ## ------------------------------------- 
  }
