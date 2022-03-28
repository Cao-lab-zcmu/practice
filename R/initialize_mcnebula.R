initialize_mcnebula <-
  function(
           sirius_path,
           output_path = sirius_path,
           output_file = "mcnebula_results",
           palette = unique(c(ggsci::pal_simpsons()(16), ggsci::pal_igv("default")(51))),
           palette_stat = palette,
           palette_label = colorRampPalette(c("#C6DBEFFF", "#3182BDFF", "red"))(10),
           rm_var = F
           ){
    if(rm_var == T){
      ls(envir = .GlobalEnv, pattern = "^.MCn.(.*)", all.names = T) %>%
        rm(envir = .GlobalEnv)
    }
    if(file.exists(sirius_path)==F | file.exists(output_path)==F){
      cat("File path not find.\n")
      return()
    }
    if(file.exists(paste0(sirius_path, "/", ".format"))==F){
      cat("SIRIUS project not find.\n")
      return()
    }
    .MCn.sirius <<- sirius_path
    .MCn.output <<- output_path
    .MCn.results <<- output_file
    .MCn.palette <<- palette
    .MCn.palette_stat <<- palette_stat
    .MCn.palette_label <<- palette_label
    dir.create(paste0(.MCn.output, "/", .MCn.results))
    cat("MCnebula project has initialized at ->", .MCn.output, "\n")
  }
