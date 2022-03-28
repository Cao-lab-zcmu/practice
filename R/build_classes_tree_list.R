build_classes_tree_list <-
  function(
           class_index="canopus.tsv", path=.MCn.sirius
           ){
    data <- read_tsv(paste0(path, "/", class_index))
    ## ---------------------------------------------------------------------- 
    ## separate each levels of classes into sub-list
    root <- data[which(data$parentId==""), ]
    list <- list()
    n = 1
    list[[n]] <- root %>% dplyr::as_tibble()
    df <- data[data$parentId %in% root$id, ]
    ## ---------------------------------------------------------------------- 
    while(nrow(df) > 0){
      n = n + 1
      list[[n]] <- df %>% dplyr::as_tibble()
      df <- data[data$parentId %in% df$id, ]
    }
    .MCn.class_tree_list <<- list
    cat("INFO: Classification Index in.MCn.sirius project --->", class_index, "\nA total of 11 levels. These classes (upon ClassyFire and CANOPUS) were separated into sub-lists.
        Use following arguments to get some specific classes:
        .MCn.class_tree_list[[3]] >>> superclass
        .MCn.class_tree_list[[4]] >>> class
        .MCn.class_tree_list[[5]] >>> subclass
        .MCn.class_tree_list[[6]] >>> level 5 \n")
}
