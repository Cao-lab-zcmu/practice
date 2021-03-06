annotate_child_nebulae <-
  function(
           nebula_name,
           compound_class_list = .MCn.nebula_class,
           write_output = T,
           output = paste0(.MCn.output, "/", .MCn.results),
           layout = "fr",
           height = "auto",
           width = "auto",
           plot_nodes_id = T,
           plot_structure = T,
           plot_ppcp = T,
           ratio_df = NULL,
           merge_image = T,
           return_plot = F,
           ...
           ){
    cat("[INFO] MCnebula run: annotate_child_nebulae\n")
    ## ------------------------------------------------------------------------
    ## all nodes in graph
    nodes <- dplyr::filter(.MCn.nebula_index, name == nebula_name)$".id"
    ## get top compound class (nodes_color data)
    ## as well as, collate metadata
    metadata <- lapply(compound_class_list, head, n = 1) %>%
      data.table::rbindlist(idcol = T) %>% # as data.frame
      dplyr::filter(.id %in% nodes) %>% # filter via nodes
      dplyr::select(.id, name) %>%
      dplyr::rename(vis_class = name)
    ## push environment name into parent.env, let some data could be catch in sub-environment via 'get' function
    assign("envir_meta", environment(), envir = parent.env(environment()))
    ## gather data for annotation (nebula_name, hierarchy)
    hierarchy <- dplyr::filter(.MCn.nebula_index, name == nebula_name) %>%
      head(n = 1)
    anno = c(nebula_index = nebula_name, hierarchy = hierarchy$hierarchy)
    ## ------------------------------------------------------------------------
    ## set a environment to store layout data
    envir_layout <- new.env() 
    ## set to remove nodes or not (set to 0, remove)
    if(plot_ppcp == T | plot_structure == T){
      remove_nodes = 0
    }else{
      remove_nodes = NULL
    }
    ## plot origin network (child network, with legend)
    p <- grid_child_nebula(
                           .MCn.child_graph_list[[nebula_name]],
                           anno = anno,
                           print_into = F,
                           layout = layout,
                           ## save layout data in this environment
                           save_layout_df = envir_layout,
                           ## remove origin nodes
                           remove_nodes = remove_nodes, 
                           ...)
    ## ---------------------------------------------------------------------- 
    ## whether plot pie diagram
    if(is.null(ratio_df) == F){
      if(is.data.frame(ratio_df) == T){
        plot_ratio = T
      }else{
        cat("is.data.frame(ratio_df) == F\n")
        plot_ratio = F
      }
    }else{
      plot_ratio = F
    }
    ## ------------------------------------------------------------------------
    ## tmp dir
    tmp_dir <- paste0(output, "/", "tmp")
    if(file.exists(tmp_dir) == F){
      dir.create(tmp_dir)
    }
    ## add annotation ---------------------------------------------------------
    ## nodes id
    if(plot_nodes_id == T & (plot_ppcp == F)){
      p <- p + ggraph::geom_node_text(aes(label = name), size = 1)
    }
    ## add annotation ---------------------------------------------------------
    ## require ChemmineOB and ChemmineR
    with_structure <- 0
    if(requireNamespace("ChemmineOB", quietly = T)){
      ## structure
      tmp_stru <- paste0(tmp_dir, "/", "structure")
      if(file.exists(tmp_stru) == F){
        dir.create(tmp_stru)
      }
      if(plot_structure == T){
        with_structure <- 1
        batch_mode_structure(metadata = metadata, tmp_stru = tmp_stru)
      }
    }
    ## add annotation ---------------------------------------------------------
    ## re draw nodes with or without ppcp bar
    tmp_ppcp <- paste0(tmp_dir, "/", "ppcp")
    if(file.exists(tmp_ppcp) == F){
      dir.create(tmp_ppcp)
    }
    if(plot_ppcp == T | plot_structure == T | plot_ratio == T){
      batch_mode_nodes(
                       metadata = metadata,
                       tmp_ppcp = tmp_ppcp,
                       with_structure = with_structure,
                       plot_ppcp = plot_ppcp,
                       plot_ratio = plot_ratio,
                       ratio_df = ratio_df,
                       ...)
    }
    ## ------------------------------------------------------------------------
    ## merge image
    if(merge_image == T){
      if(requireNamespace("ggimage", quietly = T) &
         requireNamespace("gridExtra", quietly = T)){
        ## remove legend of size
        p <- p + ggplot2::guides(size = "none")
        merge_image(p, envir_layout$layout_n, tmp_ppcp)
      }
    }
    ## ------------------------------------------------------------------------
    ## write_output ## estimate width
    if(write_output == T){
      if(height == "auto" | width == "auto"){
        ## estimate width upon legend number of 'fill'
        n = length(unique(metadata$vis_class))
        height = 8
        width = ifelse(n <= 17, 10, ## 'class' less than 17
                       ifelse(n <= 34, 12.5,
                              ifelse(n <= 51, 15, 18)))
      }
      ## output
      ggplot2::ggsave(p, file = paste0(output, "/", nebula_name, "_graph.svg"),
             width = width, height = height)
    }
    cat("[INFO] MCnebula Job Done: annotate_child_nebulae\n")
    if(return_plot == T){
      return(p)
    }
  }
## function gather all subview
gather_subview <-
  function(
           subview,
           x,
           y,
           width,
           height,
           p = get("p", envir = get("envir_meta"))
           ){
    p <- p + ggimage::geom_subview(x = x, y = y, width = width, height = height,
                          subview = subview)
    assign("p", p, envir = get("envir_meta"))
    return("Done")
    ##
  }
## funtion merge image, involves nodes (may include ppcp bar), structure, and network layout (with edges)
merge_image <-
  function(
           p, ## ggplot2 object
           layout_n,
           tmp_ppcp,
           ...
           ){
    ## ---------------------------------------------------------------------- 
    ## check svg image
    df <- dplyr::select(layout_n, x, y, name, tanimotoSimilarity) %>%
      dplyr::mutate(nodes_path = paste0(tmp_ppcp, "/", name, ".svg"),
                    check_nodes = file.exists(nodes_path)) %>%
      dplyr::filter(check_nodes == T)
    cat("## read_cairo_svg:", nrow(df), "(number)\n")
    ## read svg image
    subview_list <- pbapply::pblapply(df$name, base_read_cairo,
                                      path = tmp_ppcp, 
                                      ...)
    ## ---------------------------------------------------------------------- 
    ## calculate width and height for subview, according to attributes of tanimotoSimilarity
    df <- dplyr::mutate(df,
                        width = ifelse(is.na(tanimotoSimilarity) == T, 1,
                                       1 + tanimotoSimilarity),
                        height = width)
    ## ---------------------------------------------------------------------- 
    ## as subview 
    cat("## Advance visualization: gather_subview\n")
    pbapply::pbmapply(
                      gather_subview, ## function
                      subview_list,
                      df$x,
                      df$y,
                      df$width,
                      df$height
    )
  }
