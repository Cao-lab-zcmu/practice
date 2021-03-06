generate_parent_nebula <-
  function(
           write_output = T,
           output_format = "graphml",
           output = paste0(.MCn.output, "/", .MCn.results),
           edges_file = paste0(output, "/parent_nebula/parent_nebula_edges.tsv"), # exists edges file
           edges_method = "method_formula_based_spec_compare", # or NULL
           nodes_attributes = .MCn.formula_set,
           nodes_other_attributes = .MCn.structure_set,
           edge_filter = 0.5,
           cpu_cores = 8, 
           ...
           ){
    cat("[INFO] MCnebula run: generate_parent_nebula\n")
    ## main body
    ## ---------------------------------------------------------------------- 
    ## generate edges data
    if(is.null(edges_method) == T){
      ## no edges_method
      cat("# generate_parent_nebula: no edges_uethod used\n")
      edges <- dplyr::as_tibble(cbind(".id_1" = nodes_other_attributes$".id",
                     ".id_2" = nodes_other_attributes$".id")) %>%
        dplyr::mutate(dotproduct = 1, mass_diff = 0)
    }else if(edges_method == "method_formula_based_spec_compare"){
      ## with edges_method
      if(is.null(edges_file) == F & file.exists(edges_file)){
        cat("# generate_parent_nebula: file.exists(edges_file) == T. Escape from time-consuming computation\n")
        edges <- read_tsv(edges_file) %>%
          dplyr::mutate_at(c(".id_1", ".id_2"), as.character) %>%
          dplyr::mutate_at(c(colnames(edges)[3:4]), as.numeric)
      }else{
        edges = method_formula_based_spec_compare(edge_filter = edge_filter, cpu_cores = cpu_cores, ...)
      }
    }
    ## ---------------------------------------------------------------------- 
    ## generate nodes data
    nodes <- nodes_attributes
    ## additional nodes attributes
    if(is.null(nodes_other_attributes) == F){
      nodes <- merge(nodes, nodes_other_attributes, by = ".id", all.x = T, sort = T) %>%
        ## rename the column name, otherwise the column will be choosed as key column in igraph
        dplyr::rename(compound_name = name)
    }
    ## ---------------------------------------------------------------------- 
    ## graph
    parent_nebula <- igraph::graph_from_data_frame(edges, directed = T, vertices = nodes)
    if(write_output == T){
      dir = paste0(output, "/", "parent_nebula")
      if(file.exists(dir) == F){
        dir.create(dir)
      }
      write_graph(parent_nebula,
                  file = paste0(dir, "/", "parent_nebula.", output_format),
                  format = output_format)
      write_tsv(edges, paste0(dir, "/", "parent_nebula_edges.tsv"))
      write_tsv(nodes, paste0(dir, "/", "parent_nebula_nodes.tsv"))
    }
    ## ---------------------------------------------------------------------- 
    ## set as global var for next stage
    .MCn.parent_graph <<- parent_nebula
    .MCn.parent_nodes <<- nodes %>% as_tibble()
    .MCn.parent_edges <<- edges %>% as_tibble()
    cat("[INFO] MCnebula Job Done: generate_parent_nebula\n")
  }
