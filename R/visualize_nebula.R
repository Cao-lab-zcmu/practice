visualize_nebula <-
  function(
           nebula_name,
           ...
           ){
    p <- annotate_child_nebulae(
                                nebula_name = nebula_name,
                                write_output = F,
                                plot_structure = F,
                                plot_ppcp = F,
                                merge_image = F,
                                return_plot = T,
                                ...)
    return(p)
  }
