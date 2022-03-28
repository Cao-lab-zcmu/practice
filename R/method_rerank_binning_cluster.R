method_rerank_binning_cluster <-
  function(
           apset,
           cluster_cutoff = seq(0.9, 0.1, by = -0.1),
           least_size = 3
           ){
    ## cluster via function of ChemmineR
    cat("## method_rerank_binning_cluster: ChemmineR::cmp.cluster\n")
    meta_rank <- ChemmineR::cmp.cluster(db = apset, cutoff = cluster_cutoff) %>%
      dplyr::select(ids, starts_with("CLSZ_")) %>%
      ## convert into long table
      reshape2::melt(id.var = "ids", variable.name = "cutoff", value.name = "size") %>%
      ## get 'cutoff' and as.numeric
      dplyr::mutate(cutoff = as.numeric(gsub("CLSZ_", "", cutoff))) %>%
      ## get '.id' and 'structure_rank'
      tidyr::separate(col = "ids", into = c(".id", "structure_rank"), sep = "_", remove = T) %>%
      dplyr::mutate(structure_rank = as.numeric(structure_rank)) %>%
      dplyr::arrange(desc(cutoff), desc(size), structure_rank) %>%
      ## at least, the size of cluster reach 'least_size', contribute to re-rank
      dplyr::filter(!(cutoff >= min(cluster_cutoff) & size <= least_size)) %>%
      ## for each .id, only the top 1 (according to 'cutoff', 'size', 'structure_rank', sequentialy) retain
      dplyr::distinct(.id, .keep_all = T)
    return(meta_rank)
  }
