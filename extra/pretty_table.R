pretty_table <- 
  function(
           df,
           title = "compounds summary",
           subtitle = "LC-MS",
           footnote = "Compounds summary",
           filename = "tmp.html",
           path = ".",
           return_gt = T,
           font = "Times",
           shorter_name = T
           ){
    title = paste0("**", Hmisc::capitalize(title), "**")
    subtitle = paste0("**", Hmisc::capitalize(subtitle), "**")
    colnames(df) <- Hmisc::capitalize(colnames(df))
    t <- gt(df) %>%
      opt_table_font(font=list(font)) %>%
      tab_header(title = md(title),
                 subtitle = md(subtitle)) %>%
      opt_align_table_header(align = "left") %>%
      tab_footnote(footnote = footnote,
                   locations = cells_title(groups = c("title"))) %>%
      opt_table_lines(extent = c("none")) %>%
      cols_align(align = "left",
                 columns = everything()) %>%
      tab_style(style = cell_borders(sides = c("top", "bottom"),
                                     color = "black",
                                     weight = px(1.5),
                                     style = "solid"),
                locations = cells_column_labels()) %>%
      tab_style(style = cell_text(v_align="top"),
                locations = cells_column_labels(columns = everything())) %>%
      tab_style(style = cell_borders(sides = c("bottom"),
                                     color = "black",
                                     weight = px(1.5),
                                     style = "solid"),
                locations = cells_body(columns=everything(),
                                       rows=nrow(df))) %>%
      tab_style(style = cell_text(v_align="top"),
                locations = cells_body(columns = everything()))
    if(shorter_name == T){
      t <- t %>%
        cols_width(Name ~ px(300))
    }
    gtsave(t, filename, path)
    if(return_gt == T)
      return(t)
  }
