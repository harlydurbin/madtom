usage_summary <-
  function(process, keyword) {
    list.files(
      path = here::here(glue::glue("log/psrecord/{process}")),
      pattern = keyword,
      full.names = TRUE
    ) %>%
      set_names(nm = (basename(.) %>%
                        tools::file_path_sans_ext())) %>%
      map_dfr(~ read_table2(
        .x,
        skip = 1,
        col_names = c("time", "cpu_percent", "real_mb", "virtual_mb")
      ), .id = "file") %>%
      mutate(rule = keyword) %>%
      summarise(max_time_minutes = max(time),
                max_cpu_percent = max(cpu_percent),
                max_mb = max(real_mb)) %>%
      mutate(max_time_minutes = max_time_minutes / 60,
             max_gb = max_mb * 0.001)
    
  }

usage_facets <-
  function(process, keyword, sample_num = NULL) {
    logs <-
      list.files(
        path = here::here(glue::glue("log/psrecord/{process}")),
        pattern = keyword,
        full.names = TRUE
      ) %>%
      set_names(nm = (basename(.) %>%
                        tools::file_path_sans_ext())) %>%
      map_dfr(~ read_table2(
        .x,
        skip = 1,
        col_names = c("time", "cpu_percent", "real_mb", "virtual_mb")
      ), .id = "file") %>%
      mutate(rule = keyword) 
    
    logs <- 
      if(!is.null(sample_num)){
        logs %>% 
          filter(file %in% sample(file, sample_num))
      } else logs
    
    logs %>% 
      ggplot(aes(x = time, y = real_mb)) +
      geom_line() +
      facet_wrap( ~ file)
    
  }