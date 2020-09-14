## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(readr)
library(tidyr)
library(magrittr)

# First, generated list of contigs

# `grep -nr ">" /storage/htc/deckerlab/FROM_MUG01/MUG01_N/lwwvd/madtom/assembly/madtom_assembly_170124/CA/10-gapclose/backup_genome.scf.fasta &> output/call_genotypes/interval_list/all_contigs.txt`

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
contigs <-
  read_table2(here::here("output/call_genotypes/interval_list/all_contigs.txt"), col_names = "contig") %>% 
  separate(contig, into = c("start", "contig_name"), sep = ":>") %>% 
  mutate(
    start = as.integer(start),
    stop = lead(start)-1,
    # total number of lines to calculate last contig size
    stop = if_else(is.na(stop), 12892536, stop),
    size = stop - start) %>% 
  select(contig_name, everything()) 


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
contigs %>% 
  arrange(size)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

assigned <-
  tribble(~ contig_name, ~ group_number,
          "initialize", 1)

repeat {
  
  n <- max(assigned$group_number) + 1
  
  choose <-
    contigs[!contigs$contig_name %in% assigned$contig_name,] %>%
    # Arrange by position
    arrange(start) %>%
    mutate(cumsum = cumsum(size)) %>%
    filter(200000 >= cumsum) %>%
    mutate(group_number = n) %>%
    select(contig_name, group_number, start)
  
  choose <-
    # Limit number of contigs
    if (length(choose$contig_name) > 500) {
      choose %>%
        # Arrange by position
        arrange(start) %>% 
        slice(1:500) %>% 
        select(-start)
      
    } else{
      choose %>% 
        select(-start)
    }
  
  print(glue::glue("Assigned group {n}"))
  
  assigned %<>%
    bind_rows(choose)
  
  if (length(contigs$contig_name) == length(assigned$contig_name) - 1) {
    print("Finished")
    
    assigned %<>%
      filter(contig_name != "initialize") %>% 
      mutate(group_number = group_number - 1)
    
    break
  }
  
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
assigned %>% 
  group_by(group_number) %>% 
  tally(sort = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
contigs %>% 
  arrange(desc(size))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
assigned %>% 
  group_by(group_number) %>% 
  group_map(~ write_delim(.x, here::here(glue::glue("output/call_genotypes/interval_list/{.y}.interval_list")), col_names = FALSE))

assigned %>% 
  distinct(group_number) %>% 
  write_delim(here::here("output/call_genotypes/interval_list/all_intervals.txt"), col_names = FALSE)

assigned %>% 
  distinct(group_number) %>% 
  mutate(path = glue::glue("output/joint_genotyping/{group_number}/MDTM.combined.{group_number}.filter.SNP.vcf.gz")) %>% 
  select(path) %>% 
  write_delim(here::here("output/joint_genotyping/merge_files.list"), col_names = FALSE) 