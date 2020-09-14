library(readr)
library(dplyr)
library(stringr)

# https://steinbaugh.com/posts/illumina-sequencer.html
bind_cols(
  read_table2(here::here("output/fastq_filenames.txt"), col_names = "file"),
  read_table2(here::here("output/readnames.txt"), col_names = "header")
) %>% 
  mutate(sample = str_extract(file, "(?<=MDTM\\.)[[:digit:]]+(?=\\.)"),
         instrument = case_when(
           str_detect(header, "ACXX:") ~ "HiSeq High-Output v3",
           str_detect(header, "BCXX:") ~ "HiSeq v1.5",
         )) %>% 
  write_csv(here::here("output/fastqheader_instrument.csv"), na = "")