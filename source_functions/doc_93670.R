#' ---
#' title: "Summarize coverage for reference individual"
#' author: "Harly Durbin"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(purrr)
library(dplyr)
library(purrr)
library(readr)

#' 
#' # Notes & questions
#' 
#' * _summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases
#' * _statistics: coverage histograms (# locus with X coverage), aggregated over all bases
#' 
#' module load java/openjdk/java-1.8.0-openjdk
#' module load gatk/gatk-3.5
#' srun -p BioCompute,htc4,hpc5 --account animalsci --mem 120G -c 18 -t 24:00:00 java -Djava.io.tmpdir=temp -XX:ParallelGCThreads=6 -Xmx24g -jar /cluster/software/gatk/gatk-3.5/GenomeAnalysisTK.jar -nt 12 -R /storage/htc/deckerlab/FROM_MUG01/MUG01_N/lwwvd/madtom/assembly/madtom_assembly_170124/CA/10-gapclose/backup_genome.scf.fasta -I output/call_genotypes/93670/93670.realigned.bam -T DepthOfCoverage -o output/call_genotypes/93670/93670.wg.realigned.bam.coverage -omitBaseOutput -omitIntervals --omitLocusTable
#' 
#' # Setup 
#' 
## -----------------------------------------------------------------------------------------
interval_doc <-
  list.files(here::here("output/call_genotypes/93670"),
             pattern = "sample_summary") %>% 
  purrr::set_names() %>% 
  purrr::map_dfr(~ read_table2(here::here(glue::glue("output/call_genotypes/93670/{.x}")),
                               na = "N/A"),
                 .id = "interval") %>% 
  mutate(interval = stringr::str_extract(interval, "(?<=93670\\.)[[:digit:]]{1,3}(?=\\.)")) %>% 
  filter(sample_id == "MDTM93670") %>% 
  select(-sample_id)

#' 
## -----------------------------------------------------------------------------------------
head(interval_doc)

#' 
#' # Summarize
#' 
## -----------------------------------------------------------------------------------------
interval_doc %>% 
  summarise(mean(mean))

#' 
