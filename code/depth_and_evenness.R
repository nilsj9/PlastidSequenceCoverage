#!/usr/bin/R
# "R script to conduct statistical analyses among several plastid genomes"

#### AUTHOR AFFILIATION ####
# ---------------------------------------------------------------------------- #
#contributor  = "Nils Jenke"
#email        = "nilsj24@zedat.fu-berlin.de"
#version      = "2021.02.04.2300"
#Rversion     = R version 4.0.2 (2020-06-22)

#### LIBRARIES ####
# ---------------------------------------------------------------------------- #
if (!require("pacman")) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(RCircos,
               genbankr,
               devtools,
               tidyverse,
               vroom)
pacman::p_install_gh("nilsj9/PACVr")
library(PACVr)

#### FUNCTIONS ####
# ---------------------------------------------------------------------------- #
evennessScore <- function(coverage) {
  coverage_mean <- round(mean(coverage))
  D2 <- coverage[coverage <= coverage_mean]
  E <- 1 - (length(D2) - sum(D2) / coverage_mean) / length(coverage)
  return(E)
}

effSize_conv <- function(x) {
  return(2 * x / (sqrt(1 - x ^ 2)))
}


#### SAMPLE LIST ####
# ---------------------------------------------------------------------------- #
sample_file <- list.files(
  path = "./",
  recursive = TRUE,
  pattern = "^samples_list.csv",
  full.names = TRUE
)

if (length(file.exists(sample_file)) == 0) {
  acquire_data()
  data <- read.csv("samples_list.csv")
}

sample_list <- read.csv(sample_file[1], header = TRUE)

#### WORKING DIR ####
# ---------------------------------------------------------------------------- #
w_dir <- "/PATH/TO/WORKING/DIRECTORY/"

#### INPUT FILES ####
# ---------------------------------------------------------------------------- #
# coverage files
bed_files <- list.files(
  path = paste0(w_dir,"/samples"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = ".bed"
)
# meta data files
metadata_files <- list.files(
  path = paste0(w_dir,"/samples"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "metadata.csv"
)
# annotation files
gb_files <- list.files(
  path = paste0(w_dir,"/samples"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "[0-9].gb"
)

### READ META DATA ###
# ---------------------------------------------------------------------------- #
# lineage information
sample_list$Lineage <-
  sapply(gb_files, function(x)
    gsub(".*Embryophyta ", "", paste(
      trimws(parseGenBank(x)$SOURCE$lineage), collapse = " "
    )))

meta_data <- purrr::map_df(metadata_files, ~ vroom::vroom(
  .x,
  delim = ",",
  col_types = c(
    Number_N = "i",
    Mismatches = "i",
    avgLength = "i",
    Model = "c",
    `Assembly Method` = "c",
    `Sequencing Technology` = "c"
  )
))
meta_data <- cbind(sample_list, meta_data)

### READ COVERAGE DATA ###
# ---------------------------------------------------------------------------- #
# coverage information
coverage_data <- purrr::map(bed_files, ~ vroom::vroom(
  .x,
  col_types = c(
    Chromosome = "c",
    chromStart = "i",
    chromEnd = "i",
    coverage = "i",
    lowCoverage = "c",
    gene = "c"
  )
))
# remove regions < 250
coverage_data <-
  lapply(coverage_data, function(x)
    x[-which(x$chromEnd - x$chromStart < 249), ])

# Assign names to list elements
names(coverage_data) <- unlist(unname(lapply(sample_list$Accession,
                                             function(x)
                                               c(
                                                 paste(x, "_genes", sep = ""),
                                                 paste(x, "_noncoding", sep = ""),
                                                 paste(x, "_regions", sep = "")
                                               ))))

# Split list by coverage type
coverage_regions  <-
  coverage_data[grep("regions", names(coverage_data))]
coverage_genes <- coverage_data[grep("genes", names(coverage_data))]
coverage_noncoding <-
  coverage_data[grep("noncoding", names(coverage_data))]

# Add coverage depth data
# plastid partition
regions <- lapply(coverage_regions, function(x)
  x %>%
    tidyr::separate_rows(Chromosome) %>%
    dplyr::group_by(Chromosome) %>%
    dplyr::summarise(
      lowCoverage = sum(lowCoverage == "*", na.rm = TRUE),
      .groups = "drop"
    ))
regions <- lapply(regions, function(x)
  x %>%
    tidyr::spread(., Chromosome, lowCoverage))
regions <- select(dplyr::bind_rows(regions), IRa, IRb, LSC, SSC)
meta_data <- dplyr::bind_cols(meta_data, regions)

# coding/non-coding partition
meta_data$coding <-
  unlist(unname(lapply(coverage_genes, function(x)
    sum(
      !is.na(x$lowCoverage == "*")
    ))))
meta_data$noncoding <-
  unlist(unname(lapply(coverage_noncoding, function(x)
    sum(!is.na(x$lowCoverage == "*")))))

meta_data <- rename(meta_data, Ns = Number_N)

### CALCULATE E-SCORE ###
# ---------------------------------------------------------------------------- #
# Calculate coverage Escore
meta_data$Escore <- unname(sapply(coverage_regions,
                                  function(x)
                                    evennessScore(x$coverage)))

### SORT DATAFRAME ###
# ---------------------------------------------------------------------------- #
# sort columns
meta_data <- meta_data %>%
  dplyr::select(
    Samples,
    Accession,
    SRA,
    Escore,
    LSC,
    IRb,
    SSC,
    IRa,
    coding,
    noncoding,
    Ns,
    Mismatches,
    Model,
    avgLength,
    `Assembly Method`
  )


#### DATASET PRE-PROCESSING ####
# ---------------------------------------------------------------------------- #
# resolve spelling mistake
meta_data[which(meta_data$Samples == "Eragrostis_tef"), ]$`Assembly Method` <-
  "Geneious"

# unify assembly methods
nameList <- c(
  "Velvet",
  "Consed",
  "Newbler",
  "Spades",
  "CLC",
  "SOAP",
  "Mira",
  "Ray",
  "Geneious",
  "Obitools",
  "Flash",
  "PriceTI",
  "Spades",
  "Yasra",
  "GS",
  "abyss",
  "Allpath"
)
assignList <- c(
  "Velvet",
  "Consed",
  "Newbler Assembler",
  "SPAdes",
  "CLC Assembly",
  "SOAPdenovo",
  "MIRA",
  "Ray",
  "Geneious",
  "OBITools",
  "FLASH",
  "PriceTI",
  "SPAdes",
  "YASRA",
  "GS De novo assembler",
  "ABySS",
  "ALLPATHS-LG"
)

for (i in 1:length(nameList)) {
  meta_data$`Assembly Method`[grepl(nameList[i],
                                    meta_data$`Assembly Method`,
                                    ignore.case = TRUE)] <-
    assignList[i]
}

vroom::vroom_write(meta_data,"./stat_ready_metadata.csv", delim = ",")