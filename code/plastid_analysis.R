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
               vroom,
               ggpubr,
               xtable,
               rstatix,
               reporttools)
pacman::p_install_gh("nilsj9/PACVr")
library(PACVr)

#### WORKING DIRECTORY ####
# ---------------------------------------------------------------------------- #
if (any(commandArgs() == "--interactive") == TRUE) {
  if (!require(rstudioapi))
    suppressWarnings(install.packages('rstudioapi'))
  library(rstudioapi)
  file_path <- rstudioapi::getSourceEditorContext()$path
} else{
  return(dirname(sys.frame(1)$ofile))
}

setwd(dirname(file_path))

#### FUNCTIONS ####
# ---------------------------------------------------------------------------- #
acquire_data <- function() {
  raw_data <-
    read.table(
      "https://raw.githubusercontent.com/chloroExtractorTeam/benchmark/master/data/real/real_datasets.tsv",
      sep = "\t",
      header = FALSE
    )
  raw_data <- raw_data[c("V1", "V2", "V4")]
  colnames(raw_data) <- c("SRA", "Sample", "Accession")
  raw_data$Sample <- gsub(" ", "_", raw_data$Sample)
  
  ind <- rep(TRUE, length(raw_data$Accession))
  for (i in seq(1, nrow(raw_data))) {
    tryCatch({
      genbankr::readGenBank(genbankr::GBAccession(raw_data$Accession[i]))
    },
    error = function(e) {
      ind[i] <<- FALSE
    })
  }
  data <- raw_data[ind,]
  data <- data[order(data$Sample),]
  row.names(data) <- 1:nrow(data)
  write.csv(data,
            "genbankr_samples_list.csv",
            quote = FALSE,
            row.names = FALSE)
  return(data)
}

# Thanks to Oexle (https://www.nature.com/articles/jhg201621)
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

#### INPUT FILES ####
# ---------------------------------------------------------------------------- #
# coverage files
bed_files <- list.files(
  path = "./samples",
  full.names = TRUE,
  recursive = TRUE,
  pattern = ".bed"
)
# meta data files
metadata_files <- list.files(
  path = "./samples",
  full.names = TRUE,
  recursive = TRUE,
  pattern = "metadata.csv"
)
# annotation files
gb_files <- list.files(
  path = "./samples",
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


#### OUTPUT FILES ####
# ---------------------------------------------------------------------------- #

### META DATA ###
library(xtable)
meta_table <-
  select(
    meta_data,
    Accession,
    Model,
    `Assembly Method`,
    avgLength,
    Escore,
    LSC,
    IRb,
    SSC,
    IRa,
    coding,
    noncoding,
    Ns,
    Mismatches
  )
meta_table$`Assembly Method` <-
  replace(meta_table$`Assembly Method`, which(is.na(meta_table$`Assembly Method`)), "n.s.")
meta_table$avgLength <-
  ifelse(
    meta_table$avgLength <= 300,
    paste0(meta_table$avgLength, " (s)"),
    paste0(meta_table$avgLength, " (m)")
  )
meta_table$Model <- gsub("Illumina ", "", meta_table$Model)
meta_table$Model <- gsub("Genome Analyzer", "GA", meta_table$Model)
colnames(meta_table) <-
  c(
    "Sample",
    "Sequencing method",
    "Asy. software",
    "Avg. read length (group)",
    "E-score",
    "LSC",
    "IRb",
    "SSC",
    "IRa",
    "Cod.",
    "Non-cod.",
    "Ns",
    "Mismatches"
  )

meta_table_tex <- xtable(meta_table)
dir.create("./images")
print(
  meta_table_tex,
  file = "./images/sample_table.tex",
  compress = FALSE,
  include.rownames = FALSE,
  NA.string = getOption("xtable.NA.string", "-")
)

### LINEAGE DATA ###
sample_table <-
  select(sample_list, Samples, Accession, SRA, Lineage)
colnames(sample_table) <-
  c("Species", "NCBI Nucleotide", "NCBI SRA", "order,family,genus")
samples_xtab <- xtable(sample_list)
print(
  samples_xtab,
  file = "./images/sample_lineage.tex",
  compress = FALSE,
  include.rownames = FALSE
)

### DESCRIPTIVE STATISTICS ###
descr_stats <-
  meta_data %>% select(Escore, LSC, IRb, SSC, IRa, coding, noncoding, Ns, Mismatches)
colnames(descr_stats) <-
  c("E-score",
    "LSC",
    "IRb",
    "SSC",
    "IRa",
    "coding",
    "noncoding",
    "Ns",
    "mismatches")
descr_stats <-
  rstatix::get_summary_stats(descr_stats,
                             show = c("n", "min", "q1", "median", "mean", "q3", "max", "sd", "iqr"))
descr_stats <- descr_stats[c(2, 5, 4, 9, 3, 1, 7, 8, 6),]
descr_stats["NA"] <- nrow(meta_data) - descr_stats$n
descr_stats_xtab <- xtable(descr_stats)
print(
  descr_stats_xtab,
  file = "./images/descr_stats.tex",
  compress = FALSE,
  include.rownames = FALSE
)

nom_stats <-
  meta_data %>% select(Model, avgLength, `Assembly Method`)
nom_stats$avgLength <-
  ifelse(nom_stats$avgLength <= 300, "small", "medium")
nom_stats$Model <- gsub("Illumina ", "", nom_stats$Model)
nom_stats <-
  as.data.frame(apply(nom_stats, 2, function(x)
    replace(x, x %in% names(
      which(table(x) < 5) == TRUE
    ), NA)))
nom_stats$Escore <- meta_data$Escore
for (i in colnames(nom_stats)[-ncol(nom_stats)]) {
  tmp <- freq_table(nom_stats, i, na.rm = FALSE)
  tmp$cumsum <- cumsum(tmp$prop)
  tmp <- nom_stats %>%
    select(i, Escore) %>%
    group_by_at(i) %>%
    summarise(E_mean = mean(Escore), E_sd = sd(Escore)) %>%
    bind_cols(tmp[2:4], .)
  print(tmp)
  print(sum(tmp$n))
}

# Transform assembly quality by inverse
meta_data$Mismatches <- 1 / (meta_data$Mismatches + 1)
meta_data$Ns <- 1 / (meta_data$Ns + 1)


### CORRELATION ANALYSIS ###
# ---------------------------------------------------------------------------- #
# Number of N's correlated with Escore
sub_Ns <- meta_data[, c("Ns", "Escore")]
Ns <- ggscatter(
  sub_Ns,
  x = "Ns",
  y = "Escore",
  add = "reg.line",
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE,
  xlab = "Assembly quality (#Ns)",
  ylab = "E-score",
  xlim = c(0, 1.01)
) +
  stat_cor(
    method = "spearman",
    p.accuracy = 0.001,
    r.accuracy = 0.01,
    label.y = 1.025,
    label.x = 0.2
  ) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
Ns <- Ns +
  font("axis.title", size = 10) +
  font("axis.text", size = 10)

# Number of Mismatches correlated with Escore
sub_MMs <- meta_data[, c("Mismatches", "Escore")]
sub_MMs <- sub_MMs[!is.na(sub_MMs$Mismatches),] # remove NA's
rownames(sub_MMs) <- 1:nrow(sub_MMs)
MMs <- ggscatter(
  sub_MMs,
  x = "Mismatches",
  y = "Escore",
  add = "reg.line",
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE,
  xlab = "Assembly quality (#Mismatches)",
  ylab = "E-score",
  xlim = c(0, 1.01)
) +
  stat_cor(
    method = "spearman",
    p.accuracy = 0.001,
    r.accuracy = 0.01,
    label.y = 1.025,
    label.x = 0.2
  ) +
  theme(plot.margin = unit(c(0.05, 0, 0, 0), "cm"))

MMs <- MMs +
  font("axis.title", size = 10) +
  font("axis.text", size = 10)

### CORRELATION PLOT  ###
# ---------------------------------------------------------------------------- #
ggarrange(
  Ns,
  MMs,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  label.y = 1.02,
  label.x = -0.02
)
ggsave(
  filename = "SpearmanRankCorr.pdf",
  path = "./images",
  device = 'pdf',
  dpi = 700,
  width = 15.49779,
  height = 7.654888,
  #7.25,
  units = ("cm")
)

### VARIANCE ANALYSIS ###
# ---------------------------------------------------------------------------- #
# Is the number of regions with significantly reduced coverage significantly
# different in coding versus non-coding regions?
wilcox_df <-
  tidyr::gather(meta_data[, c("coding", "noncoding")],
                "coding",
                "noncoding",
                key = "sequences",
                value = "low_coverage")
coding_noncoding <- ggpubr::ggboxplot(
  wilcox_df,
  x = "sequences",
  y = "low_coverage",
  color = "sequences",
  add = "jitter",
  bxp.errorbar = TRUE,
  legend = "none",
  ylab = "#inst. of sig. red. cov. depth"
) +
  stat_compare_means(label.y = 100, label.x = 1.25) +
  rremove("xlab") +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
        axis.title.y = element_text(hjust = 0.1))
ggpar(coding_noncoding, ylim = c(0, 100))
coding_noncoding <- coding_noncoding +
  font("axis.text", size = 10) +
  font("ylab", size = 10)

wilcox.test(wilcox_df$low_coverage ~ wilcox_df$sequences)
car::leveneTest(wilcox_df$low_coverage, wilcox_df$sequences)

# ---------------------------------------------------------------------------- #

# Is the number of regions with significantly reduced coverage significantly
# different across the four partitions of the plastid genome (i.e., LSC, IR, SSC)?

kruskal_df <-
  tidyr::gather(meta_data[, c("LSC", "IRb", "SSC", "IRa")],
                "LSC",
                "IRb",
                "SSC",
                "IRa",
                key = "regions",
                value = "low_coverage")
plastid_partition <- ggpubr::ggboxplot(
  kruskal_df,
  x = "regions",
  y = "low_coverage",
  color = "regions",
  add = "jitter",
  bxp.errorbar = TRUE,
  legend = "none",
  ylab = "#inst. of sig. red. cov. depth"
) +
  stat_compare_means(label.y = 100, label.x = 2) +
  rremove("xlab") +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"),
        axis.title.y = element_text(hjust = 0.1))
ggpar(plastid_partition, ylim = c(0, 100))
plastid_partition <- plastid_partition +
  font("axis.text", size = 10) +
  font("ylab", size = 10)

kruskal.test(kruskal_df$low_coverage ~ kruskal_df$regions)

car::leveneTest(kruskal_df$low_coverage, kruskal_df$regions)
test <-
  pairwise.wilcox.test(kruskal_df$low_coverage,
                       kruskal_df$regions,
                       p.adjust.method = "BH")

# ---------------------------------------------------------------------------- #

# Is evenness of coverage significantly different across different assembly software tools?
kruskal_df2 <- meta_data[, c("Assembly Method", "Escore")]
colnames(kruskal_df2) <- c("Assembly", "Escore")
kruskal_df2 <-
  kruskal_df2[!is.na(kruskal_df2$Assembly),] # remove na
kruskal_df2 <- kruskal_df2 %>%
  group_by(Assembly) %>%
  filter(n() >= 5) %>% # remove sample size <= 5
  ungroup(.)
assembly <- ggpubr::ggboxplot(
  kruskal_df2,
  x = "Assembly",
  y = "Escore",
  color = "Assembly",
  add = "jitter",
  legend = "none",
  bxp.errorbar = TRUE,
  ylab = "E-score"
) +
  stat_compare_means(label.y = 1, label.x = 2) +
  rotate_x_text(20) +
  rremove("xlab") +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
ggpar(assembly, ylim = c(0, 1.0))
assembly <- assembly +
  font("axis.text", size = 10) +
  font("ylab", size = 10)

kruskal.test(kruskal_df2$Escore ~ kruskal_df2$Assembly)

car::leveneTest(kruskal_df2$Escore, kruskal_df2$Assembly)
pairwise.wilcox.test(kruskal_df2$Escore,
                     kruskal_df2$Assembly,
                     p.adjust.method = "none")

# ---------------------------------------------------------------------------- #

# Is evenness of coverage significantly different across different sequencing forms?
kruskal_df3 <- meta_data[, c("Model", "Escore")] # remove nan
kruskal_df3$Model <- gsub("Illumina ", "", kruskal_df3$Model)
kruskal_df3$Model <-
  gsub("Genome Analyzer ", "GA ", kruskal_df3$Model)
kruskal_df3 <- kruskal_df3 %>%
  group_by(Model) %>%
  filter(n() >= 5) %>% # remove sample size <= 5
  ungroup(.)
sequencing <- ggpubr::ggboxplot(
  kruskal_df3,
  x = "Model",
  y = "Escore",
  color = "Model",
  add = "jitter",
  legend = "none",
  bxp.errorbar = TRUE,
  ylab = "E-score"
) +
  stat_compare_means(label.y = 1, label.x = 2.75) +
  rotate_x_text(25) +
  rremove("xlab") +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
ggpar(sequencing, ylim = c(0, 1.0))
sequencing <- sequencing +
  font("axis.text", size = 10) +
  font("ylab", size = 10)

kruskal.test(kruskal_df3$Escore ~ kruskal_df3$Model)

car::leveneTest(kruskal_df3$Escore, kruskal_df3$Model)
test <-
  pairwise.wilcox.test(kruskal_df3$Escore, kruskal_df3$Model, p.adjust.method = "BH")

# ---------------------------------------------------------------------------- #

# Is evenness of coverage significantly different across different read lengths?
wilcox_df2 <- meta_data[, c("avgLength", "Escore")]
wilcox_df2$Readlength <-
  ifelse(wilcox_df2$avgLength >= 100 &
           wilcox_df2$avgLength <= 300,
         "short",
         "medium")
read_length <- ggpubr::ggboxplot(
  wilcox_df2,
  x = "Readlength",
  y = "Escore",
  color = "Readlength",
  legend = "none",
  add = "jitter",
  bxp.errorbar = TRUE,
  ylab = "E-score"
) +
  stat_compare_means(label.y = 1, label.x = 1) +
  rremove("xlab") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
ggpar(read_length, ylim = c(0, 1.00))
read_length <- read_length +
  font("axis.text", size = 10) +
  font("ylab", size = 10)

wilcox.test(wilcox_df2$Escore ~ wilcox_df2$Readlength)

### VARIANCE ANALYSIS PLOT ###
# ---------------------------------------------------------------------------- #
ggarrange(
  plastid_partition,
  coding_noncoding,
  sequencing,
  assembly,
  read_length,
  nrow = 3,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E"),
  label.y = 1.02,
  label.x = -0.02,
  heights = c(1, 1.1156625968355375019018071372623, 0.99)
)
ggsave(
  filename = "AmongSampleVarianceAnalysis.pdf",
  path = "./images",
  device = 'pdf',
  dpi = 700,
  width = 15.49779,
  height = 20.75,
  units = ("cm")
)

### OUTLIER PLOT ###
# ---------------------------------------------------------------------------- #

test <- ggpubr::ggboxplot(
  meta_data,
  x = NULL,
  y = "Escore",
  add = "jitter",
  bxp.errorbar = TRUE,
  label = "Accession",
  font.label = list(size = 10, color = "red"),
  label.select = list(top.down = 13),
  repel = TRUE,
  xlab = "Samples"
) +
  rremove("x.text") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  font("y.text", size = 10) +
  font("xylab", size = 10)
test$layers[[3]]$aes_params$size <- 0.5
ggsave(
  filename = "Outlier.pdf",
  path = "./images",
  device = 'pdf',
  dpi = 700,
  width = 7.748895,
  height = 7.700765,
  units = ("cm")
)

### RESULT TABLE ###
# KRUSKAL WALLIS
res <-
  rstatix::kruskal_test(kruskal_df[!is.na(kruskal_df$low_coverage), ], low_coverage ~ regions) %>%
  bind_cols(., select(
    rstatix::kruskal_effsize(kruskal_df, low_coverage ~ regions),
    effsize
  )) %>%
  bind_rows(., bind_cols(
    rstatix::kruskal_test(kruskal_df3[!is.na(kruskal_df3$Escore), ], Escore ~ Model),
    select(rstatix::kruskal_effsize(kruskal_df3, Escore ~ Model), effsize)
  )) %>%
  bind_rows(., bind_cols(
    rstatix::kruskal_test(kruskal_df2[!is.na(kruskal_df2$Escore), ], Escore ~ Assembly),
    select(
      rstatix::kruskal_effsize(kruskal_df2, Escore ~ Assembly),
      effsize
    )
  ))

rstatix::kruskal_test(kruskal_df3[!is.na(kruskal_df3$Escore), ], Escore ~ Model)
rstatix::kruskal_effsize(kruskal_df3, Escore ~ Model)

# WILCOXON RANK SUM
res <-
  rstatix::pairwise_wilcox_test(kruskal_df[!is.na(kruskal_df$low_coverage), ], low_coverage ~ regions, p.adjust.method = "BH") %>%
  bind_cols(., select(
    rstatix::wilcox_effsize(kruskal_df, low_coverage ~ regions),
    effsize
  )) %>%
  bind_rows(., bind_cols(
    rstatix::wilcox_test(wilcox_df[!is.na(wilcox_df$low_coverage), ], low_coverage ~ sequences),
    (select(
      rstatix::wilcox_effsize(wilcox_df, low_coverage ~ sequences),
      effsize
    ))
  )) %>%
  bind_rows(., bind_cols(
    rstatix::pairwise_wilcox_test(kruskal_df3[!is.na(kruskal_df3$Escore), ], Escore ~ Model, p.adjust.method = "BH"),
    (select(
      rstatix::wilcox_effsize(kruskal_df3, Escore ~ Model),
      effsize
    ))
  )) %>%
  bind_rows(., bind_cols(
    rstatix::pairwise_wilcox_test(kruskal_df2[!is.na(kruskal_df2$Escore), ], Escore ~ Assembly, p.adjust.method = "BH"),
    (select(
      rstatix::wilcox_effsize(kruskal_df2, Escore ~ Assembly),
      effsize
    ))
  )) %>%
  bind_rows(., bind_cols(
    rstatix::wilcox_test(wilcox_df2[!is.na(wilcox_df2$Escore), ], Escore ~ Readlength),
    (select(
      rstatix::wilcox_effsize(wilcox_df2, Escore ~ Readlength),
      effsize
    ))
  ))
res$effsize <- round(effSize_conv(res$effsize), 3)
res <- rename(res, variable = .y.)
res <- select(res, -p.adj.signif)



res_xtab <- xtable(res, digits = 3)
print(
  res_xtab,
  file = "./images/results.tex",
  compress = FALSE,
  include.rownames = FALSE
)


rstatix::pairwise_wilcox_test(kruskal_df,
                              low_coverage ~ regions,
                              p.adjust.method = "BH",
                              detailed = TRUE)



rstatix::wilcox_test(wilcox_df2, Escore ~ Readlength)
