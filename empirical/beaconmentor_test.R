
# was getting error running devtools::install_github(), so used this approach instead
# download.file("https://github.com/lmiratrix/clusterRCT/archive/refs/heads/main.zip", destfile = "clusterRCT.zip")
# unzip("clusterRCT.zip")
# install.packages("clusterRCT-main", repos = NULL, type = "source")

# load libraries ----
# devtools::install_github("lmiratrix/clusterRCT")
library(tidyverse)
library(haven)
library(clusterRCT)
library(tools)
library(writexl)



# load cleaned data ----

f <- "/data/share/pse/IES CRT/sd7/clean study datasets/beaconmentor.RDS"
assign(file_path_sans_ext(basename(f)),
       readRDS(f))


# Run compare_methods() ----

## Beacon Mentor ----

# outcomes_beaconmentor <- c('S8TGPA', 'S8TOTCR', 'SRPRST2')
# covars_beaconmentor <- c('RADEV', 'S8TEACHID', 'S8FULTMEC', 'ZFemale', 'RATIME', 'ZAGELT24', 'PTTOP75C', 'PTBTM25C')


# example with no covariates
cm_bm_nocovar <- compare_methods(S8TOTCR ~ racode_b | S8CLUSTE | RABLOCK,
                                 data = beaconmentor,
                                 include_method_characteristics = F)

# example with one covariate
cm_bm_onecovar <- compare_methods(S8TOTCR ~ racode_b | S8CLUSTE | RABLOCK,
                                 data = beaconmentor,
                                 control_formula = ~ ZFemale,
                                 include_method_characteristics = F)

