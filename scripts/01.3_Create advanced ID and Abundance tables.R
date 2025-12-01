## Create advanced ID and Abundance tables
## Author: Emily H
## Created: August 6, 2025

## Last edited: August 6, 2025

# Clear memory
rm(list=ls())
gc()

#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("abmi.themes")
#install.packages("Matrix")

library(readxl)
library(abmi.themes)
library(tidyverse)
library(Matrix)

setwd("~/GitHub/New aquatic invertebrate analyses")

#### import and clean CABIN/ABMI site data ####
#read in fine data
SppData.raw <- read_rds("output/cleaned CABIN fine data.rds")
str(SppData.raw)

##split into species/genus

#read in UMOS data
umosData.raw <- read_rds("output/cleaned CABIN UMOS data.rds")
str(umosData.raw)

#read in coarse data (including MCN data)
SppData.coarse <- read_rds("output/cleaned CABIN coarse data.rds") %>%
    rename(Coleoptera.adults = "Coleoptera (adults)",
           Coleoptera.larvae = "Coleoptera (larvae)",
           Other.Diptera = "Other Diptera",
           Other.Laevicaudata = "Other: Laevicaudata")
str(SppData.coarse)

#### Coarse group summary ####
coarse.summary <- SppData.coarse %>%
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 3) %>%
    dplyr::select(!c(Site, Year, SiteYear, MC_counted)) %>%
    group_by(Protocol) %>%
    group_split() %>%
    map_df(~ {
        taxa_data <- select(.x, where(is.numeric))  # assumes taxa are in numeric columns

        taxa_counts <- colSums(taxa_data > 0, na.rm = TRUE)  # number of non-zero occurrences per taxon

        tibble(
            Protocol = unique(.x$Protocol),
            total_taxa = sum(taxa_counts > 0),
            rare_taxa = sum(taxa_counts > 0 & taxa_counts < 3),
            unique_taxa = sum(taxa_counts == 1)
        )
    })
