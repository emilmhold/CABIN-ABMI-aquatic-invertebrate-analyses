## Comparison of collection times between CABIN and ABMI protocols
## Author: Emily H
## Created: July 30, 2025
## Last edited: July 31, 2025

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

#### import and clean sampling data ####
sampling <- read.csv("data/CABIN avg sampling depth and time data.csv") %>%
    mutate(ABMI.Site = str_remove(Site, "CABIN-")) %>%
    mutate(ABMI.Site = str_remove(ABMI.Site, "W")) %>% #remove extra characters from Site column
    mutate(ABMI.SiteYear = paste0(ABMI.Site, "_", Year), .after = Year) %>% #create identifying column
    mutate(Sample_Depth_Mean = as.numeric(str_remove(Sample_Depth_Mean, "M"))) %>% #remove "M" from values and format as numeric
    mutate(TOC = as.numeric(str_remove(TOC, "MIN"))) %>% #remove "MIN" from values and format as numeric
    mutate(Protocol = if_else(str_starts(Site, "CABIN"), "CABIN", "ABMI")) #create protocol column
str(sampling)

##impute values for sites where no data was collected
sampling[is.na(sampling$Sample_Depth_Mean) & sampling$Protocol == "CABIN", "Sample_Depth_Mean"] <- mean(sampling[sampling$Protocol == "CABIN" ,"Sample_Depth_Mean"], na.rm = TRUE)
sampling[is.na(sampling$Sample_Depth_Mean) & sampling$Protocol == "ABMI", "Sample_Depth_Mean"] <- mean(sampling[sampling$Protocol == "ABMI" ,"Sample_Depth_Mean"], na.rm = TRUE)
sampling[is.na(sampling$TOC) & sampling$Protocol == "CABIN", "TOC"] <- mean(sampling[sampling$Protocol == "CABIN" ,"TOC"], na.rm = TRUE)
sampling[is.na(sampling$TOC) & sampling$Protocol == "ABMI", "TOC"] <- mean(sampling[sampling$Protocol == "ABMI" ,"TOC"], na.rm = TRUE)

##check for duplicates
sampling %>%
    dplyr::count(ABMI.SiteYear, Protocol) %>%
    dplyr::filter(n > 1) #no duplicates

#### calculate means ####
sampling.summary <- sampling %>%
    group_by(Protocol) %>%
    summarise(mean.sample.depth = mean(Sample_Depth_Mean, na.rm = TRUE),
              se.sample.depth = sd(Sample_Depth_Mean, na.rm = TRUE)/sqrt(length(Sample_Depth_Mean)),
              mean.time = mean(TOC, na.rm = TRUE),
              se.time = sd(TOC, na.rm = TRUE)/sqrt(length(TOC)))

#### Wilcox test (because TOC is not normal) ####
##check for normality
hist(sampling$TOC) #not normal

## reshape data into wide format
sampling.for.t.test <- sampling %>%
    dplyr::select(ABMI.SiteYear, Protocol, TOC) %>%
    pivot_wider(names_from = Protocol, values_from = TOC)
#run test
wilcox.test(sampling.for.t.test$ABMI, sampling.for.t.test$CABIN, paired = TRUE)
