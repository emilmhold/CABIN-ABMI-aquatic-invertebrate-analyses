## Clean and prep environmental data to re-run CABIN-ABMI comparisons
## Author: Emily H
## Created: July 15, 2025
## Last edited: August 21, 2025

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
## load species data
adj.count <- read_excel("output/CABIN Ambiguous parent resolved all taxa_APTC.abundant.xlsx") %>%
    rename(SiteYear = '...1') %>%
    mutate(Site = sub("_.*", "", SiteYear), .before = 1) %>%
    mutate(Year = sub(".*_", "", SiteYear), .before = 2) %>%
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 3) %>%
    filter(SiteYear != "W478-2_2018") %>% ## remove sites with no observations
    filter(SiteYear != "W638-2_2018") %>%
    filter(SiteYear != "CABIN-W478-2_2018") %>% ## remove sites corresponding cabin rows
    filter(SiteYear != "CABIN-W638-2_2018")

## isolate base site to select relevant enviro data
cabin_sites <- adj.count %>%
    mutate(base_SiteYear = SiteYear %>% sub("^CABIN-", "", .)) %>%
    mutate(ABMI.SiteYear = base_SiteYear %>% sub("W","",.)) %>%
    separate(base_SiteYear, into = c("Site", "Year"), sep = "_") %>%
   dplyr::select(Site, Year, SiteYear, ABMI.SiteYear, Protocol)

#### import and clean environmental data ####
## read in physiochemistry data
physiochemistry <- read.csv("data/A_W04_Water_Physiochemistry_AllYear_AB.csv") %>%
    mutate(ABMI.SiteYear = paste0(ABMI.Site, "_", Year)) %>% #create identifying column
    right_join(cabin_sites, by = "ABMI.SiteYear", relationship = "many-to-many") %>% #join with cabin_sites df
    mutate(Depth..metres. = as.numeric(Depth..metres.),
        Temperature..Degrees.Celsius. = as.numeric(Temperature..Degrees.Celsius.),
        Dissolved.Oxygen..milligrams.Litre. = as.numeric(Dissolved.Oxygen..milligrams.Litre.),
        Salinity..parts.per.thousand. = as.numeric(Salinity..parts.per.thousand.),
        pH = as.numeric(pH),
        Dissolved.Organic.Carbon..milligrams.Litre. = as.numeric(Dissolved.Organic.Carbon..milligrams.Litre.)) %>% # make sure variables are numeric
    group_by(ABMI.SiteYear) %>%
    summarise(Sample_Depth_Mean = mean(Depth..metres., na.rm = TRUE),
              Deep_Samples = sum(Location == "Deepest Point", na.rm = TRUE),
              Temp_Mean = mean(Temperature..Degrees.Celsius., na.rm = TRUE),
              DO_Mean = mean(Dissolved.Oxygen..milligrams.Litre., na.rm = TRUE),
              Sal_Mean = mean(Salinity..parts.per.thousand., na.rm = TRUE),
              pH_Mean = mean(pH, na.rm = TRUE),
              DOC = mean(Dissolved.Organic.Carbon..milligrams.Litre., na.rm = TRUE),
              Max_Depth = max(Depth..metres., na.rm = TRUE),
              Conductivity_Mean = mean(Salinity..parts.per.thousand., na.rm = TRUE)) %>% #summarise to match Brandon's values
    dplyr::select("ABMI.SiteYear", "Sample_Depth_Mean", "Deep_Samples",
                                           "Temp_Mean", "DO_Mean",
                                           "Sal_Mean", "pH_Mean", "DOC", "Max_Depth", "Conductivity_Mean")%>%
    mutate(Year = sub(".*_", "", ABMI.SiteYear), .before = 1) %>% #add year column
    mutate(Site = paste0("W", sub("_.*", "", ABMI.SiteYear)), .before = 1) %>% #add site column
    mutate(SiteYear = paste0(Site, "_", Year))  %>%
    mutate(Site = sub("_.*", "", SiteYear),
           Year = sub(".*_", "", SiteYear),
           IsDuplicate = grepl("-D", Site),
           BaseSite = gsub("-D", "", Site),
           BaseSiteYear = paste0(BaseSite, "_", Year))

# Join base site values to duplicate rows
physiochemistry_imputed <- physiochemistry %>%
    left_join(
        physiochemistry %>%
            filter(!IsDuplicate) %>%
            select(-SiteYear, -IsDuplicate, -Site, -Year, -BaseSite) %>%
            rename_with(~ paste0(.x, "_base"), -BaseSiteYear),
        by = c("BaseSiteYear" = "BaseSiteYear")
    ) %>%
    mutate(across(
        c(Sample_Depth_Mean, Deep_Samples, Temp_Mean, DO_Mean,
          Sal_Mean, pH_Mean, DOC, Max_Depth, Conductivity_Mean),
        ~ if_else(IsDuplicate, get(paste0(cur_column(), "_base")), .x)
    )) %>%
    select(-ends_with("_base"), -IsDuplicate, -BaseSite, -BaseSiteYear)

# Restore rownames
physiochemistry_imputed <- physiochemistry_imputed %>%
    column_to_rownames("SiteYear")

##assess if data are missing
apply(physiochemistry_imputed, 2, function(x) table(is.na(x)))

## replace missing values with mean values
physiochemistry_imputed[is.na(physiochemistry_imputed$Temp_Mean), "Temp_Mean"] <- mean(physiochemistry_imputed$Temp_Mean, na.rm = TRUE)
physiochemistry_imputed[is.na(physiochemistry_imputed$DO_Mean), "DO_Mean"] <- mean(physiochemistry_imputed$DO_Mean, na.rm = TRUE)
physiochemistry_imputed[is.na(physiochemistry_imputed$Sal_Mean), "Sal_Mean"] <- mean(physiochemistry_imputed$Sal_Mean, na.rm = TRUE)
physiochemistry_imputed[is.na(physiochemistry_imputed$pH_Mean), "pH_Mean"] <- mean(physiochemistry_imputed$pH_Mean, na.rm = TRUE)
physiochemistry_imputed[is.na(physiochemistry_imputed$Conductivity_Mean), "Conductivity_Mean"] <- mean(physiochemistry_imputed$Conductivity_Mean, na.rm = TRUE)

#### import and clean HFI data ####
load("data/abmi-terrestrial-sites_2021HFI.Rdata") #loads in d.wide.qa
load("data/abmi-wetland-sites-simplified_2007_2023.Rdata")

# Helper function to convert a dgCMatrix to a tidy data.frame
sparse_to_df <- function(mat, matrix_name = "value") {
    triplet <- as(mat, "dgTMatrix")
    data.frame(
        ID = rownames(mat)[triplet@i + 1],
        Variable = colnames(mat)[triplet@j + 1],
        Value = triplet@x,
        Source = matrix_name
    )
}
veg.current.df <- as.data.frame(d.wide.250$veg.current)
str(veg.current.df)

HFI <- veg.current.df %>%
    mutate(total.cover = rowSums(.)) %>% #calculate total cover
    mutate(total.HF = rowSums(across(c("Industrial",
        "Mine",
        "Rural",
        "EnSoftLin",
        "HardLin",
        "TrSoftLin",
        "EnSeismic",
        "Urban",
        "Wellsites")))) %>% #calculate total human footprint
    mutate(Human_Footprint = total.HF/total.cover) %>%
    dplyr::select(Human_Footprint) %>%
    rownames_to_column(var = "ABMI.SiteYear") %>%
    mutate(ABMI.SiteYear = ABMI.SiteYear %>% sub("W","",.)) %>% #remove "W" prefix to facilitate joining
    inner_join(cabin_sites, by = "ABMI.SiteYear") %>% #select CABIN sites
    dplyr::select(ABMI.SiteYear, Human_Footprint)
str(HFI)

#### import and wetland area data ####
open.water <- read.csv("data/CABIN wetland areas.csv") %>%
    rename(Site = ABMISite,
           open.water.area = Area..m2.) %>%
    dplyr::select(Site, open.water.area)
str(open.water)

#### join environmental data ####
site.information <- physiochemistry_imputed %>%
    full_join(HFI, by = "ABMI.SiteYear") %>%
    full_join(open.water, by = "Site")

#### duplicate environmental data for CABIN sites (because the physiochemistry sampling is the same for both) ####
cabin.site.information <- site.information %>%
    mutate(Site = paste0("CABIN-",Site))
##append to bottom of df
site.information <- rbind(site.information, cabin.site.information) %>%
    mutate(SiteYear = paste0(Site, "_", Year), .after = Year) %>%
    dplyr::select(!ABMI.SiteYear)

#### export data ####
write_rds(site.information, "output/cleaned site environmental data.rds")

#### make site summary tables ####
## import lat/long data
cabin.coordinates <- read_csv("data/A_W01A_Physical_Characteristics_AllYear_AB.csv") %>%
    dplyr::select("ABMI Site", "Public Latitude", "Public Longitude") %>%
    rename(ABMI.Site = "ABMI Site",
           Public.Latitude = "Public Latitude",
           Public.Longitude = "Public Longitude") %>%
    mutate(Site = paste0("W", ABMI.Site)) %>%
    dplyr::select(Site, Public.Latitude, Public.Longitude) %>%
    distinct(Site, .keep_all = TRUE) %>%
    right_join(cabin_sites, by = "Site") %>% #join with cabin_sites df
    dplyr::select(SiteYear, Protocol, Public.Latitude, Public.Longitude)

##join dfs
site.summary.table <- cabin.coordinates %>%
    left_join(site.information, by = "SiteYear", relationship = "many-to-many") %>% #join with cabin_sites df
    filter(Protocol != "CABIN") %>%
    dplyr::select(Site, Year, Public.Latitude, Public.Longitude, open.water.area, Sample_Depth_Mean, Max_Depth, Temp_Mean, DO_Mean, pH_Mean, Conductivity_Mean, Human_Footprint) %>%
    distinct(Site, .keep_all = TRUE) %>%
    filter(!grepl("-D$", Site)) #remove duplicate samples
str(site.summary.table)

##export wetland characteristics table
wetland.characteristics.table <- site.summary.table %>%
    rename('Site Code' = Site,
           'Public Latitude' = Public.Latitude,
           'Public Longitude' = Public.Longitude,
           'Open Water Area (ha)' = open.water.area,
           'Year Sampled' = Year,
           'Average Depth (m)' = Sample_Depth_Mean,
           'Maximum Depth (m)' = Max_Depth,
           'Temperature (C)' = Temp_Mean,
           'Dissolved Oxygen (mg/L)' = DO_Mean,
           pH = pH_Mean,
           'Specific Conductivity (uS/cm)' = Conductivity_Mean,
           'Human Footprint (%)' = Human_Footprint) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(wetland.characteristics.table, "output/CABIN Site general wetland characteristics.csv", row.names = FALSE)

##export mean characteristics table
wetland.means.table <- site.summary.table %>%
    summarize(mean.area = mean(open.water.area),
              se.mean.area = sd(open.water.area)/sqrt(length(open.water.area)),
        mean.depth = mean(Sample_Depth_Mean),
              se.mean.depth = sd(Sample_Depth_Mean)/sqrt(length(Sample_Depth_Mean)),
              mean.max.depth = mean(Max_Depth),
              se.max.depth = sd(Max_Depth)/sqrt(length(Max_Depth)),
              mean.temp = mean(Temp_Mean),
              se.temp = sd(Temp_Mean)/sqrt(length(Temp_Mean)),
              mean.do = mean(DO_Mean),
              se.do = sd(DO_Mean)/sqrt(length(DO_Mean)),
              mean.ph = mean(pH_Mean),
              se.ph = sd(pH_Mean)/sqrt(length(pH_Mean)),
              mean.conductivity = mean(Conductivity_Mean),
              se.conductivity = sd(Conductivity_Mean)/sqrt(length(Conductivity_Mean)),
              mean.hf = mean(Human_Footprint),
              se.hf = sd(Human_Footprint)/sqrt(length(Human_Footprint))
              ) %>%
    mutate(across(where(is.numeric), round, 2)) %>%
mutate('Open Water Area (ha)' = paste(mean.area, "±", se.mean.area),
    'Average Depth (m)' = paste(mean.depth, "±", se.mean.depth),
           'Maximum Depth (m)' = paste(mean.max.depth, "±", se.max.depth),
           'Temperature (C)' = paste(mean.temp, "±", se.temp),
           'Dissolved Oxygen (mg/L)' = paste(mean.do, "±", se.do),
           pH = paste(mean.ph, "±", se.ph),
           'Specific Conductivity (uS/cm)' = paste(mean.conductivity, "±", se.conductivity),
           'Human Footprint (%)' = paste(mean.hf, "±", se.hf)
           ) %>%
    dplyr:: select('Open Water Area (ha)', 'Average Depth (m)', 'Maximum Depth (m)', 'Temperature (C)', 'Dissolved Oxygen (mg/L)',
                   pH, 'Specific Conductivity (uS/cm)', 'Human Footprint (%)')
write.csv(wetland.means.table, "output/mean abmi wetland characteristics table.csv", row.names = FALSE)

## summary descriptions
physiochemistry.raw <- read.csv("data/A_W04_Water_Physiochemistry_AllYear_AB.csv") %>%
    mutate(ABMI.SiteYear = paste0(ABMI.Site, "_", Year)) %>% #create identifying column
    right_join(cabin_sites, by = "ABMI.SiteYear", relationship = "many-to-many") %>% #join with cabin_sites df
    mutate(Depth..metres. = as.numeric(Depth..metres.))
str(physiochemistry.raw)

max(physiochemistry.raw$Depth..metres., na.rm = TRUE)
min(physiochemistry.raw$Depth..metres., na.rm = TRUE)

