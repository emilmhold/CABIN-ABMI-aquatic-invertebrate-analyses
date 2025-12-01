#
# Title: Make map for aquatic invertebrate addendum report
# Created: November 10, 2025
# Last Updated by Emily: November 10, 2025
# Authors: Emily Holden
# Objective: recreate a map from the initial report showing the sites sampled in this work
# and their distribution.
#

#########
# Notes #
#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The original 2018 data included the following datasets:

# Clear memory
rm(list=ls())
gc()

# Install necessary packages. Packages only need to be installed once.
#install.packages("remotes")
#remotes::install_github("mabecker89/abmi.themes")
#install.packages("tidyverse")
#install.packages("usdm")
#install.packages("vegan")
#install.packages("readr")
#install.packages("readxl")
#install.packages("ggspatial")

# Load libraries into memory
library(abmi.themes)
library(tidyverse)
library(usdm)
library(vegan)
library(readr)
library(readxl)
library(sf)
library(ggspatial)

####################
# Site information #
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CHANGE DIRECTORY
setwd("~/GitHub/New aquatic invertebrate analyses")

#### load data ####
# Load species data and isolate CABIN sites
cabin.sites <- read_excel("output/CABIN Ambiguous parent resolved all taxa_APTC.abundant.xlsx") %>%
    rename(SiteYear = '...1') %>%
    mutate(Site = sub("_.*", "", SiteYear), .before = 1) %>%
    mutate(Year = sub(".*_", "", SiteYear), .before = 2) %>%
    filter(SiteYear != "W478-2_2018") %>% ## remove sites with no observations
    filter(SiteYear != "W638-2_2018") %>%
    filter(SiteYear != "CABIN-W478-2_2018") %>% ## remove sites corresponding cabin rows
    filter(SiteYear != "CABIN-W638-2_2018") %>%
    filter(!str_starts(Site, "CABIN")) %>% # remove CABIN sites (b/c they're duplicates of ABMI sites)
    dplyr::select(Site, Year) %>%
    mutate(Site = str_remove(Site, "^W")) # remove W from the start of site names.

## count duplicate sampling sites
cabin.sites %>%
    filter(str_ends(Site, "-D")) %>%
    summarise(n = n())

##load provincial boundary
load(file = "data/provincial-boundary.Rdata")

##load public site template
load(file = "data/site-template.Rdata")

#### make map ####
# --- Ensure matching site names ---
occ.template <- occ.template %>%
    mutate(Site = as.character(trimws(Site)))

cabin.sites <- cabin.sites %>%
    mutate(Site = as.character(trimws(Site)))

# --- Identify CABIN vs other sites ---
occ.template <- occ.template %>%
    mutate(
        SiteType = ifelse(Site %in% cabin.sites$Site, "CABIN site", "ABMI site")
    )
table(occ.template$SiteType)

# create a variable that identifies CABIN vs non-CABIN sites
occ.template <- occ.template %>%
    dplyr::mutate(
        SiteType = ifelse(Site %in% cabin.sites$Site, "CABIN site", "ABMI site")
    )

map <- ggplot() +
    geom_sf(data = province.shapefile, aes(fill = NRNAME), color = "black") +
    scale_fill_manual(
        values = setNames(province.shapefile$Color, province.shapefile$NRNAME)
    ) +
    geom_sf(
        data = occ.template,
        aes(color = SiteType, size = SiteType)
    ) +
    scale_color_manual(
        values = c("CABIN site" = "blue", "ABMI site" = "black"),
        name = "Site Type"
    ) +
    scale_size_manual(
        values = c("CABIN site" = 1.5, "ABMI site" = 0.5),
        name = "Site Type"
    ) +
    annotation_scale(
        location = "tr",
        width_hint = 0.25,
        text_cex = 0.8,
        bar_cols = c("grey60", "white"),
        line_width = 0.8,
        height = unit(0.25, "cm"),
        unit_category = "metric",
        text_family = "sans"
    ) +
    annotation_north_arrow(
        location = "tr",
        which_north = "true",
        pad_x = unit(0, "in"),
        pad_y = unit(0.6, "in"),
        style = north_arrow_fancy_orienteering
    ) +
    theme_minimal() +
    labs(fill = "Natural Region") +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
    )
map

## save
ggsave( filename = "output/sampling map.jpeg",
        plot = map,
        height = 6,          # in inches
        width = 8,          # in inches
        dpi = 300,           # high resolution
        units = "in",
        quality = 100)

#### create table tallying sites in each NR ####
NR.summary.table <- occ.template %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(province.shapefile)) %>%  # convert sites to sf
    st_join(province.shapefile %>% select(NRNAME)) %>%                                   # spatial join to NR polygons
    st_drop_geometry() %>%                                                               # remove geometry for table operations
    filter(!is.na(NRNAME)) %>%                                                           # remove sites outside polygons
    group_by(NRNAME, SiteType) %>%                                                       # group by NR and SiteType
    summarise(SiteCount = n(), .groups = "drop")                                         # count sites
print(NR.summary.table)
