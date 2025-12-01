#
# Title: ABMI Invertebrate Species Groups Analysis
# Created: October 1, 2025
# Last Updated by Emily: November 3, 2025
# Authors: Brandon Allen and Emily Holden
# Objective: Perform a series of analyses that match and expand on the Hanisch et al (2020) manuscript comparing the ABMI and CABIN protocols
# Keywords: Notes, Initialization, Site information, Species level
#

#########
# Notes #
#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The original 2018 data included the following datasets:
# Adj.abund - At each site, there is subsampling performed which results in variable raw counts. Adjusted counts are a result of taking the raw numbers and adjusted to 100 boxes.
# Relative abund - Relative abundance from the raw sorting data. Will want to consider which measure we use for each analysis so we treat the abundance/composition patterns appropriately.
# CPUE-t - Total number in the sample divided by the number of minutes spent collecting in the field.
# CPUE-b - Total number in the sample divided by the number of bottles collected at the site. This is equivalent to the sample weight (each bottle weights the same).
# We should try to use the raw counts (need to know the number of counted cells). Create into a new data frame.
# There are two different subsamples (all species vs Chiro). Make sure to account for this.

# We will run analyses for both adjusted abundance and relative abundance
# Add analysis of figure 3 found in the Hannich paper.
# We could also consider a species indicator status measure (which protocol captures which species)

##################
# Initialization #
##################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# Load libraries into memory
library(abmi.themes)
library(tidyverse)
library(usdm)
library(vegan)
library(readr)
library(readxl)

# Define function for making multi-plots in ggplot
# Define functions
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

####################
# Site information #
####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load site data
# CHANGE DIRECTORY
setwd("~/GitHub/New aquatic invertebrate analyses")

## read in environmental data
site.information <- read_rds("output/cleaned site environmental data.rds")
str(site.information)

################
# Missing data #
################

# Assess degree of missing data
apply(site.information, 2, function(x) table(is.na(x)))
    ##there are NaNs for one site because they were listed as "DNC" in the original data file
site.information <- site.information %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
## replace missing values with mean values
site.information[is.na(site.information$Human_Footprint), "Human_Footprint"] <- mean(site.information$Human_Footprint, na.rm = TRUE)

site.information$Human_Footprint <- unname(site.information$Human_Footprint) #remove names so vif will run
site.information <- as.data.frame(site.information) %>% #vifstep doesn't like tibbles
    distinct(SiteYear, .keep_all = TRUE)  %>% # drop duplicate observations
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 2)
str(site.information)

# Check there are no NA values
table(is.na(site.information))

#################################
# Correlation between variables #
#################################

# We are going to assess correlations using Variance Inflation Factor and pearson correlation coefficients.

# VIF approach (threshold of 5)
vifstep(site.information[, -c(1:4)], th = 5) # Will change the variable set

# This approach identifies 19 variables as colinear.
# suggests the following variable set:
# Variables      VIF
# 1       MSP 2.452654
# 2      DD18 3.479116
# 3      bFFP 2.837100
# 4       PAS 3.025385
# 5       EMT 1.629182
# 6       MAR 1.999870
# 7        RH 2.158307


# This variable set seems reasonable
    ## Em's note: I'm going to use the variables Brandon initially included for comparibility.
cor.matrix <- cor(x = site.information[, c( "Temp_Mean", "DO_Mean",
                                           "Sal_Mean", "pH_Mean", "DOC", "Max_Depth",
                                            "Human_Footprint")], method = "pearson")
                                           ##"Open_Water","TOC", "Sample_Depth_Mean", "Deep_Samples", "Hab_Complex", "DP" ## note: these are variables I don't have access to in 2022/2023

cor.matrix

# Subset the information to the final data set
site.information <- site.information[, c("Site", "SiteYear", "Temp_Mean", "DO_Mean",
                                          "Sal_Mean", "pH_Mean", "DOC", "Max_Depth",
                                           "Human_Footprint"
                                         #"TOC", "Sample_Depth_Mean","Deep_Samples","Hab_Complex", "DP", "Open_Water",
                                         )]
#rownames(site.information) <- paste(site.information$Site, site.information$Protocol, sep = "_")

rm(cor.matrix)

# We are going to leave the variables in their current state initially
# Check for normality of the remaining variables

 for (x in 3:9) {

             hist(site.information[, x], main = colnames(site.information)[x])
             hist(log(site.information[, x]), main = paste0("Log ", colnames(site.information)[x]))

 }
## log transformation makes most variables pretty normally distributed.
#

#################
#### Convert species to coarse groups ####
#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load species data
adj.count.spp <- read_excel("output/CABIN Ambiguous parent resolved all taxa_APTC.abundant.xlsx") %>%
    rename(SiteYear = '...1') %>%
    mutate(Site = sub("_.*", "", SiteYear), .before = 1) %>%
    mutate(Year = sub(".*_", "", SiteYear), .before = 2) %>%
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 3) %>%
    filter(SiteYear != "W478-2_2018") %>% ## remove sites with no observations
    filter(SiteYear != "W638-2_2018") %>%
    filter(SiteYear != "CABIN-W478-2_2018") %>% ## remove sites corresponding cabin rows
    filter(SiteYear != "CABIN-W638-2_2018")

taxa <- read_csv("data/Look-up Table_Invert Taxonomy_2025-07-28.csv") %>%
    dplyr::select(!Row)

# --- Step 1: Keep metadata ---
metadata <- adj.count.spp[, 1:4]

# --- Step 2: Map Analysis_Name to Coarse_Group ---
species_cols <- colnames(adj.count.spp)[5:262]
name_map <- setNames(taxa$Coarse_Group, taxa$Analysis_Name)

# Replace column names using the map (only if they exist in taxa)
new_names <- ifelse(species_cols %in% names(name_map),
                    name_map[species_cols],
                    species_cols)

colnames(adj.count.spp)[5:262] <- new_names

# --- Step 3: Collapse duplicate columns by prefix before "." ---
species_data <- adj.count.spp[, 5:ncol(adj.count.spp)]

# Now extract prefixes (these should be Coarse_Group names!)
prefixes <- sub("\\..*$", "", colnames(species_data))

# Collapse columns with the same prefix
species_summary <- as.data.frame(
    lapply(split.default(species_data, prefixes),
           function(df) rowSums(df, na.rm = TRUE))
)

# --- Step 4: Combine metadata + summarized species data ---
adj.count <- cbind(metadata, species_summary)

# #### find the number of duplicate sites ####
# adj.count %>%
#     filter(str_ends(Site, "-D")) %>%
#     nrow() #48 site

#### counts of inverts identified by protocol ####
##read in coarse data
str(adj.count)
invert.counts <- adj.count %>%
    rowwise() %>%  # sum across specimen columns for each row
    mutate(row_total = sum(c_across(5:35), na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(Protocol) %>%
    summarise(Total.Specimens = sum(row_total, na.rm = TRUE))

# Previous species data contained five data frames
# adj.count: Adjusted counts of species (columns = species, rows = sites)
# rel.abund: Relative abundance of species (columns = species, rows = sites)
# cpue.t: Count per unit effort (time) (columns = species, rows = sites)
# cpue.b: Count per unit effort (bottles) (columns = species, rows = sites)
# taxa: the short_code and long_codes for each species/group/family (columns) in the species data

###################
# Adjusted Counts #
###################

data.type <- "adjusted-counts"
data.in <- adj.count


# Align the species data with the site data
rownames(data.in) <- data.in$SiteYear
#data.in <- data.in[rownames(site.information), ] # Match site order

#
# Species Accumulation Curves
#

# ABMI
abmi.sac <- specaccum(comm = data.in[data.in$Protocol == "ABMI", -c(1,2,3,4)], method = "random")

# CABIN
cabin.sac <- specaccum(comm = data.in[data.in$Protocol == "CABIN", -c(1,2,3,4)], method = "random")

# Visualize
spp.accum <- data.frame(Protocol = c(rep("ABMI", length(abmi.sac$sites)), rep("CABIN", length(cabin.sac$sites))),
                        Site = c(abmi.sac$sites, cabin.sac$sites),
                        Richness = c(abmi.sac$richness, cabin.sac$richness),
                        SD = c(abmi.sac$sd, cabin.sac$sd))
spp.accum$Upper <- spp.accum$Richness + spp.accum$SD
spp.accum$Lower <- spp.accum$Richness - spp.accum$SD

png(filename = paste0("figures/invertebrate-protocol-analyses/coarse group/species-accumulation-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

ggplot(data = spp.accum) +
    geom_pointrange(aes(x = Site, y = Richness, ymin = Lower, ymax = Upper, color = Protocol, fill = Protocol)) +
    scale_color_manual(values = abmi_pal("main")(2)) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = paste("Correlation = ", round(cor(spp.accum[spp.accum$Protocol == "ABMI", "Richness"], spp.accum[spp.accum$Protocol == "CABIN", "Richness"]), 3))) +
    #ylim(c(0,300)) +
    theme_bw()

dev.off()

# ##find the number of species observed in ABMI vs. CABIN
# cabin.spp <- data.in %>%
#     filter(Protocol == "CABIN") %>%
#     select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
# print(ncol(cabin.spp)-4) #199 species
#
# abmi.spp <- data.in %>%
#     filter(Protocol == "ABMI") %>%
#     select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
# print(ncol(abmi.spp)-4) #212 species
#
# General diversity characteristics
#

site.richness <- data.frame(Site = data.in$Site,
                            Protocol = data.in$Protocol,
                            Simpson = diversity(x = data.in[, -c(1:4)], index = "simpson", MARGIN = 1),
                            Shannon = diversity(x = data.in[, -c(1:4)], index = "shannon", MARGIN = 1),
                            Richness = specnumber(x = data.in[, -c(1:4)], MARGIN = 1))


richness.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Richness"], y = site.richness[site.richness$Protocol == "CABIN", "Richness"], paired = TRUE)
simpson.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Simpson"], y = site.richness[site.richness$Protocol == "CABIN", "Simpson"], paired = TRUE)
shannon.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Shannon"], y = site.richness[site.richness$Protocol == "CABIN", "Shannon"], paired = TRUE)

# Name change
#site.richness$Protocol <- gsub("ABMI", "ABMI", site.richness$Protocol)
#site.richness$Protocol <- gsub("CABIN", "CABIN", site.richness$Protocol)

png(filename = paste0("figures/invertebrate-protocol-analyses/coarse group/simpson-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Simpson, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Simpson Diversity; P-value = ", round(simpson.results$p.value, 4))) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = "Simpson Diversity; P-value = 0.007") +
    theme_bw()


site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Simpson")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Simpson")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = paste0("Simpson Diversity; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()

multiplot(fig.1, fig.2, cols = 2)

dev.off()

png(filename = paste0("figures/invertebrate-protocol-analyses/coarse group/richness-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Richness, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Richness; P-value = ", round(richness.results$p.value, 4))) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = "Richness; P-value = 0.032") +
    theme_bw()

site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Richness")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Richness")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = paste0("Richness; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()


multiplot(fig.1, fig.2, cols = 2)

dev.off()

png(filename = paste0("figures/invertebrate-protocol-analyses/coarse group/shannon-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Shannon, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Shannon Diversity; P-value = ", round(shannon.results$p.value, 2))) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = "Shannon Diversity; P-value = 0.009") +
    theme_bw()

site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Shannon")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Shannon")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "Taxonomically resolved data for all years",
         subtitle = paste0("Shannon Diversity; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()

multiplot(fig.1, fig.2, cols = 2)

dev.off()

#### Table 8: Abundance by taxa (Wilcox test) ####

##add base_SiteYear column to identify pairs for wilcox.test
wilcox.data <- data.in %>%
    mutate(base_SiteYear = SiteYear %>% sub("^CABIN-", "", .), .after = 1) %>%
    dplyr::select(!c(Site, Year))
rownames(wilcox.data) <- NULL
wilcox.data <- wilcox.data %>% column_to_rownames(var = "SiteYear")

site.abundance <- NULL
sign.results <- NULL

# for (x in 3:ncol(data.in)) {
#
#     temp.abundance <- data.frame(Protocol = data.in$Protocol,
#                                  Species = colnames(data.in)[x],
#                                  Abundance = data.in[, x])
#
#     # Wilcoxon Test
#     wilcox.results <- wilcox.test(temp.abundance$Abundance, temp.abundance$Protocol, paired = TRUE)
#     if(wilcox.results$p.value < 0.05) {
#
#         site.abundance <- rbind(site.abundance, temp.abundance)
#         sign.results <- rbind(sign.results, data.frame(Species = colnames(data.in)[x],
#                                                        Significance = data.in$p.value))
#
#     }
#
#     rm(temp.abundance, wilcox.results)
#
# }

site.abundance <- data.frame()
sign.results <- data.frame()

# Loop through each species column
for (x in 3:ncol(wilcox.data)) {

    species.name <- colnames(wilcox.data)[x]

    # Create long-format data for this species
    temp.abundance <- wilcox.data %>%
        select(base_SiteYear, Protocol, !!sym(species.name)) %>%
        rename(Abundance = !!sym(species.name))

    # Convert to wide format by Protocol
    temp.wide <- temp.abundance %>%
        pivot_wider(names_from = Protocol, values_from = Abundance)

    # Keep only rows with both CABIN and ABMI values
    temp.wide <- temp.wide %>% drop_na(CABIN, ABMI)

    # Skip if too few pairs
    if (nrow(temp.wide) < 2) next

    # Use non-formula method for paired test
    test <- wilcox.test(temp.wide$CABIN, temp.wide$ABMI, paired = TRUE)
    print(paste("Wilcoxon p for", species.name, ":", test$p.value))

    # If significant, store results
    if (!is.na(test$p.value) && test$p.value < 0.05) {
        temp.abundance <- temp.abundance %>%
            mutate(Species = species.name)

        site.abundance <- bind_rows(site.abundance, temp.abundance)

        sign.results <- bind_rows(sign.results,
                                  data.frame(Species = species.name,
                                             P_value = test$p.value))
    }
}



for(spp.id in unique(site.abundance$Species)) {

    temp.abundance <- site.abundance[site.abundance$Species == spp.id, ]
    #long.code <- taxa[taxa$Analysis_Name == spp.id, "long_code"]
    #temp.abundance$Species <- gsub(spp.id, long.code, temp.abundance$Species)

    # Name change
    # temp.abundance$Protocol <- gsub("ABMI", "ABMI", temp.abundance$Protocol)
    # temp.abundance$Protocol <- gsub("CABIN", "CABIN", temp.abundance$Protocol)

    png(filename = paste0("figures/invertebrate-protocol-analyses/coarse group/species-abundance", spp.id, "-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
        width = 1200,
        height = 1200,
        res = 300)

    print(ggplot(data = temp.abundance) +
              geom_boxplot(mapping = aes(x = Species, y = Abundance, fill = Protocol), show.legend = TRUE) +
              scale_fill_manual(values = abmi_pal("main")(2)) +
              ggtitle(paste0("Wilcoxon = ", round(sign.results[sign.results$Species == spp.id, "P_value"], 3))) +
              theme_bw())

    dev.off()
    message("Finished plotting for species: ", spp.id)

}

# Species by site abundance
# Only for species with sufficient detections using each protocol
# site.abundance <- NULL
#
# for (x in 5:ncol(data.in)) {
#
#     temp.abundance <- data.frame(Protocol = data.in$Protocol,
#                                  Species = colnames(data.in)[x],
#                                  Abundance = data.in[, x])
#     site.abundance <- rbind(site.abundance, temp.abundance)
#
#     rm(temp.abundance)
#
# }

# Initialize output with consistent column names
site.abundance <- data.frame(Protocol = character(),
                             Species = character(),
                             Abundance = numeric(),
                             stringsAsFactors = FALSE)

for (x in 5:ncol(data.in)) {

    temp.abundance <- data.frame(
        Protocol = data.in$Protocol,
        Species = rep(colnames(data.in)[x], nrow(data.in)),
        Abundance = data.in[[x]]
    )

    site.abundance <- rbind(site.abundance, temp.abundance)

    rm(temp.abundance)
}

##summarize site abundance data
site.abundance.summary <- site.abundance %>%
    group_by(Protocol, Species) %>%
    summarize(Mean = mean(Abundance, na.rm = TRUE),
              SD = sd(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Protocol, values_from = c(Mean, SD), names_sep ="_")

#join with Wilcox results
wilcox.results <- sign.results %>%
    left_join(site.abundance.summary, by = "Species")

#reformat to export
table8 <- wilcox.results %>%
    mutate(ABMI = paste0(round(Mean_ABMI, 1), "(", round(SD_ABMI,1), ")")) %>%
    mutate(CABIN = paste0(round(Mean_CABIN,1), "(", round(SD_CABIN,1), ")")) %>%
    mutate(P_value = round(P_value, 3)) %>%
    dplyr::rename("P-value" = P_value) %>%
    dplyr::select(Species, CABIN, ABMI, "P-value")
write_csv(table8, "tables/coarse group/wilcox analyses for species detections by protocol.csv")

#### Table 9: Model species responses to abiotic data - need abiotic data to do this ####
# --- Ensure output folder exists ---
dir.create("figures/species/models", recursive = TRUE, showWarnings = FALSE)

# --- Initialize storage lists/tables ---
all.coef <- list()
all.cor <- data.frame(Species = character(), Correlation = numeric(), stringsAsFactors = FALSE)

# --- Add protocol column to site information ---
site.information <- site.information %>%
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 2)

# --- Loop over species ---
for (spp.id in unique(site.abundance$Species)) {

    # --- Step 1: Check species exists in data.in ---
    if (!(spp.id %in% colnames(data.in))) {
        warning(paste("Species", spp.id, "not found in data.in. Skipping."))
        next
    }

    # --- Step 2: Extract species column ---
    spp.col <- data.frame(
        SiteYear = rownames(data.in),
        Spp = data.in[[spp.id]]
    )

    # --- Step 3: Merge with site information and standardize environment variables ---
    site.information.scaled <- site.information
    env.cols <- 4:10  # adjust to your numeric environmental columns
    site.information.scaled[, env.cols] <- scale(site.information.scaled[, env.cols])

    abmi.data <- site.information.scaled %>%
        left_join(spp.col, by = "SiteYear")

    # --- Step 4: Skip if all abundance values are NA or zero ---
    if (all(is.na(abmi.data$Spp)) || sum(abmi.data$Spp > 0, na.rm = TRUE) <= 10) {
        next
    }

    # --- Step 5: Subset by protocol ---
    if (!"Protocol" %in% colnames(abmi.data)) {
        stop("Protocol column not found in site.information.scaled!")
    }

    cabin.data <- abmi.data %>% filter(Protocol == "CABIN")
    abmi.data  <- abmi.data %>% filter(Protocol == "ABMI")

    # Skip if not enough nonzero observations in either
    if (sum(abmi.data$Spp > 0) <= 10 || sum(cabin.data$Spp > 0) <= 10) {
        next
    }

    # --- Step 6: Abundance correlation ---
    site.cor <- data.frame(
        CTA = abmi.data$Spp,
        TSA = cabin.data$Spp
    )

    cor.value <- if (nrow(site.cor) > 1) cor(site.cor$CTA, site.cor$TSA) else NA
    all.cor <- rbind(all.cor, data.frame(Species = spp.id, Correlation = cor.value))

    fig.1 <- ggplot(site.cor) +
        geom_point(aes(x = CTA, y = TSA), color = abmi_pal("main")(1)) +
        geom_abline(intercept = 0, slope = 1) +
        ggtitle(paste0(spp.id, " Correlation = ", round(cor.value, 2))) +
        theme_bw()

    # --- Step 7: Modeling ---
    abmi.data$Spp <- scale(log(abmi.data$Spp + 0.0001))
    cabin.data$Spp <- scale(log(cabin.data$Spp + 0.0001))

    spp.model <- list(
        CTA = glm(Spp ~ Temp_Mean + DO_Mean + Sal_Mean + pH_Mean + DOC + Max_Depth + Human_Footprint,
                  family = "gaussian", data = abmi.data),
        TSA = glm(Spp ~ Temp_Mean + DO_Mean + Sal_Mean + pH_Mean + DOC + Max_Depth + Human_Footprint,
                  family = "gaussian", data = cabin.data)
    )

    # --- Step 8: Extract coefficients ---
    coef.store <- do.call(rbind, lapply(names(spp.model), function(model.name) {
        temp <- summary(spp.model[[model.name]])$coefficients
        data.frame(
            Protocol = model.name,
            Variable = rownames(temp),
            Estimate = temp[, 1],
            StdErr = temp[, 2],
            Pvalue = temp[, 4],
            stringsAsFactors = FALSE
        )
    }))

    # Remove intercept
    coef.store <- coef.store[coef.store$Variable != "(Intercept)", ]

    # Remove extreme StdErrs
    coef.store[max(abs(coef.store$Estimate)) * 10 < abs(coef.store$StdErr), c("Estimate", "StdErr")] <- NA

    # Add plotting columns
    coef.store$LowInner <- coef.store$Estimate - coef.store$StdErr
    coef.store$HighInner <- coef.store$Estimate + coef.store$StdErr
    coef.store$Coefficient <- factor(coef.store$Variable,
                                     levels = c("Deep_Samples", "Hab_Complex", "Temp_Mean", "DO_Mean", "Sal_Mean", "pH_Mean",
                                                "DP", "DOC", "Max_Depth", "Open_Water", "Human_Footprint"))

    # --- Step 9: Save coefficient table ---
    all.coef[[spp.id]] <- coef.store
    write.csv(coef.store,
              file = paste0("figures/species/models/", spp.id, "-coefficients_", Sys.Date(), ".csv"),
              row.names = FALSE)

    # --- Step 10: Coefficient plot ---
    fig.2 <- ggplot(coef.store, aes(x = Coefficient, y = Estimate, color = Protocol, shape = Protocol)) +
        geom_pointrange(aes(ymin = LowInner, ymax = HighInner),
                        position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        coord_flip() +
        theme_bw() +
        ggtitle(paste("CTA vs TSA:", spp.id))

    # --- Step 11: Save figures ---
    png(filename = paste0("figures/species/models/",
                          spp.id, "-model-", data.type, "-species-level_", Sys.Date(), ".png"),
        height = 2400, width = 4800, res = 300)

    grid.arrange(fig.1, fig.2, ncol = 2)
    dev.off()
}

# --- Step 12: Save all correlations ---
write.csv(all.cor,
          file = paste0("figures/species/models/species_correlations_", Sys.Date(), ".csv"),
          row.names = FALSE)

#### find number of species with overlapping coefficient estimates ####
# Initialize a table to store results
variables <- unique(unlist(lapply(all.coef, function(df) df$Variable)))
overlap.count <- data.frame(Variable = variables, SpeciesWithOverlap = 0, stringsAsFactors = FALSE)

# Loop over each variable
for (var in variables) {

    count <- 0

    for (spp.id in names(all.coef)) {
        df <- all.coef[[spp.id]]

        # Extract ABMI and CABIN rows for this variable
        abmi.row <- df %>% filter(Variable == var & Protocol == "CTA")
        cabin.row <- df %>% filter(Variable == var & Protocol == "TSA")

        # Skip if either row is missing
        if (nrow(abmi.row) == 0 || nrow(cabin.row) == 0) next

        # Compute intervals
        abmi.int <- c(abmi.row$Estimate - abmi.row$StdErr, abmi.row$Estimate + abmi.row$StdErr)
        cabin.int <- c(cabin.row$Estimate - cabin.row$StdErr, cabin.row$Estimate + cabin.row$StdErr)

        # Check if intervals overlap
        if (!(abmi.int[2] < cabin.int[1] || cabin.int[2] < abmi.int[1])) {
            count <- count + 1
        }
    }

    overlap.count$SpeciesWithOverlap[overlap.count$Variable == var] <- count
}

# View the table
overlap.count

##write table
write.csv(overlap.count, "tables/coarse group/overlap_counts_per_variable.csv", row.names = FALSE)


######################
# Relative Abundance #
######################
#### Multivariate stats ####

data.type <- "relative-abundance"
# data.in <- rel.abund
data.in <- decostand(adj.count[,5:35], method = "total") #convert to relative abundance
data.in <- cbind(adj.count[, 1:4], data.in)
sum(is.na(data.in))
summary(data.in)

# Align the species data with the site data
rownames(data.in) <- paste(data.in$Site, data.in$Protocol, sep = "_")
#data.in <- data.in[rownames(site.information), ] # Match site order

#
# SIMPER
#

# The SIMPER function takes a community matrix (data.in[, -c(1,2)]) and a grouping factor (Protocol)
# and converts the matrix into a brays-curtis dissimilarity. Columns in the community matrix (e.g., species, familes, etc)
# are listed in order of highest to lowest contribution.

simper.results <- simper(comm = data.in[, -c(1:4)], group = data.in$Protocol,
                         permutations = 1000,
                         trace = TRUE)

# Summary table of results for all columns in the community matrix.
# ChiroUA, OligoUA, AmphiUA, ChaobUA, GastrUA, and EphemUA are the top group
summary(simper.results)
write.csv(summary(simper.results)$ABMI_CABIN,
          file = paste0("output/simper-analysis-adjusted-counts-coarse-group-level_2025-10-02.csv", row.names = TRUE))


# #
# # RDA analyses
# #
#
# # The RDA analysis will help us assess which environmental covariates are associated with the communities collected
# # using both the ABMI and CABIN methods.
# # Consider transformations to the count based on the type of data.
# # Look into adding a conditioning call so we can pull out the protocol effect # Feb 24th, 2021 Note Look into this more.
#
# rda.results <- rda(X = data.in[, -c(1,2)], Y = site.information[, -c(1:2)], z = site.information$Protocol, scale = TRUE)
# rda.results <- summary(rda.results)
#
# # Visualize
# species.coord <- data.frame(RDA1 = rda.results$species[, 1],
#                             RDA2 = rda.results$species[, 2])
#
# species.subset.coord <- species.coord[rownames(summary(simper.results)$ABMI_CABIN)[1:5], ]# When adding labels, we are only adding the most import species identified in the SIMPER analysis. Figure gets muddied with too many labeles.
#
# site.coord <- data.frame(RDA1 = rda.results$sites[, 1],
#                          RDA2 = rda.results$sites[, 2],
#                          Protocol = site.information$Protocol)
#
# biplot.coord <- data.frame(RDA1 = rda.results$biplot[, 1],
#                            RDA2 = rda.results$biplot[, 2])
#
# png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/rda-", data.type, "-species-level_", Sys.Date(), ".png"),
#     width = 2400,
#     height = 2400,
#     res = 300)
#
# print(ggplot() +
#           geom_point(data = site.coord, aes(x = RDA1, y = RDA2, color = Protocol), size = 3) +
#           scale_color_manual(values = abmi_pal("main")(2)) +
#           geom_segment(data = species.coord, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
#                        arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),
#                                      type = "closed"),linetype = 1, size = 0.6, colour = "#E8A396") +
#           geom_text(data = species.subset.coord, aes(x = RDA1, y = RDA2, label = row.names(species.subset.coord))) +
#           geom_segment(data = biplot.coord, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
#                        arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),
#                                      type = "closed"),linetype = 1, size = 0.6,colour = "#829EBC") +
#           geom_text(data = biplot.coord, aes(x = RDA1, y = RDA2, label = row.names(biplot.coord))) +
#           labs(x = paste0("RDA 1 (", format(100 *rda.results$cont[[1]][2,1], digits=4), "%)"),
#                y = paste0("RDA 2 (", format(100 *rda.results$cont[[1]][2,2], digits=4), "%)")) +
#           geom_hline(yintercept = 0, linetype = 2,size=  1) +
#           geom_vline(xintercept = 0,linetype = 2,size = 1) +
#           guides(shape=guide_legend(title=NULL,color="black"),
#                  fill=guide_legend(title=NULL))+
#           theme_bw())
#
# dev.off()

#
# NMDS
#

# The metaMDS function takes a community matrix, transforms it based on the defined distance (brays-curtis), and a defined
# number axes to reduce the dimentionality. We are trying to simplify the data, so we don't want too many axes (e.g., greater than 3)
# but we don't want there to be high stress (poor ability to simplify the data (greater than 0.2 is poor). The try options are defining how many
# iterations the function should run for. Higher values take longer to run, but will result in more stable results between runs.

stress.values <- NULL
for (axis.size in 1:10) {

    stress.values <- c(stress.values, metaMDS(comm = data.in[, -c(1:4)], distance = "bray",
                                              k = axis.size, try = 100, trymax = 500)$stress)

}

stress.values <- data.frame(k = 1:10, Stress = stress.values)

# Create stress plot and identify smallest index with value less than 0.2
ggplot(data = stress.values, aes(x = k, y = Stress, color = "#829EBC")) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = 0.2) +
    theme_bw()

optimal.k <- as.numeric(table(stress.values$Stress < 0.2)["FALSE"]) + 1

# Optimal value of K (3)
nmds.results <- metaMDS(comm = data.in[, -c(1:4)], distance = "bray",
                        k = optimal.k, try = 100, trymax = 1000)

# Lets assess the stress of the nmds using a shepard plot
# It looks like there is some scatter between the ordination and dissimilarity distances, but it isn't too bad.
stressplot(nmds.results)

# Lets visualize the first two axes of the nmds

# store the site level information
data.scores <- as.data.frame(scores(nmds.results, display = "sites"))
#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Protocol <- data.in$Protocol #  add the protocol variable

# Name change
data.scores$Protocol <- gsub("ABMI", "ABMI", data.scores$Protocol)
data.scores$Protocol <- gsub("CABIN", "CABIN", data.scores$Protocol)

# Name change
data.in$Protocol <- gsub("ABMI", "ABMI", data.in$Protocol)
data.in$Protocol <- gsub("CABIN", "CABIN", data.in$Protocol)

# Store the ellipse information
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {

    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))

}

plot.new()
ordiplot(nmds.results)
ord <- ordiellipse(nmds.results, as.factor(data.in$Protocol), display = "sites", kind ="sd", conf = 0.9, label = FALSE)
dev.off()

#Generate ellipse points
ellipse.df <- data.frame()
for(g in unique(data.scores$Protocol)){
    if(g!="" && (g %in% names(ord))){

        ellipse.df <- rbind(ellipse.df, cbind(as.data.frame(with(data.scores[data.scores$Protocol==g,],
                                                                 veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                              ,Protocol=g))
    }
}

# Name change
ellipse.df$Protocol <- gsub("ABMI", "ABMI", ellipse.df$Protocol)
ellipse.df$Protocol <- gsub("CABIN", "CABIN", ellipse.df$Protocol)

# Store the hull information
grp.a <- data.scores[data.scores$Protocol == "ABMI", ][chull(data.scores[data.scores$Protocol ==
                                                                             "ABMI", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI
grp.b <- data.scores[data.scores$Protocol == "CABIN", ][chull(data.scores[data.scores$Protocol ==
                                                                              "CABIN", c("NMDS1", "NMDS2")]), ]  # hull values for CABIN
# Name change
grp.a$Protocol <- gsub("ABMI", "ABMI", grp.a$Protocol)
grp.b$Protocol <- gsub("CABIN", "CABIN", grp.b$Protocol)

hull.data <- rbind(grp.a, grp.b)  # combine the hull data

# Store the species level information
species.scores <- as.data.frame(scores(nmds.results, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  # look at the data
species.scores <- species.scores[rownames(summary(simper.results)$ABMI_CABIN)[1:10], ]

png(filename = "figures/invertebrate-protocol-analyses/coarse group/nmds-adjusted-counts-coarse-group-level_2025-10-01.png",
    width = 2400,
    height = 2400,
    res = 300)

# Looks like there is a lot of overlap between the two protocols
print(ggplot() +
          #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Protocol, group = Protocol), alpha = 0.30) + # add the convex hulls
          geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species), alpha = 0.5) +  # add the species labels
          geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, shape = Protocol, colour = Protocol), size = 4) + # add the point markers
          geom_path(data = ellipse.df, aes(x = NMDS1, y = NMDS2, group = Protocol, colour = Protocol)) +
          scale_color_manual(values = abmi_pal("main")(2)) +
          scale_fill_manual(values = abmi_pal("main")(2)) +
          labs(title = "Taxonomically resolved data for all years") +
          coord_equal() +
          theme_bw())

dev.off()

# Permanova
site.information <- site.information %>%
    mutate(Protocol = if_else(str_starts(Site, "CABIN"),"CABIN", "ABMI")) %>%
    filter(SiteYear %in% data.in$SiteYear) %>%
    distinct(SiteYear, .keep_all = TRUE)
str(site.information)
perma.results <- adonis2(data.in[, -c(1:4)] ~ Protocol, data = site.information, method = "bray", permutations = 999)
write.csv(perma.results, file = "output/coarse group/permanova results.csv", row.names = TRUE) # No difference in their dispersion
#rm(perma.results)

# Analysis of dispersion (permadis analysis) Which one is more dispersed
beta.results <- betadisper(d = vegdist(x = data.in[, -c(1:4)], method = "bray"), group = data.in$Protocol)
anova(beta.results)
write.csv(anova(beta.results), file = paste0("tables/coarse group/permadisp-analysis-", data.type, "-coarse-group-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion


#
# Procrustes test
#

# The procrustes test assesses the configuration of two matrices to maximize similarity
# I believe you can use the original matrix, but it is more informative to test on the ordination results (NMDS)
# We can take the ordination scores from the data.scores object created in the previous step

pt.visual <- procrustes(Y = data.scores[data.scores$Protocol == "ABMI", 1:2],
                        X = data.scores[data.scores$Protocol == "CABIN", 1:2], symmetric = TRUE)
pt.test <- protest(Y = data.scores[data.scores$Protocol == "ABMI", 1:2],
                   X = data.scores[data.scores$Protocol == "CABIN", 1:2], symmetric = TRUE)

# # Distance between the two plots
residuals(pt.test)
plot(residuals(pt.test))
# lm(residuals(pt.test) ~ covariates)# FIX Use only the measures that are consistent between protocols # NEED TO UPDATE

write.csv(data.frame(SS = pt.test$ss,
                     Correlation = pt.test$scale,
                     Significance = pt.test$signif), file = paste0("output/coarse group/procrustes-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = FALSE) # No difference in their dispersion

# Change plot

before.data <- data.frame(pt.visual$X)
colnames(before.data) <- colnames(pt.visual$X)
after.data <- data.frame(pt.visual$Yrot)
colnames(after.data) <- colnames(pt.visual$X)

png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/procrustes-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

print(ggplot() +
          geom_point(data = after.data, aes(x = NMDS1, y = NMDS2, color = "ABMI"), size = 3, show.legend = FALSE) + # add the point markers
          geom_segment(aes(x = after.data$NMDS1, y = after.data$NMDS2, xend = before.data$NMDS1, yend = before.data$NMDS2), color = abmi_pal("main")(2)[2], size = 0.5, lty = 1, arrow = arrow(length=unit(0.30,"cm"), type = "closed"), show.legend = FALSE) + # add the point markers
          scale_color_manual(values = abmi_pal("main")(1)) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          xlab("Axis 1") +
          ylab("Axis 2") +
          coord_equal() +
          theme_bw())

dev.off()

#### 2018 data only ####
## filter data sets
adj.count <- adj.count %>%
    filter(Year == "2018")
site.information <- site.information %>%
    filter(str_ends(SiteYear, "2018"))
###################
# Adjusted Counts #
###################

data.type <- "adjusted-counts"
data.in <- adj.count

# Align the species data with the site data
rownames(data.in) <- data.in$SiteYear
#data.in <- data.in[rownames(site.information), ] # Match site order

#
# Species Accumulation Curves
#

# ABMI
abmi.sac <- specaccum(comm = data.in[data.in$Protocol == "ABMI", -c(1,2,3,4)], method = "random")

# CABIN
cabin.sac <- specaccum(comm = data.in[data.in$Protocol == "CABIN", -c(1,2,3,4)], method = "random")

# Visualize
spp.accum <- data.frame(Protocol = c(rep("ABMI", length(abmi.sac$sites)), rep("CABIN", length(cabin.sac$sites))),
                        Site = c(abmi.sac$sites, cabin.sac$sites),
                        Richness = c(abmi.sac$richness, cabin.sac$richness),
                        SD = c(abmi.sac$sd, cabin.sac$sd))
spp.accum$Upper <- spp.accum$Richness + spp.accum$SD
spp.accum$Lower <- spp.accum$Richness - spp.accum$SD

png(filename = paste0("figures/coarse group/2018 data/species-accumulation-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

ggplot(data = spp.accum) +
    geom_pointrange(aes(x = Site, y = Richness, ymin = Lower, ymax = Upper, color = Protocol, fill = Protocol)) +
    scale_color_manual(values = abmi_pal("main")(2)) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "2018 data",
         subtitle = paste("Correlation = ", round(cor(spp.accum[spp.accum$Protocol == "ABMI", "Richness"], spp.accum[spp.accum$Protocol == "CABIN", "Richness"]), 3))) +
    #ylim(c(0,300)) +
    theme_bw()

dev.off()

# ##find the number of species observed in ABMI vs. CABIN
# cabin.spp <- data.in %>%
#     filter(Protocol == "CABIN") %>%
#     select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
# print(ncol(cabin.spp)-4) #199 species
#
# abmi.spp <- data.in %>%
#     filter(Protocol == "ABMI") %>%
#     select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
# print(ncol(abmi.spp)-4) #212 species
#
# General diversity characteristics
#

site.richness <- data.frame(Site = data.in$Site,
                            Protocol = data.in$Protocol,
                            Simpson = diversity(x = data.in[, -c(1:4)], index = "simpson", MARGIN = 1),
                            Shannon = diversity(x = data.in[, -c(1:4)], index = "shannon", MARGIN = 1),
                            Richness = specnumber(x = data.in[, -c(1:4)], MARGIN = 1))


richness.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Richness"], y = site.richness[site.richness$Protocol == "CABIN", "Richness"], paired = TRUE)
simpson.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Simpson"], y = site.richness[site.richness$Protocol == "CABIN", "Simpson"], paired = TRUE)
shannon.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Shannon"], y = site.richness[site.richness$Protocol == "CABIN", "Shannon"], paired = TRUE)

# Name change
#site.richness$Protocol <- gsub("ABMI", "ABMI", site.richness$Protocol)
#site.richness$Protocol <- gsub("CABIN", "CABIN", site.richness$Protocol)

png(filename = paste0("figures/coarse group/2018 data/simpson-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Simpson, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Simpson Diversity; P-value = ", round(simpson.results$p.value, 4))) +
    labs(title = "2018 data",
         subtitle = "Simpson Diversity; P-value < 0.001") +
    theme_bw()


site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Simpson")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Simpson")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "2018 data",
         subtitle = paste0("Simpson Diversity; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()

multiplot(fig.1, fig.2, cols = 2)

dev.off()

png(filename = paste0("figures/coarse group/2018 data/richness-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Richness, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Richness; P-value = ", round(richness.results$p.value, 4))) +
    labs(title = "2018 data",
         subtitle = "Richness; P-value = 0.016") +
    theme_bw()

site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Richness")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Richness")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "2018 data",
         subtitle = paste0("Richness; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()


multiplot(fig.1, fig.2, cols = 2)

dev.off()

png(filename = paste0("figures/coarse group/2018 data/shannon-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Shannon, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Shannon Diversity; P-value = ", round(shannon.results$p.value, 2))) +
    labs(title = "2018 data",
         subtitle = "Shannon Diversity; P-value < 0.001") +
    theme_bw()

site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Shannon")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Shannon")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "2018 data",
         subtitle = paste0("Shannon Diversity; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()

multiplot(fig.1, fig.2, cols = 2)

dev.off()

######################
# Relative Abundance #
######################
#### Multivariate stats ####

data.type <- "relative-abundance"
# data.in <- rel.abund
data.in <- decostand(adj.count[,5:35], method = "total") #convert to relative abundance
data.in <- cbind(adj.count[, 1:4], data.in)
sum(is.na(data.in))
summary(data.in)

# Align the species data with the site data
rownames(data.in) <- paste(data.in$Site, data.in$Protocol, sep = "_")
#data.in <- data.in[rownames(site.information), ] # Match site order

# #
# # SIMPER
# #
#
# # The SIMPER function takes a community matrix (data.in[, -c(1,2)]) and a grouping factor (Protocol)
# # and converts the matrix into a brays-curtis dissimilarity. Columns in the community matrix (e.g., species, familes, etc)
# # are listed in order of highest to lowest contribution.
#
# simper.results <- simper(comm = data.in[, -c(1:4)], group = data.in$Protocol,
#                          permutations = 1000,
#                          trace = TRUE)
#
# # Summary table of results for all columns in the community matrix.
# # ChiroUA, OligoUA, AmphiUA, ChaobUA, GastrUA, and EphemUA are the top group
# summary(simper.results)
# write.csv(summary(simper.results)$ABMI_CABIN,
#           file = paste0("output/simper-analysis-adjusted-counts-coarse-group-level_2025-10-02.csv", row.names = TRUE))
#
#
# # #
# # # RDA analyses
# # #
# #
# # # The RDA analysis will help us assess which environmental covariates are associated with the communities collected
# # # using both the ABMI and CABIN methods.
# # # Consider transformations to the count based on the type of data.
# # # Look into adding a conditioning call so we can pull out the protocol effect # Feb 24th, 2021 Note Look into this more.
# #
# # rda.results <- rda(X = data.in[, -c(1,2)], Y = site.information[, -c(1:2)], z = site.information$Protocol, scale = TRUE)
# # rda.results <- summary(rda.results)
# #
# # # Visualize
# # species.coord <- data.frame(RDA1 = rda.results$species[, 1],
# #                             RDA2 = rda.results$species[, 2])
# #
# # species.subset.coord <- species.coord[rownames(summary(simper.results)$ABMI_CABIN)[1:5], ]# When adding labels, we are only adding the most import species identified in the SIMPER analysis. Figure gets muddied with too many labeles.
# #
# # site.coord <- data.frame(RDA1 = rda.results$sites[, 1],
# #                          RDA2 = rda.results$sites[, 2],
# #                          Protocol = site.information$Protocol)
# #
# # biplot.coord <- data.frame(RDA1 = rda.results$biplot[, 1],
# #                            RDA2 = rda.results$biplot[, 2])
# #
# # png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/rda-", data.type, "-species-level_", Sys.Date(), ".png"),
# #     width = 2400,
# #     height = 2400,
# #     res = 300)
# #
# # print(ggplot() +
# #           geom_point(data = site.coord, aes(x = RDA1, y = RDA2, color = Protocol), size = 3) +
# #           scale_color_manual(values = abmi_pal("main")(2)) +
# #           geom_segment(data = species.coord, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
# #                        arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),
# #                                      type = "closed"),linetype = 1, size = 0.6, colour = "#E8A396") +
# #           geom_text(data = species.subset.coord, aes(x = RDA1, y = RDA2, label = row.names(species.subset.coord))) +
# #           geom_segment(data = biplot.coord, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
# #                        arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),
# #                                      type = "closed"),linetype = 1, size = 0.6,colour = "#829EBC") +
# #           geom_text(data = biplot.coord, aes(x = RDA1, y = RDA2, label = row.names(biplot.coord))) +
# #           labs(x = paste0("RDA 1 (", format(100 *rda.results$cont[[1]][2,1], digits=4), "%)"),
# #                y = paste0("RDA 2 (", format(100 *rda.results$cont[[1]][2,2], digits=4), "%)")) +
# #           geom_hline(yintercept = 0, linetype = 2,size=  1) +
# #           geom_vline(xintercept = 0,linetype = 2,size = 1) +
# #           guides(shape=guide_legend(title=NULL,color="black"),
# #                  fill=guide_legend(title=NULL))+
# #           theme_bw())
# #
# # dev.off()

#
# NMDS
#

# The metaMDS function takes a community matrix, transforms it based on the defined distance (brays-curtis), and a defined
# number axes to reduce the dimentionality. We are trying to simplify the data, so we don't want too many axes (e.g., greater than 3)
# but we don't want there to be high stress (poor ability to simplify the data (greater than 0.2 is poor). The try options are defining how many
# iterations the function should run for. Higher values take longer to run, but will result in more stable results between runs.

stress.values <- NULL
for (axis.size in 1:10) {

    stress.values <- c(stress.values, metaMDS(comm = data.in[, -c(1:4)], distance = "bray",
                                              k = axis.size, try = 100, trymax = 500)$stress)

}

stress.values <- data.frame(k = 1:10, Stress = stress.values)

# Create stress plot and identify smallest index with value less than 0.2
ggplot(data = stress.values, aes(x = k, y = Stress, color = "#829EBC")) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = 0.2) +
    theme_bw()

optimal.k <- as.numeric(table(stress.values$Stress < 0.2)["FALSE"]) + 1

# Optimal value of K (3)
nmds.results <- metaMDS(comm = data.in[, -c(1:4)], distance = "bray",
                        k = optimal.k, try = 100, trymax = 1000)

# Lets assess the stress of the nmds using a shepard plot
# It looks like there is some scatter between the ordination and dissimilarity distances, but it isn't too bad.
stressplot(nmds.results)

# Lets visualize the first two axes of the nmds

# store the site level information
data.scores <- as.data.frame(scores(nmds.results, display = "sites"))
#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Protocol <- data.in$Protocol #  add the protocol variable

# Name change
data.scores$Protocol <- gsub("ABMI", "ABMI", data.scores$Protocol)
data.scores$Protocol <- gsub("CABIN", "CABIN", data.scores$Protocol)

# Name change
data.in$Protocol <- gsub("ABMI", "ABMI", data.in$Protocol)
data.in$Protocol <- gsub("CABIN", "CABIN", data.in$Protocol)

# Store the ellipse information
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {

    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))

}

plot.new()
ordiplot(nmds.results)
ord <- ordiellipse(nmds.results, as.factor(data.in$Protocol), display = "sites", kind ="sd", conf = 0.9, label = FALSE)
dev.off()

#Generate ellipse points
ellipse.df <- data.frame()
for(g in unique(data.scores$Protocol)){
    if(g!="" && (g %in% names(ord))){

        ellipse.df <- rbind(ellipse.df, cbind(as.data.frame(with(data.scores[data.scores$Protocol==g,],
                                                                 veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                              ,Protocol=g))
    }
}

# Name change
ellipse.df$Protocol <- gsub("ABMI", "ABMI", ellipse.df$Protocol)
ellipse.df$Protocol <- gsub("CABIN", "CABIN", ellipse.df$Protocol)

# Store the hull information
grp.a <- data.scores[data.scores$Protocol == "ABMI", ][chull(data.scores[data.scores$Protocol ==
                                                                             "ABMI", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI
grp.b <- data.scores[data.scores$Protocol == "CABIN", ][chull(data.scores[data.scores$Protocol ==
                                                                              "CABIN", c("NMDS1", "NMDS2")]), ]  # hull values for CABIN
# Name change
grp.a$Protocol <- gsub("ABMI", "ABMI", grp.a$Protocol)
grp.b$Protocol <- gsub("CABIN", "CABIN", grp.b$Protocol)

hull.data <- rbind(grp.a, grp.b)  # combine the hull data

# Store the species level information
species.scores <- as.data.frame(scores(nmds.results, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  # look at the data
species.scores <- species.scores[rownames(summary(simper.results)$ABMI_CABIN)[1:10], ]

png(filename = "figures/coarse group/2018 data/nmds-adjusted-counts-coarse-group-level_2025-10-01.png",
    width = 2400,
    height = 2400,
    res = 300)

# Looks like there is a lot of overlap between the two protocols
print(ggplot() +
          #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Protocol, group = Protocol), alpha = 0.30) + # add the convex hulls
          geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species), alpha = 0.5) +  # add the species labels
          geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, shape = Protocol, colour = Protocol), size = 4) + # add the point markers
          geom_path(data = ellipse.df, aes(x = NMDS1, y = NMDS2, group = Protocol, colour = Protocol)) +
          scale_color_manual(values = abmi_pal("main")(2)) +
          scale_fill_manual(values = abmi_pal("main")(2)) +
          labs(title = "2018 data") +
          coord_equal() +
          theme_bw())

dev.off()

# Permanova
site.information <- site.information %>%
    mutate(Protocol = if_else(str_starts(Site, "CABIN"),"CABIN", "ABMI")) %>%
    filter(SiteYear %in% data.in$SiteYear) %>%
    distinct(SiteYear, .keep_all = TRUE)
str(site.information)
perma.results <- adonis2(data.in[, -c(1:4)] ~ Protocol, data = site.information, method = "bray", permutations = 999)
write.csv(perma.results, file = "output/coarse group/2018 permanova results.csv", row.names = TRUE) # No difference in their dispersion
#rm(perma.results)

# Analysis of dispersion (permadis analysis) Which one is more dispersed
beta.results <- betadisper(d = vegdist(x = data.in[, -c(1:4)], method = "bray"), group = data.in$Protocol)
anova(beta.results)
write.csv(anova(beta.results), file = paste0("tables/coarse group/2018 permadisp-analysis-", data.type, "-coarse-group-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion


#
# Procrustes test
#

# The procrustes test assesses the configuration of two matrices to maximize similarity
# I believe you can use the original matrix, but it is more informative to test on the ordination results (NMDS)
# We can take the ordination scores from the data.scores object created in the previous step

pt.visual <- procrustes(Y = data.scores[data.scores$Protocol == "ABMI", 1:2],
                        X = data.scores[data.scores$Protocol == "CABIN", 1:2], symmetric = TRUE)
pt.test <- protest(Y = data.scores[data.scores$Protocol == "ABMI", 1:2],
                   X = data.scores[data.scores$Protocol == "CABIN", 1:2], symmetric = TRUE)

# # Distance between the two plots
residuals(pt.test)
plot(residuals(pt.test))
# lm(residuals(pt.test) ~ covariates)# FIX Use only the measures that are consistent between protocols # NEED TO UPDATE

write.csv(data.frame(SS = pt.test$ss,
                     Correlation = pt.test$scale,
                     Significance = pt.test$signif), file = paste0("output/coarse group/2018-procrustes-analysis-", data.type, "-coarse-group-level_", Sys.Date(), ".csv"), row.names = FALSE) # No difference in their dispersion

# Change plot

before.data <- data.frame(pt.visual$X)
colnames(before.data) <- colnames(pt.visual$X)
after.data <- data.frame(pt.visual$Yrot)
colnames(after.data) <- colnames(pt.visual$X)

png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/procrustes-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

print(ggplot() +
          geom_point(data = after.data, aes(x = NMDS1, y = NMDS2, color = "ABMI"), size = 3, show.legend = FALSE) + # add the point markers
          geom_segment(aes(x = after.data$NMDS1, y = after.data$NMDS2, xend = before.data$NMDS1, yend = before.data$NMDS2), color = abmi_pal("main")(2)[2], size = 0.5, lty = 1, arrow = arrow(length=unit(0.30,"cm"), type = "closed"), show.legend = FALSE) + # add the point markers
          scale_color_manual(values = abmi_pal("main")(1)) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          xlab("Axis 1") +
          ylab("Axis 2") +
          coord_equal() +
          theme_bw())

dev.off()

#### 2023 data only ####
## filter data sets
adj.count <- cbind(metadata, species_summary)
adj.count <- adj.count %>%
    filter(Year != "2018")

## read in environmental data
site.information <- read_rds("output/cleaned site environmental data.rds")

# Assess degree of missing data
apply(site.information, 2, function(x) table(is.na(x)))
##there are NaNs for one site because they were listed as "DNC" in the original data file
site.information <- site.information %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
## replace missing values with mean values
site.information[is.na(site.information$Human_Footprint), "Human_Footprint"] <- mean(site.information$Human_Footprint, na.rm = TRUE)

site.information$Human_Footprint <- unname(site.information$Human_Footprint) #remove names so vif will run
site.information <- as.data.frame(site.information) %>% #vifstep doesn't like tibbles
    distinct(SiteYear, .keep_all = TRUE)  %>% # drop duplicate observations
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 2)
str(site.information)
site.information <- site.information %>%
    filter(str_ends(SiteYear, "2018"))

###################
# Adjusted Counts #
###################

data.type <- "adjusted-counts"
data.in <- adj.count

# Align the species data with the site data
rownames(data.in) <- data.in$SiteYear
#data.in <- data.in[rownames(site.information), ] # Match site order

#
# Species Accumulation Curves
#

# ABMI
abmi.sac <- specaccum(comm = data.in[data.in$Protocol == "ABMI", -c(1,2,3,4)], method = "random")

# CABIN
cabin.sac <- specaccum(comm = data.in[data.in$Protocol == "CABIN", -c(1,2,3,4)], method = "random")

# Visualize
spp.accum <- data.frame(Protocol = c(rep("ABMI", length(abmi.sac$sites)), rep("CABIN", length(cabin.sac$sites))),
                        Site = c(abmi.sac$sites, cabin.sac$sites),
                        Richness = c(abmi.sac$richness, cabin.sac$richness),
                        SD = c(abmi.sac$sd, cabin.sac$sd))
spp.accum$Upper <- spp.accum$Richness + spp.accum$SD
spp.accum$Lower <- spp.accum$Richness - spp.accum$SD

png(filename = paste0("figures/coarse group/2022 and 2023 data/species-accumulation-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

ggplot(data = spp.accum) +
    geom_pointrange(aes(x = Site, y = Richness, ymin = Lower, ymax = Upper, color = Protocol, fill = Protocol)) +
    scale_color_manual(values = abmi_pal("main")(2)) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "2022 and 2023 data",
         subtitle = paste("Correlation = ", round(cor(spp.accum[spp.accum$Protocol == "ABMI", "Richness"], spp.accum[spp.accum$Protocol == "CABIN", "Richness"]), 3))) +
    #ylim(c(0,300)) +
    theme_bw()

dev.off()

# ##find the number of species observed in ABMI vs. CABIN
# cabin.spp <- data.in %>%
#     filter(Protocol == "CABIN") %>%
#     select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
# print(ncol(cabin.spp)-4) #199 species
#
# abmi.spp <- data.in %>%
#     filter(Protocol == "ABMI") %>%
#     select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
# print(ncol(abmi.spp)-4) #212 species
#
# General diversity characteristics
#

site.richness <- data.frame(Site = data.in$Site,
                            Protocol = data.in$Protocol,
                            Simpson = diversity(x = data.in[, -c(1:4)], index = "simpson", MARGIN = 1),
                            Shannon = diversity(x = data.in[, -c(1:4)], index = "shannon", MARGIN = 1),
                            Richness = specnumber(x = data.in[, -c(1:4)], MARGIN = 1))


richness.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Richness"], y = site.richness[site.richness$Protocol == "CABIN", "Richness"], paired = TRUE)
simpson.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Simpson"], y = site.richness[site.richness$Protocol == "CABIN", "Simpson"], paired = TRUE)
shannon.results <- t.test(x = site.richness[site.richness$Protocol == "ABMI", "Shannon"], y = site.richness[site.richness$Protocol == "CABIN", "Shannon"], paired = TRUE)

# Name change
#site.richness$Protocol <- gsub("ABMI", "ABMI", site.richness$Protocol)
#site.richness$Protocol <- gsub("CABIN", "CABIN", site.richness$Protocol)

png(filename = paste0("figures/coarse group/2022 and 2023 data/simpson-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Simpson, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Simpson Diversity; P-value = ", round(simpson.results$p.value, 4))) +
    labs(title = "2022 and 2023 data",
         subtitle = paste0("Simpson Diversity; P-value = ", round(simpson.results$p.value, 3))) +
    theme_bw()


site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Simpson")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Simpson")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "2022 and 2023 data",
         subtitle = paste0("Simpson Diversity; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()

multiplot(fig.1, fig.2, cols = 2)

dev.off()

png(filename = paste0("figures/coarse group/2022 and 2023 data/richness-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Richness, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Richness; P-value = ", round(richness.results$p.value, 4))) +
    labs(title = "2022 and 2023 data",
         subtitle = paste0("Richness; P-value = ", round(richness.results$p.value, 3))) +
         theme_bw()

site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Richness")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Richness")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "2022 and 2023 data",
         subtitle = paste0("Richness; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()


multiplot(fig.1, fig.2, cols = 2)

dev.off()

png(filename = paste0("figures/coarse group/2022 and 2023 data/shannon-", data.type, "-coarse-group-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Shannon, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    #ggtitle(paste0("Shannon Diversity; P-value = ", round(shannon.results$p.value, 2))) +
    labs(title = "2022 and 2023 data",
         subtitle = paste0("Shannon Diversity; P-value = ", round(shannon.results$p.value, 3))) +
    theme_bw()

site.cor <- data.frame(ABMI = site.richness[site.richness$Protocol == "ABMI", c("Shannon")],
                       CABIN = site.richness[site.richness$Protocol == "CABIN", c("Shannon")])

fig.2 <- ggplot(data = site.cor) +
    geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
    scale_color_manual(values = abmi_pal("main")(1)) +
    #geom_abline(intercept = 0, slope = 1) +
    labs(title = "2022 and 2023 data",
         subtitle = paste0("Shannon Diversity; Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
    theme_bw()

multiplot(fig.1, fig.2, cols = 2)

dev.off()

######################
# Relative Abundance #
######################
#### Multivariate stats ####

data.type <- "relative-abundance"
# data.in <- rel.abund
data.in <- decostand(adj.count[,5:35], method = "total") #convert to relative abundance
data.in <- cbind(adj.count[, 1:4], data.in)
sum(is.na(data.in))
summary(data.in)

# Align the species data with the site data
rownames(data.in) <- paste(data.in$Site, data.in$Protocol, sep = "_")
#data.in <- data.in[rownames(site.information), ] # Match site order

# #
# # SIMPER
# #
#
# # The SIMPER function takes a community matrix (data.in[, -c(1,2)]) and a grouping factor (Protocol)
# # and converts the matrix into a brays-curtis dissimilarity. Columns in the community matrix (e.g., species, familes, etc)
# # are listed in order of highest to lowest contribution.
#
# simper.results <- simper(comm = data.in[, -c(1:4)], group = data.in$Protocol,
#                          permutations = 1000,
#                          trace = TRUE)
#
# # Summary table of results for all columns in the community matrix.
# # ChiroUA, OligoUA, AmphiUA, ChaobUA, GastrUA, and EphemUA are the top group
# summary(simper.results)
# write.csv(summary(simper.results)$ABMI_CABIN,
#           file = paste0("output/simper-analysis-adjusted-counts-coarse-group-level_2025-10-02.csv", row.names = TRUE))
#
#
# # #
# # # RDA analyses
# # #
# #
# # # The RDA analysis will help us assess which environmental covariates are associated with the communities collected
# # # using both the ABMI and CABIN methods.
# # # Consider transformations to the count based on the type of data.
# # # Look into adding a conditioning call so we can pull out the protocol effect # Feb 24th, 2021 Note Look into this more.
# #
# # rda.results <- rda(X = data.in[, -c(1,2)], Y = site.information[, -c(1:2)], z = site.information$Protocol, scale = TRUE)
# # rda.results <- summary(rda.results)
# #
# # # Visualize
# # species.coord <- data.frame(RDA1 = rda.results$species[, 1],
# #                             RDA2 = rda.results$species[, 2])
# #
# # species.subset.coord <- species.coord[rownames(summary(simper.results)$ABMI_CABIN)[1:5], ]# When adding labels, we are only adding the most import species identified in the SIMPER analysis. Figure gets muddied with too many labeles.
# #
# # site.coord <- data.frame(RDA1 = rda.results$sites[, 1],
# #                          RDA2 = rda.results$sites[, 2],
# #                          Protocol = site.information$Protocol)
# #
# # biplot.coord <- data.frame(RDA1 = rda.results$biplot[, 1],
# #                            RDA2 = rda.results$biplot[, 2])
# #
# # png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/rda-", data.type, "-species-level_", Sys.Date(), ".png"),
# #     width = 2400,
# #     height = 2400,
# #     res = 300)
# #
# # print(ggplot() +
# #           geom_point(data = site.coord, aes(x = RDA1, y = RDA2, color = Protocol), size = 3) +
# #           scale_color_manual(values = abmi_pal("main")(2)) +
# #           geom_segment(data = species.coord, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
# #                        arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),
# #                                      type = "closed"),linetype = 1, size = 0.6, colour = "#E8A396") +
# #           geom_text(data = species.subset.coord, aes(x = RDA1, y = RDA2, label = row.names(species.subset.coord))) +
# #           geom_segment(data = biplot.coord, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
# #                        arrow = arrow(angle = 22.5,length = unit(0.35,"cm"),
# #                                      type = "closed"),linetype = 1, size = 0.6,colour = "#829EBC") +
# #           geom_text(data = biplot.coord, aes(x = RDA1, y = RDA2, label = row.names(biplot.coord))) +
# #           labs(x = paste0("RDA 1 (", format(100 *rda.results$cont[[1]][2,1], digits=4), "%)"),
# #                y = paste0("RDA 2 (", format(100 *rda.results$cont[[1]][2,2], digits=4), "%)")) +
# #           geom_hline(yintercept = 0, linetype = 2,size=  1) +
# #           geom_vline(xintercept = 0,linetype = 2,size = 1) +
# #           guides(shape=guide_legend(title=NULL,color="black"),
# #                  fill=guide_legend(title=NULL))+
# #           theme_bw())
# #
# # dev.off()

#
# NMDS
#

# The metaMDS function takes a community matrix, transforms it based on the defined distance (brays-curtis), and a defined
# number axes to reduce the dimentionality. We are trying to simplify the data, so we don't want too many axes (e.g., greater than 3)
# but we don't want there to be high stress (poor ability to simplify the data (greater than 0.2 is poor). The try options are defining how many
# iterations the function should run for. Higher values take longer to run, but will result in more stable results between runs.

stress.values <- NULL
for (axis.size in 1:10) {

    stress.values <- c(stress.values, metaMDS(comm = data.in[, -c(1:4)], distance = "bray",
                                              k = axis.size, try = 100, trymax = 500)$stress)

}

stress.values <- data.frame(k = 1:10, Stress = stress.values)

# Create stress plot and identify smallest index with value less than 0.2
ggplot(data = stress.values, aes(x = k, y = Stress, color = "#829EBC")) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = 0.2) +
    theme_bw()

optimal.k <- as.numeric(table(stress.values$Stress < 0.2)["FALSE"]) + 1

# Optimal value of K (3)
nmds.results <- metaMDS(comm = data.in[, -c(1:4)], distance = "bray",
                        k = optimal.k, try = 100, trymax = 1000)

# Lets assess the stress of the nmds using a shepard plot
# It looks like there is some scatter between the ordination and dissimilarity distances, but it isn't too bad.
stressplot(nmds.results)

# Lets visualize the first two axes of the nmds

# store the site level information
data.scores <- as.data.frame(scores(nmds.results, display = "sites"))
#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Protocol <- data.in$Protocol #  add the protocol variable

# Name change
data.scores$Protocol <- gsub("ABMI", "ABMI", data.scores$Protocol)
data.scores$Protocol <- gsub("CABIN", "CABIN", data.scores$Protocol)

# Name change
data.in$Protocol <- gsub("ABMI", "ABMI", data.in$Protocol)
data.in$Protocol <- gsub("CABIN", "CABIN", data.in$Protocol)

# Store the ellipse information
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {

    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))

}

plot.new()
ordiplot(nmds.results)
ord <- ordiellipse(nmds.results, as.factor(data.in$Protocol), display = "sites", kind ="sd", conf = 0.9, label = FALSE)
dev.off()

#Generate ellipse points
ellipse.df <- data.frame()
for(g in unique(data.scores$Protocol)){
    if(g!="" && (g %in% names(ord))){

        ellipse.df <- rbind(ellipse.df, cbind(as.data.frame(with(data.scores[data.scores$Protocol==g,],
                                                                 veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                              ,Protocol=g))
    }
}

# Name change
ellipse.df$Protocol <- gsub("ABMI", "ABMI", ellipse.df$Protocol)
ellipse.df$Protocol <- gsub("CABIN", "CABIN", ellipse.df$Protocol)

# Store the hull information
grp.a <- data.scores[data.scores$Protocol == "ABMI", ][chull(data.scores[data.scores$Protocol ==
                                                                             "ABMI", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI
grp.b <- data.scores[data.scores$Protocol == "CABIN", ][chull(data.scores[data.scores$Protocol ==
                                                                              "CABIN", c("NMDS1", "NMDS2")]), ]  # hull values for CABIN
# Name change
grp.a$Protocol <- gsub("ABMI", "ABMI", grp.a$Protocol)
grp.b$Protocol <- gsub("CABIN", "CABIN", grp.b$Protocol)

hull.data <- rbind(grp.a, grp.b)  # combine the hull data

# Store the species level information
species.scores <- as.data.frame(scores(nmds.results, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  # look at the data
species.scores <- species.scores[rownames(summary(simper.results)$ABMI_CABIN)[1:10], ]

png(filename = "figures/coarse group/2022 and 2023 data/nmds-adjusted-counts-coarse-group-level_2025-10-01.png",
    width = 2400,
    height = 2400,
    res = 300)

# Looks like there is a lot of overlap between the two protocols
print(ggplot() +
          #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Protocol, group = Protocol), alpha = 0.30) + # add the convex hulls
          geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species), alpha = 0.5) +  # add the species labels
          geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, shape = Protocol, colour = Protocol), size = 4) + # add the point markers
          geom_path(data = ellipse.df, aes(x = NMDS1, y = NMDS2, group = Protocol, colour = Protocol)) +
          scale_color_manual(values = abmi_pal("main")(2)) +
          scale_fill_manual(values = abmi_pal("main")(2)) +
          labs(title = "2022 and 2023 data") +
          coord_equal() +
          theme_bw())

dev.off()

# Permanova
site.information <- site.information %>%
    #mutate(Protocol = if_else(str_starts(Site, "CABIN"),"CABIN", "ABMI")) %>%
    filter(SiteYear %in% data.in$SiteYear) %>%
    distinct(SiteYear, .keep_all = TRUE)
str(site.information)
perma.results <- adonis2(data.in[, -c(1:4)] ~ Protocol, data = site.information, method = "bray", permutations = 999)
write.csv(perma.results, file = "output/coarse group/2022 and 2023 permanova results.csv", row.names = TRUE) # No difference in their dispersion
#rm(perma.results)

# Analysis of dispersion (permadis analysis) Which one is more dispersed
beta.results <- betadisper(d = vegdist(x = data.in[, -c(1:4)], method = "bray"), group = data.in$Protocol)
anova(beta.results)
write.csv(anova(beta.results), file = paste0("tables/coarse group/2022 and 2023 permadisp-analysis-", data.type, "-coarse-group-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion


#
# Procrustes test
#

# The procrustes test assesses the configuration of two matrices to maximize similarity
# I believe you can use the original matrix, but it is more informative to test on the ordination results (NMDS)
# We can take the ordination scores from the data.scores object created in the previous step

pt.visual <- procrustes(Y = data.scores[data.scores$Protocol == "ABMI", 1:2],
                        X = data.scores[data.scores$Protocol == "CABIN", 1:2], symmetric = TRUE)
pt.test <- protest(Y = data.scores[data.scores$Protocol == "ABMI", 1:2],
                   X = data.scores[data.scores$Protocol == "CABIN", 1:2], symmetric = TRUE)

# # Distance between the two plots
residuals(pt.test)
plot(residuals(pt.test))
# lm(residuals(pt.test) ~ covariates)# FIX Use only the measures that are consistent between protocols # NEED TO UPDATE

write.csv(data.frame(SS = pt.test$ss,
                     Correlation = pt.test$scale,
                     Significance = pt.test$signif), file = paste0("output/coarse group/2022 and 2023-procrustes-analysis-", data.type, "-coarse-group-level_", Sys.Date(), ".csv"), row.names = FALSE) # No difference in their dispersion

# Change plot

before.data <- data.frame(pt.visual$X)
colnames(before.data) <- colnames(pt.visual$X)
after.data <- data.frame(pt.visual$Yrot)
colnames(after.data) <- colnames(pt.visual$X)

png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/procrustes-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

print(ggplot() +
          geom_point(data = after.data, aes(x = NMDS1, y = NMDS2, color = "ABMI"), size = 3, show.legend = FALSE) + # add the point markers
          geom_segment(aes(x = after.data$NMDS1, y = after.data$NMDS2, xend = before.data$NMDS1, yend = before.data$NMDS2), color = abmi_pal("main")(2)[2], size = 0.5, lty = 1, arrow = arrow(length=unit(0.30,"cm"), type = "closed"), show.legend = FALSE) + # add the point markers
          scale_color_manual(values = abmi_pal("main")(1)) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          xlab("Axis 1") +
          ylab("Axis 2") +
          coord_equal() +
          theme_bw())

dev.off()
