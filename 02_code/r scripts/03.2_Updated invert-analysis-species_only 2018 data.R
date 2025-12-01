#
# Title: ABMI Invertebrate Species Groups Analysis for 2018 data
# Created: February 5th, 2021
# Last Updated Emily: November 3, 2025
# Author: Brandon Allen
# Objective: Perform a series of analyses that match and expand on the Hanisch et al (2020) manuscript comparing the ABMI and CABIN protocols
# Keywords: Notes, Initialization, Site information, Species level
#

#########
# Notes #
#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
#install.packages("ggplot2")
#install.packages("usdm")
#install.packages("vegan")

# Load libraries into memory
library(abmi.themes)
library(usdm)
library(vegan)
library(readr)
library(tidyverse)
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
site.information <- read_rds("output/cleaned wetland and climate info.rds")

################
# Missing data #
################

# Assess degree of missing data
apply(site.information, 2, function(x) table(is.na(x)))

# 27 sites are missing lat/long info
# Remaining sites are complete.
## Em's note: I think it's ok if this info isn't there for now.

# Check there are no NA values
table(is.na(site.information))

#################################
# Correlation between variables #
#################################

# We are going to assess correlations using Variance Inflation Factor and pearson correlation coefficients.

# VIF approach (threshold of 5)
vifstep(site.information[, -c(1:7)], th = 5) # Will change the variable set

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
# cor.matrix <- cor(x = site.information[, c("TOC", "Sample_Depth_Mean", "Deep_Samples",
#                                            "Hab_Complex", "Temp_Mean", "DO_Mean",
#                                            "Sal_Mean", "pH_Mean", "DP", "DOC", "Max_Depth",
#                                            "Open_Water", "Human_Footprint")], method = "pearson")
#
# cor.matrix

# Subset the information to the final data set
# site.information <- site.information[, c("Site", "Protocol", "TOC", "Sample_Depth_Mean", "Deep_Samples",
#                                          "Hab_Complex", "Temp_Mean", "DO_Mean",
#                                          "Sal_Mean", "pH_Mean", "DP", "DOC", "Max_Depth",
#                                          "Open_Water", "Human_Footprint")]
# rownames(site.information) <- paste(site.information$Site, site.information$Protocol, sep = "_")
#
# rm(cor.matrix)

# We are going to leave the variables in their current state initially
# Check for normality of the remaining variables

# for (x in 3:14) {
#
#             hist(site.information[, x], main = colnames(site.information)[x])
#             hist(log(site.information[, x]), main = paste0("Log ", colnames(site.information)[x]))
#
# }

#################
# Species level #
#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load species data
adj.count <- read_excel("output/CABIN Ambiguous parent resolved all taxa_APTC.abundant.xlsx") %>%
    rename(SiteYear = '...1') %>%
    mutate(Site = sub("_.*", "", SiteYear), .before = 1) %>%
    mutate(Year = sub(".*_", "", SiteYear), .before = 2) %>%
    mutate(Protocol = if_else(startsWith(Site, "CABIN"), "CABIN", "ABMI"), .after = 3) %>%
    filter(SiteYear != "W478-2_2018") %>% ## remove sites with no observations
    filter(SiteYear != "W638-2_2018") %>%
    filter(SiteYear != "CABIN-W478-2_2018") %>% ## remove sites corresponding cabin rows
    filter(SiteYear != "CABIN-W638-2_2018") %>%
    filter(Year == "2018") #filter to only include 2018 data

taxa <- read_csv("data/Look-up Table_Invert Taxonomy_2025-07-28.csv") %>%
    dplyr::select(!Row)

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
max(abmi.sac$richness)

# CABIN
cabin.sac <- specaccum(comm = data.in[data.in$Protocol == "CABIN", -c(1,2,3,4)], method = "random")
max(cabin.sac$richness)

# Visualize
spp.accum <- data.frame(Protocol = c(rep("ABMI", length(abmi.sac$sites)), rep("CABIN", length(cabin.sac$sites))),
                        Site = c(abmi.sac$sites, cabin.sac$sites),
                        Richness = c(abmi.sac$richness, cabin.sac$richness),
                        SD = c(abmi.sac$sd, cabin.sac$sd))
spp.accum$Upper <- spp.accum$Richness + spp.accum$SD
spp.accum$Lower <- spp.accum$Richness - spp.accum$SD

png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/species-accumulation-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

ggplot(data = spp.accum) +
    geom_pointrange(aes(x = Site, y = Richness, ymin = Lower, ymax = Upper, color = Protocol, fill = Protocol)) +
    scale_color_manual(values = abmi_pal("main")(2)) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "2018 data",
        subtitle = paste("Correlation = ", round(cor(spp.accum[spp.accum$Protocol == "ABMI", "Richness"], spp.accum[spp.accum$Protocol == "CABIN", "Richness"]), 3))) +
    ylim(c(0,300)) +
    theme_bw()

dev.off()

##find the number of species observed in ABMI vs. CABIN
cabin.spp <- data.in %>%
    filter(Protocol == "CABIN") %>%
    select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
str(cabin.spp) #60-4=56 spp

abmi.spp <- data.in %>%
    filter(Protocol == "ABMI") %>%
    select(1:4, where(~ !(is.numeric(.) && sum(., na.rm = TRUE) == 0))) #drop columns where no individuals are observed.
str(abmi.spp) #156-4 = 152 species
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

png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/simpson-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Simpson, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "2018 data",
         #subtitle = paste0("Simpson Diversity; P-value = ", round(simpson.results$p.value, 2))
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

png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/richness-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Richness, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "2018 data",
         #subtitle = paste0("Richness; P-value = ", round(richness.results$p.value, 2))
         subtitle = "Richness; P-value < 0.001") +
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

png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/shannon-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 1200,
    res = 300)

fig.1 <- ggplot(data = site.richness) +
    geom_boxplot(mapping = aes(x = Protocol, y = Shannon, fill = Protocol), show.legend = FALSE) +
    scale_fill_manual(values = abmi_pal("main")(2)) +
    labs(title = "2018 data",
         #subtitle = paste0("Shannon Diversity; P-value = ", round(shannon.results$p.value, 2))
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

# Abundance by taxa
#### Em's note: this is not working. Need to troubleshoot more. ####

##add base_SiteYear column to identify pairs for wilcox.test
wilcox.data <- data.in %>%
    mutate(base_SiteYear = SiteYear %>% sub("^CABIN-", "", .), .after = 1) %>%
    dplyr::select(!c(Site, Year)) %>%
    column_to_rownames(var = "SiteYear")

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
    long.code <- taxa[taxa$Analysis_Name == spp.id, "long_code"]
    temp.abundance$Species <- gsub(spp.id, long.code, temp.abundance$Species)

    # Name change
    temp.abundance$Protocol <- gsub("ABMI", "ABMI", temp.abundance$Protocol)
    temp.abundance$Protocol <- gsub("CABIN", "CABIN", temp.abundance$Protocol)

    png(filename = paste0("figures/invertebrate-protocol-analyses/species-abundance/", spp.id, "-", data.type, "-species-level_", Sys.Date(), ".png"),
        width = 1200,
        height = 1200,
        res = 300)

    print(ggplot(data = temp.abundance) +
              geom_boxplot(mapping = aes(x = Species, y = Abundance, fill = Protocol), show.legend = TRUE) +
              scale_fill_manual(values = abmi_pal("main")(2)) +
              ggtitle(paste0("Wilcoxon = ", round(sign.results[sign.results$Species == spp.id, "Significance"], 3))) +
              theme_bw())

    dev.off()

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
#### this is also not working yet ####

for(spp.id in unique(site.abundance$Species)) {

    # Create dataset
    site.cor <- data.frame(ABMI = site.abundance[site.abundance$Protocol == "ABMI" & site.abundance$Species == spp.id, "Abundance"],
                           CABIN = site.abundance[site.abundance$Protocol == "CABIN" & site.abundance$Species == spp.id, "Abundance"])

    abmi.test <- as.numeric(table(site.cor$ABMI > 0)["TRUE"]) > 10
    if(is.na(abmi.test)) {abmi.test <- FALSE}

    cabin.test <- as.numeric(table(site.cor$CABIN > 0)["TRUE"]) > 10
    if(is.na(cabin.test)) {cabin.test <- FALSE}

    if(abmi.test == TRUE & cabin.test == TRUE) {

        # Abundance Correlation

        fig.1 <- ggplot(data = site.cor) +
            geom_point(mapping = aes(x = ABMI, y = CABIN, color = "ABMI"), show.legend = FALSE) +
            scale_color_manual(values = abmi_pal("main")(1)) +
            geom_abline(intercept = 0, slope = 1) +
            ggtitle(paste0(spp.id, " Correlation = ", round(cor(site.cor$ABMI, site.cor$CABIN), 2))) +
            theme_bw()


        # Standardize environment variables

        site.information.scaled <- site.information

        for (x in c(5:15)) {

            site.information.scaled[, x] <- as.numeric(scale(site.information.scaled[, x]))

        }

        abmi.data <- cbind.data.frame(site.information.scaled, data.frame(Spp = data.in[, c(spp.id)]))
        cabin.data <- abmi.data[abmi.data$Protocol == "CABIN", ]
        abmi.data <- abmi.data[abmi.data$Protocol == "ABMI", ]

        # Modeling

        spp.model <- list()

        abmi.data$Spp <- as.numeric(scale(log(abmi.data$Spp + 0.0001)))

        spp.model[["ABMI"]] <- glm(Spp ~ Deep_Samples + Hab_Complex +
                                       Temp_Mean + DO_Mean + Sal_Mean + pH_Mean + DP + DOC +
                                       Max_Depth + Open_Water +
                                       Human_Footprint, family = "gaussian", data = abmi.data)

        cabin.data$Spp <- as.numeric(scale(log(cabin.data$Spp + 0.0001)))

        spp.model[["CABIN"]] <- glm(Spp ~ Deep_Samples + Hab_Complex +
                                        Temp_Mean + DO_Mean + Sal_Mean + pH_Mean + DP + DOC +
                                        Max_Depth + Open_Water +
                                        Human_Footprint, family = "gaussian", data = cabin.data)

        # Store Coefficients

        model.name <- names(spp.model)
        coef.store <- NULL

        for (x in 1:length(model.name)) {

            temp.results <- summary(spp.model[[x]])$coefficients
            temp.results <- data.frame(Protocol = rep(model.name[x], nrow(temp.results)),
                                       Variable = rownames(temp.results),
                                       Estimate = temp.results[, 1],
                                       StdErr = temp.results[, 2],
                                       Pvalue = temp.results[, 4])

            coef.store <- rbind(coef.store, temp.results)

            rm(temp.results)

        }

        coef.store <- coef.store[, c(1, 2, 3, 4, 5)]

        # Coefficient visualization

        coef.store <- coef.store[!(coef.store$Variable == "(Intercept)"), ]

        # If StdErr is 10x greater than max/min remove
        coef.store[max(abs(coef.store$Estimate))*10 < abs(coef.store$StdErr), c("Estimate", "StdErr")] <- NA

        # Rename variables for legibility
        coef.store$Variable <- as.character(coef.store$Variable)

        # Sort factor order
        coef.store["Coefficient"] <- factor(coef.store$Variable, levels = c("Deep_Samples", "Hab_Complex", "Temp_Mean", "DO_Mean", "Sal_Mean", "pH_Mean",
                                                                            "DP", "DOC", "Max_Depth", "Open_Water", "Human_Footprint"))

        coef.store["LowInner"] <- coef.store$Estimate - coef.store$StdErr
        coef.store["HighInner"] <- coef.store$Estimate + coef.store$StdErr


        fig.2 <- ggplot(coef.store, aes(colour = Protocol, shape = Protocol)) +
            ggtitle(paste0("ABMI vs CABIN: ", spp.id)) +
            xlab("Variables") +
            ylab("Standardized Estimate") +
            geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
            geom_linerange(aes(x = Coefficient, ymin = LowInner,  ymax = HighInner),lwd = 1, position = position_dodge(width = 1/2)) +
            scale_color_manual(values = abmi_pal("main")(2)) +
            geom_pointrange(aes(x = Coefficient, y = Estimate, ymin = LowInner,
                                ymax = HighInner),
                            lwd = 1/2, position = position_dodge(width = 1/2)) +
            coord_flip() +
            theme_bw() +
            theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
                  axis.text = element_text(size = rel(1), colour = "grey30"))

        png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/models/", spp.id, "-model-", data.type, "-species-level_", Sys.Date(), ".png"),
            height = 2400,
            width = 4800,
            res = 300)

        multiplot(fig.1, fig.2, cols = 2)

        dev.off()

    }

}

#
# SIMPER
#

# The SIMPER function takes a community matrix (data.in[, -c(1,2)]) and a grouping factor (Protocol)
# and converts the matrix into a brays-curtis dissimilarity. Columns in the community matrix (e.g., species, familes, etc)
# are listed in order of highest to lowest contribution.

simper.results <- simper(comm = data.in[, -c(1,2,3,4)], group = data.in$Protocol,
                         permutations = 1000,
                         trace = TRUE)

# Summary table of results for all columns in the community matrix.
# ChiroUA, OligoUA, AmphiUA, ChaobUA, GastrUA, and EphemUA are the top group
write.csv(summary(simper.results)$ABMI_CABIN,
          file = paste0("output/simper-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = TRUE)

#
# RDA analyses
#

# The RDA analysis will help us assess which environmental covariates are associated with the communities collected
# using both the ABMI and CABIN methods.
# Consider transformations to the count based on the type of data.
# Look into adding a conditioning call so we can pull out the protocol effect # Feb 24th, 2021 Note Look into this more.
# RDA does not appear to add more than the regression models.

# # ABMI RDA
# abmi.data.in <- data.in[data.in$Protocol == "ABMI", ]
# abmi.site.information <- site.information[rownames(abmi.data.in), ]
#
# rda.results <- rda(X = abmi.data.in[, -c(1,2,3,4)], Y = abmi.site.information[, -c(1:3)], scale = TRUE)
# rda.results <- summary(rda.results)
#
# format(100 *rda.results$cont[[1]][2,1], digits=4)
# format(100 *rda.results$cont[[1]][2,2], digits=4)
#
# # CABIN RDA
# cabin.data.in <- data.in[data.in$Protocol == "CABIN", ]
# cabin.site.information <- site.information[rownames(cabin.data.in), ]
#
# rda.results <- rda(X = cabin.data.in[, -c(1,2,3,4)], Y = cabin.site.information[, -c(1:3)], scale = TRUE)
# rda.results <- summary(rda.results)
#
# format(100 *rda.results$cont[[1]][2,1], digits=4)
# format(100 *rda.results$cont[[1]][2,2], digits=4)
#
# # Original
# rda.results <- rda(X = data.in[, -c(1,2,3,4)], Y = site.information[, -c(1:3)], z = site.information$Protocol, scale = TRUE)
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
# png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/rda-", data.type, "-species-level_", Sys.Date(), ".png"),
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

# Based on this grouping, only ChaobUA is associated with HF. The remaining species are all associated with native conditions.

#
# NMDS
#

# The metaMDS function takes a community matrix, transforms it based on the defined distance (brays-curtis), and a defined
# number axes to reduce the dimentionality. We are trying to simplify the data, so we don't want too many axes (e.g., greater than 3)
# but we don't want there to be high stress (poor ability to simplify the data (greater than 0.2 is poor). The try options are defining how many
# iterations the function should run for. Higher values take longer to run, but will result in more stable results between runs.

stress.values <- NULL
for (axis.size in 1:10) {

    stress.values <- c(stress.values, metaMDS(comm = data.in[, -c(1,2,3,4)], distance = "bray",
                                              k = axis.size, try = 100, trymax = 500)$stress)

}

stress.values <- data.frame(k = 1:10, Stress = stress.values)

optimal.k <- as.numeric(table(stress.values$Stress < 0.2)["FALSE"]) + 1

# Create stress plot and identify smallest index with value less than 0.2
png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/nmds-k-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

print(ggplot(data = stress.values, aes(x = k, y = Stress, color = "#829EBC")) +
          geom_point(show.legend = FALSE) +
          geom_hline(yintercept = 0.2) +
          ggtitle(paste0("Stress = ", round(stress.values$Stress[optimal.k], 3))) +
          theme_bw())

dev.off()

# Optimal value of K (3)
nmds.results <- metaMDS(comm = data.in[, -c(1,2,3,4)], distance = "bray",
                        k = optimal.k, try = 100, trymax = 1000)

# Lets assess the stress of the nmds using a shepard plot
# It looks like there is some scatter between the ordination and dissimilarity distances, but it isn't too bad.
png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/stressplot-best-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

stressplot(nmds.results)

dev.off()

# Lets visualize the first two axes of the nmds

# store the site level information
data.scores <- as.data.frame(scores(nmds.results, display = "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.in)  # create a column of site names, from the rownames of data.scores
data.scores$Protocol <- data.in$Protocol #  add the protocol variable

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

# Store the hull information
grp.a <- data.scores[data.scores$Protocol == "ABMI", ][chull(data.scores[data.scores$Protocol ==
                                                                             "ABMI", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI
grp.b <- data.scores[data.scores$Protocol == "CABIN", ][chull(data.scores[data.scores$Protocol ==
                                                                              "CABIN", c("NMDS1", "NMDS2")]), ]  # hull values for CABIN

hull.data <- rbind(grp.a, grp.b)  # combine the hull data

# Store the species level information
species.scores <- as.data.frame(scores(nmds.results, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  # look at the data
species.scores <- species.scores[rownames(summary(simper.results)$ABMI_CABIN)[1:10], ]

png(filename = paste0("figures/invertebrate-protocol-analyses/2018 species/nmds-", data.type, "-2018-species-level_", Sys.Date(), ".png"),
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
perma.results <- adonis2(data.in[, -c(1,2,3,4)] ~ data.in$Protocol, data = site.information, method = "bray", permutations = 999)
write.csv(perma.results, file = "output/2018 data only/permanova results.csv", row.names = TRUE) # No difference in their dispersion
rm(perma.results)

# Analysis of dispersion (permadis analysis) Which one is more dispersed
beta.results <- betadisper(d = vegdist(x = data.in[, -c(1,2,3,4)], method = "bray"), group = data.in$Protocol)
anova(beta.results)
write.csv(anova(beta.results), file = "output/2018 data only/permadisp results.csv", row.names = TRUE) # No difference in their dispersion

# But we can test for the similarity between them using a procruste test

#
# Procrustes test
#

# The procrustes test assesses the configuration of two matrices to maximize similarity
# I believe you can use the original matrix, but it is more informative to test on the ordination results (NMDS)
# We can take the ordination scores from the data.scores object created in the previous step
##remove extra sites from CABIN df
# data.scores <- data.scores %>%
#     filter(site != "CABIN-W478-2_2018") %>%
#     filter(site != "CABIN-W638-2_2018")

pt.visual <- procrustes(Y = data.scores[data.scores$Protocol == "ABMI", 1:3],
                        X = data.scores[data.scores$Protocol == "CABIN", 1:3], symmetric = TRUE)
pt.test <- protest(Y = data.scores[data.scores$Protocol == "ABMI", 1:3],
                   X = data.scores[data.scores$Protocol == "CABIN", 1:3], symmetric = TRUE)

# # Distance between the two plots
residuals(pt.test)
plot(residuals(pt.test))
# lm(residuals(pt.test) ~ covariates)# FIX Use only the measures that are consistent between protocols # NEED TO UPDATE

write.csv(data.frame(SS = pt.test$ss,
                     Correlation = pt.test$scale,
                     Significance = pt.test$signif), file = paste0("tables/procrustes-analysis-", data.type, "-2018-species-level_", Sys.Date(), ".csv"), row.names = FALSE) # No difference in their dispersion


# Change plot

before.data <- data.frame(pt.visual$X)
colnames(before.data) <- colnames(pt.visual$X)
after.data <- data.frame(pt.visual$Yrot)
colnames(after.data) <- colnames(pt.visual$X)

png(filename = paste0("figures/ordination/procrustes-", data.type, "-species-level_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

print(ggplot() +
          geom_point(data = after.data, aes(x = NMDS1, y = NMDS2, color = "ABMI"), size = 3, show.legend = FALSE) + # add the point markers
          geom_segment(aes(x = after.data$NMDS1, y = after.data$NMDS2, xend = before.data$NMDS1, yend = before.data$NMDS2), color = abmi_pal("main")(2)[2], linewidth = 0.5, lty = 1, arrow = arrow(length=unit(0.30,"cm"), type = "closed"), show.legend = FALSE) + # add the point markers
          scale_color_manual(values = abmi_pal("main")(1)) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          xlab("Axis 1") +
          ylab("Axis 2") +
          coord_equal() +
          theme_bw())

dev.off()

#
#### Shallow vs. Deep - can't complete right now ####
#

# # store the site level information
# data.scores <- as.data.frame(scores(nmds.results))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
# data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
# data.scores$Group <- paste(data.in$Protocol, ifelse(site.information$Deep_Samples < 4, "Shallow",
#                                                      ifelse(site.information$Deep_Samples < 6, "Medium", "Deep")), sep = "_") #  add the protocol variable
#
# # Store the ellipse information
# veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {
#
#     theta <- (0:npoints) * 2 * pi/npoints
#     Circle <- cbind(cos(theta), sin(theta))
#     t(center + scale * t(Circle %*% chol(cov)))
#
# }
#
#
# plot.new()
# ordiplot(nmds.results)
# ord <- ordiellipse(nmds.results, as.factor(data.scores$Group), display = "sites", kind ="sd", conf = 0.9, label = FALSE)
# dev.off()
#
# # Remove CABIN_Deep
# data.scores <- data.scores[data.scores$Group != "CABIN_Deep", ]
#
# #Generate ellipse points
# ellipse.df <- data.frame()
# for(g in unique(data.scores$Group)){
#     if(g!="" && (g %in% names(ord))){
#
#         ellipse.df <- rbind(ellipse.df, cbind(as.data.frame(with(data.scores[data.scores$Group==g,],
#                                                                  veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
#                                               ,Group=g))
#     }
# }
#
#
# # Store the hull information
# grp.a <- data.scores[data.scores$Group == "ABMI_Deep", ][chull(data.scores[data.scores$Group ==
#                                                                              "ABMI_Deep", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI_Deep
# grp.b <- data.scores[data.scores$Group == "ABMI_Medium", ][chull(data.scores[data.scores$Group ==
#                                                                               "ABMI_Medium", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI_Medium
# grp.c <- data.scores[data.scores$Group == "ABMI_Shallow", ][chull(data.scores[data.scores$Group ==
#                                                                            "ABMI_Shallow", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI_Shallow
# grp.d <- data.scores[data.scores$Group == "CABIN_Shallow", ][chull(data.scores[data.scores$Group ==
#                                                                            "CABIN_Shallow", c("NMDS1", "NMDS2")]), ]  # hull values for CABIN_Shallow
#
# hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  # combine the hull data
#
# png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/nmds-", data.type, "-species-level-shallow-vs-deep_", Sys.Date(), ".png"),
#     width = 2400,
#     height = 2400,
#     res = 300)
#
# # Looks like there is a lot of overlap between the two protocols
# print(ggplot() +
#           #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Group, group = Group), alpha = 0.30) + # add the convex hulls
#           geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, shape = Group, colour = Group), size = 4) +
#           geom_path(data = ellipse.df, aes(x = NMDS1, y = NMDS2, group = Group, colour = Group)) + # add the point markers# add the point markers
#           scale_color_manual(values = abmi_pal("main")(4)) +
#           scale_fill_manual(values = abmi_pal("main")(4)) +
#           coord_equal() +
#           theme_bw())
#
# dev.off()

# PERMANOVA
# Organize data
perm.site <- site.information
perm.site$Group <- paste(perm.site$Protocol, ifelse(perm.site$Deep_Samples < 4, "Shallow",
                                                    ifelse(perm.site$Deep_Samples < 6, "Medium", "Deep")), sep = "_")
# Remove CABIN_Deep
perm.site <- perm.site[perm.site$Group != "CABIN_Deep", ]
perm.species <- data.in[rownames(perm.site), ]

perma.results <- adonis2(perm.species[, -c(1,2,3,4)] ~ Group, data = perm.site, method = "bray", permutations = 999)
write.csv(perma.results, file = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/tables/species/permanova-depth-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion
rm(perma.results)

# #
# # Group bootstrap
# #
#
# # Optimal value of K (3)
# nmds.results <- metaMDS(comm = data.in[, -c(1,2,3)], distance = "bray",
#                         k = optimal.k, try = 100, trymax = 1000)
#
# nmds.boot <- NULL
#
# abmi.boot <- data.in[data.in$Protocol == "ABMI", ]
# cabin.boot <- data.in[data.in$Protocol == "CABIN", ]
#
# for (boot.iter in 1:100) {
#
#     boot.comm <- rbind(abmi.boot[sample(1:nrow(abmi.boot), nrow(abmi.boot), TRUE), ],
#                        cabin.boot[sample(1:nrow(cabin.boot), nrow(cabin.boot), TRUE), ])
#     nmds.results <- metaMDS(comm = boot.comm[, -c(1,2,3)], distance = "bray",
#                             k = optimal.k, try = 100, trymax = 1000)
#     nmds.results <- as.data.frame(scores(nmds.results))
#
#
# }

######################
#### Relative Abundance #####
######################

#data.type <- "relative-abundance"
# data.in <- rel.abund
data.in <- decostand(adj.count[,5:262], method = "total") #convert to relative abundance
data.in <- cbind(adj.count[, 1:4], data.in)
sum(is.na(data.in))
summary(data.in)
str(data.in)

# Align the species data with the site data
rownames(data.in) <- paste(data.in$Site, data.in$Protocol, sep = "_")
#data.in <- data.in[rownames(site.information), ] # Match site order
data.in[, 5:ncol(data.in)][is.na(data.in[, 5:ncol(data.in)])] <- 0 #gets rid of NAs
summary(data.in[, -c(1:4)])

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
          file = "output/2018 data only/species-simper-analysis.csv", row.names = TRUE)


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

stress.values <- numeric(0)  # ensures it's a clean numeric vector
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

png(filename = "figures/invertebrate-protocol-analyses/2018 species/nmds-adjusted-counts-2018-species-level_2025-09-19.png",
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
## fix up site.information so permanova will run
site.information <- site.information %>%
    mutate(Protocol = if_else(str_starts(Site, "CABIN"),"CABIN", "ABMI")) %>%
    filter(SiteYear %in% data.in$SiteYear) %>%
    distinct(SiteYear, .keep_all = TRUE)
str(site.information)
perma.results <- adonis2(data.in[, -c(1:4)] ~ Protocol, data = site.information, method = "bray", permutations = 999)
write.csv(perma.results, file = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/tables/species/permanova-protocol-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion
rm(perma.results)

# Analysis of dispersion (permadis analysis) Which one is more dispersed
beta.results <- betadisper(d = vegdist(x = data.in[, -c(1:4)], method = "bray"), group = data.in$Protocol)
anova(beta.results)
write.csv(anova(beta.results), file = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/tables/species/permadisp-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion


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
                     Significance = pt.test$signif), file = paste0("output/2018 data only/procrustes-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = FALSE) # No difference in their dispersion

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

#
# Shallow vs. Deep
#

# store the site level information
data.scores <- as.data.frame(scores(nmds.results))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Group <- paste(data.in$Protocol, ifelse(site.information$Deep_Samples < 4, "Shallow",
                                                    ifelse(site.information$Deep_Samples < 6, "Medium", "Deep")), sep = "_") #  add the protocol variable

# Store the ellipse information
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {

    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))

}

plot.new()
ordiplot(nmds.results)
ord <- ordiellipse(nmds.results, as.factor(data.scores$Group), display = "sites", kind ="sd", conf = 0.9, label = FALSE)
dev.off()

# Remove CABIN_Deep
data.scores <- data.scores[data.scores$Group != "TSA_Deep", ]

#Generate ellipse points
ellipse.df <- data.frame()
for(g in unique(data.scores$Group)){
    if(g!="" && (g %in% names(ord))){

        ellipse.df <- rbind(ellipse.df, cbind(as.data.frame(with(data.scores[data.scores$Group==g,],
                                                                 veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                              ,Group=g))
    }
}


# Store the hull information
grp.a <- data.scores[data.scores$Group == "CTA_Deep", ][chull(data.scores[data.scores$Group ==
                                                                              "CTA_Deep", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI_Deep
grp.b <- data.scores[data.scores$Group == "CTA_Medium", ][chull(data.scores[data.scores$Group ==
                                                                                "CTA_Medium", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI_Medium
grp.c <- data.scores[data.scores$Group == "CTA_Shallow", ][chull(data.scores[data.scores$Group ==
                                                                                 "CTA_Shallow", c("NMDS1", "NMDS2")]), ]  # hull values for ABMI_Shallow
grp.d <- data.scores[data.scores$Group == "TSA_Shallow", ][chull(data.scores[data.scores$Group ==
                                                                                 "TSA_Shallow", c("NMDS1", "NMDS2")]), ]  # hull values for CABIN_Shallow

hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  # combine the hull data

png(filename = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/figures/species/ordination/nmds-", data.type, "-species-level-shallow-vs-deep_", Sys.Date(), ".png"),
    width = 2400,
    height = 2400,
    res = 300)

# Looks like there is a lot of overlap between the two protocols
print(ggplot() +
          #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Group, group = Group), alpha = 0.30) + # add the convex hulls
          geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, shape = Group, colour = Group), size = 4) + # add the point markers
          geom_path(data = ellipse.df, aes(x = NMDS1, y = NMDS2, group = Group, colour = Group)) + # add the point markers# add the point markers
          scale_color_manual(values = abmi_pal("main")(4)) +
          scale_fill_manual(values = abmi_pal("main")(4)) +
          coord_equal() +
          theme_bw())

dev.off()

# PERMANOVA
# Organize data
perm.site <- site.information
perm.site$Group <- paste(perm.site$Protocol, ifelse(perm.site$Deep_Samples < 4, "Shallow",
                                                    ifelse(perm.site$Deep_Samples < 6, "Medium", "Deep")), sep = "_")
# Remove CABIN_Deep
perm.site <- perm.site[perm.site$Group != "CABIN_Deep", ]
perm.species <- data.in[rownames(perm.site), ]

perma.results <- adonis2(perm.species[, -c(1,2)] ~ Group, data = perm.site, method = "bray", permutations = 999)
write.csv(perma.results, file = paste0("D:/ABMI-covid-19/general-requests/RobH/invertebrate-protocol-analyses_2021/tables/species/permanova-depth-analysis-", data.type, "-species-level_", Sys.Date(), ".csv"), row.names = TRUE) # No difference in their dispersion
rm(perma.results)
