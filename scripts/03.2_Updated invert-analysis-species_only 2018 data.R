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

##add base_SiteYear column to identify pairs for wilcox.test
wilcox.data <- data.in %>%
    mutate(base_SiteYear = SiteYear %>% sub("^CABIN-", "", .), .after = 1) %>%
    dplyr::select(!c(Site, Year)) %>%
    column_to_rownames(var = "SiteYear")

site.abundance <- NULL
sign.results <- NULL

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