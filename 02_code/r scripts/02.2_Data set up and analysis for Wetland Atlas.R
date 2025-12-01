## Data setup and analysis final
## Author: Ermias T. Azeria
## Last edited by Emily: September 18, 2025
#Load raw data
#Look up table
#Site information

#install.packages("openxlsx")
#install.packages("abmi.themes")
#install.packages("tidyverse")
#install.packages("backports")

library(openxlsx)
library(abmi.themes)
library(tidyverse)
library(backports)

setwd("~/GitHub/New aquatic invertebrate analyses")

natcol <- scale_fill_abmi(palette = "main")$palette(6)

### this script requires functions from "Taxonomic ambiguous assignment abmi_clean.R"
source("scripts/Taxonomic ambiguous assignment abmi_clean.R")
##Em's note, there will be an error message:"Error in eval(ei, envir) :
#object 'fix.sites' not found"
## that's ok - the fix.sites function was a fragment in Ermias' script that
#I couldn't figure out where it came from

#First step: USE "Analysis_Grouping" to aggregate "Analysis_Name";
#but keep taxon names with "Distribute to Children" in the Analysis_Grouping as is in raw data
#handle this in LOOK up table, add Taxon_Analysis

#Invert_lookup <- xlsx::read.xlsx("C:/ABMI/Wetland Invertebrate/Taxonomc assignment processing/Data/Look-up Table_Invert Taxonomy_2023-09-28.xlsx", 2)
# ROOTout <- "C:/ABMI/Wetland Invertebrate/Taxonomc assignment processing/Results/"
# ROOTdat <- "C:/ABMI/Wetland Invertebrate/Taxonomc assignment processing/Data/"

Invert_lookup <-  read.csv("Data/Look-up Table_Invert Taxonomy_2025-07-28.csv")
#Add a column with analysis name for ambiguous parents designated as "distribute to children" in analysis grouping
Invert_lookup$Taxon_Analysis <- Invert_lookup$Analysis_Grouping
parent.id <- which(Invert_lookup$Taxon_Analysis=="Distribute.to.Children")
Invert_lookup$Taxon_Analysis[parent.id ] <- Invert_lookup$Analysis_Name[parent.id]

#Get fine species data, and use lookup table to first aggregate them to analysis-grouping
#use first step taxon grouping to aggregate data
SppData.raw <- read.csv("~/GitHub/invertebrate-protocol-analyses_2021/data/base/Aquatic invertebrate wide_all fine id_UMOS excluded_2024-03-07.csv", row.names = 1)
umosData <- read.csv("~/GitHub/invertebrate-protocol-analyses_2021/data/base/Aquatic invertebrate wide_targeted UMOS_2024-03-07.csv", row.names = 1)
#GET MCN INFO from COARSE DATA
SppData.coarse <- read.csv("~/GitHub/invertebrate-protocol-analyses_2021/data/base/Aquatic invertebrate Coarse aggregate by sites_2024-03-07.csv", row.names = 1)
#site strata
WetlandInfo <- get(load("~/GitHub/invertebrate-protocol-analyses_2021/Taxonomic assignment processing/Data/wetland-climate_2023 updated.RData"))

#Omit non-applicable wetlands if any with missing info,
cabin <- which(startsWith(rownames(SppData.raw  ), "CABIN"))
head(SppData.raw [cabin, ])

SppData <- SppData.raw #[-cabin, ]
SppData$Site <-  SppData$Year <- SppData$SiteYear <- NULL
mnames <- match(colnames(SppData), Invert_lookup$Analysis_Name) #check
colnames(SppData)[is.na(mnames)] #check for missing names


#Match and aggregate to taxon group (analysis grouping)
mnames <- match(colnames(SppData), Invert_lookup$Analysis_Name)
taxon_gp <- Invert_lookup$Taxon_Analysis[mnames]
SppData.in  <- as.data.frame(as.matrix(mefa4::groupSums(as.matrix(SppData) ,2, by= as.character(taxon_gp), na.rm=TRUE)))
#site info
misssites <- setdiff(rownames(SppData.in ), rownames(WetlandInfo))
misssites

#Auxiliary data
cabin <- which(startsWith(rownames(umosData ), "CABIN"))
head(umosData [cabin, ])
umosData <- umosData[-cabin, ]
umosData$Site <-  umosData$Year <- umosData$SiteYear <- NULL
mnames <- match(colnames(umosData), Invert_lookup$Analysis_Name)
colnames(umosData)[is.na(mnames)]


mnames <- match(colnames(umosData), Invert_lookup$Analysis_Name)
colnames(umosData)[is.na(mnames)]


taxon_gp <- Invert_lookup$Taxon_Analysis[mnames]
umosData.in  <- as.data.frame(as.matrix(mefa4::groupSums(as.matrix(umosData) ,2, by= as.character(taxon_gp), na.rm=TRUE)))
setdiff(rownames(umosData.in ), rownames(WetlandInfo))
umosData.in$Do.not.analyze<- NULL

#RUN ANALYSIS

Taxon_Lookup <- Invert_lookup

#Match main data and auxiliary (UMOS) data
mmdata <- fn_matchdata (data=SppData.in, aux.data=umosData.in)
SppData.in <- mmdata$data
umosData.in<- mmdata$aux.data



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FIRST ADJUST FOR MARCHANT CELLS (balanced)
#Alternative is to run the analysis first then adjust for proportion of MCN after processing data

SppData.bal <- SppData.in
setdiff(rownames( SppData.bal) , rownames( SppData.coarse ))#wetlands with no coarse data
setdiff(rownames( SppData.coarse ),rownames( SppData.bal)  )#wetlands with no fine data

commsites.b <- intersect(rownames(SppData.bal), rownames(SppData.coarse))
SppData.bal <- SppData.bal[ commsites.b ,  ]
SppData.coarse.bal <- SppData.coarse[ commsites.b ,  ]
SppData.bal  <- (SppData.bal /SppData.coarse.bal$MC_counted) *100 #adjust to 100 cells
SppData.bal  <- round(SppData.bal)

umosData.bal <- umosData.in[ commsites.b ,  ]
NR.site.b <- as.character(WetlandInfo[  commsites.b  , "NR"])
NSR.site.b <- as.character(WetlandInfo[  commsites.b , "NSR"])

#### run algorithm ####
## Em's note: I've commented out the DPAC code because we are going with the
## APTC-f approach

#WITH STRATA
#most frequent/most abundant
#APTC


resAPTC.frq.bal <- aptc_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
                                     strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="frequent", reg.prop=F)
# resAPTC.abu.bal  <- aptc_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
#                                       strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=F)



write.xlsx(resAPTC.frq.bal , rowNames=TRUE, file = "output/MCN adjusted abundance_Ambiguous parent resolved all taxa_APTC.frequent.xlsx")
# write.xlsx(resAPTC.abu.bal, rowNames=TRUE, file = paste0(ROOTout, "MCN adjusted abundance/Ambiguous parent resolved all taxa_APTC.abundant.xlsx"))


# #DPAC
# resDPAC.frq.bal <- dpac_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
#                                      strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="frequent", reg.prop=F)
# resDPAC.abu.bal  <- dpac_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
#                                       strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=F)
#
#
#
# write.xlsx(resDPAC.frq.bal , rowNames=TRUE, file = paste0(ROOTout, "MCN adjusted abundance/Ambiguous parent resolved all taxa_DPAC.frequent.xlsx"))
# write.xlsx(resDPAC.abu.bal, rowNames=TRUE, file = paste0(ROOTout, "MCN adjusted abundance/Ambiguous parent resolved all taxa_DPAC.abundant.xlsx"))



#PROPORTIONAL DATA
#APTC
resAPTC.frqp.bal  <- aptc_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
                                       strata.primary=NSR.site.b,strata.secondary=NR.site.b, method.region="frequent", reg.prop=T)
# resAPTC.abup.bal  <- aptc_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
#                                        strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=T)


write.xlsx(resAPTC.frqp.bal , rowNames=TRUE, file = "output/MCN adjusted abundance_Ambiguous parent resolved all taxa_APTC.Propfrequency.xlsx")
# write.xlsx(resAPTC.abup.bal, rowNames=TRUE, file = paste0(ROOTout, "MCN adjusted abundance/Ambiguous parent resolved all taxa_APTC.Propabundance.xlsx"))



# #DPAC
# resDPAC.frqp.bal  <- dpac_abmi (data=SppData.bal,aux.data=umosData.bal, taxalookup=Taxon_Lookup,
#                                        strata.primary=NSR.site.b,strata.secondary=NR.site.b, method.region="frequent", reg.prop=T)
# resDPAC.abup.bal  <- dpac_abmi (data=SppData.bal, aux.data=umosData.bal,taxalookup=Taxon_Lookup,
#                                        strata.primary=NSR.site.b,strata.secondary=NR.site.b,  method.region="abundant", reg.prop=T)
#
#
# write.xlsx(resDPAC.frqp.bal , rowNames=TRUE, file = paste0(ROOTout, "MCN adjusted abundance/Ambiguous parent resolved all taxa_DPAC.Propfrequency.xlsx"))
# write.xlsx(resDPAC.abup.bal, rowNames=TRUE, file = paste0(ROOTout, "MCN adjusted abundance/Ambiguous parent resolved all taxa_DPAC.Propabundance.xlsx"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Biodiversity metrics
#Most frequent or most abundant in region
APTC.frq_bdiv.adj <- fn.commdiv (resAPTC.frq.bal)
# APTC.abu_bdiv.adj <- fn.commdiv (resAPTC.abu.bal)
# DPAC.frq_bdiv.adj <- fn.commdiv (resDPAC.frq.bal)
# DPAC.abu_bdiv.adj <- fn.commdiv (resDPAC.abu.bal)

#proportional to regional abundance /frequency

APTC.Propfrq_bdiv.adj <- fn.commdiv (resAPTC.frqp.bal)
# APTC.Propabu_bdiv.adj <- fn.commdiv (resAPTC.abup.bal)
# DPAC.Propfrq_bdiv.adj <- fn.commdiv (resDPAC.frqp.bal)
# DPAC.Propabu_bdiv.adj <- fn.commdiv (resDPAC.abup.bal)



#Export diversity data
BDIV.out_MCN.adj <- list(APTC.frq=APTC.frq_bdiv.adj , #APTC.abu=APTC.abu_bdiv.adj, DPAC.frq=DPAC.frq_bdiv.adj , DPAC.abu=DPAC.abu_bdiv.adj  ,
                         APTC.Propfrq=APTC.Propfrq_bdiv.adj) #APTC.Propabu=APTC.Propabu_bdiv.adj) #DPAC.Propfrq=DPAC.Propfrq_bdiv.adj , DPAC.Propabu=DPAC.Propabu_bdiv.adj )
write.xlsx(BDIV.out_MCN.adj , rowNames=TRUE, file = "output/Biodiversity metrics for resolved taxonomy all taxa_Marchant count adjusted.xlsx")

save(BDIV.out_MCN.adj ,  file = "output/Biodiversity metrics for resolved taxonomy all taxa_Marchant count adjusted.RData")

#### make summary table for wetland atlas update ####
wetland.atlas.summary.data <- BDIV.out_MCN.adj$APTC.frq %>%
    dplyr::select(Total.abu, Taxa.rich, ShanonDiv)%>%
    rownames_to_column(var = "SiteYear") %>%
    left_join( WetlandInfo %>% dplyr::select(SiteYear, NR, NSR),
        by = "SiteYear") %>%
    group_by(NR) %>%
    summarise(mean.abundance = mean(Total.abu),
              se.abundance = sd(Total.abu)/sqrt(length(Total.abu)),
              mean.richness = mean(Taxa.rich),
              se.richness = sd(Taxa.rich)/sqrt(length(Taxa.rich)),
              mean.shannon.div = mean(ShanonDiv),
              se.shannon.div = sd(ShanonDiv)/sqrt(length(ShanonDiv)),
              max.shannon.div = max(ShanonDiv),
              min.shannon.div = min(ShanonDiv))
## export
write.csv(wetland.atlas.summary.data,"output/Aquatic invert summary data for wetland atlas figures_09-15-25.csv", row.names = FALSE)
#### PDF - Most frequent or abundant ####
#~~~

for (j in 1:length( BDIV.out_MCN.adj )){
  datIN <- BDIV.out_MCN.adj[[j]]
  dataname <- names( BDIV.out_MCN.adj)[j]
  #  dataname <- deparse(substitute(DPAC.frq_bdiv.adj))
  # dataname <- gsub ("_bdiv", "",dataname, fixed=T )
  ###################


  pdf.out <- paste0(ROOTout,"Summary of wetland invertebrate biodiversity metrics by natural region resolved by ",dataname, "_MCN adjusted.pdf")
  pdf(pdf.out , width=6, height=5)
  for (i in 1:ncol(datIN)){
    varin<- 1
    agg <- aggregate(datIN[ , i] ,list(NR.site.b), mean)
    rownames(agg)<-agg$Group.1
    agg$Group.1 <- NULL
    agg<- t(agg)
    agg <- agg[ , c ("Rocky Mountain","Foothills" , "Boreal" ,"Canadian Shield" , "Parkland"  , "Grassland")]
    op<- par(oma=c(0,0,0,0),mar=c(6,4,2,2))
    x <- barplot(agg,srt=45, col=natcol[1:6] , xaxt="n", ylab="Average", main=colnames(datIN)[i])
    labs <- c ("Rocky Mountain","Foothills" , "Boreal" ,"Canadian Shield" , "Parkland"  , "Grassland")

    text(x, par("usr")[3]-(0.03*max(agg)),   srt = 45, adj = 1, xpd = TRUE,
         labels =  labs, cex = 1, col=natcol , font=2)
    box(bty="l",col="#2d2a2a")
  }
  dev.off()
  par(op)

}








