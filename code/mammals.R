#small mammal analysis
rm(list=ls())

#devtools::install_github('NEON-biodiversity/Ostats')
library(tidyverse)
library(spData)
library(sf)
library(lubridate)
library(plyr)
library(Ostats)
library(ggplot2)
library(iNEXT)

#note that I am setting awd here in the shared files from the G-drive, not using the path in the Rproject
setwd("G:/Shared drives/MacrosystemsBiodiversity/data/organism/L0/neon_downloads")

#read in neon mammal data from g-drive
load("DP1.10072.001.Rdata")

#make tibbles
mam_dat<-lapply(neondata, as_tibble)
head(mam_dat)

#variable descriptions
#vars<-mam_dat$variables_10072

#this is the tibble with the itv data (similar data sheet to what Q used for the ecography paper)
itv_mammal_data<-mam_dat$mam_pertrapnight
head(itv_mammal_data)

# Get rid of all non-target taxa and combine different columns that were used for individual IDs in different contexts.

#this one does not remove juveniles doubles etc
#mammal_data <- mutate(itv_mammal_data, 
                     # individual =  as.character(tagID), na.rm=TRUE,
                     # year = year(collectDate), logweight=log(weight)) 


mammal_data2 <- itv_mammal_data%>%
                mutate(year = year(collectDate), logweight=log(weight))%>% #make a year column and a log weight column
                drop_na(tagID, scientificName)%>% # get rid of empty trap (i.e. NA tags) and species designation na
                filter(lifeStage=="adult")%>% #take only adults
                filter( !grepl('sp\\.', scientificName))%>%#remove identification with sp. (e.g.,Peromyscus sp.)
                group_by(tagID) %>% #group individuals by tag (i.e., recaptures) 
                filter(collectDate==min(collectDate))%>% #for recaptures, take the earliest record
                ungroup()%>%
                mutate(tax_Site = paste(taxonID, siteID, sep = "_"))#make taxa by site designation



####calculating richness as a covaraite
#note that this is calculating richness after removing unidentified sp but before removing site/taxa combos with <5 inividuals--talk about with group
# generate vectors of abundances by species for each site
mammaltables <- mammal_data2  %>% 
  group_by(siteID) %>% 
  do(t = table(.$taxonID))

# Name the list of vectors
mamx <- lapply(mammaltables$t, as.numeric)
names(mamx) <- mammaltables$siteID

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

#Estimtes for each site )species rich, shannon, simpson)
richness_estimators$AsyEst

# Calculate Chao1 richness estimator and combine all richness estimators by site. (not sure what this does of if it used?)
#chao <- function(x) {
 # xcomm <- table(x$taxonID)
 # S_obs <- length(xcomm)
  #f1 <- sum(xcomm == 1)
  #f2 <- sum(xcomm == 2)
  #return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
#}


#just grab richness by site
asymptotic_richness <- richness_estimators$AsyEst %>% 
                       filter(Diversity == 'Species richness') 

#name site id "siteID" to match other data frame
asymptotic_richness$siteID <- factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])#make siteID column to left join

#join species richness to data frame
mammal_data3<-mammal_data2%>%
      left_join(., asymptotic_richness, by= "siteID")

####remove species that have fewer than 5 individuals in a given site (we can change this threshold, talk to group)
#also note, we can recalculate site richness after removing these species...
#make list of species/site combos that have <5 individuals  
abund_filt<-mammal_data3%>%
            dplyr::count(tax_Site) %>%
            filter(n>4)

#remove those species/site combos and filter sites with <1 species 
high_abun_mam<-mammal_data3[mammal_data3$tax_Site %in% abund_filt$tax_Site, ]%>%     
      filter(Observed >1)







# 6. Calculate overlap statistics and null effect sizes -------------------



#prep data for OSTATs----

#  filter for only 3 sites to test then select the relevant columns required by the function site, id, svl.
# Use the mutate function to add a new column named "log_SVL" to log-transform the measurements.

o_stat_mam <- high_abun_mam %>%
  #need to filter for 3 sites if you want the mean plots to work 
  select(siteID, scientificName, logweight) %>%
  filter(!is.na(logweight)) 



#select the env columns from the matrix to use with ostats output
mam_env <- high_abun_mam %>%
           select(siteID, elevation ,nlcdClass,decimalLatitude,decimalLongitude )

# Group the svl data by siteID and taxonID and look at the summary 
o_stat_mam  %>%
  group_by(siteID, scientificName) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(o_stat_mam)


####run Ostats function: copied from vignette####

overlap_mam<- Ostats(traits = as.matrix(o_stat_mam[,'logweight']),
                 sp = factor(o_stat_mam$scientificName),
                 plots = factor(o_stat_mam$siteID),
                 data_type = "linear",
                 nperm=1)

head(overlap_mam)

#make ostats a data frame

ostats_output<-as.data.frame(overlap_mam)

#give Ostats output a site id (SITE2) column from the current rownames and join to env data

sal_output<-ostats_output%>%
  mutate(SITE2 = row.names(ostats_output))%>%#give Ostats output a site id column from the current rownames
  left_join(.,unique(o_env), by = "SITE2") #join site env data to ostats_output

#need code here to save out OSTATS
#write.csv(sal_output,"outputs/overlap_5_14.csv")



####Analyze ostats output####

svl_overlap<-sal_output #output from above if you don't call it in
#read.csv("outputs/ostats_outputv1.csv")#all data with only 1 species sites removed

svl_overlap2<-na.omit(svl_overlap)#remove rows with NA


#run some models...
mod<-lm(overlaps_norm~as.numeric(Latitude)+count+ Elevation+BIO1, data=svl_overlap2)
summary(mod)
plot(mod)
plot(svl_overlap2$count, svl_overlap2$overlaps_norm)

#Plot univariate relationships
ggplot(svl_overlap2, aes(x=BIO1, y=overlaps_norm)) + 
  geom_point()+
  geom_smooth(method=lm)




####calculating richness####

# Get rid of poorly identified individuals (those marked sp.) because they should not count for our sampling . . . they are spurious singletons.
# Then generate vectors of abundances by species for each site
#note that this is calculated on al
mammaltables <- high_abun_mam %>% 
  group_by(siteID) %>% 
  do(t = table(.$taxonID))


# Name the list of vectors
mamx <- lapply(mammaltables$t, as.numeric)
names(mamx) <- mammaltables$siteID

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

#Estimtes for each site )species rich, shannon, simpson)
richness_estimators$AsyEst

# Calculate Chao1 richness estimator and combine all richness estimators by site. (not sure what this does of if it used?)
chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}


#just grab richness by site
asymptotic_richness <- richness_estimators$AsyEst %>% 
  filter(Diversity == 'Species richness') 
asymptotic_richness$siteID <- factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])#make siteID column to left join



