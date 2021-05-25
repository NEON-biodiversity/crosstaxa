#small mammal analysis


#devtools::install_github('NEON-biodiversity/Ostats')
library(tidyverse)
library(spData)
library(sf)
library(lubridate)
library(plyr)
library(Ostats)
library(ggplot2)

#note that I am setting awd here in the shared files from the G-drive, not using the path in the Rproject
setwd("G:/Shared drives/MacrosystemsBiodiversity/data/organism/L0/neon")

#read in neon mammal data from g-drive
load("DP1.10072.001.Rdata")

#make tibbles
mam_dat<-lapply(neondata, as_tibble)
head(mam_dat)
vars=mam_dat$variables_10072[individualCode]

#this is the tibble with the itv data (similar data sheet to what Q used for the ecography paper)
itv_mammal_data<-mam_dat$mam_pertrapnight
head(itv_mammal_data)

# Get rid of all non-target taxa and combine different columns that were used for individual IDs in different contexts.
mammal_data <- mutate(itv_mammal_data, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  filter(order == 'Rodentia', taxonProtocolCategory == 'target')