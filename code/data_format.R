#Libraries

library(dplyr)
library(tidyr)
library(plyr)

# Load data from Read et al. (2018) from Figshare web archive
dat <- read.csv('https://ndownloader.figshare.com/files/9167548')


# Keep only sites "HARV" and "JORN" using the filter function, and select the relevant columns required  # by the function.
# Use the mutate function to add a new column named "log_weight" to log-transform the measurements.

dat <- dat %>%
  filter(siteID %in% c('HARV','JORN')) %>%
  select(siteID, taxonID, weight) %>%
  filter(!is.na(weight)) %>%
  mutate(log_weight = log10(weight))

# Group the data by siteID and taxonID and look at the summary 
dat %>%
  group_by(siteID, taxonID) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(dat)


<<<<<<< HEAD
#salamander data----

#individual svl measurements
=======
#salamander data

>>>>>>> 02581cbde9e545b42f77b6c8ed5ececc58516443
sal_SVL<-read.csv("./data/MasterSVLData.csv")

head(sal_SVL)

<<<<<<< HEAD
#sire data
=======
>>>>>>> 02581cbde9e545b42f77b6c8ed5ececc58516443
sal_site<-read.csv("./data/sitedata.csv")

head(sal_site)

<<<<<<< HEAD
#rename latitude and longitude column names in sal_site to match sal_SVL
sal_site<-dplyr::rename(sal_site,
              Latitude = Lat2,
              Longitude = Long2)

#make a new column for a concatination of lat and long called lat_long 
sal_site <- sal_site%>% mutate(
            lat_long = paste(Latitude, Longitude, sep = "_"))

#look for duplicates
length(unique(sal_site$lat_long))
sal_site$lat_long[duplicated(sal_site$lat_long)]

#filter out all duplicate rows (i.e., both pairs of the duplicate)
singletons <- names(which(table(sal_site$lat_long) == 1))
final_site<-sal_site[sal_site$lat_long %in% singletons, ]



=======
rename(sal_site$Latitude2 = Latitude), Longitude2 = Longitude)

dplyr::rename(sal_site,
  Latitude = Latitude2, # better naming conventions
  Longitude = Longitude2)
>>>>>>> 02581cbde9e545b42f77b6c8ed5ececc58516443
