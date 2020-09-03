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


#salamander data

sal_SVL<-read.csv("./data/MasterSVLData.csv")

head(sal_SVL)

sal_site<-read.csv("./data/sitedata.csv")

head(sal_site)

rename(sal_site$Latitude2 = Latitude), Longitude2 = Longitude)

dplyr::rename(sal_site,
  Latitude = Latitude2, # better naming conventions
  Longitude = Longitude2)
