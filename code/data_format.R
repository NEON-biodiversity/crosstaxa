#Libraries

library(dplyr)
library(tidyr)
library(plyr)




#salamander data----

#individual svl measurements

#salamander data


sal_SVL<-read.csv("./data/MasterSVLData.csv")

head(sal_SVL)


#sire data


sal_site<-read.csv("./data/sitedata.csv")

head(sal_site)


#prep data for OSTATs----
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

#rename latitude and longitude column names in sal_site to match sal_SVL
sal_site<-dplyr::rename(sal_site,
              Latitude = Lat2,
              Longitude = Long2)

#make a new column from a concatination of lat and long called lat_long for the site data set
sal_site <- sal_site%>% mutate(
            lat_long = paste(Latitude, Longitude, sep = "_"))

#make a new column from a concatination of lat and long called lat_long for the svl data set
sal_SVL <- sal_SVL%>% mutate(
  lat_long = paste(Latitude, Longitude, sep = "_"))

#number of lat/longs in SVL data = 3907
length(unique(sal_SVL$lat_long))

#number of unique lat/longs in site data = 4475
length(unique(sal_site$lat_long))

# duplicated lat/longs in site data = 87 (some may be triplicate)
sal_site$lat_long[duplicated(sal_site$lat_long)]

#filter out all duplicate rows (i.e., both pairs of the duplicate) from the site data
singletons <- names(which(table(sal_site$lat_long) == 1))
final_site<-sal_site[sal_site$lat_long %in% singletons, ]

#filter out all duplicate rows (i.e., both pairs of the duplicate) from the svl data
final_svl<-sal_SVL[sal_SVL$lat_long %in% singletons, ]
head(final_site)


#joing site data and svl data
svl_site<-left_join(final_svl, final_site[,-(1:2),], by = "lat_long")  
head(svl_site)
