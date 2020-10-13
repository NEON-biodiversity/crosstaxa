#Libraries

library(dplyr)
library(tidyr)
library(plyr)
library(Ostats)
devtools::install_github('NEON-biodiversity/Ostats')



#salamander data----

#individual svl measurements

#salamander data


sal_SVL<-read.csv("./data/MasterSVLData.csv")

head(sal_SVL)


#site data


sal_site<-read.csv("./data/sitedata.csv")

head(sal_site)




#rename latitude and longitude column names in sal_site to match sal_SVL
sal_site<-dplyr::rename(sal_site,
              Latitude = Lat2,
              Longitude = Long2)


#make a new column from a concatenation of lat and long called lat_long for the svl data set
sal_SVL <- sal_SVL%>% mutate(
           lat_long = paste(Latitude, Longitude, sep = "_"))

#make a new column from a concatenation of lat and long called lat_long for the site data set
sal_site <- sal_site%>% mutate(
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


#works to get species counts per site and filter out sites with < 2 species
svl_site_filt<-ddply(svl_site, .(SITE), mutate, count = length(unique(ID)))%>%
  filter(count >1)

head(svl_site_filt)




#prep data for OSTATs----
# Load data from Read et al. (2018) from Figshare web archive
#dat <- read.csv('https://ndownloader.figshare.com/files/9167548')


#  filter for only 3 sites to test then select the relevant columns required by the function site, id, svl.
# Use the mutate function to add a new column named "log_SVL" to log-transform the measurements.

dat <- svl_site_filt %>%
  #filter(SITE %in% c('14','83', '1473'))%>% 
  select(SITE, ID, SVL) %>%
  filter(!is.na(SVL)) %>%
  mutate(log_SVL = log10(SVL))

# Group the data by siteID and taxonID and look at the summary 
dat %>%
  group_by(SITE, ID) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(dat)


#run Ostats function: copied from vingette

Ostats_example <- Ostats(traits = as.matrix(dat[,'log_SVL']),
                         sp = factor(dat$ID),
                         plots = factor(dat$SITE),
                         data_type = "linear")

#make ostats a data frame

ostats_output<-as.data.frame(Ostats_example)

#give Ostats output a site id column from the current rownames

final_output<-ostats_output%>%
              mutate(SITE = as.integer(row.names(ostats_output)))%>%#give Ostats output a site id column from the current rownames
              left_join(.,final_site, by = "SITE") #join site data to ostats_output
  
#missing elevation grrrrrrrrrrrrrrrrrr

mod<-lm(overlaps_norm~elevation, data=final_output)

####Work on plotting

sites2use<-c('14','83', '1473')
Ostats_plot(indiv_dat = dat, plots = dat$SITE, sp = dat$ID, trait = dat$log_SVL, overlap_dat = Ostats_example, sites2use = sites2use, name_x = 'SVL (log-transformed)', means=T)
?Ostats_plot
