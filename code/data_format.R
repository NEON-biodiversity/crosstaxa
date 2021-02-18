#Libraries

library(dplyr)
library(tidyr)
library(plyr)
library(Ostats)
library(ggplot2)
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

#number of lat/longs in SVL data = 3907
length(unique(sal_SVL$lat_long))

#make a new column from a concatenation of lat and long called lat_long for the site data set
sal_site <- sal_site%>% mutate(
            lat_long = paste(Latitude, Longitude, sep = "_"))



#number of unique lat/longs in site data = 4475 - number of toal lat/longs (4562) = 87
length(sal_site$lat_long)-length(unique(sal_site$lat_long))






# duplicated lat/longs in site data = 80 
length(names(which(table(sal_site$lat_long) > 1))) 



#filter out all duplicate rows (i.e., both pairs of the duplicate) from the site data
singletons <- names(which(table(sal_site$lat_long) == 1))
final_site<-sal_site[sal_site$lat_long %in% singletons, ]

length(names(which(table(final_site$lat_long) > 1)))#show that there are no duplicate lat_longs in "final_site"


#filter out all duplicate rows (i.e., both pairs of the duplicate) from the svl data
final_svl<-sal_SVL[sal_SVL$lat_long %in% final_site$lat_long, ]
head(final_svl)

length(unique(final_svl$lat_long))




#join site data and svl data, this has all the env data attached to single measurements
svl_site<-left_join(final_svl, final_site[,-(1:2),], by = "lat_long")  
head(svl_site)

length(unique(svl_site$lat_long))


#merge site and elevation to make a unique identifer for each site because some sites have the same site # but different elevations
svl_site_merged <- svl_site%>% mutate(
  SITE2 = paste(SITE, Elevation, sep = "_"))

length(unique(svl_site_merged$SITE2))




#works to get species counts per site and filter out sites with < 2 species
svl_site_filt<-ddply(svl_site_merged, .(SITE2), mutate, count = length(unique(ID)))%>%
  filter(count >1)

head(svl_site_filt)

length(duplicated(svl_site_filt$SITE))

#try to get an env data set where there is one row per site. this will be used for post-ostats analysis.
site_vars<-svl_site_filt%>%
          select (-c(USNM, SVL,ID))%>%#####problem here is that there are other vars that diffeer e.g.count. do count after? and drop others
          distinct(.)
  
length(unique(site_vars$SITE2)) #####problem here is 
length(site_env$SITE)  

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


####run Ostats function: copied from vignette####

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
 #need code here to save out OSTATS
  
####Analyze ostats output####

overlap<-read.csv("outputs/ostats_outputv1.csv")
overlap2<-na.omit(overlap)#remove rows with NA

mod<-lm(overlaps_norm~BIO1+Elevation+Latitude, data=overlap2)
summary(mod)
plot(mod)
plot(overlap2$Richness, overlap2$overlaps_norm)


ggplot(overlap2, aes(x=BIO1, y=overlaps_norm)) + 
  geom_point()+
  geom_smooth(method=lm)



####Work on plotting

sites2use<-c('14','83', '1473')
Ostats_plot(indiv_dat = dat, plots = dat$SITE, sp = dat$ID, trait = dat$log_SVL, overlap_dat = Ostats_example, sites2use = sites2use, name_x = 'SVL (log-transformed)', means=T)
?Ostats_plot
