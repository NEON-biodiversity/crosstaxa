#Libraries
#devtools::install_github('NEON-biodiversity/Ostats')
library(dplyr)
library(tidyr)
library(plyr)
library(Ostats)
library(ggplot2)




#salamander data----

#individual svl measurements

sal_SVL<-read.csv("./data/MasterSVLData.csv")

head(sal_SVL)


#site data

sal_site<-read.csv("./data/sitedata.csv")

head(sal_site)


#clean data----

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

#find duplicate lat/longs (i.e. sites had two plots and they were given the same coarse lat/long)
length(names(which(table(sal_site$lat_long) > 1))) # duplicated lat/longs in site data = 80 


#filter out all duplicate rows (i.e., both pairs of the duplicate) from the site data; Dean said it is fine to do this
singletons <- names(which(table(sal_site$lat_long) == 1))
final_site<-sal_site[sal_site$lat_long %in% singletons, ]

#show that there are no duplicate lat_longs in "final_site"
length(names(which(table(final_site$lat_long) > 1)))


#filter out all duplicate rows (i.e., both pairs of the duplicate) from the svl data
final_svl<-sal_SVL[sal_SVL$lat_long %in% final_site$lat_long, ]
head(final_svl)

length(unique(final_svl$lat_long))


#join site data and svl data, this has all the env data attached to single measurements
svl_site<-left_join(final_svl, final_site[,-(1:2),], by = "lat_long")  
head(svl_site)

length(unique(svl_site$lat_long))


#merge site and elevation to make a unique identifier for each site because some sites have the same site # but different elevations
svl_site_merged <- svl_site%>% mutate(
  SITE2 = paste(SITE, Elevation, sep = "_"))



#get species counts per site and filter out sites with < 2 species
svl_site_filt<-ddply(svl_site_merged, .(SITE2), mutate, count = length(unique(ID)))%>%
  filter(count >1)

head(svl_site_filt)


#Calculate species counts for each site and take species with 5 or more individuals

hi_abund<-svl_site_filt %>%
  dplyr::count(SITE2, ID) %>%
  filter(n >4)

#take svl data for only those sites where all species have >4 individuals
svl_site_input<-svl_site_filt[svl_site_filt$SITE2 %in% hi_abund$SITE2, ]
#filter(count >1) #in case we want to keep sites where we remove species, we can filter count again here


####dont think i need this anymore (check)
#try to get an env data set where there is one row per site. this will be used for post-ostats analysis.
#site_vars<-svl_site_filt%>%
          #select (-c(USNM, SVL,ID))%>%#####problem here is that there are other vars that differ e.g.count. do count after? and drop others
          #distinct(.)
  
#length(unique(site_vars$SITE2)) #####problem here is 
#length(site_env$SITE)  

#prep data for OSTATs----

#  filter for only 3 sites to test then select the relevant columns required by the function site, id, svl.
# Use the mutate function to add a new column named "log_SVL" to log-transform the measurements.

o_data <- svl_site_filt %>%
  #filter(SITE %in% c('14','83', '1473'))%>% 
  select(SITE2, ID, SVL) %>%
  filter(!is.na(SVL)) %>%
  mutate(log_SVL = log10(SVL))


#select the env colums from the matrix to use with ostats output
o_env <- svl_site_filt %>%
  select(-SITE, -USNM ,-ID, -SVL)

# Group the svl data by siteID and taxonID and look at the summary 
o_data %>%
  group_by(SITE2, ID) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(o_data)


####run Ostats function: copied from vignette####

overlap<- Ostats(traits = as.matrix(o_data[,'log_SVL']),
                         sp = factor(o_data$ID),
                         plots = factor(o_data$SITE2),
                         data_type = "linear")

#make ostats a data frame

ostats_output<-as.data.frame(overlap)

#give Ostats output a site id (SITE2) column from the current rownames and join to env data

final_output<-ostats_output%>%
              mutate(SITE2 = row.names(ostats_output))%>%#give Ostats output a site id column from the current rownames
              left_join(.,unique(o_env), by = "SITE2") #join site env data to ostats_output

#need code here to save out OSTATS
  
####Analyze ostats output####

svl_overlap<-final_output
#read.csv("outputs/ostats_outputv1.csv")#all data with only 1 species sites removed
svl_overlap2<-na.omit(svl_overlap)#remove rows with NA

mod<-lm(overlaps_norm~as.numeric(Latitude)+count+ Elevation+BIO1, data=svl_overlap2)
summary(mod)
plot(mod)
plot(svl_overlap2$count, svl_overlap2$overlaps_norm)


ggplot(svl_overlap2, aes(x=BIO1, y=overlaps_norm)) + 
  geom_point()+
  geom_smooth(method=lm)



####Work on plotting

sites2use<-c('14','83', '1473')
Ostats_plot(indiv_dat = dat, plots = dat$SITE, sp = dat$ID, trait = dat$log_SVL, overlap_dat = Ostats_example, sites2use = sites2use, name_x = 'SVL (log-transformed)', means=T)
?Ostats_plot
