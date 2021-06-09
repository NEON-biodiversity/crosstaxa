#Salamander Analysis
rm(list=ls())

#Libraries
#devtools::install_github('NEON-biodiversity/Ostats')
library(dplyr)
library(tidyr)
library(plyr)
library(Ostats)
library(ggplot2)




#salamander data----
#note that I am setting a wd here in the shared files from the G-drive, not using the path in the Rproject
setwd("G:/Shared drives/MacrosystemsBiodiversity/data/organism/L0/non_neon_observations/Adams_salamander")

#individual svl measurements
sal_SVL<-read.csv("MasterSVLData.csv")
head(sal_SVL)


#site data

sal_site<-read.csv("sitedata.csv")
head(sal_site)


#clean data----

#rename latitude and longitude column names in sal_site to match sal_SVL
sal_site<-dplyr::rename(sal_site,
              Latitude = Lat2,
              Longitude = Long2)


#make a new column from a concatenation of lat and long called lat_long as a unique site identifier for the svl data set
sal_SVL <- sal_SVL%>% mutate(
           lat_long = paste(Latitude, Longitude, sep = "_"))

#number of unique lat/longs in SVL data = 3907
length(unique(sal_SVL$lat_long))

#make a new column from a concatenation of lat and long called lat_long for the site data set
sal_site <- sal_site%>% mutate(
            lat_long = paste(Latitude, Longitude, sep = "_"))

#find the number of duplicate lat/longs (i.e. some sites had two plots and they were given the same coarse lat/long)
length(names(which(table(sal_site$lat_long) > 1))) # duplicated lat/longs in site data = 80 


#filter out all duplicate rows (i.e., both pairs of the duplicate) from the site data; Dean Adams said it is fine to do this
singletons <- names(which(table(sal_site$lat_long) == 1))
final_site<-sal_site[sal_site$lat_long %in% singletons, ]

#show that there are no duplicate lat_longs in "final_site"
length(names(which(table(final_site$lat_long) > 1)))


#filter out all duplicate rows (i.e., both pairs of the duplicate) from the svl data
final_svl<-sal_SVL[sal_SVL$lat_long %in% final_site$lat_long, ]
head(final_svl)

length(unique(final_svl$lat_long))


#join site data and svl data, this has all the env data attached to individual measurements
svl_site<-left_join(final_svl, final_site[,-(1:2),], by = "lat_long")  
head(svl_site)

length(unique(svl_site$lat_long))


#merge site and elevation to make a unique identifier for each site because some sites have the same site # but different elevations
svl_site_merged <- svl_site%>% mutate(
  SITE2 = paste(SITE, Elevation, sep = "_"))


####calculating richness as a covariate
#note that this is calculating richness after removing unidentified sp but before removing site/taxa combos with <5 individuals --talk about with group
# generate vectors of abundances by species for each site
salTables <- svl_site_merged  %>% 
  group_by( SITE2) %>% 
  do(t = table(.$ID))

# Name the list of vectors
mamx <- lapply(salTables$t, as.numeric)
names(mamx) <- salTables$SITE2

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

#Estimtes for each site (species rich, shannon, simpson)
richness_estimators$AsyEst


#just grab richness by site (or grab other diversity indicators?)
asymptotic_richness <- richness_estimators$AsyEst %>% 
  filter(Diversity == 'Species richness') 

#name site id "siteID" to match other data frame
asymptotic_richness$SITE2 <- factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])#make siteID column to left join

#join species richness to data frame
svl_wRich<-svl_site_merged%>%
  left_join(., asymptotic_richness, by= "SITE2")



#remove sites with < 2 species
svl_site_filt<-svl_wRich%>%
  filter(Observed >1)

head(svl_site_filt)


#Figure out which species/site combinations  with greater than 5 individuals 
hi_abund<-svl_site_filt %>%
          dplyr::count(SITE2, ID) %>%
          filter(n >4)
#filter(n >4) #will allow us to keep sites where some species have 5 or more individual
# however, we would have to filter out this with only 1 ( filter(count >1))species left, discuss with group

#take svl data for only those sites where all species have >4 individuals
svl_site_input<-svl_site_filt[svl_site_filt$SITE2 %in% hi_abund$SITE2, ]%>%
                filter(Observed >1)

#this leaves us with 323 valid sites
length(unique(svl_site_input$SITE2))
length(unique(hi_abund$SITE2))

####don't think i need this anymore (check)
#try to get an env data set where there is one row per site. this will be used for post-ostats analysis.
#site_vars<-svl_site_filt%>%
          #select (-c(USNM, SVL,ID))%>%#####problem here is that there are other vars that differ e.g.count. do count after? and drop others
          #distinct(.)
  
#length(unique(site_vars$SITE2)) #####problem here is 
#length(site_env$SITE)  

#prep data for OSTATs----

#  filter for only 3 sites to test then select the relevant columns required by the function site, id, svl.
# Use the mutate function to add a new column named "log_SVL" to log-transform the measurements.

o_data <- svl_site_input %>%
  #filter(SITE %in% c('14','83', '3045'))%>% #need to filter for 3 sites if you want the mean plots to work 
  select(SITE2, ID, SVL) %>%
  filter(!is.na(SVL)) %>%
  mutate(log_SVL = log10(SVL))



#select the env columns from the matrix to use with ostats output
o_env <- svl_site_input %>%
  select(-SITE, -USNM ,-ID, -SVL)

# Group the svl data by siteID and taxonID and look at the summary 
o_data %>%
  group_by(SITE2, ID) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(o_data)


####run Ostats function: copied from vignette####

overlap_sal<- Ostats(traits = as.matrix(o_data[,'log_SVL']),
                         sp = factor(o_data$ID),
                         plots = factor(o_data$SITE2),
                         data_type = "linear"
                         nperm=1)

#make ostats a data frame

ostats_output<-as.data.frame(overlap_sal)

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
mod<-lm(overlaps_norm~+Observed, data=svl_overlap2)
summary(mod)
plot(mod)
plot( svl_overlap2$overlaps_norm, svl_overlap2$Observed)

#Plot univariate relationships
ggplot(svl_overlap2, aes(x=BIO1, y=overlaps_norm)) + 
  geom_point()+
  geom_smooth(method=lm)



####Work on plotting
#get inputs for the plot function
sites2use<-c(unique(o_data$SITE2))
#sites2use<-c("4450_2090" ,"4476_2900" ,"4477_2600")
plots <- o_data$SITE2
sp <- o_data$ID
traits <- o_data$log_SVL

Ostats_plot(plots = plots, sp = sp, traits = traits,
            overlap_dat = overlap,
            use_plots = sites2use, means = TRUE)

