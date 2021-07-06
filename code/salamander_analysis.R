#Salamander Analysis
rm(list=ls())

#Libraries
#devtools::install_github('NEON-biodiversity/Ostats')
library(tidyverse)
library(dplyr)
library(tidyr)
library(plyr)
library(Ostats)
library(ggplot2)
library(iNEXT)
library(piecewiseSEM)
library(effects)





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
sal_SVL <- sal_SVL%>% 
           mutate(lat_long = paste(Latitude, Longitude, sep = "_"))

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
svl_site_merged <- svl_site%>% 
                   mutate(SITE2 = paste(SITE, Elevation, sep = "_"))%>% 
                   mutate(tax_Site = paste(ID, SITE2, sep = "_"))#unique identifier for species/site combos


####calculating richness as a covariate
#note that this is calculating richness after removing unidentified sp but before removing site/taxa combos with <5 individuals --talk about with group
# generate vectors of abundances by species for each site
salTables <- svl_site_merged  %>% 
             group_by(SITE2) %>% 
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

#name site id "SITE2" to match other data frame
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
          dplyr::count(tax_Site) %>%
          filter(n >4)

#filter(n >4) #will allow us to keep sites where some species have 5 or more individual
# however, we would have to filter out this with only 1 ( filter(count >1))species left, discuss with group

#take svl data for only those sites where all species have >4 individuals
svl_site_input<-svl_site_filt[svl_site_filt$tax_Site %in% hi_abund$tax_Site, ]

#this leaves us with 883 valid sites
length(unique(svl_site_input$SITE2))



####recalculating richness as a covariate after removing species/site combos with less than 5 individuals (i.e., species richness of actual species we have itv for)

salTables2 <- svl_site_input  %>% 
  group_by(SITE2) %>% 
  do(t = table(.$ID))

# Name the list of vectors
mamx <- lapply(salTables2$t, as.numeric)
names(mamx) <- salTables2$SITE2

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

#Estimtes for each site (species rich, shannon, simpson)
richness_estimators$AsyEst


#just grab richness by site (or grab other diversity indicators?)
asymptotic_richness2 <- richness_estimators$AsyEst %>% 
                       filter(Diversity == 'Species richness') 

#check out new distribution
hist(asymptotic_richness2$Observed)

#name site id "SITE2" to match other data frame
asymptotic_richness2$SITE2 <- factor(asymptotic_richness2$Site, levels=asymptotic_richness2$Site[order(asymptotic_richness2$Observed)])#make siteID column to left join

#join species richness to data frame
svl_wRich2<-svl_site_input%>%
            left_join(., asymptotic_richness2, by= "SITE2")%>%
            filter(Observed.y>1)

#prep data for OSTATs----

#  filter for only 3 sites to test then select the relevant columns required by the function site, id, svl.
# Use the mutate function to add a new column named "log_SVL" to log-transform the measurements.

o_data <- svl_wRich2 %>%
         #filter(SITE %in% c('14','83', '3045'))%>% #need to filter for 3 sites if you want the mean plots to work 
         dplyr::select(SITE2, ID, SVL) %>%
         filter(!is.na(SVL)) %>%
         mutate(log_SVL = log10(SVL))



#select the env columns from the matrix to use with ostats output (cold probably get rid omany more in this data set...)
o_env <- svl_site_input %>%
         dplyr::select(-SITE, -USNM ,-ID, -SVL,-tax_Site)

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
                         nperm=1)

#make ostats a data frame

ostats_output<-as.data.frame(overlap_sal)

#give Ostats output a site id (SITE2) column from the current rownames and join to env data

sal_output<-ostats_output%>%
              mutate(SITE2 = row.names(ostats_output))%>%#give Ostats output a site id column from the current rownames
              left_join(.,unique(o_env), by = "SITE2") #join site env data to ostats_output


####bioclim for these sites####

library(raster)
library(sp)

r <- getData("worldclim",var="bio",res=10)

r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")


coords <- sal_output %>%
         dplyr::select(.,Latitude, Longitude)%>%
         mutate_if(is.character, as.numeric)
        

values <- extract(r,coords[,2:1])

clim <- cbind.data.frame(coords,values,sal_output$SITE2)
colnames(clim)[5]<-"SITE2"



#join new climte data with sal_out
sal_results<-clim%>%
            dplyr::select(-Latitude, -Longitude)%>%
            left_join(.,sal_output, by = "SITE2")%>%
            mutate(Temp = Temp / 10)



#need code here to save out OSTATS
#write.csv(sal_output,"outputs/overlap_5_14.csv")



####Analyze ostats output####


svl_overlap<-sal_results #output from above if you don't call it in
#svl_overlap<-read.csv("outputs/ostats_outputv1.csv")#all data with only 1 species sites removed




#rename variables so it is easier...
overlap2 = dplyr::rename(svl_overlap, Overlap = overlaps_norm,
                    Richness1 = Observed.x,Richness2 = Observed.y, Precipitation=Prec, Temperature=Temp)

#view data ordered by any variable
dplyr::select(overlap2,SITE2, Overlap, Richness2)%>%
  arrange(.,Richness2)


#run some exploratory models...
mod1<-lm(Richness1~Overlap, data=overlap2)
mod2<-lm(Richness~Temperature, data=overlap2)
mod3<-lm(Richness~Precipitation, data=overlap2)
mod4<-lm(Richness~Temperature*Precipitation, data=overlap2)
mod5<-lm(Overlap~Temperature, data=overlap2)
mod6<-lm(Overlap~Precipitation, data=overlap2)
mod7<-lm(Overlap~Temperature+Precipitation, data=overlap2)

#look at models
summary(mod1)
plot(mod1)
car::vif(mod6)
cor(mam_output$logweight,mam_output$field_mean_annual_temperature_C)

#Plot univariate relationships
ggplot(overlap2, aes(x=Overlap, y=Richness1)) + 
geom_point()+
geom_smooth(method=lm)+
xlab("Overlap")+
ylab ("Richness")


####Work on plotting
#get inputs for the plot function
#sites2use<-c(unique(o_data$SITE2))
sites2use<-c("1473_3350" ,"6398_3600" ,"3041_3900")
plots <- o_data$SITE2
sp <- o_data$ID
traits <- o_data$log_SVL

Ostats_plot(plots = plots, sp = sp, traits = traits,
            overlap_dat = overlap_sal,
            use_plots = sites2use, means = TRUE)






####piecewise SEM####
#models

#predicting species richness (Observed)
rich<-lm(Richness1~Overlap+Temperature+Precipitation, data=overlap2)
plot(rich)
summary(rich)
car::vif(rich)

#predicting niche overlap
niche_overlap<-lm(Overlap~Temperature+Precipitation, data=overlap2)
plot(niche_overlap)
summary(niche_overlap)
car::vif(niche_overlap)

####Path model using psem
model1<-psem(rich,niche_overlap)
summary(model1, .progressBar = F)
AIC(model1)

#save out coefficients table
mod1_coefs<-coefs(model1)


#write.csv(mod1_coefs, file = "results/mpd_model.csv", quote = FALSE, row.names = F)

plot(model1)

AIC(model1)
plot(model1, node_attrs = list(
  shape = "rectangle", color = "black",
  fillcolor = "orange", x = 3, y=1:12))

