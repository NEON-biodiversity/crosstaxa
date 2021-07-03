#clean MAPS BIRD DATA
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

#loading R data file for maps data

#note that I am setting a wd here in the shared files from the G-drive, not using the path in the Rproject
setwd("G:/Shared drives/MacrosystemsBiodiversity/data/organism/L0/non_neon_observations/MAPS_birds")

load("MAPSexport.Rdata")

#make a tibble
dat<-lapply(maps.data, as_tibble)

rm(maps.data)#why remove maps data because its o big and will slow things down? Ask Daijiang

head(dat)

#this is the individual level data
dat$band#this is the individual level data
dat$stations#this is the station data



#############################################################################
#data formatting for OSTATS 
levels(dat$band$AGE)

#works to get species counts per site (combining across years for now) and filter out sites with < 2 species
#bird_site_filt<-ddply(dat$band, .(STATION), mutate,  count = length(unique(SPEC)))%>%#get species counts per station (over time)
                      #mutate(Spec_Stat = paste(SPEC, STATION, sep = "_"))%>% #Spec_Stat is the station/species combo 
                      #filter(count >1)

#head(bird_site_filt)


bird_site_filt<-dat$band%>%
                filter(AGE %in% c(1,5,6,7,8 ))%>%#TAKE ONLY ADULT BIRDS (i.e., age = 1,5,6,7,8);seems like thats all of them
                mutate(year=year(DATE))%>%#make a year column
                mutate(Spec_Stat = paste(SPEC, STATION, sep = "_"))#Spec_Stat is the station/species combo 
              

#most species rich site
BSOL<-dat$band%>%
      filter(STATION== "BSOL")
BSOL%>%
  count("SPEC")



###calculating richness as a covariate
#note that this is calculating richness before removing site/taxa combos with <5 individuals --talk about with group
# generate vectors of abundances by species for each site
birdtables <- bird_site_filt  %>% 
  group_by(STATION) %>% 
  do(t = table(.$SPEC))

# Name the list of vectors
mamx <- lapply(birdtables$t, as.numeric)
names(mamx) <- birdtables$STATION

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500))

#Estimtes for each site (species rich, shannon, simpson)
richness_estimators$AsyEst


#just grab richness by site (or grab other diverstiy indicators?)
asymptotic_richness <- richness_estimators$AsyEst %>% 
                       filter(Diversity == 'Species richness') 

#name site id "STATION" to match other data frame
asymptotic_richness$STATION <- factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])#make siteID column to left join

#join species richness to data frame
bird_dat<-bird_site_filt%>%
          left_join(., asymptotic_richness, by= "STATION")


####remove species that have fewer than 5 individuals in a given site (we can change this threshold, talk to group)
#also note, we can recalculate site richness after removing these species...

#make list of species/site combos that have <5 individuals  
abund_filt<-bird_dat%>%
            dplyr::count(Spec_Stat) %>%
            filter(n>4)

#remove those species/site combos and filter sites with <1 species 
high_abun_birds<-bird_dat[bird_dat$Spec_Stat %in% abund_filt$Spec_Stat, ]%>%     
                 filter(Observed >1)

write.csv(high_abun_birds,"ben_birds_forQ.csv")

####calculate species richness again with species/site combos with less that 5 removed####
# generate vectors of abundances by species for each site
birdtables2 <- high_abun_birds %>% 
               group_by(STATION) %>% 
               do(t = table(.$SPEC))

# Name the list of vectors
mamx <- lapply(birdtables2$t, as.numeric)
names(mamx) <- birdtables2$STATION

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

#Estimtes for each site (species rich, shannon, simpson)
richness_estimators$AsyEst


#just grab richness by site (or grab other diversity indicators?)
asymptotic_richness2 <- richness_estimators$AsyEst %>% 
  filter(Diversity == 'Species richness') 

#name site id "STATION" to match other data frame
asymptotic_richness2$STATION <- factor(asymptotic_richness2$Site, levels=asymptotic_richness2$Site[order(asymptotic_richness2$Observed)])#make siteID column to left join

#join species richness to data frame
#bird_dat<-high_abun_birds%>% this needs anew name for later that will go through ostats...
  left_join(., asymptotic_richness2, by= "STATION")

#prep data for OSTATs----

#subset number of stations to run in reasonable time... 
sub_station<-c("VINS","PATT","FTGI")

dat_in <- high_abun_birds %>%
          filter(STATION %in% sub_station)%>% 
          select(STATION, SPEC, WEIGHT) %>%
          filter(!is.na(WEIGHT)) %>%
          mutate(log_WEIGHT = log10(WEIGHT))



# Group the data by Station and Species and look at the summary (why? ask quentin or Daijiang)
dat_in %>%
  group_by(STATION, SPEC) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(dat_in)


####run Ostats function: copied from vignette####

Ostats_example3 <- Ostats(traits = as.matrix(dat_in[,'log_WEIGHT']),
                         sp = factor(dat_in$SPEC),
                         plots = factor(dat_in$STATION),
                         output = "mean",
                         nperm = 1)
?Ostats
Ostats_example3


head(Ostats_example3$overlaps_norm)

#make ostats a data frame
ostats_output<-as.data.frame(Ostats_example3)
colnames(ostats_output)

#give Ostats output a site id column from the current rownames and join to env data
bird_output<-ostats_output%>%
            mutate(STATION= row.names(ostats_output))%>%#give Ostats output a site id column from the current rownames
            left_join(., asymptotic_richness, by= "STATION")%>% #join site env data to ostats_output
            drop_na(log_WEIGHT)

#need code here to save out OSTATS
#write.csv(bird_output,"bird_overlap_6_15.csv")


####Analyze ostats output####

#output from above if you don't call it in
bird_output<-read.csv("bird_overlap_6_15.csv", row=1)#all data with only 1 species sites removed

#join species richness 2 to data frame
test<-bird_output%>%
      left_join(., asymptotic_richness2, by= "STATION")

range<-select(test,STATION, log_WEIGHT, Observed.y)%>%
       arrange(.,Observed.y)

#run some models...
mod<-lm(Observed.y~log_WEIGHT, data=test)
summary(mod)
plot(mod)
car::vif(mod)

#long names for vars in model (cut/paste)
field_mean_canopy_height_m++field_mean_annual_temperature_C+field_mean_annual_precipitation_mm

#Plot univariate relationships
ggplot(test, aes(x=log_WEIGHT, y=log(Observed.x)) )+ 
  geom_point()+
  geom_smooth(method=lm)+
  #geom_smooth(method= "loess")+
  xlab("Overlap")+
  ylab ("Richness")











#ostats plots for three sites/stations

#inputs for "Ostats_plot" function

sites2use<-c("VINS", "PATT", "FTGI") #pick specific sites or
#sites2use<-unique(dat_in$STATION)# if you filter above then use this
plots <- dat_in$STATION
sp <- dat_in$SPEC
traits <- dat_in$log_WEIGHT

#plot distributions and means
Ostats_plot(plots = plots, sp = sp, traits = traits,
            overlap_dat = Ostats_example3,
            use_plots = sites2use, means = TRUE)



#sites of interest

#most species rich site
BSOL<-dat$band%>%
  filter(STATION== "VINS")
BSOL%>%
  count("SPEC")


#############################################################################3
#alternative way to plot with inputs included
Ostats_plot(plots = dat_in$STATION, sp = dat_in$SPEC, traits = dat_in$log_WEIGHT,
            overlap_dat = Ostats_example,
            use_plots = sites2use, means = TRUE)


#make ostats a data frame for further analysis
ostats_bird_output<-as.data.frame(Ostats_example)

#make a data frame of site richness
site_richness<-dat_in %>% 
  distinct(STATION, count)

#give Ostats output a site id column from the current rownames

#****need to add in other env. vars with code below

final_output<-ostats_bird_output%>%
  mutate(STATION= row.names(ostats_output))%>%#give Ostats output a site id column from the current rownames
  left_join(.,site_richness, by = "STATION") #join site data to ostats_output

#rename(final_output, SITE = STATION)

mod<-lm(overlaps_norm~count, data=final_output)
summary(mod)

plot(final_output$count,final_output$overlaps_norm)
#need code here to save out OSTATS






#ostats plots

#inputs for "Ostats_plot" function
#sites2use<-c("0004", "0005", "0006")
sites2use<-unique(dat_in$STATION)
plots <- dat_in$STATION
sp <- dat_in$SPEC
traits <- dat_in$log_WEIGHT

#plot distributions and means
Ostats_plot(plots = plots, sp = sp, traits = traits,
            overlap_dat = Ostats_example,
            use_plots = sites2use, means = TRUE)

#alternative way to plot with inputs included
Ostats_plot(plots = dat_in$STATION, sp = dat_in$SPEC, traits = dat_in$log_WEIGHT,
            overlap_dat = Ostats_example,
            use_plots = sites2use, means = TRUE)

