#small mammal analysis
rm(list=ls())

#packages
#devtools::install_github('NEON-biodiversity/Ostats')
library(tidyverse)
library(spData)
library(sf)
library(lubridate)
library(plyr)
library(Ostats)
library(ggplot2)
library(iNEXT)
library(piecewiseSEM)


#note that I am setting a wd here in the shared files from the G-drive, not using the path in the Rproject
setwd("G:/Shared drives/MacrosystemsBiodiversity/data/organism/L0/neon_downloads")

#read in neon mammal data from g-drive
load("DP1.10072.001.Rdata")

#read in the NEON taxonomy data for all species from g-drive and select relevant columns
tax<-read.csv("../neon_taxa/OS_TAXON_SMALL_MAMMAL-20200129T161511.csv")%>%
     select(taxonID, acceptedTaxonID,vernacularName,taxonProtocolCategory,taxonRank,order,family,subfamily,tribe,genus)

#filter taxa list keeping targeted rodents with species designations (Note that T&E species do not have sp designations)
tax_reduced<-tax%>%
             filter(order == 'Rodentia',taxonProtocolCategory == 'target', taxonRank%in% c('species','subspecies'))#keeping target, rodents with species designation

#read in NEON site environmental data
site_env<-read.csv("C:/Users/bbaiser/Documents/GitHub/ITV/crosstaxa/data/NEON_Field_Site_Metadata_20210226_0_mod.csv")%>%
          dplyr::rename (.,siteID = field_site_id)

  
###WHAT TO DO WITH SUBSPECIES? HERE WE KEEP THEM


#make tibbles
mam_dat<-lapply(neondata, as_tibble)
head(mam_dat)

#variable descriptions from NEON, i.e. metadata
#vars<-mam_dat$variables_10072

#this is the tibble with the itv data 
itv_mammal_data<-mam_dat$mam_pertrapnight
colnames(itv_mammal_data)

####clean data#####

mammal_dataT <- itv_mammal_data%>%
                filter(taxonID%in%tax_reduced$taxonID)%>% #take taxa from neon taxa filtered list
                mutate(year = year(collectDate), logweight=log(weight))%>% #CREATE YEAR COLUMN and log weight column
                filter(lifeStage=="adult")%>% #Keep only adults
                group_by(tagID) %>% #group individuals by tag (i.e., recaptures) 
                filter(collectDate==min(collectDate))%>% #for recaptures, take the earliest record
                ungroup()%>%
                mutate(tax_Site = paste(taxonID, siteID, sep = "_")) #make a unique species-by-site identifier called "tax_Site"

          
           
####calculating richness as a covariate
#note that this is calculating richness after removing unidentified sp but before removing site/taxa combos with <5 individuals --talk about with group
# generate vectors of abundances by species for each site
mammaltables <- mammal_dataT  %>% 
                group_by(siteID) %>% 
                do(t = table(.$taxonID))

# Name the list of vectors
mamx <- lapply(mammaltables$t, as.numeric)
names(mamx) <- mammaltables$siteID

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

#Estimtes for each site (species rich, shannon, simpson)
richness_estimators$AsyEst

# Calculate Chao1 richness estimator and combine all richness estimators by site. (not sure what this does of if it used?)
#chao <- function(x) {
 # xcomm <- table(x$taxonID)
 # S_obs <- length(xcomm)
  #f1 <- sum(xcomm == 1)
  #f2 <- sum(xcomm == 2)
  #return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
#}


#just grab richness by site (or grab other diversity indicators?)
asymptotic_richness <- richness_estimators$AsyEst %>% 
                       filter(Diversity == 'Species richness') 

#name site id "siteID" to match other data frame
asymptotic_richness$siteID <- factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])#make siteID column to left join

#join species richness to data frame
mammal_data<-mammal_dataT%>%
              left_join(., asymptotic_richness, by= "siteID")


####remove species that have fewer than 5 individuals in a given site (we can change this threshold, talk to group)
#also note, we can recalculate site richness after removing these species...

#make list of species/site combos that have <5 individuals  
abund_filt<-mammal_data%>%
            dplyr::count(tax_Site) %>%
            filter(n>4)

#remove those species/site combos and filter sites with <1 species 
high_abun_mam<-mammal_data[mammal_data$tax_Site %in% abund_filt$tax_Site, ]%>%     
               filter(Observed >1)



# 6. Calculate overlap statistics and null effect sizes -------------------



#prep data for OSTATs----

#subset number of stations to run in reasonable time... 
#sub_Site<-c("GRSM", "SCBI", "JORN")
sub_Site<-c("GUAN", "DEJU")#sites with one species to remove

o_stat_mam <- high_abun_mam %>%
              filter(!siteID %in% sub_Site)%>% #remove sites with one species because it messes up plotting
              select(siteID, scientificName, logweight) %>%
              filter(!is.na(logweight)) 

#write.csv(o_stat_mam, "../../L1/cross_taxa_ITV/mammal/o_stat_mam.csv")

#select the env columns from the matrix to use with ostats output (need site level vars here...)
site_rich <- high_abun_mam %>%#take the site richness and site id to join to env below
             select(siteID,Observed)

#join site richness to other env. vars
env<-site_env[site_env$siteID %in% o_stat_mam$siteID, ]%>%
     left_join(.,unique(site_rich), by = "siteID")
     



# Group the mam data by siteID and taxonID and look at the summary 
o_stat_mam  %>%
  group_by(siteID, scientificName) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(o_stat_mam)


####run Ostats function: copied from vignette####

overlap_mam<- Ostats(traits = as.matrix(o_stat_mam[,'logweight']),
                 sp = factor(o_stat_mam$scientificName),
                 plots = factor(o_stat_mam$siteID),
                 #density_args=list(bw=.05),
                 nperm=1)



#write.csv(overlap_mam,"../../L1/cross_taxa_ITV/mammal/overlap_mam.csv")

#make ostats a data frame
ostats_output<-as.data.frame(overlap_mam)
colnames(ostats_output)

#give Ostats output a site id column from the current rownames and join to env data
mam_output<-ostats_output%>%
            mutate(siteID= row.names(ostats_output))%>%#give Ostats output a site id column from the current rownames
            left_join(.,env, by = "siteID")#%>% #join site env data to ostats_output
           # drop_na()




#need code here to save out OSTATS
#write.csv(mam_output,"../../L1/cross_taxa_ITV/mammal/overlap_7_3.csv")


####Analyze ostats output####

#output from above if you don't call it in
#read.csv("outputs/ostats_outputv1.csv")#all data with only 1 species sites removed
select(mam_output,siteID, logweight, Observed)%>%
        arrange(.,logweight)

select(mam_output,siteID, field_mean_annual_precipitation_mm,field_mean_annual_temperature_C,logweight, Observed)%>%
  arrange(.,Observed)


#run some exploratory models...
mod1<-lm(Observed~logweight, data=mam_output)
mod2<-lm(logweight~field_mean_annual_temperature_C, data=mam_output)
mod3<-lm(Observed~field_mean_annual_temperature_C, data=mam_output)
mod4<-lm(logweight~field_mean_annual_precipitation_mm*field_mean_annual_temperature_C, data=mam_output)
mod5<-lm(Observed~field_mean_annual_precipitation_mm, data=mam_output)
mod6<-lm(Observed~field_mean_annual_precipitation_mm*field_mean_annual_temperature_C+logweight, data=mam_output)
mod6<-lm(Observed~field_mean_annual_precipitation_mm+field_mean_annual_temperature_C+logweight, data=mam_output)
summary(mod1)
summary(mod6)
plot(mod6)
car::vif(mod6)
cor(mam_output$logweight,mam_output$field_mean_annual_temperature_C)

#long names for vars in model (cut/paste)
field_mean_canopy_height_m++field_mean_annual_temperature_C+field_mean_annual_precipitation_mm


#Plot univariate relationships
ggplot(mam_output, aes(x=logweight, y=Observed)) + 
  geom_point()+
  geom_smooth(method=lm)+
  #geom_smooth(method= "loess")+
  xlab("Overlap")+
  ylab ("Richness")

str()
plot(Overlap, Richness)

#ostats plots

#inputs for "Ostats_plot" function
sites2use<-c("GRSM", "SCBI", "JORN")#even if you set these sites, if there is an NA from other sites for overlap you get an error
#sites2use<-unique(dat_in$STATION)
plots <- o_stat_mam$siteID
sp <- o_stat_mam$scientificName
traits <- o_stat_mam$logweight

#plot distributions and means
Ostats_plot(plots = plots, sp = sp, traits = traits,
            overlap_dat = overlap_mam,
            use_plots = sites2use, means = TRUE)



####piecewise SEM
#models
#predicting speciers richness (Observed)
rich<-lm(Observed~logweight+field_mean_annual_precipitation_mm*field_mean_annual_temperature_C, data=mam_output)
plot(rich)
summary(rich)
car::vif(rich)
#field_mean_annual_precipitation_mm+field_mean_annual_temperature_C+
niche_overlap<-lm(logweight~field_mean_annual_precipitation_mm*field_mean_annual_temperature_C, data=mam_output)
plot(niche_overlap)
car::vif(rich)

####Path model using psem
model1<-psem(rich,niche_overlap)
summary(model1, .progressBar = F)

#save out coefficents table
mod1_coefs<-coefs(model1)
write.csv(mod1_coefs, file = "results/mpd_model.csv", quote = FALSE, row.names = F)

plot(model1)

AIC(model1)
plot(model1, node_attrs = list(
  shape = "rectangle", color = "black",
  fillcolor = "orange", x = 3, y=1:12))

