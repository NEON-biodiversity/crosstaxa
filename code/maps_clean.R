#clean MAPS BIRD DATA

library(tidyverse)
library(spData)
library(sf)
library(lubridate)


#loading R data file for maps data
load("./data/MAPSexport.Rdata")

#make a tibble
dat = lapply(maps.data, as_tibble)

rm(maps.data)#why remove maps data because its o big and will slow things down?

head(dat)
dat$band

n_distinct(dat$band$SPEC)# calculate the species richness

#data year range (1989-2017)
range(dat$band$DATE)

#data weight range (4.1g-353.4g) and histogram
range(dat$band$WEIGHT)
hist(dat$band$WEIGHT)

#rename year columns and make yr2 (year number in data set) and log mass columns
dat$band = mutate(dat$band, mass_g_log10 = log10(WEIGHT),
                  yr = year(DATE),
                  yr2 = yr - min(yr))

#histogram of log weight
hist(dat$band$mass_g_log10)

#I am guessing this is individual plots for each species? look into it

# pdf("mass.pdf", width = 8, height = 5, onefile = T)
# for(i in unique(dat$band$SPEC)){
#   cat(i, "\t")
#   d = filter(dat$band, SPEC == i)
#   if(nrow(d) > 1){
#     m = lm(WEIGHT ~ yr, d)
#     m2 = coef(summary(m))
#     plot(x = d$yr, y = d$WEIGHT, col = d$SEX,
#          main = paste0(i, ": slope = ", round(m2[2, "Estimate"], 5), ", p = ", round(m2[2, "Pr(>|t|)"], 4)))
#     abline(m)
#   } else {
#     plot(x = d$yr, y = d$WEIGHT, main = i)
#   }
# }
# dev.off()


#change column names to lowercase
names(dat$band) = tolower(names(dat$band))#change to lowercase

# Identify species with < 10 observations (ten individual observations across the entire data set) needs to be 10 per site (or per year)
sp_rm = group_by(dat$band, spec) %>% 
  tally() %>% 
  filter(n < 10) %>% 
  pull(spec) %>% as.character()

dat$band = filter(dat$band, !spec %in% sp_rm)

n_distinct(dat$band$spec) # 267 species




# remove outliers function (why trim=.1, why .25 and 2 *mean?)
rm_outlier = function(d, trim_pct = 0.1){
  t_mean = mean(d$weight, na.rm = TRUE, trim = trim_pct)
  min_mass = t_mean * 0.25
  max_mass = t_mean * 2.0
  d2 = filter(d, weight >= min_mass & weight <= max_mass)
  cat(nrow(d) - nrow(d2), "\t")
  d2
}


#remove the outliers
dat_weight = group_by(dat$band, spec) %>% 
  do(rm_outlier(.)) %>% ungroup()

# pdf("mass_2.pdf", width = 8, height = 5, onefile = T)
# for(i in unique(dat_weight$spec)){
#   cat(i, "\t")
#   d = filter(dat_weight, spec == i)
#   if(nrow(d) > 1){
#     m = lm(weight ~ yr, d)
#     m2 = coef(summary(m))
#     plot(x = d$yr, y = d$weight, col = d$sex,
#          main = paste0(i, ": slope = ", round(m2[2, "Estimate"], 5), ", p = ", round(m2[2, "Pr(>|t|)"], 4)))
#     abline(m)
#   } else {
#     plot(x = d$yr, y = d$weight, main = i)
#   }
# }
# dev.off()

# library(lme4)
# mod_1 = lmer(mass_g_log10 ~ yr2 + (1|spec) + (1|loc/station) + (0 + yr|spec), data = dat_weight)
# summary(mod_1)


#summary table for each station including # of years sampled, avg. mass, species richness
summ_loc = group_by(dat_weight, loc, station) %>% 
  summarise(n_yr = n_distinct(year(date)),
            n_sp = n_distinct(spec),
            ave_mass = mean(weight, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(rename(dat$stations, loc = LOC, station = STATION), by = c("loc", "station"))#join to station data

?year
#break apart function
a= group_by(dat_weight, loc, station) 
b=summarise(a,n_yr = n_distinct(year(date))),
            n_sp = n_distinct(spec),
            ave_mass = mean(weight, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(rename(dat$stations, loc = LOC, station = STATION), by = c("loc", "station"))#join to station data


n_sp = n_distinct(dat_weight$spec)



#histogram of years sampled across locations
hist(summ_loc$n_yr)


#maps of locations
north_ame = spData::world %>% 
  filter(continent == "North America")
north_ame_p = ggplot() +
  geom_sf(data = north_ame)


#map colored by sampling years
north_ame_p +
  geom_point(data = summ_loc, aes(x = DECLNG, y = DECLAT, color = n_yr))


#map colored by sampling years
north_ame_p +
  geom_point(data = summ_loc, aes(x = DECLNG, y = DECLAT, color = n_sp))


#############################################################################
#############################################################################
#data formatting for OSTATS (FROM data_format.R script)

#works to get species counts per site (combining across years for now) and filter out sites with < 2 species
bird_site_filt<-ddply(dat$band, .(station), mutate,  count = length(unique(spec)))%>%
  filter(count >1)

##stopped here, ready to run ostats


#prep data for OSTATs----
# Load data from Read et al. (2018) from Figshare web archive
#dat <- read.csv('https://ndownloader.figshare.com/files/9167548')


#  filter for only 3 sites to test then select the relevant columns required by the function site, id, svl.
# Use the mutate function to add a new column named "log_SVL" to log-transform the measurements.

dat <- bird_site_filt %>%
  filter(SITE %in% c('10','505', '1468'))%>% 
  select(SITE, ID, SVL) %>%
  filter(!is.na(SVL)) %>%
  mutate(log_SVL = log10(SVL))

# Group the data by siteID and taxonID and look at the summary 
dat %>%
  group_by(SITE, ID) %>%
  slice(1)

#look at data that is input for OSTATS functions
head(dat)


