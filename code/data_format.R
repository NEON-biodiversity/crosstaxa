# Load data from Read et al. (2018) from Figshare web archive
dat <- read.csv('https://ndownloader.figshare.com/files/9167548')
# Keep only sites "HARV" and "JORN" using the filter function, and select the relevant columns required  # by the function.
# Use the mutate function to add a new column named "log_weight" to log-transform the measurements.
library(dplyr)
dat <- dat %>%
  filter(siteID %in% c('HARV','JORN')) %>%
  select(siteID, taxonID, weight) %>%
  filter(!is.na(weight)) %>%
  mutate(log_weight = log10(weight))

# Group the data by siteID and taxonID and look at the summary 
dat %>%
  group_by(siteID, taxonID) %>%
  slice(1)

head(dat)