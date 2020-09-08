#clean MAPS BIRD DATA

library(tidyverse)
library(spData)
library(sf)
library(lubridate)
HHH

load("MAPSexport.Rdata")

dat = lapply(maps.data, as_tibble)
rm(maps.data)
dat$band
n_distinct(dat$band$SPEC)
range(dat$band$DATE)
range(dat$band$WEIGHT)
hist(dat$band$WEIGHT)
dat$band = mutate(dat$band, mass_g_log10 = log10(WEIGHT),
                  yr = year(DATE),
                  yr2 = yr - min(yr))
hist(dat$band$mass_g_log10)

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

# remove species with < 10 observations
names(dat$band) = tolower(names(dat$band))
sp_rm = group_by(dat$band, spec) %>% 
  tally() %>% 
  filter(n < 10) %>% 
  pull(spec) %>% as.character()
dat$band = filter(dat$band, !spec %in% sp_rm)
n_distinct(dat$band$spec) # 267 species
# remove outliers
rm_outlier = function(d, trim_pct = 0.1){
  t_mean = mean(d$weight, na.rm = TRUE, trim = trim_pct)
  min_mass = t_mean * 0.25
  max_mass = t_mean * 2.0
  d2 = filter(d, weight >= min_mass & weight <= max_mass)
  cat(nrow(d) - nrow(d2), "\t")
  d2
}

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

summ_loc = group_by(dat_weight, loc, station) %>% 
  summarise(n_yr = n_distinct(year(date)),
            n_sp = n_distinct(spec),
            ave_mess = mean(weight, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(rename(dat$stations, loc = LOC, station = STATION), by = c("loc", "station"))

north_ame = spData::world %>% 
  filter(continent == "North America")
north_ame_p = ggplot() +
  geom_sf(data = north_ame)

hist(summ_loc$n_yr)

north_ame_p +
  geom_point(data = summ_loc, aes(x = DECLNG, y = DECLAT, color = n_yr))

north_ame_p +
  geom_point(data = summ_loc, aes(x = DECLNG, y = DECLAT, color = n_sp))





