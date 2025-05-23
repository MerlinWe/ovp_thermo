

# Activate package libraries 
library(patchwork)
library(feather)
library(mgcv)
library(car)
library(ggeffects)
library(MASS)
library(emmeans)
library(gamm4)
library(effects)
library(stats)
library(lme4)
library(lmerTest)
library(performance)
library(tidyverse)


setwd("/Users/serpent/Documents/VHL/OVP/Code/Analysis")
source("/Users/serpent/Documents/VHL/OVP/Code/Analysis/SEM_(MWE)/functions.R")
dat <- read_csv("/Users/serpent/Documents/VHL/OVP/Data/ovp_data_26_03_25.csv")

#summer day


summer_day <- dat %>% prep_ovp("Summer", "day") # ok

summer_day <- summer_day %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

top_BT_summer_day <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI+
    day_season + day_season_sq + mean_activity_percent,
  data = summer_day,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_BT_summer_day) #-568.3372

Summer_day_ACT_BT <- plot(effect(top_BT_summer_day, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Summer_day_ACT_BT)

# summmer_night
summer_night <- dat %>% prep_ovp("Summer", "night") 
summer_night <- summer_night %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))


top_BT_summer_night <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI + 
     day_season_sq + mean_activity_percent,
  data = summer_night,
  random = ~1+mean_activity_percent | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_BT_summer_night) #-204.7824
Summer_night_THI_HR <- plot(effect(top_BT_summer_night, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Summer_night_THI_HR)



#winter_day
winter_day <- dat %>% prep_ovp("Winter", "day")

top_BT_winter_day <- nlme::lme(
  fixed = mean_BT_smooth ~  phase_mean_THI +
    day_season + mean_activity_percent,
  data = winter_day,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_BT_winter_day)# -175.4768
Winter_day_ACT_BT <- plot(effect(top_BT_winter_day, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Winter_day_ACT_BT)

#winter_night
winter_night <- dat %>% prep_ovp("Winter", "night")

top_BT_winter_night <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI +
    day_season + mean_activity_percent,
  data = winter_night,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_BT_winter_night) # -285.4662

Winter_night_ACT_BT <- plot(effect(top_BT_winter_night, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Winter_night_ACT_BT)


#spring day
spring_day <- dat %>% prep_ovp("Spring", "day")

top_BT_spring_day <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI +
    day_season  + mean_activity_percent,
  data = spring_day,
  random = ~1+ day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_BT_spring_day)  #-786.7148 

Spring_day_ACT_BT <- plot(effect(top_BT_spring_day , term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Spring_day_ACT_BT)

#spring night
spring_night <- dat %>% prep_ovp("Spring", "night")
spring_night <- spring_night %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))


top_BT_spring_night <- nlme::lme(
  fixed = mean_BT_smooth ~ day_season_sq + phase_mean_THI +
    mean_activity_percent + prev_phase_mean_THI,
  data = spring_night,
  random = ~1 + day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML")

AIC(top_BT_spring_night ) #-443.0661

Spring_night_ACT_BT <- plot(effect(top_BT_spring_night , term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Spring_night_ACT_BT)

#Autumn day

autumn_day <- dat %>% prep_ovp("Fall", "day") 

top_BT_autumn_day <- nlme::lme(
  fixed = mean_BT_smooth ~ season_year + phase_mean_THI+
    mean_activity_percent,
  data = autumn_day,
  random = ~1+ prev_phase_mean_THI  | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_BT_autumn_day) #-384.4325

Autumn_day_ACT_BT <- plot(effect(top_BT_autumn_day, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Autumn_day_ACT_BT)


#Autumn night

autumn_night <- dat %>% prep_ovp("Fall", "night") 

top_BT_autumn_night <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI + mean_activity_percent+
    prev_phase_mean_THI ,
  data = autumn_night,
  random = ~1+ day_season  | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude)#,
#control = lmeControl(opt = "optim"))

AIC(top_BT_autumn_night ) #-369.9326
Autumn_night_THI_HR <- plot(effect(top_BT_autumn_night, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Autumn_night_THI_HR)

#combined plot

library(tidyverse)
library(nlme)
library(effects)
library(ggplot2)

# Get effects for summer night and summer day
Summer_night_ACT_BT <- as.data.frame(effect(top_BT_summer_night, term = "mean_activity_percent"))
Summer_night_ACT_BT$time_period <- "Summer Night"

Summer_day_ACT_BT <- as.data.frame(effect(top_BT_summer_day, term = "mean_activity_percent"))
Summer_day_ACT_BT$time_period <- "Summer Day"

Winter_day_ACT_BT <- as.data.frame(effect(top_BT_winter_day, term = "mean_activity_percent"))
Winter_day_ACT_BT$time_period <- "Winter Day"

Winter_night_ACT_BT <- as.data.frame(effect(top_BT_winter_night, term = "mean_activity_percent"))
Winter_night_ACT_BT$time_period <- "Winter Night"

Spring_day_ACT_BT <- as.data.frame(effect(top_BT_spring_day, term = "mean_activity_percent"))
Spring_day_ACT_BT$time_period <- "Spring Day"

Spring_night_ACT_BT <- as.data.frame(effect(top_BT_spring_night, term = "mean_activity_percent"))
Spring_night_ACT_BT$time_period <- "Spring Night"

Autumn_day_ACT_BT <- as.data.frame(effect(top_BT_autumn_day, term = "mean_activity_percent"))
Autumn_day_ACT_BT$time_period <- "Autumn Day"

Autumn_night_ACT_BT <- as.data.frame(effect(top_BT_autumn_night, term = "mean_activity_percent"))
Autumn_night_ACT_BT$time_period <- "Autumn Night"


# Combine data
combined_data_ACT_BT <- bind_rows(Summer_night_ACT_BT, Summer_day_ACT_BT, Winter_day_ACT_BT,Winter_night_ACT_BT,
                           Spring_day_ACT_BT,Spring_night_ACT_BT, Autumn_day_ACT_BT, Autumn_night_ACT_BT )
