rm(list = ls()); gc() # make sure environment is clean 
export = TRUE # Export?

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

# Set WD and get source code functions 
setwd("C:/Users/kuipe_003/OneDrive - Van Hall Larenstein/project OVP WTN en MeW/model THI Merlin")
#setwd("/Users/serpent/Documents/VHL/OVP/Code/Analysis")
#source(url("https://raw.githubusercontent.com/MerlinWe/ovp_thermo/main/Analysis/SEM_(MWE)/functions.R"))
source("C:/Users/kuipe_003/OneDrive - Van Hall Larenstein/project OVP WTN en MeW/model THI Merlin/SEM_(MWE)/functions.R")

# Read data
dat <- read_csv("ovp_data_26_03_25.csv")

#summer day
summer_day <- dat %>% prep_ovp("Summer", "day") # ok
top_HR_summer_day <- nlme::lme(
  fixed = mean_heartrate ~ phase_mean_THI + 
    day_season + mean_activity_percent,
  data = summer_day,
  random = ~1 + mean_activity_percent | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML" )
AIC(top_HR_summer_day ) #3859.477

Summer_day_ACT_HR <- plot(effect(top_HR_summer_day, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Summer_day_ACT_HR)

# summmer_night
summer_night <- dat %>% prep_ovp("Summer", "night") 
top_HR_summer_night <- nlme::lme(
  fixed = mean_heartrate ~ phase_mean_THI+ 
    day_season +  + weight + mean_activity_percent,
  data = summer_night,
  random = ~1 + mean_activity_percent | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC (top_HR_summer_night) #2587.597

Summer_night_ACT_HR <- plot(effect(top_HR_summer_night, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Summer_night_ACT_HR)



#winter_day
winter_day <- dat %>% prep_ovp("Winter", "day")
winter_day <- winter_day %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))


top_HR_winter_day <- nlme::lme(
  fixed = mean_heartrate ~ season_year + phase_mean_THI+
    + day_season+ day_season_sq+ weight + mean_activity_percent,
  data = winter_day,
  random = ~1| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_HR_winter_day) #2060.731

Winter_day_ACT_HR <- plot(effect(top_HR_winter_day, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Winter_day_ACT_HR)

#winter_night
winter_night <- dat %>% prep_ovp("Winter", "night")
winter_night <- winter_night %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    prev_phase_max_THI_sq = prev_phase_max_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

top_HR_winter_night <- nlme::lme(
  fixed = mean_heartrate ~ season_year + phase_mean_THI +
    day_season+ day_season_sq+ weight+ mean_activity_percent,
  data = winter_night,
  random = ~1+phase_mean_THI| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_HR_winter_night) #2340.535

Winter_night_ACT_HR <- plot(effect(top_HR_winter_night, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Winter_night_ACT_HR)


#spring day
spring_day <- dat %>% prep_ovp("Spring", "day")
spring_day <- spring_day %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

top_HR_spring_day <- nlme::lme(
  fixed = mean_heartrate ~  phase_mean_THI +
    + day_season+ day_season_sq+ mean_activity_percent,
  data = spring_day,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_HR_spring_day) #3701.214

Spring_day_ACT_HR <- plot(effect(top_HR_spring_day, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Spring_day_ACT_HR)

#spring night
spring_night <- dat %>% prep_ovp("Spring", "night")


top_HR_spring_night <- nlme::lme(
  fixed = mean_heartrate ~ phase_mean_THI +
    + day_season+ mean_activity_percent,
  data = spring_night,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_HR_spring_night ) # 3103.69

Spring_night_ACT_HR <- plot(effect(top_HR_spring_night, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Spring_night_ACT_HR)

#Autumn day

autumn_day <- dat %>% prep_ovp("Fall", "day") 
autumn_day <- autumn_day %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))


top_HR_autumn_day <- nlme::lme(
  fixed = mean_heartrate ~  phase_mean_THI +
    + day_season+ day_season_sq+ poly(mean_activity_percent,2),
  data = autumn_day,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_HR_autumn_day)  # AIC= 1949.236 
plot(allEffects(top_HR_autumn_day))

Autumn_day_ACT_HR <- plot(effect(top_HR_autumn_day, term="mean_activity_percent"),confint=list(style="none"),
                          colors="blue",  rug=F, 
                          main="")
plot(Autumn_day_ACT_HR)

#Autumn night

autumn_night <- dat %>% prep_ovp("Fall", "night") 
autumn_night <- autumn_night %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    prev_phase_max_THI_sq = prev_phase_max_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

top_HR_autumn_night <- nlme::lme(
  fixed = mean_heartrate ~  phase_mean_THI  +
    poly(mean_activity_percent,2),
  data = autumn_night,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

AIC(top_HR_autumn_night ) #1883.462  lower so I will remove day_season
Autumn_night_ACT_HR <- plot(effect(top_HR_autumn_night, term="mean_activity_percent"),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")
plot(Autumn_night_ACT_HR)

Autumn_night_ACT_HR <- plot(effect(top_HR_autumn_night, term="mean_activity_percent",
                                   xlevels = list(mean_activity_percent = seq(0, 23, length.out = 100)) ),confint=list(style="none"),
                            colors="blue",  rug=F, 
                            main="")





#combined plot

library(tidyverse)
library(nlme)
library(effects)
library(ggplot2)

# Get effects for summer night and summer day
Summer_night_ACT_HR <- as.data.frame(effect(top_HR_summer_night, term = "mean_activity_percent"))
Summer_night_ACT_HR$time_period <- "Summer Night"

Summer_day_ACT_HR <- as.data.frame(effect(top_HR_summer_day, term = "mean_activity_percent"))
Summer_day_ACT_HR$time_period <- "Summer Day"

Winter_day_ACT_HR <- as.data.frame(effect(top_HR_winter_day, term = "mean_activity_percent"))
Winter_day_ACT_HR$time_period <- "Winter Day"

Winter_night_ACT_HR <- as.data.frame(effect(top_HR_winter_night, term = "mean_activity_percent"))
Winter_night_ACT_HR$time_period <- "Winter Night"

Spring_day_ACT_HR <- as.data.frame(effect(top_HR_spring_day, term = "mean_activity_percent"))
Spring_day_ACT_HR$time_period <- "Spring Day"

Spring_night_ACT_HR <- as.data.frame(effect(top_HR_spring_night, term = "mean_activity_percent"))
Spring_night_ACT_HR$time_period <- "Spring Night"

Autumn_day_ACT_HR <- as.data.frame(effect(top_HR_autumn_day, term = "mean_activity_percent",xlevels = list(mean_activity_percent = seq(0, 22, length.out = 100))))
Autumn_day_ACT_HR$time_period <- "Autumn Day"

Autumn_night_ACT_HR <- as.data.frame(effect(top_HR_autumn_night, term = "mean_activity_percent",xlevels = list(mean_activity_percent = seq(0, 23, length.out = 100))))
Autumn_night_ACT_HR$time_period <- "Autumn Night"


# Combine data
combined_data <- bind_rows(Summer_night_ACT_HR, Summer_day_ACT_HR, Winter_day_ACT_HR,Winter_night_ACT_HR,
                           Spring_day_ACT_HR,Spring_night_ACT_HR, Autumn_day_ACT_HR, Autumn_night_ACT_HR )
# Plot
ggplot(combined_data, aes(x = mean_activity_percent, y = fit, color = time_period)) +
  geom_line(size = 1.2)+
  scale_color_manual(values = c("Summer Night" = "yellow3", "Summer Day" = "yellow1","Winter Day" = "grey","Winter Night" = "black",
                                "Spring Day" = "pink1", "Spring Night" = "pink3", "Autumn Day" = "brown2", "Autumn Night" = "brown4" )) +
  scale_fill_manual(values = c("Summer Night" = "yellow3", "Summer Day" = "yellow1","Winter Day" = "grey","Winter Night" = "black",
                               "Spring Day" = "pink1","Spring Night" = "pink3",  "Autumn Day" = "brown2","Autumn Night" = "brown4")) +
  labs(
    x = "mean_activity_percent",
    y = "Mean Heart Rate",
    color = "Time Period",
    fill = "Time Period"
  ) +
  theme_minimal()


