##################################################################################
############################ OVP analysis: spring day ############################
##################################################################################

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

# Get spring day subset 
spring <- dat %>% prep_ovp("Spring", "day") # ok

##### ----------- Exploration ---------- #####

colSums(is.na(spring)) # no missing values 

# Check response variables 
plot1 <- check_response_var(spring, "mean_BT_smooth", "Mean BT Smooth")
plot2 <- check_response_var(spring, "mean_activity_percent", "Mean Activity Percent")
plot3 <- check_response_var(spring, "mean_heartrate", "Mean Heartrate")
(plot1 | plot2 | plot3)

rm(plot1, plot2, plot3) # clean environment

## Note: we disregard WCF (wind chill factor) due to high collinearity 
response_vars <- c("mean_BT_smooth", "mean_heartrate", "mean_activity_percent")
predictor_vars <- c("phase_mean_THI", "day_season", "weight", "mean_activity_percent","prev_phase_mean_THI")
## Check linearity
exploration_results <- explore_relationships(spring, response_vars, predictor_vars, random_effect = list(ID = ~1))
edf_summary <- exploration_results$edf_summary  # Extract EDF summary

## Check autocorrelation 
autocorr_BT <- assess_autocorrelation(spring, "mean_BT_smooth", "ID", "season_year", lags = 40)
autocorr_HR <- assess_autocorrelation(spring, "mean_heartrate", "ID", "season_year", lags = 40)
autocorr_AC <- assess_autocorrelation(spring, "mean_activity_percent", "ID", "season_year", lags = 40)

## Check random effect
random_factor_BT <- explore_random_factor(spring, "mean_BT_smooth", "ID")
random_factor_HR <- explore_random_factor(spring, "mean_heartrate", "ID")
random_factor_AC <- explore_random_factor(spring, "mean_activity_percent", "ID")

## Check Interactions
interaction_results <- list(
	"mean_BT_smooth ~ mean_activity_percent * season_year" = test_interaction(spring, "mean_BT_smooth", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ mean_activity_percent * season_year" = test_interaction(spring, "mean_heartrate", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_BT_smooth ~ phase_mean_THI * season_year" = test_interaction(spring, "mean_BT_smooth", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ phase_mean_THI * season_year" = test_interaction(spring, "mean_heartrate", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_activity_percent ~ phase_mean_THI * season_year" = test_interaction(spring, "mean_activity_percent", "phase_mean_THI", "season_year", random_effect = ~1 | ID_phase))

# Get variable descriptives
descriptive_stats <- compute_descriptive_stats(spring, c("mean_BT_smooth", "mean_activity_percent", "mean_heartrate", "phase_mean_THI","prev_phase_mean_THI"))
print(descriptive_stats)

##### ----------- Model Selection & Validation ---------- #####

# •	Test random intercept vs. random slopes.
# •	Use AIC/BIC to decide which random effects to retain.

# •	Use stepAIC() to iteratively refine the fixed effect structure.
# •	Keep only significant predictors.

# •	Check residual normality.
# •	Test autocorrelation of residuals.
# •	Assess model fit using observed vs. predicted values.

# •	Fit final models for HR, BT, ACT.
# •	Combine into a structural equation model (SEM) using piecewiseSEM.

# try random slopes for covariates (dynamic quadratic selection based on EDF)
random_slope_candidates <- c("mean_activity_percent", "phase_mean_THI","prev_phase_mean_THI", "day_season")
best_random_effects_model_HR <- test_random_effects(spring, "mean_heartrate", random_slope_candidates, edf_summary)
best_random_effects_model_BT <- test_random_effects(spring, "mean_BT_smooth", random_slope_candidates, edf_summary)
best_random_effects_model_ACT <- test_random_effects(spring, "mean_activity_percent", random_slope_candidates, edf_summary)

# Apply updated fixed effect selection
best_fixed_model_HR <- select_fixed_effects(best_random_effects_model_HR, "mean_heartrate", edf_summary)
best_fixed_model_BT <- select_fixed_effects(best_random_effects_model_BT, "mean_BT_smooth", edf_summary)
best_fixed_model_ACT <- select_fixed_effects(best_random_effects_model_ACT, "mean_activity_percent", edf_summary)

# Validate top models
validate_HR <- validate_model(best_fixed_model_HR, "mean_heartrate")
validate_BT <- validate_model(best_fixed_model_BT, "mean_BT_smooth")
validate_ACT <- validate_model(best_fixed_model_ACT, "mean_activity_percent")

# Assess model fit
fit_HR <- assess_model_fit(best_fixed_model_HR, "mean_heartrate")
fit_BT <- assess_model_fit(best_fixed_model_BT, "mean_BT_smooth")
fit_ACT <- assess_model_fit(best_fixed_model_ACT, "mean_activity_percent")

spring <- spring %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

glimpse(spring)

# check top models 
best_fixed_model_HR
formula(best_fixed_model_HR)
AIC(best_fixed_model_HR) #3703.336
summary(best_fixed_model_HR)
# Refit top models using the transformed quadratic terms

top_HR1 <- nlme::lme(
  fixed = mean_heartrate ~ season_year + phase_mean_THI+ phase_mean_THI_sq +
    + day_season+ day_season_sq+ mean_activity_percent+ mean_activity_percent_sq,
  data = spring,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR1)
AIC(top_HR1)  #3721.815
#hetzelfde met poly() functie
top_HR2 <- nlme::lme(
  fixed = mean_heartrate ~ season_year + poly(phase_mean_THI, 2) + poly(day_season,2) + 
    poly(mean_activity_percent, 2) ,
  data = spring,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR2)
AIC(top_HR2) # 3721.815

#ik snap niet waarom hier niet hetzelfde als AIC(best_fixed_model_HR) #3703.336 uitkomt.

step_modelHR <- stepAIC(top_HR1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")
summary(step_modelHR) #phase_mean_THI_sq and mean_activity_percent_sq and season_year are removed.

summary(step_modelHR )
AIC(step_modelHR)  #3701.214
top_HR<-step_modelHR

best_fixed_model_BT
formula(best_fixed_model_BT)
AIC(best_fixed_model_BT) #-784.9139

top_BT1 <- nlme::lme(
  fixed = mean_BT_smooth ~ season_year + phase_mean_THI + phase_mean_THI_sq+
    day_season + day_season_sq + mean_activity_percent,
  data = spring,
  random = ~1+ day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

summary(top_BT1)
AIC(top_BT1) # -784.5665

step_modelBT <- stepAIC(top_BT1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")
summary(step_modelBT) #day_season, phase_mean_THI_sq season_year are removed.
#does the model improve if i remove the square term phase_mean_THI_sq ?

summary(step_modelBT )
AIC(step_modelBT)  #-787.5615

#can i replace the quadratic term with a linear term (day_season)

top_BT2 <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI +
    day_season  + mean_activity_percent,
  data = spring,
  random = ~1+ day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

summary(top_BT2)
AIC(top_BT2) # -786.7148 #difference <1 pt so i prefer the linear function for day_season


top_BT<-top_BT2


best_fixed_model_ACT
formula(best_fixed_model_ACT)
AIC(best_fixed_model_ACT) #3513.627

top_ACT1 <- nlme::lme(
  fixed = mean_activity_percent ~ season_year + phase_mean_THI + day_season,
  data = spring,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_ACT1 )#3513.627

step_modelACT <- stepAIC(top_ACT1, scope = list(lower = ~ phase_mean_THI, upper = ~ .),
                        direction = "both")
summary(step_modelACT) #season_year is removed
AIC(step_modelACT) #3511.669
top_ACT<-step_modelACT

summary(top_ACT)
summary(top_BT)
summary(top_HR)
# Create the piecewise SEM
psem_model <- piecewiseSEM::psem(
  top_HR,
  top_BT,
  top_ACT)


summary(psem_model)

psem_model2<-psem(top_HR,top_BT,mean_heartrate%~~%mean_BT_smooth,top_ACT)
summary(psem_model2)
