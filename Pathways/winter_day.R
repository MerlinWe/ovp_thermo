##################################################################################
############################ OVP analysis: winter day ############################
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

# Get winter day subset 
winter <- dat %>% prep_ovp("Winter", "day") # ok
table(winter$ID_phase)
##### ----------- Exploration ---------- #####

colSums(is.na(winter)) # no missing values 

# Check response variables 
plot1 <- check_response_var(winter, "mean_BT_smooth", "Mean BT Smooth")
plot2 <- check_response_var(winter, "mean_activity_percent", "Mean Activity Percent")
plot3 <- check_response_var(winter, "mean_heartrate", "Mean Heartrate")
(plot1 | plot2 | plot3)

rm(plot1, plot2, plot3) # clean environment

## Note: we disregard WCF (wind chill factor) due to high collinearity 
response_vars <- c("mean_BT_smooth", "mean_heartrate", "mean_activity_percent")
predictor_vars <- c("phase_mean_THI", "day_season", "weight", "mean_activity_percent","prev_phase_mean_THI")
## Check linearity
exploration_results <- explore_relationships(winter, response_vars, predictor_vars, random_effect = list(ID = ~1))
edf_summary <- exploration_results$edf_summary  # Extract EDF summary

## Check autocorrelation 
autocorr_BT <- assess_autocorrelation(winter, "mean_BT_smooth", "ID", "season_year", lags = 40)
autocorr_HR <- assess_autocorrelation(winter, "mean_heartrate", "ID", "season_year", lags = 40)
autocorr_AC <- assess_autocorrelation(winter, "mean_activity_percent", "ID", "season_year", lags = 40)

## Check random effect
random_factor_BT <- explore_random_factor(winter, "mean_BT_smooth", "ID")
random_factor_HR <- explore_random_factor(winter, "mean_heartrate", "ID")
random_factor_AC <- explore_random_factor(winter, "mean_activity_percent", "ID")

## Check Interactions
interaction_results <- list(
	"mean_BT_smooth ~ mean_activity_percent * season_year" = test_interaction(winter, "mean_BT_smooth", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ mean_activity_percent * season_year" = test_interaction(winter, "mean_heartrate", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_BT_smooth ~ phase_mean_THI * season_year" = test_interaction(winter, "mean_BT_smooth", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ phase_mean_THI * season_year" = test_interaction(winter, "mean_heartrate", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_activity_percent ~ phase_mean_THI * season_year" = test_interaction(winter, "mean_activity_percent", "phase_mean_THI", "season_year", random_effect = ~1 | ID_phase))

# Get variable descriptives

descriptive_stats <- compute_descriptive_stats(winter, c("mean_BT_smooth", "mean_activity_percent", "mean_heartrate", "phase_mean_THI","prev_phase_mean_THI"))
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
best_random_effects_model_HR <- test_random_effects(winter, "mean_heartrate", random_slope_candidates, edf_summary)
best_random_effects_model_BT <- test_random_effects(winter, "mean_BT_smooth", random_slope_candidates, edf_summary)
best_random_effects_model_ACT <- test_random_effects(winter, "mean_activity_percent", random_slope_candidates, edf_summary)

# Apply updated fixed effect selection
best_fixed_model_HR <- select_fixed_effects(best_random_effects_model_HR, "mean_heartrate", edf_summary)
best_fixed_model_BT <- select_fixed_effects(best_random_effects_model_BT, "mean_BT_smooth", edf_summary)
best_fixed_model_ACT <- select_fixed_effects(best_random_effects_model_ACT, "mean_activity_percent", edf_summary)
#best_fixed_model_ACT gives error message
AIC(best_fixed_model_HR)#2086.142
# Validate top models
validate_HR <- validate_model(best_fixed_model_HR, "mean_heartrate")
validate_BT <- validate_model(best_fixed_model_BT, "mean_BT_smooth")
validate_ACT <- validate_model(best_fixed_model_ACT, "mean_activity_percent")

# Assess model fit
fit_HR <- assess_model_fit(best_fixed_model_HR, "mean_heartrate")
fit_BT <- assess_model_fit(best_fixed_model_BT, "mean_BT_smooth")
fit_ACT <- assess_model_fit(best_fixed_model_ACT, "mean_activity_percent")

# check top models 
top_HR <- best_fixed_model_HR
formula(top_HR)

top_BT <- best_fixed_model_BT
formula(top_BT)

top_ACT <- best_fixed_model_ACT
formula(top_ACT)

winter <- winter %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

glimpse(winter)

# Refit top models using the transformed quadratic terms
best_fixed_model_HR
formula(best_fixed_model_HR)
AIC(best_fixed_model_HR) #2060.126
summary(best_fixed_model_HR)

top_HR1 <- nlme::lme(
  fixed = mean_heartrate ~ season_year + phase_mean_THI+phase_mean_THI_sq+
    + day_season+ day_season_sq+ weight + mean_activity_percent,
  data = winter,
  random = ~1| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR1)
AIC(top_HR1) #2060.158


step_modelHR <- stepAIC(top_HR1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")
summary(step_modelHR) #nothing is removed.

#does the model improve if i remove the square term phase_mean_THI_sq ?
top_HR2 <- nlme::lme(
  fixed = mean_heartrate ~ season_year + phase_mean_THI+
    + day_season+ day_season_sq+ weight + mean_activity_percent,
  data = winter,
  random = ~1| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR2)
AIC(top_HR2) #2060.731  difference 0,6 pts i prefer the simpler model

#remove season_year?
top_HR3 <- nlme::lme(
  fixed = mean_heartrate ~ phase_mean_THI+
    + day_season+ day_season_sq+ weight + mean_activity_percent,
  data = winter,
  random = ~1| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR3)
AIC(top_HR3) # 2062.769 2 pts higher so i won't remove season_year 

top_HR<- top_HR2
summary(top_HR)
AIC(top_HR)

best_fixed_model_BT
formula(best_fixed_model_BT)
AIC(best_fixed_model_BT) #-174.4174

top_BT1 <- nlme::lme(
  fixed = mean_BT_smooth ~ season_year + phase_mean_THI + phase_mean_THI_sq+
    day_season + mean_activity_percent,
  data = winter,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

summary(top_BT1)
AIC(top_BT1) # -174.4174

step_modelBT <- stepAIC(top_BT1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")
summary(step_modelBT) #phase_mean_THI_sq is removed.
AIC(step_modelBT) #-175.7164

#does removing season_yesr (N.S) lead to a better fit?
top_BT2 <- nlme::lme(
  fixed = mean_BT_smooth ~  phase_mean_THI +
    day_season + mean_activity_percent,
  data = winter,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

summary(top_BT2)
AIC(top_BT2) # -175.4768  0.3 pts higher so I will remove season_year
top_BT<-top_BT2
summary(top_BT)
AIC(top_BT)

best_fixed_model_ACT
formula(best_fixed_model_ACT)
AIC(best_fixed_model_ACT) #2018.07

top_ACT1 <- nlme::lme(
  fixed = mean_activity_percent ~ season_year + phase_mean_THI +day_season+weight,
  data = winter,
  random = ~1+day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_ACT1 )# 2018.07

step_modelACT <- stepAIC(top_ACT1, scope = list(lower = ~ phase_mean_THI , upper = ~ .),
                        direction = "both")
summary(step_modelACT) #season_year is removed.
AIC(step_modelACT) #2017.04



#remove day_season (N.S p=0.1630)?
top_ACT2 <- nlme::lme(
  fixed = mean_activity_percent ~ phase_mean_THI +weight,
  data = winter,
  random = ~1 + day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_ACT2)
AIC(top_ACT2) # 2017.083 0.4 pts higher so i will remove day_season.

#remove weight (N.S p=0.2066)?
top_ACT3 <- nlme::lme(
  fixed = mean_activity_percent ~ phase_mean_THI,
  data = winter,
  random = ~1 + day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_ACT3)
AIC(top_ACT3) # 2016.568 lowerr so i will remove day_season.

top_ACT <- top_ACT3

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


# Compute standardized estimates manually

# Get standard deviations of numerical predictors
sd_phase_mean_THI <- sd(winter$phase_mean_THI, na.rm = TRUE)
sd_weight <- sd(winter$weight, na.rm = TRUE)
sd_mean_activity_percent <- sd(winter$mean_activity_percent, na.rm = TRUE)
sd_mean_heartrate <- sd(winter$mean_heartrate, na.rm = TRUE)
sd_day_season <- sd(winter$day_season, na.rm = TRUE)
sd_day_season_sq <- sd(winter$day_season_sq, na.rm = TRUE)

# Extract unstandardized estimates from model output
beta_phase_mean_THI <- 0.0366 
beta_weight <- -0.0931
beta_mean_activity_percent <- 0.2471
beta_day_season <- -0.4546
beta_day_season_sq <- 0.0029

# Compute standardized estimates
std_beta_phase_mean_THI <- beta_phase_mean_THI * (sd_phase_mean_THI / sd_mean_heartrate)
std_beta_weight <- beta_weight * (sd_weight / sd_mean_heartrate)
std_beta_mean_activity_percent <- beta_mean_activity_percent * (sd_mean_activity_percent / sd_mean_heartrate)
std_beta_day_season  <- beta_day_season  * (sd_day_season  / sd_mean_heartrate)
std_beta_day_season_sq <- beta_day_season_sq  * (sd_day_season_sq  / sd_mean_heartrate)


# Print results
std_beta_phase_mean_THI
std_beta_weight
std_beta_mean_activity_percent
std_beta_day_season
std_beta_day_season_sq


#REML topmodels
top_HR_REML <- nlme::lme(
  fixed = mean_heartrate ~ season_year + phase_mean_THI+
    + day_season+ day_season_sq+ weight + mean_activity_percent,
  data = winter,
  random = ~1| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "REML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR_REML)
AIC(top_HR_REML) #2060.731  difference 0,6 pts i prefer the simpler model

top_BT_REML <- nlme::lme(
  fixed = mean_BT_smooth ~  phase_mean_THI +
    day_season + mean_activity_percent,
  data = winter,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "REML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

top_ACT_REML <- nlme::lme(
  fixed = mean_activity_percent ~ phase_mean_THI,
  data = winter,
  random = ~1 + day_season | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "REML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))

psem_model_REML<-psem(top_HR_REML,top_BT_REML,mean_heartrate%~~%mean_BT_smooth,top_ACT_REML)
summary(psem_model_REML)

