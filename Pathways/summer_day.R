##################################################################################
############################ OVP analysis: summer day ############################
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
library(nlme)
library(piecewiseSEM)
library(tidyverse)

# Set WD and get source code functions 
setwd("C:/Users/kuipe_003/OneDrive - Van Hall Larenstein/project OVP WTN en MeW/model THI Merlin")
#setwd("/Users/serpent/Documents/VHL/OVP/Code/Analysis")
#source(url("https://raw.githubusercontent.com/MerlinWe/ovp_thermo/main/Analysis/SEM_(MWE)/functions.R"))
source("C:/Users/kuipe_003/OneDrive - Van Hall Larenstein/project OVP WTN en MeW/model THI Merlin/SEM_(MWE)/functions.R")

# Read data
dat <- read_csv("ovp_data_26_03_25.csv")

# Get summer day subset 
summer <- dat %>% prep_ovp("Summer", "day") # ok

##### ----------- Exploration ---------- #####

colSums(is.na(summer)) # no missing values 

# Check response variables 
plot1 <- check_response_var(summer, "mean_BT_smooth", "Mean BT Smooth")
plot2 <- check_response_var(summer, "mean_activity_percent", "Mean Activity Percent")
plot3 <- check_response_var(summer, "mean_heartrate", "Mean Heartrate")
(plot1 | plot2 | plot3)

rm(plot1, plot2, plot3) # clean environment

## Note: we disregard WCF (wind chill factor) due to high collinearity 
response_vars <- c("mean_BT_smooth", "mean_heartrate", "mean_activity_percent")
predictor_vars <- c("phase_mean_THI", "day_season", "weight", "mean_activity_percent","prev_phase_mean_THI")
cor(summer$phase_mean_THI,summer$prev_phase_mean_THI) #0.22
## Check linearity
exploration_results <- explore_relationships(summer, response_vars, predictor_vars, random_effect = list(ID = ~1))
edf_summary <- exploration_results$edf_summary  # Extract EDF summary

## Check autocorrelation 
autocorr_BT <- assess_autocorrelation(summer, "mean_BT_smooth", "ID", "season_year", lags = 40)
autocorr_HR <- assess_autocorrelation(summer, "mean_heartrate", "ID", "season_year", lags = 40)
autocorr_AC <- assess_autocorrelation(summer, "mean_activity_percent", "ID", "season_year", lags = 40)

## Check random effect
random_factor_BT <- explore_random_factor(summer, "mean_BT_smooth", "ID")
random_factor_HR <- explore_random_factor(summer, "mean_heartrate", "ID")
random_factor_AC <- explore_random_factor(summer, "mean_activity_percent", "ID")

## Check Interactions
interaction_results <- list(
	"mean_BT_smooth ~ mean_activity_percent * season_year" = test_interaction(summer, "mean_BT_smooth", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ mean_activity_percent * season_year" = test_interaction(summer, "mean_heartrate", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_BT_smooth ~ phase_mean_THI * season_year" = test_interaction(summer, "mean_BT_smooth", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ phase_mean_THI * season_year" = test_interaction(summer, "mean_heartrate", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_activity_percent ~ phase_mean_THI * season_year" = test_interaction(summer, "mean_activity_percent", "phase_mean_THI", "season_year", random_effect = ~1 | ID_phase))

# Get variable descriptives
descriptive_stats <- compute_descriptive_stats(summer, c("mean_BT_smooth", "mean_activity_percent", "mean_heartrate", "phase_mean_THI","prev_phase_mean_THI"))
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
random_slope_candidates <- c("mean_activity_percent", "phase_mean_THI", "day_season")
best_random_effects_model_HR <- test_random_effects(summer, "mean_heartrate", random_slope_candidates, edf_summary)
best_random_effects_model_BT <- test_random_effects(summer, "mean_BT_smooth", random_slope_candidates, edf_summary)
best_random_effects_model_ACT <- test_random_effects(summer, "mean_activity_percent", random_slope_candidates, edf_summary)

# Apply updated fixed effect selection
best_fixed_model_HR <- select_fixed_effects(best_random_effects_model_HR, "mean_heartrate", edf_summary)
best_fixed_model_BT <- select_fixed_effects(best_random_effects_model_BT, "mean_BT_smooth", edf_summary)
best_fixed_model_ACT <- select_fixed_effects(best_random_effects_model_ACT, "mean_activity_percent", edf_summary)

# Validate top models
validate_HR <- validate_model(best_fixed_model_HR, "mean_heartrate", summer, "ID")
validate_BT <- validate_model(best_fixed_model_BT, "mean_BT_smooth", summer, "ID")
validate_ACT <- validate_model(best_fixed_model_ACT, "mean_activity_percent", summer, "ID")

# Assess model fit
fit_HR <- assess_model_fit(best_fixed_model_HR, "mean_heartrate")
fit_BT <- assess_model_fit(best_fixed_model_BT, "mean_BT_smooth")
fit_ACT <- assess_model_fit(best_fixed_model_ACT, "mean_activity_percent")

# check top models 
top_HR <- best_fixed_model_HR
top_HR$call

top_BT <- best_fixed_model_BT
top_BT$call

top_ACT <- best_fixed_model_ACT
top_ACT$call

##### ----- SEM modelling for summer day ----- #####

# note: piecewiseSEM does not natively support quadratic terms. 
# We therefore linearise them before inclusion.

# Create quadratic terms explicitly in data
summer <- summer %>%
	mutate(
		phase_mean_THI_sq = phase_mean_THI^2,
		day_season_sq = day_season^2,
		weight_sq = weight^2,
		day_season = as.numeric(day_season),
		ID_phase = as.character(ID_phase),
		season_year= as.character(season_year))

glimpse(summer)

best_fixed_model_HR
formula(best_fixed_model_HR)
AIC(best_fixed_model_HR) #3860.588 

# Refit top models using the transformed quadratic terms

top_HR <- nlme::lme(
	fixed = mean_heartrate ~ season_year + phase_mean_THI+ phase_mean_THI_sq +
		day_season + day_season_sq + weight + mean_activity_percent,
	data = summer,
	random = ~1 + mean_activity_percent | ID_phase, 
	correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
	method = "ML", 
	na.action = na.exclude,
	control = lmeControl(opt = "optim"))
summary(top_HR)
anova(top_HR)
AIC(top_HR) #3860.588
#
#does the model improve if i remove the square term phase_mean_THI_sq (p=0.0960) ?
top_HR1 <- nlme::lme(
  fixed = mean_heartrate ~ season_year+phase_mean_THI + 
    day_season + day_season_sq  + weight+ mean_activity_percent ,
  data = summer,
  random = ~1 + mean_activity_percent | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML" )
summary(top_HR1)
AIC(top_HR1) #3860.8; 0,22 pt higher, i prefer the simpler model (without phase_mean_THI_sq N)

step_modelHR <- stepAIC(top_HR1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")

summary(step_modelHR) #season_year en weight zijn verwijderd

AIC(step_modelHR)#3859.097

#does the model improve if i remove the square term day_season_sq  ?
top_HR2 <- nlme::lme(
  fixed = mean_heartrate ~ phase_mean_THI + 
    day_season + mean_activity_percent,
  data = summer,
  random = ~1 + mean_activity_percent | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML" )
summary(top_HR2)
AIC(top_HR2) #3859.477; 0,4 pt higher, i prefer the simpler model (day_season_sq N.S.)

top_HR<-top_HR2

top_BT <- nlme::lme(
	fixed = mean_BT_smooth ~ season_year + phase_mean_THI + phase_mean_THI_sq+
		day_season + day_season_sq + mean_activity_percent,
	data = summer,
	random = ~1 | ID_phase, 
	correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
	method = "ML", 
	na.action = na.exclude,
	control = lmeControl(opt = "optim"))

summary(top_BT)
AIC(top_BT) #-566.6704

step_modelBT <- stepAIC(top_BT, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")

summary(step_modelBT) #phase_mean_THI_sq and season_year are removed
AIC(step_modelBT) #-568.3372

#does the model improve if i remove day_season_sq 
top_BT1 <- nlme::lme(
  fixed = mean_BT_smooth ~  phase_mean_THI + 
    day_season + mean_activity_percent,
  data = summer,
  random = ~1 | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_BT1 ) #-560.9586  much higher!! step_modelBT is topmodel
top_BT<- step_modelBT
summary(top_BT)
#

top_ACT <- nlme::lme(
	fixed = mean_activity_percent ~ season_year+phase_mean_THI + weight+ prev_phase_mean_THI,
	data = summer,
	random = ~1 + phase_mean_THI | ID_phase, 
	correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
	method = "ML", 
	na.action = na.exclude,
	control = lmeControl(opt = "optim"))
summary(top_ACT)
step_modelACT <- stepAIC(top_ACT, scope = list(lower = ~ phase_mean_THI, upper = ~ .),
                        direction = "both")
summary(step_modelACT) #weight is verwijderd
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


# Compute standardized estimates manually

# Compute standard deviations of continuous predictors
sd_phase_mean_THI <- sd(summer$phase_mean_THI, na.rm = TRUE)
sd_prev_phase_mean_THI <- sd(summer$prev_phase_mean_THI, na.rm = TRUE)
sd_mean_activity_percent <- sd(summer$mean_activity_percent, na.rm = TRUE)

# Extract unstandardized estimates from model output
beta_phase_mean_THI <- 0.0493
beta_prev_phase_mean_THI <- 0.0487

# Compute standardized estimates
std_beta_phase_mean_THI <- beta_phase_mean_THI * (sd_phase_mean_THI / sd_mean_activity_percent)
std_beta_prev_phase_mean_THI <- beta_prev_phase_mean_THI * (sd_prev_phase_mean_THI / sd_mean_activity_percent)

# Print results
std_beta_phase_mean_THI
std_beta_prev_phase_mean_THI

path_testa<-psem(modsel1HRtop1,modsel1BTtop1,HR%~~%BT,modsel1ACTtop)


class(top_HR)
class(top_BT)
class(top_ACT)
conflicts()

R.version.string
packageVersion("piecewiseSEM")
packageVersion("nlme")

cor(summer %>% select(phase_mean_THI, phase_max_THI, prev_phase_mean_THI, prev_phase_max_THI))

