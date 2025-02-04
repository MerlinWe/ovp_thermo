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
library(tidyverse)

# Set WD and get source code functions 
setwd("/Users/serpent/Documents/VHL/OVP/Code/Analysis")
source("SEM (MWE)/functions.R")

# Read data
dat <- read_csv("/Users/serpent/Documents/VHL/OVP/Data/ovp_data_10_12_24.csv")

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
predictor_vars <- c("phase_mean_CT", "day_season", "weight")

## Check linearity 
lin_results <- explore_relationships(summer, response_vars, predictor_vars, random_effect = list(ID = ~1))

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
	"mean_BT_smooth ~ phase_mean_CT * season_year" = test_interaction(summer, "mean_BT_smooth", "phase_mean_CT", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ phase_mean_CT * season_year" = test_interaction(summer, "mean_heartrate", "phase_mean_CT", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_activity_percent ~ phase_mean_CT * season_year" = test_interaction(summer, "mean_activity_percent", "phase_mean_CT", "season_year", random_effect = ~1 | ID_phase)
	)

# Get variable descriptives
descriptive_stats <- compute_descriptive_stats(summer, c("mean_BT_smooth", "mean_activity_percent", "mean_heartrate", "phase_mean_CT"))
print(descriptive_stats)

##### ----------- SEM Model Selection & Validation ---------- #####

# •	Test random intercept vs. random slopes.
# •	Use AIC/BIC to decide which random effects to retain.

# •	Use stepAIC() to iteratively refine the fixed effect structure.
# •	Keep only significant predictors.

# •	Check residual normality.
# •	Test autocorrelation of residuals.
# •	Assess model fit using observed vs. predicted values.

# •	Fit final models for HR, BT, ACT.
# •	Combine into a structural equation model (SEM) using piecewiseSEM.

# try random slopes for covariates 
random_slope_candidates <- c("mean_activity_percent", "phase_mean_CT", "day_season")
best_random_effects_model_HR <- test_random_effects(summer, "mean_heartrate", random_slope_candidates)
best_random_effects_model_BT <- test_random_effects(summer, "mean_BT_smooth", random_slope_candidates)
best_random_effects_model_ACT <- test_random_effects(summer, "mean_activity_percent", random_slope_candidates)

# Apply to models
best_fixed_model_HR <- select_fixed_effects(best_random_effects_model_HR, "mean_heartrate")
best_fixed_model_BT <- select_fixed_effects(best_random_effects_model_BT, "mean_BT_smooth")
best_fixed_model_ACT <- select_fixed_effects(best_random_effects_model_ACT, "mean_activity_percent")

# Validate top models 
validate_HR <- validate_model(best_fixed_model_HR, "mean_heartrate")
validate_BT <- validate_model(best_fixed_model_BT, "mean_BT_smooth")
validate_ACT <- validate_model(best_fixed_model_ACT, "mean_activity_percent")

# Asses top model fit
fit_HR <- assess_model_fit(best_fixed_model_HR, "mean_heartrate")
fit_BT <- assess_model_fit(best_fixed_model_BT, "mean_BT_smooth")
fit_ACT <- assess_model_fit(best_fixed_model_ACT, "mean_activity_percent")

