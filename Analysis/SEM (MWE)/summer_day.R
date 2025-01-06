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

# ----- Exploration -----

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
results <- explore_relationships(summer, response_vars, predictor_vars, random_effect = list(ID = ~1))

## Check autocorrelation 
autocorr_BT <- assess_autocorrelation(summer, "mean_BT_smooth", "ID", "season_year", lags = 40)
autocorr_HR <- assess_autocorrelation(summer, "mean_heartrate", "ID", "season_year", lags = 40)
autocorr_AC <- assess_autocorrelation(summer, "mean_activity_percent", "ID", "season_year", lags = 40)

## Check random effect
random_factor_BT <- explore_random_factor(summer, "mean_BT_smooth", "ID")
random_factor_HR <- explore_random_factor(summer, "mean_heartrate", "ID")
random_factor_AC <- explore_random_factor(summer, "mean_activity_percent", "ID")

## Check Interactions
