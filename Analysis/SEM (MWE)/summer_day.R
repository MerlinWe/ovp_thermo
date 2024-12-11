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

