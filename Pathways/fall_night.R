##################################################################################
############################ OVP analysis: fall day ############################
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
fall <- dat %>% prep_ovp("Fall", "night") # ok

##### ----------- Exploration ---------- #####
table(fall$ID_phase)
id_phase_counts <- table(fall$ID_phase)

# Filter and display only non-zero values
non_zero_id_phase <- id_phase_counts[id_phase_counts > 0]

print(non_zero_id_phase)


colSums(is.na(fall)) # no missing values 

# Check response variables 
plot1 <- check_response_var(fall, "mean_BT_smooth", "Mean BT Smooth")
plot2 <- check_response_var(fall, "mean_activity_percent", "Mean Activity Percent")
plot3 <- check_response_var(fall, "mean_heartrate", "Mean Heartrate")
(plot1 | plot2 | plot3)

rm(plot1, plot2, plot3) # clean environment

## Note: we disregard WCF (wind chill factor) due to high collinearity 
response_vars <- c("mean_BT_smooth", "mean_heartrate", "mean_activity_percent")
predictor_vars <- c("phase_mean_THI", "day_season", "weight", "mean_activity_percent","prev_phase_mean_THI")
## Check linearity
exploration_results <- explore_relationships(fall, response_vars, predictor_vars, random_effect = list(ID = ~1))
edf_summary <- exploration_results$edf_summary  # Extract EDF summary

## Check autocorrelation 
autocorr_BT <- assess_autocorrelation(fall, "mean_BT_smooth", "ID", "season_year", lags = 40)
autocorr_HR <- assess_autocorrelation(fall, "mean_heartrate", "ID", "season_year", lags = 40)
autocorr_AC <- assess_autocorrelation(fall, "mean_activity_percent", "ID", "season_year", lags = 40)

## Check random effect
random_factor_BT <- explore_random_factor(fall, "mean_BT_smooth", "ID")
random_factor_HR <- explore_random_factor(fall, "mean_heartrate", "ID")
random_factor_AC <- explore_random_factor(fall, "mean_activity_percent", "ID")

## Check Interactions
interaction_results <- list(
	"mean_BT_smooth ~ mean_activity_percent * season_year" = test_interaction(fall, "mean_BT_smooth", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ mean_activity_percent * season_year" = test_interaction(fall, "mean_heartrate", "mean_activity_percent", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_BT_smooth ~ phase_mean_THI * season_year" = test_interaction(fall, "mean_BT_smooth", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ phase_mean_THI * season_year" = test_interaction(fall, "mean_heartrate", "phase_mean_THI", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_activity_percent ~ phase_mean_THI * season_year" = test_interaction(fall, "mean_activity_percent", "phase_mean_THI", "season_year", random_effect = ~1 | ID_phase))

# Get variable descriptives
descriptive_stats <- compute_descriptive_stats(fall, c("mean_BT_smooth", "mean_activity_percent", "mean_heartrate", "phase_mean_THI","prev_phase_mean_THI"))
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
best_random_effects_model_HR <- test_random_effects(fall, "mean_heartrate", random_slope_candidates, edf_summary)
best_random_effects_model_BT <- test_random_effects(fall, "mean_BT_smooth", random_slope_candidates, edf_summary)
best_random_effects_model_ACT <- test_random_effects(fall, "mean_activity_percent", random_slope_candidates, edf_summary)

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

fall <- fall %>%
  mutate(
    phase_mean_THI_sq = phase_mean_THI^2,
    prev_phase_max_THI_sq = prev_phase_max_THI^2,
    mean_activity_percent_sq = mean_activity_percent^2,
    day_season_sq = day_season^2,
    weight_sq = weight^2,
    day_season = as.numeric(day_season),
    ID_phase = as.character(ID_phase),
    season_year= as.character(season_year))

glimpse(fall)

# check top models 
#previous stepAIC gave errors
top_HR1 <- nlme::lme(
  fixed = mean_heartrate ~ season_year + weight + phase_mean_THI+ phase_mean_THI_sq +
    + day_season+ day_season_sq+ mean_activity_percent+ mean_activity_percent_sq+
    prev_phase_mean_THI+ prev_phase_max_THI_sq,
  data = fall,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
 summary(top_HR1)
AIC(top_HR1)  # 1898.409

step_modelHR <- stepAIC(top_HR1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both")
summary(step_modelHR) #season_year +day_season_sq  +day_season  + weight
#prev_phase_mean_THI + prev_phase_max_THI_sq + phase_mean_THI_sq+ removed.
#steAIC gave error  but worked until:
#phase_mean_THI + day_season + mean_activity_percent + mean_activity_percent_sq

top_HR2 <- nlme::lme(
  fixed = mean_heartrate ~  day_season+ phase_mean_THI  +
    mean_activity_percent+ mean_activity_percent_sq,
  data = fall,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR2)
AIC(top_HR2 ) #1884.552


#remove day_season (p=0.9316)
top_HR2 <- nlme::lme(
  fixed = mean_heartrate ~  phase_mean_THI  +
    mean_activity_percent+ mean_activity_percent_sq,
  data = fall,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
summary(top_HR2)
AIC(top_HR2 ) #1883.462  lower so I will remove day_season

top_HR<-top_HR2


#step AIC() gave errors for BT model

top_BT1 <- nlme::lme(
  fixed = mean_BT_smooth ~ season_year + weight + phase_mean_THI+ phase_mean_THI_sq +
    + day_season+ day_season_sq+ mean_activity_percent+ mean_activity_percent_sq+
    prev_phase_mean_THI+ prev_phase_max_THI_sq,
  data = fall,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude)#,
 # control = lmeControl(opt = "optim"))
summary(top_BT1)
anova(top_BT1)
AIC(top_BT1) #-358.7657

#control optie uitgeschakeld anders werkt stepAIC niet
step_modelBT <- stepAIC(top_BT1, scope = list(lower = ~ phase_mean_THI + mean_activity_percent, upper = ~ .),
                        direction = "both") #gives error.
summary(step_modelBT) # season_year+weight+  phase_mean_THI_sq + mean_activity_percent_sq  
#  +day_season_sq     are removed.
AIC(step_modelBT) # -369.9787

top_BT2 <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI + mean_activity_percent + prev_phase_mean_THI + 
    prev_phase_max_THI_sq,
  data = fall,
  random = ~1+day_season| ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude)

summary(top_BT2)
anova(top_BT2)
AIC(top_BT2) #-369.9787

#remove prev_phase_max_THI_sq N.S.
top_BT3 <- nlme::lme(
  fixed = mean_BT_smooth ~ phase_mean_THI + mean_activity_percent+
  prev_phase_mean_THI ,
  data = fall,
  random = ~1+ day_season  | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude)#,
  #control = lmeControl(opt = "optim"))

summary(top_BT3)
AIC(top_BT3) # -369.9326 0,04 pts higher so i will remove prev_phase_max_THI_sq 
top_BT<-top_BT3


best_fixed_model_ACT
formula(best_fixed_model_ACT)
AIC(best_fixed_model_ACT) #1758.36

top_ACT1 <- nlme::lme(
  fixed = mean_activity_percent ~ season_year + phase_mean_THI + day_season
  + day_season_sq + prev_phase_mean_THI,
  data = fall,
  random = ~1 + phase_mean_THI  | ID_phase, 
  correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE), 
  method = "ML", 
  na.action = na.exclude,
  control = lmeControl(opt = "optim"))
AIC(top_ACT1 )#1759.197  
summary(top_ACT1)
step_modelACT <- stepAIC(top_ACT1, scope = list(lower = ~ phase_mean_THI, upper = ~ .),
                         direction = "both")
summary(step_modelACT) #season_year is removed
AIC(step_modelACT) #1757.525  much lower so i will remove season_year
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

# Load necessary package
library(dplyr)

# Standard deviation of the response variables
sd_heartrate <- sd(fall$mean_heartrate, na.rm = TRUE)
sd_BT_smooth <- sd(fall$mean_BT_smooth, na.rm = TRUE)
sd_activity <- sd(fall$mean_activity_percent, na.rm = TRUE)

# Standard deviation of predictor variables
sd_THI <- sd(fall$phase_mean_THI, na.rm = TRUE)
sd_THI_sq <- sd(fall$phase_mean_THI^2, na.rm = TRUE)
sd_day_season <- sd(fall$day_season, na.rm = TRUE)
sd_day_season_sq <- sd(fall$day_season^2, na.rm = TRUE)
sd_activity_sq <- sd(fall$mean_activity_percent^2, na.rm = TRUE)

# Load necessary package
library(dplyr)


# Estimates from the table
estimates <- list(
  mean_BT_smooth_phase_mean_THI = -0.0012,
  mean_BT_smooth_mean_activity_percent = -0.0014,
  mean_activity_percent_phase_mean_THI = 0.0457,
  mean_activity_percent_day_season = -0.0165
)

# Compute missing standardized estimates
std_estimates <- data.frame(
  Predictor = c("phase_mean_THI (BT_smooth)", "mean_activity_percent (BT_smooth)",
                "phase_mean_THI (activity_percent)", "day_season (activity_percent)"),
  Estimate = unlist(estimates),
  Std.Estimate = c(
    (-0.0012 * sd_THI / sd_BT_smooth),
    (-0.0014 * sd_activity / sd_BT_smooth),
    (0.0457 * sd_THI / sd_activity),
    (-0.0165 * sd_day_season / sd_activity)
  )
)

# Print results
print(std_estimates)

# Compute standardized estimates for not missing values
std_estimates <- data.frame(
  Predictor = c("phase_mean_THI", "phase_mean_THI_sq", "day_season", "day_season_sq", "mean_activity_percent", "mean_activity_percent_sq"),
  Estimate = c(-0.5305, 0.0046, 0.1766, -0.0022, 0.9293, -0.0285),
  Std.Estimate = c(
    (-0.5305 * sd_THI / sd_heartrate),
    (0.0046 * sd_THI_sq / sd_heartrate),
    (0.1766 * sd_day_season / sd_heartrate),
    (-0.0022 * sd_day_season_sq / sd_heartrate),
    (0.9293 * sd_activity / sd_heartrate),
    (-0.0285 * sd_activity_sq / sd_heartrate)
  )
)
# Print results
print(std_estimates)




