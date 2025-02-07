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
predictor_vars <- c("phase_mean_CT", "day_season", "weight", "mean_activity_percent")

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
	"mean_BT_smooth ~ phase_mean_CT * season_year" = test_interaction(summer, "mean_BT_smooth", "phase_mean_CT", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_heartrate ~ phase_mean_CT * season_year" = test_interaction(summer, "mean_heartrate", "phase_mean_CT", "season_year", random_effect = ~1 + mean_activity_percent | ID_phase),
	"mean_activity_percent ~ phase_mean_CT * season_year" = test_interaction(summer, "mean_activity_percent", "phase_mean_CT", "season_year", random_effect = ~1 | ID_phase))

# Get variable descriptives
descriptive_stats <- compute_descriptive_stats(summer, c("mean_BT_smooth", "mean_activity_percent", "mean_heartrate", "phase_mean_CT"))
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
random_slope_candidates <- c("mean_activity_percent", "phase_mean_CT", "day_season")
best_random_effects_model_HR <- test_random_effects(summer, "mean_heartrate", random_slope_candidates, edf_summary)
best_random_effects_model_BT <- test_random_effects(summer, "mean_BT_smooth", random_slope_candidates, edf_summary)
best_random_effects_model_ACT <- test_random_effects(summer, "mean_activity_percent", random_slope_candidates, edf_summary)

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

# check top models 
top_HR <- best_fixed_model_HR
top_HR$call

top_BT <- best_fixed_model_BT
top_BT$call

top_ACT <- best_fixed_model_ACT
top_ACT$call





# !!!!!! Experimental from here on ... 




##### ----- SEM modelling for summer day ----- #####

# note: piecewiseSEM does not natively support quadratic terms. 
# We therefore linearise them before inclusion.

# Create quadratic terms explicitly in data
summer <- summer %>%
	mutate(
		phase_mean_CT_sq = phase_mean_CT^2,
		day_season_sq = day_season^2,
		weight_sq = weight^2)

# Refit top models using the transformed quadratic terms

top_HR <- lme(
	fixed = mean_heartrate ~ season_year + phase_mean_CT + phase_mean_CT_sq + 
		day_season + day_season_sq + weight + mean_activity_percent,
	data = summer,
	random = ~1 + mean_activity_percent | ID_phase, 
	correlation = corGaus(form = ~ as.numeric(day_season) | ID_phase, nugget = TRUE), 
	method = "ML", 
	na.action = na.exclude,
	control = lmeControl(opt = "optim"))

top_BT <- lme(
	fixed = mean_BT_smooth ~ season_year + phase_mean_CT + 
		day_season + day_season_sq + mean_activity_percent,
	data = summer,
	random = ~1 | ID_phase, 
	correlation = corGaus(form = ~ as.numeric(day_season) | ID_phase, nugget = TRUE), 
	method = "ML", 
	na.action = na.exclude,
	control = lmeControl(opt = "optim")
)


top_ACT <- lme(
	fixed = mean_activity_percent ~ season_year + phase_mean_CT + 
		weight + weight_sq,
	data = summer,
	random = ~1 + phase_mean_CT | ID_phase, 
	correlation = corGaus(form = ~ as.numeric(day_season) | ID_phase, nugget = TRUE), 
	method = "ML", 
	na.action = na.exclude,
	control = lmeControl(opt = "optim")
)

# Create the piecewise SEM
psem_model <- psem(
	top_HR,
	top_BT,
	top_ACT)

# **Check model summary**
summary(psem_model)

# **Plot the SEM**
plot(psem_model)


### SEM fails... trying a bayesian sem?

library(brms)
library(bayestestR)

# Prepare dataset (center and scale all continuous variables)
summer <- summer %>%
	mutate(
		phase_mean_CT = scale(phase_mean_CT, center = TRUE, scale = TRUE),
		phase_mean_CT_sq = phase_mean_CT^2,
		day_season = scale(day_season, center = TRUE, scale = TRUE),
		day_season_sq = day_season^2,
		weight = scale(weight, center = TRUE, scale = TRUE),
		weight_sq = weight^2,
		mean_activity_percent = scale(mean_activity_percent, center = TRUE, scale = TRUE)
	)

# Bayesian model for mean_heartrate
bayes_HR <- brm(
	mean_heartrate ~ season_year + phase_mean_CT + phase_mean_CT_sq + 
		day_season + day_season_sq + weight + mean_activity_percent + 
		(1 + mean_activity_percent | ID_phase),  # Random slope for repeated measures
	data = summer,
	family = gaussian(),
	prior = c(set_prior("normal(0,1)", class = "b")),  # Regularizing prior
	chains = 4, cores = 4, iter = 4000, warmup = 1000, 
	control = list(adapt_delta = 0.95)
)

bayes_BT <- brm(
	mean_BT_smooth ~ season_year + phase_mean_CT + 
		day_season + day_season_sq + mean_activity_percent + 
		(1 | ID_phase),  # Random intercept
	data = summer,
	family = gaussian(),
	prior = c(set_prior("normal(0,1)", class = "b")),
	chains = 4, cores = 4, iter = 4000, warmup = 1000, 
	control = list(adapt_delta = 0.95)
)

bayes_ACT <- brm(
	mean_activity_percent ~ season_year + phase_mean_CT + 
		weight + weight_sq + 
		(1 + phase_mean_CT | ID_phase),  # Random slope for repeated measures
	data = summer,
	family = gaussian(),
	prior = c(set_prior("normal(0,1)", class = "b")),
	chains = 4, cores = 4, iter = 4000, warmup = 1000, 
	control = list(adapt_delta = 0.95)
)

bayes_SEM <- bf(mean_heartrate ~ season_year + phase_mean_CT + phase_mean_CT_sq + 
									day_season + day_season_sq + weight + mean_activity_percent + 
									(1 + mean_activity_percent | ID_phase)) +
	bf(mean_BT_smooth ~ season_year + phase_mean_CT + 
		 	day_season + day_season_sq + mean_activity_percent + 
		 	(1 | ID_phase)) +
	bf(mean_activity_percent ~ season_year + phase_mean_CT + 
		 	weight + weight_sq + 
		 	(1 + phase_mean_CT | ID_phase))

# Fit the full Bayesian SEM
bayes_model <- brm(
	formula = bayes_SEM,
	data = summer,
	family = gaussian(),
	prior = c(set_prior("normal(0, 1)", class = "b")),
	chains = 4, cores = 4, iter = 8000, warmup = 2000,  # Increase iterations
	control = list(adapt_delta = 0.99, max_treedepth = 15)  # More adaptive sampling
)

summary(bayes_model)
pp_check(bayes_model)
bayestestR::describe_posterior(bayes_model)

# Extract posterior means & 95% credible intervals
coef_df <- as.data.frame(fixef(bayes_model))
coef_df <- coef_df %>%
	mutate(Variable = rownames(coef_df)) %>%
	select(Variable, Estimate = Estimate)

# Store values in a named list for substitution
coef_values <- list(
	HR_BT = round(coef_df$Estimate[coef_df$Variable == "meanheartrate_mean_BT_smooth"], 2),
	HR_ACT = round(coef_df$Estimate[coef_df$Variable == "meanheartrate_mean_activity_percent"], 2),
	BT_ACT = round(coef_df$Estimate[coef_df$Variable == "meanBTsmooth_mean_activity_percent"], 2),
	PHASE_HR = round(coef_df$Estimate[coef_df$Variable == "meanheartrate_phase_mean_CT"], 2),
	PHASE_BT = round(coef_df$Estimate[coef_df$Variable == "meanBTsmooth_phase_mean_CT"], 2),
	PHASE_ACT = round(coef_df$Estimate[coef_df$Variable == "meanactivitypercent_phase_mean_CT"], 2),
	WEIGHT_ACT = round(coef_df$Estimate[coef_df$Variable == "meanactivitypercent_weight"], 2),
	DAY_HR = round(coef_df$Estimate[coef_df$Variable == "meanheartrate_day_season"], 2),
	DAY_BT = round(coef_df$Estimate[coef_df$Variable == "meanBTsmooth_day_season"], 2)
)

# Replace missing coefficients with default values (optional)
coef_values <- lapply(coef_values, function(x) ifelse(is.na(x), "NA", x))

library(DiagrammeR)
library(rsvg)

# Save as SVG
grViz_output <- grViz(sprintf("
digraph BayesianSEM {
  node [shape = box, style = filled, fillcolor = lightblue]

  Mean_Heartrate -> Mean_BT_Smooth [label = '%s']
  Mean_Heartrate -> Mean_Activity_Percent [label = '%s']
  Mean_BT_Smooth -> Mean_Activity_Percent [label = '%s']
  Phase_Mean_CT -> Mean_Heartrate [label = '%s']
  Phase_Mean_CT -> Mean_BT_Smooth [label = '%s']
  Phase_Mean_CT -> Mean_Activity_Percent [label = '%s']
  Weight -> Mean_Activity_Percent [label = '%s']
  Day_Season -> Mean_Heartrate [label = '%s']
  Day_Season -> Mean_BT_Smooth [label = '%s']
}
",
coef_values$HR_BT,
coef_values$HR_ACT,
coef_values$BT_ACT,
coef_values$PHASE_HR,
coef_values$PHASE_BT,
coef_values$PHASE_ACT,
coef_values$WEIGHT_ACT,
coef_values$DAY_HR,
coef_values$DAY_BT
))

# Convert to SVG and save
grViz_output %>%
	export_svg() %>%
	charToRaw() %>%
	rsvg_pdf("BayesianSEM.pdf")

# Open the file manually
browseURL("BayesianSEM.pdf")




