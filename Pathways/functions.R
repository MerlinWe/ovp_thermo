# Function to pre-process OVP data prior to analysis
prep_ovp <- function(data, season, phase) {
	
	data <- data %>%
		# Set variable formats
		mutate(
			season = as.factor(season),
			ID = as.factor(ID),
			phase = as.factor(phase),
			block = as.factor(block),
			year = year(date),
			
			# Redefine season_year to not fail in winter
			season_year = case_when(
				season == "Winter" & month(date) == 12 ~ paste0("Winter_", year(date)),
				season == "Winter" & month(date) %in% 1:2 ~ paste0("Winter_", year(date) - 1),
				TRUE ~ paste0(as.character(season), "_", year(date))
			),
			
			season_year = as.factor(season_year),
			ID_phase = interaction(ID, season_year)
		) %>%
		# Set factor levels
		mutate(
			phase = relevel(phase, ref = "day"),
			season = relevel(season, ref = "Fall")
		) %>%
		# Filter and calculate day_season
		dplyr::filter(season == {{season}}, phase == {{phase}}) %>%
		dplyr::group_by(season_year) %>%
		dplyr::mutate(day_season = match(date, seq(min(date), max(date), by = "days"))) %>%
		dplyr::ungroup() %>%
		# Remove sparse data
		group_by(season_year, ID) %>%
		filter(n() >= 5) %>%
		ungroup()
	
	plot <- data %>% 
		ggplot(aes(x = day_season, y = ID, group = ID, fill = ID, colour = ID)) +
		geom_point() +
		theme_minimal() +
		labs(x = "cumulative day per season", y = NULL, 
				 title = paste({{season}}, {{phase}})) +
		theme(legend.position = "none")
	
	print(plot)
	return(data)
}
	
# Function to explore response vars 
check_response_var <- function(data, variable, var_name) {
	
	# Calculate statistics
	skew <- moments::skewness(data[[variable]], na.rm = TRUE)
	kurt <- moments::kurtosis(data[[variable]], na.rm = TRUE)
	
	# Line plot with labels
	line_plot <- ggplot(data, aes(x = seq_along(data[[variable]]), 
																y = data[[variable]])) +
		geom_point() +
		geom_text(aes(label = seq_along(data[[variable]])), size = 3, vjust = -0.5) +
		theme_minimal() +
		labs(title = paste(var_name), 
				 x = "Case Number", 
				 y = var_name)
	
	# Histogram
	hist_plot <- ggplot(data, aes(x = data[[variable]])) +
		geom_histogram(bins = 30, fill = "blue", color = "white", alpha = 0.7) +
		theme_minimal() +
		labs(title = paste(paste("Skewness:", round(skew, 2)), paste("; Kurtosis:", round(kurt, 2))), 
				 x = var_name, 
				 y = "Frequency")
	
	# Combine plots
	line_plot / hist_plot
}

# Function to assess random factor
explore_random_factor <- function(data, response_vars, random_factor) {
	results <- list()  # Store results
	
	for (var in response_vars) {
		message(sprintf("Exploring random factor effects for: %s", var))
		
		# Boxplot and violin plot combined
		p <- ggplot(data, aes_string(x = random_factor, y = var)) +
			geom_boxplot(outlier.color = "red", alpha = 0.6) +
			geom_violin(fill = "lightblue", alpha = 0.4) +
			labs(
				title = sprintf("Distribution of %s by %s", var, random_factor),
				x = random_factor,
				y = var
			) +
			theme_minimal()
		
		print(p)  # Show plot
		
		# Fit a random-effects model
		formula <- as.formula(sprintf("%s ~ 1 + (1|%s)", var, random_factor))
		model <- lmer(formula, data = data)
		
		# Compute ICC
		icc_value <- performance::icc(model)
		message(sprintf("ICC for %s: %.3f", var, icc_value$ICC))
		
		# Save results
		results[[var]] <- list(
			plot = p,
			model = model,
			ICC = icc_value$ICC
		)
	}
	
	return(results)
}

# Function to check linearity and if a quadratic term is neccessary
explore_relationships <- function(data, response_vars, predictor_vars, random_effect = NULL) {
	
	results <- list()
	edf_results <- data.frame(response = character(),
														predictor = character(),
														EDF = numeric(),
														correlation = numeric(),
														stringsAsFactors = FALSE)
	
	require(mgcv)
	
	for (response in response_vars) {
		# Exclude `mean_activity_percent` as predictor when response is `mean_activity_percent`
		predictors <- if (response == "mean_activity_percent") {
			setdiff(predictor_vars, "mean_activity_percent")
		} else {
			predictor_vars
		}
		
		for (predictor in predictors) {
			cat("\n", strrep("-", 50), "\n")
			cat(sprintf("Exploring %s ~ %s\n", response, predictor))
			
			# Check correlation
			if (is.numeric(data[[response]]) && is.numeric(data[[predictor]])) {
				corr <- cor(data[[response]], data[[predictor]], use = "complete.obs")
				cat(sprintf("Correlation: %.3f\n", corr))
			} else {
				corr <- NA
			}
			
			# Fit GAM model
			if (is.numeric(data[[response]]) && is.numeric(data[[predictor]])) {
				gam_formula <- as.formula(sprintf("%s ~ s(%s)", response, predictor))
				
				if (!is.null(random_effect)) {
					gam_model <- gamm(gam_formula, random = random_effect, data = data)
				} else {
					gam_model <- gam(gam_formula, data = data)
				}
				
				edf <- summary(gam_model$gam)$edf[1]  # Extract EDF
				cat(sprintf("EDF: %.3f\n", edf))
			} else {
				gam_model <- NULL
				edf <- NA
			}
			
			# Store results
			results[[sprintf("%s_%s", response, predictor)]] <- list(
				response = response,
				predictor = predictor,
				correlation = corr,
				gam_model = gam_model,
				edf = edf
			)
			
			# Store EDF in summary dataframe
			edf_results <- rbind(edf_results, data.frame(response = response, 
																									 predictor = predictor, 
																									 EDF = edf, 
																									 correlation = corr))
		}
		cat("\n", strrep("=", 50), "\n")
	}
	
	return(list(results = results, edf_summary = edf_results))
}

# Function to assess temporal autocorrelation 
assess_autocorrelation <- function(data, variables, id_col, period_col, lags = 40) {
	results <- list()  # Store results
	
	for (var in variables) {
		message(sprintf("Assessing autocorrelation for variable: %s", var))
		
		# Loop over unique IDs and periods
		for (id in unique(data[[id_col]])) {
			for (period in unique(data[[period_col]])) {
				
				# Subset data for current ID and period
				subset_data <- data %>% filter(!!sym(id_col) == id, !!sym(period_col) == period)
				
				if (nrow(subset_data) > 0) {
					# Compute ACF and PACF
					acf_values <- acf(subset_data[[var]], lag.max = lags, plot = FALSE)
					pacf_values <- pacf(subset_data[[var]], lag.max = lags, plot = FALSE)
					
					# Create a data frame for visualization
					plot_data <- tibble(
						Lag = seq_along(acf_values$acf[-1]),
						ACF = acf_values$acf[-1],
						PACF = pacf_values$acf
					)
					
					# Plot ACF and PACF with a legend
					p <- ggplot(plot_data, aes(x = Lag)) +
						geom_line(aes(y = ACF, color = "ACF")) +
						geom_line(aes(y = PACF, color = "PACF")) +
						scale_color_manual(values = c("ACF" = "blue", "PACF" = "red")) +
						labs(
							title = sprintf("Autocorrelation: %s (ID: %s, Period: %s)", var, id, period),
							x = "Lag",
							y = "Correlation",
							color = "Type"
						) +
						theme_minimal() +
						theme(legend.position = "top")
					
					print(p)  # Show plot
					
					# Save results
					results[[paste(var, id, period, sep = "_")]] <- list(
						variable = var,
						id = id,
						period = period,
						acf = acf_values,
						pacf = pacf_values,
						plot = p
					)
				}
			}
		}
	}
	
	return(results)
}

# Function to assess random factor
explore_relationships <- function(data, response_vars, predictor_vars, random_effect = NULL) {
	
	results <- list()
	edf_results <- data.frame(response = character(),
														predictor = character(),
														EDF = numeric(),
														correlation = numeric(),
														stringsAsFactors = FALSE)
	
	require(mgcv)
	
	for (response in response_vars) {
		# Exclude `mean_activity_percent` as predictor when response is `mean_activity_percent`
		predictors <- if (response == "mean_activity_percent") {
			setdiff(predictor_vars, "mean_activity_percent")
		} else {
			predictor_vars
		}
		
		for (predictor in predictors) {
			cat("\n", strrep("-", 50), "\n")
			cat(sprintf("Exploring %s ~ %s\n", response, predictor))
			
			# Check correlation
			if (is.numeric(data[[response]]) && is.numeric(data[[predictor]])) {
				corr <- cor(data[[response]], data[[predictor]], use = "complete.obs")
				cat(sprintf("Correlation: %.3f\n", corr))
			} else {
				corr <- NA
			}
			
			# Check number of unique values of the predictor
			unique_vals <- length(unique(na.omit(data[[predictor]])))
			
			# Define the smoothing parameter `k`
			k_value <- min(10, unique_vals - 1)  # Ensure k is not greater than unique values
			
			# Fit GAM model only if there are enough unique values
			if (is.numeric(data[[response]]) && is.numeric(data[[predictor]])) {
				if (unique_vals > 5) {  # Ensure enough unique values before smoothing
					gam_formula <- as.formula(sprintf("%s ~ s(%s, k=%d)", response, predictor, k_value))
					
					if (!is.null(random_effect)) {
						gam_model <- gamm(gam_formula, random = random_effect, data = data)
						edf <- summary(gam_model$gam)$edf[1]
					} else {
						gam_model <- gam(gam_formula, data = data)
						edf <- summary(gam_model)$edf[1]
					}
					
					cat(sprintf("EDF: %.3f (k = %d)\n", edf, k_value))
					
				} else {
					# If too few unique values, use linear regression as fallback
					cat("Too few unique values for smoothing, using linear regression instead.\n")
					gam_model <- lm(as.formula(sprintf("%s ~ %s", response, predictor)), data = data)
					edf <- 1  # Linear regression has 1 degree of freedom for the predictor
				}
			} else {
				gam_model <- NULL
				edf <- NA
			}
			
			# Store results
			results[[sprintf("%s_%s", response, predictor)]] <- list(
				response = response,
				predictor = predictor,
				correlation = corr,
				gam_model = gam_model,
				edf = edf
			)
			
			# Store EDF in summary dataframe
			edf_results <- rbind(edf_results, data.frame(response = response, 
																									 predictor = predictor, 
																									 EDF = edf, 
																									 correlation = corr))
		}
		cat("\n", strrep("=", 50), "\n")
	}
	
	return(list(results = results, edf_summary = edf_results))
}

# Function to run interaction test
test_interaction <- function(data, response, predictor, moderator, random_effect = NULL) {
	
	# Dynamically exclude `mean_activity_percent` as a predictor when response is `mean_activity_percent`
	predictors <- if (response == "mean_activity_percent") {
		setdiff(predictor_vars, "mean_activity_percent")
	} else {
		predictor_vars
	}
	
	message(sprintf("Testing interaction: %s ~ %s * %s", response, predictor, moderator))
	
	# Boxplot to check overlap
	p1 <- ggplot(data, aes_string(x = moderator, y = predictor)) +
		geom_boxplot() +
		theme_minimal() +
		labs(title = sprintf("Boxplot of %s by %s", predictor, moderator))
	
	print(p1)
	
	# Scatter plot with interaction
	p2 <- ggplot(data, aes_string(x = predictor, y = response, color = moderator)) +
		geom_point() +
		geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
		theme_minimal() +
		labs(title = sprintf("Interaction of %s and %s on %s", predictor, moderator, response))
	
	print(p2)
	
	# Ensure correct correlation grouping
	if ("ID_phase" %in% colnames(data)) {
		correlation_structure <- corGaus(form = ~ as.numeric(day_season) | ID_phase, nugget = TRUE)
	} else {
		correlation_structure <- corGaus(form = ~ as.numeric(day_season) | ID.per, nugget = TRUE)
	}
	
	# Fit model
	formula <- as.formula(sprintf("%s ~ %s * %s", response, predictor, moderator))
	use_lme <- !is.null(random_effect)  # Use lme if a random effect is provided
	
	model <- tryCatch(
		{
			if (use_lme) {
				lme(formula, 
						random = random_effect, 
						method = "ML",
						control = lmeControl(opt = "optim"),  # More stable optimizer
						data = data, 
						correlation = correlation_structure, 
						na.action = na.exclude)
			} else {
				lm(formula, data = data)
			}
		},
		error = function(e) {
			message("Error fitting model. Trying without correlation structure...")
			if (use_lme) {
				return(lme(formula, random = random_effect, method = "ML", data = data, na.action = na.exclude))
			} else {
				return(lm(formula, data = data))
			}
		}
	)
	
	# Model summary
	model_summary <- summary(model)
	print(model_summary)
	
	# ANOVA
	anova_res <- anova(model)
	print(anova_res)
	
	# Extract p-value for interaction term (handling name issues)
	interaction_term <- paste0(predictor, ":", moderator)
	interaction_p <- NA
	
	# Check if interaction term exists in ANOVA output
	if (interaction_term %in% rownames(anova_res)) {
		interaction_p <- anova_res[interaction_term, "p-value"]
	} else {
		# Try alternative name formats for factors
		interaction_term_alt <- grep(":", rownames(anova_res), value = TRUE)
		if (length(interaction_term_alt) > 0) {
			interaction_p <- anova_res[interaction_term_alt, "p-value"]
		}
	}
	
	# AIC
	aic_value <- AIC(model)
	message(sprintf("AIC: %.2f", aic_value))
	
	# Extract interaction coefficient, handling different model types
	interaction_coef <- NA
	interaction_se <- NA
	
	if (use_lme) {
		# For lme models, use tTable
		coef_table <- model_summary$tTable
		if (interaction_term %in% rownames(coef_table)) {
			interaction_coef <- coef_table[interaction_term, "Value"]
			interaction_se <- coef_table[interaction_term, "Std.Error"]
		}
	} else {
		# For lm models, use coefficients
		coef_table <- coef(model)
		if (interaction_term %in% names(coef_table)) {
			interaction_coef <- coef_table[interaction_term]
			interaction_se <- summary(model)$coefficients[interaction_term, "Std. Error"]
		}
	}
	
	# Interpretation message
	message("\n--- Interpretation ---")
	if (!is.na(interaction_p)) {
		if (interaction_p < 0.05) {
			message("The interaction term is statistically significant (p < 0.05). It should likely be included in the model.")
		} else if (interaction_p < 0.1) {
			message(" The interaction term is borderline significant (0.05 < p < 0.1). Consider including it if the effect size is meaningful.")
		} else {
			message("The interaction term is **not statistically significant (p > 0.1)**. Consider removing it from the model.")
		}
	} else {
		message("The interaction term was not found in the model. Check variable names.")
	}
	
	# Check if effect size is meaningful
	if (!is.na(interaction_coef)) {
		if (abs(interaction_coef) < interaction_se * 2) {
			message("The interaction effect is small relative to its standard error. The effect may not be practically meaningful.")
		} else {
			message("The interaction effect size appears meaningful.")
		}
	}
	
	return(list(
		model = model, 
		anova = anova_res, 
		AIC = aic_value, 
		interaction_p = interaction_p,
		interaction_coef = interaction_coef,
		interaction_se = interaction_se,
		plots = list(p1, p2)
	))
}


# Function to compute descriptive statistics
compute_descriptive_stats <- function(data, variables) {
	stats <- data.frame(
		Variable = character(),
		Min = numeric(),
		Max = numeric(),
		Mean = numeric(),
		SD = numeric(),
		CV_Percent = numeric(),
		Q2_5 = numeric(),
		Q97_5 = numeric(),
		stringsAsFactors = FALSE
	)
	
	for (var in variables) {
		if (var %in% colnames(data)) {
			min_val <- min(data[[var]], na.rm = TRUE)
			max_val <- max(data[[var]], na.rm = TRUE)
			mean_val <- mean(data[[var]], na.rm = TRUE)
			sd_val <- sd(data[[var]], na.rm = TRUE)
			cv_val <- (sd_val / mean_val) * 100
			quantiles <- quantile(data[[var]], probs = c(0.025, 0.975), na.rm = TRUE)
			
			stats <- rbind(stats, data.frame(
				Variable = var,
				Min = round(min_val, 2),
				Max = round(max_val, 2),
				Mean = round(mean_val, 2),
				SD = round(sd_val, 2),
				CV_Percent = round(cv_val, 2),
				Q2_5 = round(quantiles[1], 2),
				Q97_5 = round(quantiles[2], 2)
			), make.row.names = FALSE)  # Ensure row names do not get inherited
		} else {
			message(sprintf("Warning: Variable '%s' not found in dataset.", var))
		}
	}
	
	rownames(stats) <- NULL  # Explicitly reset row names
	return(stats)
}

# Test random effects and random slopes
test_random_effects <- function(data, response, random_slope_candidates, edf_summary) {
	message(sprintf("Testing random effects for: %s", response))
	
	# Dynamically exclude `mean_activity_percent` as a predictor when response is `mean_activity_percent`
	predictors <- if (response == "mean_activity_percent") {
		setdiff(predictor_vars, "mean_activity_percent")
	} else {
		predictor_vars
	}
	
	# Identify quadratic terms based on EDF
	relevant_edfs <- edf_summary[edf_summary$response == response, ]
	quadratic_terms <- relevant_edfs$predictor[relevant_edfs$EDF > 1.5]  # Use cutoff EDF > 1.5
	
	# Construct model formula dynamically
	model_formula <- paste(response, "~ season_year", collapse = " ")
	
	# Add predictors dynamically based on EDF insights
	if ("mean_activity_percent" %in% relevant_edfs$predictor) {
		if ("mean_activity_percent" %in% quadratic_terms) {
			model_formula <- paste(model_formula, "+ poly(mean_activity_percent,2)")
		} else {
			model_formula <- paste(model_formula, "+ mean_activity_percent")
		}
	}
	
	if ("phase_mean_THI" %in% relevant_edfs$predictor) {
		if ("phase_mean_THI" %in% quadratic_terms) {
			model_formula <- paste(model_formula, "+ poly(phase_mean_THI,2)")
		} else {
			model_formula <- paste(model_formula, "+ phase_mean_THI")
		}
	}
	
	if ("day_season" %in% relevant_edfs$predictor) {
		if ("day_season" %in% quadratic_terms) {
			model_formula <- paste(model_formula, "+ poly(day_season,2)")
		} else {
			model_formula <- paste(model_formula, "+ day_season")
		}
	}
	
	model_formula <- as.formula(model_formula)
	
	# **Step 1: Fit the baseline model with only a random intercept**
	base_model <- lme(
		fixed = model_formula,
		random = ~ 1 | ID_phase,  # Random intercept only
		method = "REML",
		data = data,
		control = lmeControl(opt = "optim"),  # More stable optimizer
		correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE),
		na.action = na.exclude
	)
	
	best_model <- base_model
	best_AIC <- AIC(base_model)
	best_structure <- "~ 1 | ID_phase"  # Baseline: Random intercept only
	
	# **Step 2: Loop through candidate random slopes**
	for (var in random_slope_candidates) {
		message(sprintf("Testing random slope for: %s", var))
		
		# **Step 3: Try fitting a model with a random slope**
		model_candidate <- tryCatch(
			lme(
				fixed = model_formula,
				random = as.formula(sprintf("~ 1 + %s | ID_phase", var)),  # Adding random slope
				method = "REML",
				data = data,
				control = lmeControl(opt = "optim"),  # More stable optimizer
				correlation = corGaus(form = ~ day_season | ID_phase, nugget = TRUE),
				na.action = na.exclude
			),
			error = function(e) {
				message(sprintf("Model with random slope %s failed. Skipping...", var))
				return(NULL)
			}
		)
		
		# **Step 4: Compare AIC values**
		if (!is.null(model_candidate)) {
			model_AIC <- AIC(model_candidate)
			message(sprintf("AIC for model with %s as random slope: %.2f", var, model_AIC))
			
			# **Step 5: If AIC is lower, update best model**
			if (model_AIC < best_AIC) {
				best_model <- model_candidate
				best_AIC <- model_AIC
				best_structure <- sprintf("~ 1 + %s | ID_phase", var)
			}
		}
	}
	
	message(sprintf("Best random structure: %s (AIC: %.2f)", best_structure, best_AIC))
	return(best_model)
}

# Fixed effect structure model selection 
select_fixed_effects <- function(model, response, edf_summary) {
	message(sprintf("Performing stepwise AIC selection (but keeping SEM relevant vars in scope) for: %s", response))
	
	# Dynamically exclude `mean_activity_percent` as a predictor when response is `mean_activity_percent`
	predictors <- if (response == "mean_activity_percent") {
		setdiff(predictor_vars, "mean_activity_percent")
	} else {
		predictor_vars
	}
	
	# Extract relevant EDF values for this response variable
	relevant_edfs <- edf_summary[edf_summary$response == response, ]
	quadratic_terms <- relevant_edfs$predictor[relevant_edfs$EDF > 1.5]  
	
	# Construct full model with necessary quadratic terms
	formula_terms <- c("season_year")  # Always included
	
	for (predictor in predictors) {
		if (predictor %in% quadratic_terms) {
			formula_terms <- c(formula_terms, paste0("poly(", predictor, ",2)"))
		} else {
			formula_terms <- c(formula_terms, predictor)
		}
	}
	
	# Full model formula
	full_formula <- as.formula(sprintf("%s ~ %s", response, paste(formula_terms, collapse = " + ")))
	
	# Define the minimum model (lower bound of scope)
	min_model_terms <- c("season_year")
	
	# Ensure `mean_activity_percent` is included and match the quadratic form if necessary
	if ("mean_activity_percent" %in% quadratic_terms) {
		min_model_terms <- c(min_model_terms, "poly(mean_activity_percent,2)")
	} else {
		min_model_terms <- c(min_model_terms, "mean_activity_percent")
	}
	
	# Ensure `mean_activity_percent` is not included when it is the response variable
	if (response == "mean_activity_percent") {
		min_model_terms <- setdiff(min_model_terms, "mean_activity_percent")
	}
	
	# Ensure `phase_mean_THI` is included and match the quadratic form if necessary
	if ("phase_mean_THI" %in% quadratic_terms) {
		min_model_terms <- c(min_model_terms, "poly(phase_mean_THI,2)")
	} else {
		min_model_terms <- c(min_model_terms, "phase_mean_THI")
	}
	
	lower_formula <- as.formula(sprintf("%s ~ %s", response, paste(min_model_terms, collapse = " + ")))
	
	# Fit full model using ML for AIC-based model selection
	full_model <- lme(
		fixed = full_formula,
		random = formula(model$modelStruct$reStruct),
		data = model$data,
		method = "ML",
		control = lmeControl(opt = "optim"),  # More stable optimizer
		correlation = model$modelStruct$corStruct,
		na.action = na.exclude
	)
	
	# Perform stepwise selection with corrected scope
	step_model <- stepAIC(full_model, scope = list(lower = lower_formula, upper = full_formula), direction = "both", trace = TRUE)
	
	message("Final Model Summary:")
	print(summary(step_model))
	
	return(step_model)
}

# Validate model
validate_model <- function(model, response, data, id_var) {
	message(sprintf("Validating model for: %s", response))
	
	# Extract residuals
	model_residuals <- resid(model, type = "normalized")
	fitted_values <- fitted(model, type = "normalized")
	
	# Add residuals to dataset
	data$residuals <- model_residuals
	data$fitted_values <- fitted_values
	
	# Overall residual histogram
	hist(model_residuals, main = paste("Residuals for", response), xlab = "Residuals", col = "gray", breaks = 30)
	
	# Skewness and Kurtosis
	library(moments)
	skew_val <- skewness(model_residuals, na.rm = TRUE)
	kurt_val <- kurtosis(model_residuals, na.rm = TRUE)
	message(sprintf("Skewness: %.3f, Kurtosis: %.3f", skew_val, kurt_val))
	
	# Residual vs. Fitted plot
	plot(fitted_values, model_residuals, main = paste("Fitted vs. Residuals -", response),
			 xlab = "Fitted Values", ylab = "Residuals", pch = 19, col = "blue")
	abline(h = 0, lty = 2, col = "red")
	
	# Get list of unique individuals
	unique_ids <- unique(data[[id_var]])
	num_ids <- length(unique_ids)
	
	#### ACF Plots ####
	rows <- ceiling(sqrt(num_ids))  
	cols <- ceiling(num_ids / rows)  
	par(mfrow = c(rows, cols))
	
	for (i in unique_ids) {
		id_data <- subset(data, data[[id_var]] == i)
		
		if (nrow(id_data) > 10) {
			acf(id_data$residuals, main = paste("ACF -", response, "ID:", i), lag.max = 40)
		} else {
			message(sprintf("Skipping ACF for ID %s due to insufficient data.", i))
		}
	}
	
	par(mfrow = c(1, 1))  # Reset layout
	
	#### PACF Plots ####
	par(mfrow = c(rows, cols))
	
	for (i in unique_ids) {
		id_data <- subset(data, data[[id_var]] == i)
		
		if (nrow(id_data) > 10) {
			pacf(id_data$residuals, main = paste("PACF -", response, "ID:", i), lag.max = 40)
		} else {
			message(sprintf("Skipping PACF for ID %s due to insufficient data.", i))
		}
	}
	
	par(mfrow = c(1, 1))  # Reset layout
	
	return(list(
		residuals = model_residuals,
		skewness = skew_val,
		kurtosis = kurt_val
	))
}

# Assess model fit 
assess_model_fit <- function(model, response) {
	message(sprintf("Assessing model fit for: %s", response))
	
	# Extract Fitted vs. Observed
	fitted_values <- fitted(model, type = "normalized")
	observed_values <- model$data[[response]]
	
	# Compute R² (Squared Correlation)
	model_r2 <- cor(fitted_values, observed_values, use = "complete.obs")^2
	message(sprintf("R² (correlation squared): %.3f", model_r2))
	
	# Model AIC & BIC
	model_AIC <- AIC(model)
	model_BIC <- BIC(model)
	message(sprintf("AIC: %.2f, BIC: %.2f", model_AIC, model_BIC))
	
	# Scatter plot of Observed vs. Fitted Values
	plot(observed_values, fitted_values, 
			 main = paste("Observed vs. Fitted Values -", response),
			 xlab = "Observed", ylab = "Fitted", pch = 19, col = "blue")
	abline(a = 0, b = 1, col = "red", lty = 2)
	
	return(list(
		R2 = model_r2,
		AIC = model_AIC,
		BIC = model_BIC
	))
}
