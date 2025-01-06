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
			ID_phase = interaction(phase, ID),
			season_year = paste0(season, "_", format(date, format = "%y"))
		) %>%
		# Set factor levels
		mutate(
			phase = relevel(phase, ref="day"),
			season = relevel(season, ref = "Fall")
		) %>%
		# Get clean day_season counter
		dplyr::filter(season == {{season}}, phase == {{phase}}) %>%
		dplyr::group_by(season, year) %>%
		dplyr::mutate(day_season = match(date, seq(min(date), max(date), by = "days"))) %>%
		dplyr::ungroup() %>%
		
		# Remove data with >5 obs of IDs per season
		group_by(year, ID) %>%                   
		filter(n() >= 5) %>%                       
		ungroup() 
	
	plot <- data %>% 
		ggplot(aes(x=day_season, y=ID, group=ID, fill = ID, colour = ID))+
		geom_point() +
		theme_minimal() +
		labs(x="cumulative day per season", y = NULL, 
				 title = paste({{season}},{{phase}})) +
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

# Function to check linearity 
explore_relationships <- function(data, response_vars, predictor_vars, random_effect = NULL) {
	results <- list()  # To store results for each combination
	
	require(mgcv)
	
	for (response in response_vars) {
		for (predictor in predictor_vars) {
			cat("\n", strrep("-", 50), "\n")  # Add a separator
			cat(sprintf("Exploring %s ~ %s\n", response, predictor))
			
			# Check correlation (numeric only)
			if (is.numeric(data[[response]]) && is.numeric(data[[predictor]])) {
				corr <- cor(data[[response]], data[[predictor]], use = "complete.obs")
				cat(sprintf("Correlation: %.3f\n", corr))
			} else {
				corr <- NA
			}
			
			# Plot relationships with legends
			p <- ggplot(data, aes_string(x = predictor, y = response)) +
				geom_point(alpha = 0.6, color = "gray") +
				geom_smooth(method = "loess", se = FALSE, aes(color = "LOESS"), linetype = "dashed", formula = 'y ~ x') +
				geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE, aes(color = "GAM")) +
				scale_color_manual(values = c("LOESS" = "blue", "GAM" = "red")) +
				labs(color = "Fit Type", 
						 title = sprintf("%s ~ %s", response, predictor),
						 x = predictor, 
						 y = response) +
				theme_minimal()
			print(p)
			
			# Fit GAM (if both variables are numeric)
			if (is.numeric(data[[response]]) && is.numeric(data[[predictor]])) {
				gam_formula <- as.formula(sprintf("%s ~ s(%s)", response, predictor))
				
				if (!is.null(random_effect)) {
					gam_model <- gamm(gam_formula, random = random_effect, data = data)
				} else {
					gam_model <- gam(gam_formula, data = data)
				}
				
				edf <- summary(gam_model$gam)$edf[1]
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
				edf = edf,
				plot = p
			)
		}
		cat("\n", strrep("=", 50), "\n")  # Add a separator for the next response variable
	}
	
	return(results)
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

