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
			ID_phase = interaction(phase, ID)
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