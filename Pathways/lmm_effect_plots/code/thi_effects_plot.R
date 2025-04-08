library(gdata)
keep(combined_data_ACT, combined_data_HR, combined_data_BT, sure = TRUE)

glimpse(combined_data_ACT)
glimpse(combined_data_HR)
glimpse(combined_data_BT)

combined_data_HR <- combined_data_HR %>%
	mutate(focus_variable = "Heart Rate")

combined_data_ACT <- combined_data_ACT %>%
	mutate(focus_variable = "Activity")

combined_data_BT <- combined_data_BT %>%
	mutate(focus_variable = "Body Temperature")

combined_data_all <- bind_rows(combined_data_HR, combined_data_ACT, combined_data_BT)

combined_data_all <- combined_data_all %>%
	separate(time_period, into = c("season", "timeofday"), sep = " ") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night")),
		focus_variable = factor(focus_variable, levels = c("Heart Rate", "Activity", "Body Temperature"))
	)

ggplot(combined_data_all, aes(x = phase_mean_THI, y = fit, color = season, linetype = timeofday)) +
	geom_line(size = 1.2) +
	scale_color_manual(values = c(
		"Winter" = "grey20",
		"Spring" = "springgreen4",
		"Summer" = "goldenrod2",
		"Autumn" = "brown3"
	)) +
	scale_linetype_manual(values = c("Day" = "solid", "Night" = "dotted")) +
	facet_wrap(~focus_variable, scales = "free_y") +
	labs(
		x = "Phase Mean THI",
		y = "Fitted Value",
		color = "Season",
		linetype = "Time of Day"
	) +
	theme_bw(base_size = 14) +
	theme(
		strip.background = element_rect(fill = "white", color = NA),
		strip.text = element_text(face = "bold")
	)

ggplot(combined_data_all, aes(x = phase_mean_THI, y = fit, color = season)) +
	# Softer THI stress zones
	geom_rect(aes(xmin = -Inf, xmax = 60, ymin = -Inf, ymax = Inf), 
						fill = "grey95", alpha = 1, inherit.aes = FALSE) +
	geom_rect(aes(xmin = 60, xmax = 70, ymin = -Inf, ymax = Inf), 
						fill = "#fffcc2", alpha = 1, inherit.aes = FALSE) +
	geom_rect(aes(xmin = 70, xmax = 80, ymin = -Inf, ymax = Inf), 
						fill = "#ffe0b3", alpha = 1, inherit.aes = FALSE) +
	geom_rect(aes(xmin = 80, xmax = Inf, ymin = -Inf, ymax = Inf), 
						fill = "#ffd6d6", alpha = 1, inherit.aes = FALSE) +
	
	# THI stress labels
	annotate("text", x = 55, y = -Inf, label = "No stress", vjust = -0.5, size = 3.2, color = "black") +
	annotate("text", x = 65, y = -Inf, label = "Mild", vjust = -0.5, size = 3.2, color = "black") +
	annotate("text", x = 75, y = -Inf, label = "Moderate", vjust = -0.5, size = 3.2, color = "black") +
	annotate("text", x = 85, y = -Inf, label = "Severe", vjust = -0.5, size = 3.2, color = "black") +
	
	# Fitted lines
	geom_line(size = 1.4) +
	
	# Updated season color palette
	scale_color_manual(values = c(
		"Winter" = "black",   # green
		"Spring" = "springgreen3",   # orange
		"Summer" = "orangered3",   # purple
		"Autumn" = "goldenrod4"    # magenta
	)) +
	
	facet_grid(focus_variable ~ timeofday, scales = "free_y") +
	labs(
		x = "Mean Temperature Humidity Index (THI) per Phase",
		y = "Fitted Values",
		color = "Season"
	) +
	theme_bw(base_size = 14) +
	theme(
		strip.background = element_rect(fill = "white", color = NA),
		strip.text = element_text(face = "bold"),
		panel.grid.minor = element_blank()
	)




