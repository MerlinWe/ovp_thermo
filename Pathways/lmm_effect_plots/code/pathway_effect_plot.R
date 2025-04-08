
library(gdata)
keep(combined_data_ACT_BT, combined_data_THI_ACT, combined_data_ACT_HR, sure = TRUE)

# THI → Activity
combined_data_THI_ACT_clean <- combined_data_THI_ACT %>%
	rename(x = phase_mean_THI) %>%
	mutate(response = "Activity")

# Activity → Heart Rate
combined_data_ACT_HR_clean <- combined_data_ACT_HR %>%
	rename(x = mean_activity_percent) %>%
	mutate(response = "Heart Rate")

# Activity → Body Temperature
combined_data_ACT_BT_clean <- combined_data_ACT_BT %>%
	rename(x = mean_activity_percent) %>%
	mutate(response = "Body Temperature")

combined_pathways <- bind_rows(
	combined_data_THI_ACT_clean,
	combined_data_ACT_HR_clean,
	combined_data_ACT_BT_clean
)

combined_pathways <- combined_pathways %>%
	separate(time_period, into = c("season", "timeofday"), sep = " ") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night")),
		response = factor(response, levels = c("Activity", "Heart Rate", "Body Temperature"))
	)


ggplot(combined_pathways, aes(x = x, y = fit, color = season, linetype = timeofday)) +
	geom_line(size = 1.3) +
	facet_wrap(~response, scales = "free_y") +
	scale_color_manual(values = c(
		"Winter" = "#1b9e77", "Spring" = "#d95f02",
		"Summer" = "#7570b3", "Autumn" = "#e7298a"
	)) +
	scale_linetype_manual(values = c("Day" = "solid", "Night" = "dashed")) +
	labs(
		x = "Predictor",
		y = "Fitted Value",
		color = "Season",
		linetype = "Time of Day",
		title = "Pathway Effects of THI and Activity on Physiology"
	) +
	theme_bw(base_size = 14) +
	theme(
		strip.background = element_rect(fill = "white"),
		strip.text = element_text(face = "bold")
	)

ggplot(combined_pathways, aes(x = x, y = fit, color = season)) +
	geom_line(size = 1.3) +
	facet_grid(response ~ timeofday, scales = "free_y") +
	scale_color_manual(values = c(
		"Winter" = "#1b9e77",   # green
		"Spring" = "#d95f02",   # orange
		"Summer" = "#7570b3",   # purple
		"Autumn" = "#e7298a"    # pink
	)) +
	labs(
		x = "Predictor",
		y = "Fitted Value",
		color = "Season",
		title = "Pathway Effects: THI → Activity → HR/BT"
	) +
	theme_bw(base_size = 14) +
	theme(
		strip.background = element_rect(fill = "white"),
		strip.text = element_text(face = "bold"),
		panel.grid.minor = element_blank()
	)



combined_data_activity_paths <- bind_rows(combined_data_ACT_HR, combined_data_ACT_BT)

combined_data_activity_paths <- combined_data_activity_paths %>%
	separate(time_period, into = c("season", "timeofday"), sep = " ") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night")),
		response = factor(response, levels = c("Heart Rate", "Body Temperature"))
	)



library(patchwork)

# --------------------
# THI → Activity
# --------------------
plot_thi_to_act <- combined_data_THI_ACT %>%
	separate(time_period, into = c("season", "timeofday"), sep = " ") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = phase_mean_THI, y = fit, color = season)) +
	geom_line(size = 1.3) +
	facet_wrap(~timeofday) +
	scale_color_manual(values = c(
		"Winter" = "black", "Spring" = "springgreen4",
		"Summer" = "orangered3", "Autumn" = "goldenrod4"
	)) +
	labs(
		x = "Phase Temperature-Humidity Index (THI) per Phase",
		y = "Non-Recumbent (%)",
		color = "Season",
		linetype = "Time of Day",
		title = "Temperature-Humidity-Index → Non-Recumbent (Activity)"
	) +
	theme_bw(base_size = 13) +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		text = element_text(size = 12)
	)

# --------------------
# Activity → Heart Rate
# --------------------
plot_act_to_hr <- combined_data_ACT_HR %>%
	separate(time_period, into = c("season", "timeofday"), sep = " ") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = mean_activity_percent, y = fit, color = season)) +
	geom_line(size = 1.3) +
	facet_wrap(~timeofday) +
	scale_color_manual(values = c(
		"Winter" = "black", "Spring" = "springgreen4",
		"Summer" = "orangered3", "Autumn" = "goldenrod4"
	)) +
	labs(
		x = "Mean Non-Recumbent (%, Activity) per Phase",
		y = "Mean Heart Rate (bpm) per phase",
		color = "Season",
		linetype = "Time of Day",
		title = "Non-Recumbent (Activity) → Heart Rate"
	) +
	theme_bw(base_size = 13) +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		text = element_text(size = 12)
	)

# --------------------
# Activity → Body Temperature
# --------------------
plot_act_to_bt <- combined_data_ACT_BT %>%
	separate(time_period, into = c("season", "timeofday"), sep = " ") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = mean_activity_percent, y = fit, color = season)) +
	geom_line(size = 1.3) +
	facet_wrap(~timeofday) +
	scale_color_manual(values = c(
		"Winter" = "black", "Spring" = "springgreen4",
		"Summer" = "orangered3", "Autumn" = "goldenrod4"
	)) +
	labs(
		x = "Mean Non-Recumbent (%, Activity) per Phase",
		y = "Mean Body Temperature (°C) per phase",
		color = "Season",
		linetype = "Time of Day",
		title = "Non-Recumbent (Activity) → Body Temperature"
	) +
	theme_bw(base_size = 13) +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		text = element_text(size = 12)
	)

# --------------------
# Combine into a single figure
# --------------------
final_pathway_plot <- plot_thi_to_act / plot_act_to_hr / plot_act_to_bt +
	plot_layout(guides = "collect") & theme(legend.position = "right")

# Display the plot
final_pathway_plot

# Save the plot
ggsave(
	filename = "/Users/serpent/Desktop/pathway_plot.png",
	plot = final_pathway_plot,
	width = 260,
	height = 260,
	units = "mm",
	dpi = 300)
