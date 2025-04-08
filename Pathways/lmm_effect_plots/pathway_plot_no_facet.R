rm(list = ls()); gc()

library(patchwork)
library(tidyverse)

# Color and linetype setup
season_colors <- c(
	"Winter" = "black",
	"Spring" = "springgreen4",
	"Summer" = "orangered3",
	"Autumn" = "goldenrod4"
)

season_color_scale <- scale_color_manual(values = season_colors)
season_linetype_scale <- scale_linetype_manual(values = c("Day" = "solid", "Night" = "dashed"))

# Load data
combined_pathways <- read_csv("/Users/serpent/Documents/VHL/OVP/Code/Pathways/lmm_effect_plots/pathway_effects.csv")
combined_thi      <- read_csv("/Users/serpent/Documents/VHL/OVP/Code/Pathways/lmm_effect_plots/combined_THI.csv")

# THI background
thi_background <- list(
	annotate("rect", xmin = -Inf, xmax = 60, ymin = -Inf, ymax = Inf, fill = "grey95", alpha = 1),
	annotate("rect", xmin = 60, xmax = 70, ymin = -Inf, ymax = Inf, fill = "#fffcc2", alpha = 1),
	annotate("rect", xmin = 70, xmax = 80, ymin = -Inf, ymax = Inf, fill = "#ffe0b3", alpha = 1),
	annotate("rect", xmin = 80, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#ffd6d6", alpha = 1),
	annotate("text", x = 55, y = -Inf, label = "No stress", vjust = -0.5, size = 3.2),
	annotate("text", x = 65, y = -Inf, label = "Mild", vjust = -0.5, size = 3.2),
	annotate("text", x = 75, y = -Inf, label = "Moderate", vjust = -0.5, size = 3.2),
	annotate("text", x = 85, y = -Inf, label = "Severe", vjust = -0.5, size = 3.2)
)

# Plot 1: Activity → Heart Rate
plot_act_to_hr <- combined_pathways %>%
	filter(response == "Heart Rate") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = x, y = fit, color = season, linetype = timeofday)) +
	geom_line(size = 1.3) +
	season_color_scale +
	season_linetype_scale +
	labs(
		x = "Non-Recumbent (%, Activity)", y = "Heart Rate (bpm)",
		color = "Season", linetype = "Time of Day",
		title = "Non-Recumbent (Activity) → Heart Rate"
	) +
	theme_bw(base_size = 13) +
	theme(
		text = element_text(size = 12),
		plot.title = element_text(face = "bold")
	)

# Plot 2: Activity → Body Temperature
plot_act_to_bt <- combined_pathways %>%
	filter(response == "Body Temperature") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = x, y = fit, color = season, linetype = timeofday)) +
	geom_line(size = 1.3) +
	season_color_scale +
	season_linetype_scale +
	labs(
		x = "Non-Recumbent (%, Activity)", y = "Body Temperature (°C)",
		color = "Season", linetype = "Time of Day",
		title = "Non-Recumbent (Activity) → Body Temperature"
	) +
	theme_bw(base_size = 13) +
	theme(
		text = element_text(size = 12),
		plot.title = element_text(face = "bold")
	)

# Plot 3: THI → Heart Rate
plot_thi_hr <- combined_thi %>%
	filter(focus_variable == "Heart Rate") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = phase_mean_THI, y = fit, color = season, linetype = timeofday)) +
	thi_background +
	geom_line(size = 1.2) +
	season_color_scale +
	season_linetype_scale +
	labs(
		x = "Temperature-Humidity-Index", y = "Heart Rate (bpm)",
		title = "Temperature-Humidity-Index → Heart Rate"
	) +
	theme_bw(base_size = 14) +
	theme(
		legend.position = "none",
		panel.grid.minor = element_blank(),
		plot.title = element_text(size = 14, face = "bold"),
		axis.title = element_text(size = 12)
	)

# Plot 4: THI → Activity
plot_thi_act <- combined_thi %>%
	filter(focus_variable == "Activity") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = phase_mean_THI, y = fit, color = season, linetype = timeofday)) +
	thi_background +
	geom_line(size = 1.2) +
	season_color_scale +
	season_linetype_scale +
	labs(
		x = "Temperature-Humidity-Index", y = "Non-Recumbent (%)",
		title = "Temperature-Humidity-Index → Non-Recumbent (Activity)"
	) +
	theme_bw(base_size = 14) +
	theme(
		legend.position = "none",
		panel.grid.minor = element_blank(),
		plot.title = element_text(size = 14, face = "bold"),
		axis.title = element_text(size = 12)
	)

# Plot 5: THI → Body Temperature
plot_thi_bt <- combined_thi %>%
	filter(focus_variable == "Body Temperature") %>%
	mutate(
		season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
		timeofday = factor(timeofday, levels = c("Day", "Night"))
	) %>%
	ggplot(aes(x = phase_mean_THI, y = fit, color = season, linetype = timeofday)) +
	thi_background +
	geom_line(size = 1.2) +
	season_color_scale +
	season_linetype_scale +
	labs(
		x = "Temperature-Humidity-Index", y = "Body Temperature (°C)",
		title = "Temperature-Humidity-Index → Body Temperature"
	) +
	theme_bw(base_size = 14) +
	theme(
		legend.position = "none",
		panel.grid.minor = element_blank(),
		plot.title = element_text(size = 14, face = "bold"),
		axis.title = element_text(size = 12)
	)

# Combine plots
plot_direct <- (plot_act_to_hr / plot_act_to_bt) + 
	plot_layout(guides = "collect") &
	theme(
		legend.position = "bottom",
		legend.background = element_rect(colour = "black")
	)

plot_thi_combined <- (plot_thi_act / plot_thi_bt / plot_thi_hr) +
	plot_layout(guides = "collect", heights = c(1, 1, 1)) &
	theme(legend.position = "none")

pathways_combined <- (plot_thi_combined | plot_direct) +
	plot_layout(widths = c(1.3, 1))

# Save the figure
ggsave(
	filename = "/Users/serpent/Desktop/pathway_plot_compressed.png",
	plot = pathways_combined,
	width = 420,
	height = 290,
	units = "mm",
	dpi = 500
)
