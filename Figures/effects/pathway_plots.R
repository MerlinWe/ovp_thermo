rm(list = ls()); gc()

library(tidyverse)
library(patchwork)

# === Load & prepare data ===
pathways <- read_csv("/Users/merlin/Documents/VHL/OVP/Code/Figures/effects/pathway_effects.csv")
thi <- read_csv("/Users/merlin/Documents/VHL/OVP/Code/Figures/effects/combined_THI.csv")

sig <- read_csv("/Users/merlin/Documents/VHL/OVP/Code/Figures/effects/sig_lines_ovp.csv") %>% 
  mutate(
    season = str_to_title(season),
    response = case_when(
      y == "heartrate" ~ "Heart Rate",
      y == "bt" ~ "Body Temperature",
      y == "activity" ~ "Activity",
      TRUE ~ y
    ),
    predictor = case_when(
      x == "thi" ~ "THI",
      x == "activity" ~ "Activity",
      TRUE ~ x
    ),
    timeofday = str_to_title(phase)
  ) %>%
  dplyr:: select(-x, -y, -phase)

pathways <- pathways %>%
  left_join(sig %>%
              filter(predictor == "Activity") %>%
              dplyr::select(response, season, timeofday, sig) %>%
              rename(sig_flag = sig),
            by = c("response", "season", "timeofday")) %>%
  mutate(sig_flag = ifelse(is.na(sig_flag), FALSE, sig_flag))

thi <- thi %>%
  left_join(sig %>%
              filter(predictor == "THI") %>%
              dplyr::select(response, season, timeofday, sig) %>%
              rename(sig_flag = sig),
            by = c("focus_variable" = "response", "season", "timeofday")) %>%
  mutate(sig_flag = ifelse(is.na(sig_flag), FALSE, sig_flag))

table(pathways$sig_flag, pathways$response)
table(thi$sig_flag, thi$focus_variable)


# Add linetype
pathways <- pathways %>%
	mutate(linetype = ifelse(timeofday == "Day", "solid", "dotted")) %>%
	rename(season_name = season)

thi <- thi %>%
	mutate(linetype = ifelse(timeofday == "Day", "solid", "dotted")) %>%
	rename(season_name = season)

# === Plot styles ===
season_colors <- c(
	"Winter" = "grey25", 
	"Spring" = "springgreen3",
	"Summer" = "violetred3", 
	"Autumn" = "tan3")

color_scale <- scale_color_manual(values = season_colors)
linetype_scale <- scale_linetype_manual(values = c("Day" = "solid", "Night" = "11"))

# THI thresholds & stress labels Mid-level classification
thi_thresholds <- c(72, 79)
stress_labels <- c("I", "II", "III")

# Position stress labels centered between breaks
stress_label_x <- c(69, 75.5, 85)

# Create subdued vertical lines + labels
thi_lines <- list(
	geom_vline(xintercept = thi_thresholds, linetype = "dashed", color = "grey60", size = 0.4),
	annotate("text", x = stress_label_x, y = -Inf,
					 label = stress_labels, vjust = -0.5, size = 3, color = "grey40")
)

# === Plot functions ===
plot_pathway <- function(data, resp, xlab, ylab, title) {
  ggplot(filter(data, response == resp),
         aes(x = x, y = fit, color = season_name, linetype = timeofday, linewidth = sig_flag)) +
    geom_line() +
    color_scale + linetype_scale +
    scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.4)) +
    labs(x = xlab, y = ylab, title = title, color = NULL, linetype = NULL, linewidth = NULL) +
    theme_bw(base_line_size = 0) +
    theme(text = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 10))
}

plot_thi <- function(data, focus_var, ylab, title) {
  ggplot(filter(data, focus_variable == focus_var),
         aes(x = phase_mean_THI, y = fit, color = season_name, linetype = timeofday, linewidth = sig_flag)) +
    geom_line() +
    color_scale + linetype_scale +
    scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.4)) +
    thi_lines +
    labs(x = "Temperature-Humidity-Index", y = ylab, title = title) +
    theme_bw(base_line_size = 0) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          text = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 10))
}

# === Individual plots ===
plot_act_to_hr <- plot_pathway(pathways, "Heart Rate", "Non-Recumbent (%, Activity)", "Heart Rate (bpm)", "")
plot_act_to_bt <- plot_pathway(pathways, "Body Temperature", "Non-Recumbent (%, Activity)", "Body Temperature (°C)", "")
plot_thi_hr    <- plot_thi(thi, "Heart Rate", "Heart Rate (bpm)", "")
plot_thi_act   <- plot_thi(thi, "Activity", "Non-Recumbent (%)", "")
plot_thi_bt    <- plot_thi(thi, "Body Temperature", "Body Temperature (°C)", "")

# === Combine plots ===
plot_direct <- (plot_act_to_hr / plot_act_to_bt) + 
	plot_layout(guides = "collect") &
	guides(color = guide_legend(nrow = 1, byrow = TRUE), linetype = guide_legend(nrow = 1)) &
	theme(legend.position = "bottom")

plot_thi_combined <- (plot_thi_act / plot_thi_bt / plot_thi_hr) +
	plot_layout(guides = "collect", heights = c(1, 1, 1)) &
	theme(legend.position = "none")

# Final layout
pathways_combined <- (plot_thi_combined | plot_direct) +
	plot_layout(widths = c(1, 1))

ggsave(
	filename = "/Users/merlin/Desktop/pathway_plot_compressed.png",
	plot = pathways_combined,
	width = 240,
	height = 200,
	units = "mm",
	dpi = 500)

# Simplified stress levels
thi_legend_df <- tibble(
	x = c(68, 75.5, 84),
	label = c("I: comfort–mild discomfort", 
						"II: discomfort–alert", 
						"III: danger–emergency")
)

thi_legend_plot <- ggplot(thi_legend_df, aes(x = x, y = 1, label = label)) +
	geom_text(size = 2.8, hjust = 0.5, color = "grey40", fontface = "italic") +
	scale_x_continuous(limits = c(60, 90), expand = c(0, 0)) +
	theme_void()

ggsave(
	filename = "/Users/merlin/Desktop/thi_legend.png",
	plot = thi_legend_plot,
	width = 240,
	height = 200,
	units = "mm",
	dpi = 500)

ggsave(
	filename = "/Users/merlin/Desktop/pathway_plot_legend.png",
	plot = pathways_combined,
	width = 440,
	height = 200,
	units = "mm",
	dpi = 500)



