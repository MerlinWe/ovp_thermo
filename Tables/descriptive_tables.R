library(dplyr)
library(gt)
library(purrr)

prep_data <- dat %>%
	filter(phase %in% c("day", "night")) %>%
	mutate(
		thi_range = case_when(
			heat_stress == "comfort"   ~ "< 68",
			heat_stress == "mild discomfort"  ~ "68–71.9",
			heat_stress == "discomfort"  ~ "72–74.9",
			heat_stress == "alert"     ~ "75–78.9",
			heat_stress == "danger"    ~ "79–83.9",
			heat_stress == "emergency" ~ "≥ 84",
			TRUE ~ NA_character_),
		heat_stress = factor(heat_stress,
												 levels = c("comfort", "mild discomfort", "discomfort", "alert", "danger", "emergency")),
		thi_range = factor(thi_range,
											 levels = c("< 68", "68–71.9", "72–74.9", "75–78.9", "79–83.9", "≥ 84")),
		phase = factor(phase, levels = c("day", "night"))
	) %>%
	filter(!is.na(heat_stress), !is.na(season), !is.na(mean_BT_smooth), !is.na(mean_activity_percent), !is.na(mean_heartrate))

# Unique season × phase combinations
season_phase_combos <- prep_data %>%
	distinct(season, phase)

# Function to summarise one subset
summarise_subset <- function(df, season_i, phase_i) {
	df_sub <- df %>%
		filter(season == season_i, phase == phase_i) %>%
		group_by(thi_range, heat_stress) %>%
		summarise(
			N = n(),
			BodyTemp_mean = mean(mean_BT_smooth, na.rm = TRUE),
			BodyTemp_sd   = sd(mean_BT_smooth, na.rm = TRUE),
			
			Activity_mean = mean(mean_activity_percent, na.rm = TRUE),
			Activity_sd   = sd(mean_activity_percent, na.rm = TRUE),
			
			HR_mean       = mean(mean_heartrate, na.rm = TRUE),
			HR_sd         = sd(mean_heartrate, na.rm = TRUE),
			.groups = "drop"
		) %>%
		mutate(
			BodyTemp = ifelse(N == 1,
												sprintf("%.2f", BodyTemp_mean),
												sprintf("%.2f ± %.2f", BodyTemp_mean, BodyTemp_sd)),
			Activity = ifelse(N == 1,
												sprintf("%.2f", Activity_mean),
												sprintf("%.2f ± %.2f", Activity_mean, Activity_sd)),
			HR = ifelse(N == 1,
									sprintf("%.2f", HR_mean),
									sprintf("%.2f ± %.2f", HR_mean, HR_sd)),
			Season = season_i,
			Phase = phase_i
		) %>%
		dplyr::select(Season, Phase, thi_range, heat_stress, N, BodyTemp, Activity, HR)
	
	return(df_sub)
}

# Apply over all combinations
combined_summary <- purrr::map2_dfr(season_phase_combos$season, season_phase_combos$phase,
																		~ summarise_subset(prep_data, .x, .y))

combined_summary %>%
	arrange(Season, thi_range, Phase) %>%
	gt(groupname_col = "Season") %>%
	tab_header(
		title = md("**Descriptive Summary by Season, Phase, and Heat Stress**")
	) %>%
	cols_label(
		thi_range = "THI Range",
		heat_stress = "Heat Stress",
		Phase = "Phase",
		N = "n obs",
		BodyTemp = "Body Temp (°C)",
		Activity = "Activity (%)",
		HR = "Heart Rate (bpm)"
	) %>%
	cols_align(align = "center", columns = everything()) %>%
	tab_style(
		style = cell_fill(color = "#f7f7f7"),
		locations = cells_body(
			rows = Phase == "night"
		)
	) %>%
	tab_options(
		row_group.as_column = TRUE,
		table.font.size = px(12),
		heading.title.font.size = px(14)
	)

library(flextable)
library(officer)

# Convert the long summary table to flextable
ft <- combined_summary %>%
	arrange(Season, thi_range, Phase) %>%
	flextable() %>%
	set_header_labels(
		thi_range = "THI Range",
		heat_stress = "Heat Stress",
		Phase = "Phase",
		N = "n obs",
		BodyTemp = "Body Temp (°C)",
		Activity = "Activity (%)",
		HR = "Heart Rate (bpm)"
	) %>%
	autofit() %>%
	add_header_lines("Descriptive Summary by Season, Phase, and Heat Stress") %>%
	theme_box()

# Export to Word
doc <- read_docx() %>%
	body_add_flextable(ft) %>%
	body_add_par("", style = "Normal")

print(doc, target = "heat_stress_summary.docx")




