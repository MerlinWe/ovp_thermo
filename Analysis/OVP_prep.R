###########################################################
########## Data preparation Script: RAAK Pro OVP ##########
###########################################################

# ---------- Getting started ----------
rm(list = ls()) # clean environment 
setwd("/Users/serpent/Documents/VHL/RAAK/Data")
options(digits = 3, digits.secs = 3, scipen = 999, datatable.quiet = TRUE)

# Install package libraries if needed 
library(data.table)
library(plyr)
library(reshape2)
library(gdata) 
library(lubridate)
library(oce)
library(suncalc)
library(cowplot)
library(purrr)
library(tidyverse)

## Metadata sensors: FIWI Measurements in 3 minute intervals
## KNMI weather station 269 (Lelystad) in hourly intervals split into 20 x 3 minutes

## --- Available FIWI data --- 

# OVP 001 has very little HR & TEMP data (77)
# OVP 002 has very little HR & TEMP data
# OVP 003 has NO HR data 
# OVP 004 has HR and TEMPT data (3334)
# OVP 005 (NEW) has HR and TEMP data (1157)
# OVP 006 not on the sharepoint? 
# OVP 007 (NEW) has NO HR data 
# OVP 008 has very little HR & TEMP data (136)
# OVP 009 has HR & TEMP data (1236)
# OVO 010 has HR and TEMP data (3358) 
# OVP 011 (NEW) has little HR & TEMP data (423)
# OVP 012 has NO HR data
# OVP 013 does not exist 
# OVP 014 (NEW) has no HR data 
# OVP 015 has HR & TEMP data (22467)
# OVP 016 (NEW) has very little HR & TEMP data (16)
# OVP 017 (NEW) has very little HR & TEMP data (5)
# OVP 018 (NEW) has TEMP & HR data (1034)
# OVP 019 (NEW) has TEMP & HR data (2460)
# OVP 020 (NEW) has TEMP & HR data (1010)

########## Importing & adjusting raw data ##########

# First we create a time sequence (21.11.18 - 24.03.21) in 3 minute intervals (i.e., the shortest common denominator)
# that we can use as a margin to then stepwise merge the sensor data. 

time_sq <- seq(from = ISOdatetime(2018,11,21,0,0,0), 
               to = ISOdatetime(2021,3,23,23,57,0), 
               by = (60*3)) %>%
  as_tibble() %>%
  dplyr::rename(Date = value) %>%
  dplyr::mutate(Date = as.POSIXct(Date))

# >>>>>>>>>> Import FIWI Data <<<<<<<<<<

# Make sure WD is set to a folder containing ONLY the FIWI data files. 
setwd("/Users/serpent/Documents/VHL/RAAK/Data/FIWI") 

# We will use a custom function to process the FIWI data: 
process_fiwi <- function(file, name) {
  # We read the file and set the date format
  read.csv(file, header = TRUE) %>%
    as_tibble() %>% 
    # Make sure ID is correct and set date format 
    mutate(ID = name, Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%OS")) %>% 
    # arrange chronologically 
    arrange(Date) %>% 
    # Exclude NAs in body temperature and heart rate 
    filter(!is.na(BodyTemp), !is.na(HeartRate)) %>%
    # Filter out body temperature values below 25 
    filter(BodyTemp > 25) %>%
    # De-spike body temperature using reference values 
    mutate(BT_smooth = despike(BodyTemp, reference = "median", 
                               n = 0.5, k = 61, replace = "reference")) %>%
    # Cut off Heart Rate values below 25 AND low values if values go up again
    filter(HeartRate > 25) %>%
    group_by(ID) %>%
    arrange(Date) %>%
    mutate(Above25Later = rev(cumsum(rev(HeartRate > 25)) > 0)) %>%
    ungroup() %>%
    filter(HeartRate >= 25 | Above25Later)
}

# * Note that we de-spike data using a smoothed reference BT from a 3h running mean 
# and replace spikes (= where the difference between BT and its reference > half 
# the standard deviation) with the corresponding reference values  

# We get the file names that we pass to the function from the WD and get a vector of the corresponding names 
fiwi_files <- list.files() 
names <- paste("OVP", gsub("[^0-9]", "", fiwi_files), sep = "_")

# Now we apply the processing function and set our names 
fiwi <- map2(fiwi_files, names, ~ process_fiwi(.x, .y)) %>%
  # Keep only cows with >20 observations
  purrr::keep(~ nrow(.x) >= 20)

# Retrieve those objects as a list
fiwi <- fiwi %>%
  # Create database 
  bind_rows() %>%
  # exclude nonsense values
  mutate(BT_smooth = despike(BT_smooth, reference = "trim", min = 35.5, max = 41),  
         # Double check date format  is correct 
         Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%OS"))


# Check if values make sense
plot_grid(
  
  # Heart Rate 
  ggplot(fiwi, aes(x=Date, y = HeartRate, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Heart Rate"),
  
  # Activity
  ggplot(fiwi, aes(x=Date, y = ActivityPercent, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Activity percent"), 
  
 # Raw body temperature 
  ggplot(fiwi, aes(x=Date, y = BodyTemp, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Body temperature (RAW)"),
 
  # Clean body temperature 
  ggplot(fiwi, aes(x=Date, y = BT_smooth, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
   ggtitle("Body temperature (CLEAN)"),
 
  nrow = 2, ncol = 2)

# >>>>>>>>>> Import KNMI 269 data <<<<<<<<<<
# Set WD to a folder containing the KNMI data 
setwd("/Users/serpent/Documents/VHL/RAAK/Data/KNMI") 

knmi_dat <- read.delim("KNMI_2011_2020.txt", header = TRUE, sep = ",") %>%
  # Remove data before collaring
  slice(-(1:69144)) %>%
  # add KNMI 2021
  bind_rows(read.delim("KNMI_2021.txt", header = TRUE, sep = ",") %>%
              slice(1:1968)) %>%
  # Set intervals 
  slice(rep(1:n(), each = 20)) %>%
  # Rename columns
  dplyr::select(FH, T, SQ, DR, RH, N, R, S) %>%
  dplyr::rename(wind_spd = FH, 
                knmi_temp = T, 
                sun_dur = SQ, 
                prcp_dur = DR, 
                prcp_mm = RH, 
                cloud_cov = N, 
                prcp = R, 
                snow = S) %>%
  # Set formats 
  mutate(
    Date = time_sq$Date,
    Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%OS"),
    wind_spd = ((wind_spd * 0.1) * 3600) / 1000,
    knmi_temp = knmi_temp * 0.1,
    sun_dur = sun_dur * 0.1,
    prcp_dur = prcp_dur * 0.1) %>%
  # Set column order
  select(wind_spd, knmi_temp, sun_dur, prcp_dur, prcp_mm, cloud_cov, prcp, snow) %>%
  # Set date format 
  mutate(Date = as.POSIXct(time_sq$Date, format = "%Y-%m-%d %H:%M:%OS"),
         Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%OS"))

########## Combine data & Create database ##########

# >>>>>>>>>> Aggregate and merge data <<<<<<<<<<

# Merge fiwi and KNMI data with our time sequence margin
ovp_dat <- fiwi %>%
  # Set format 
  as_tibble() %>%
  # Merge 
  full_join(knmi_dat, by = "Date") %>%
  full_join(time_sq, by = "Date") %>%
  # Define 3 minute intervals
  mutate(Date = cut(Date, breaks = "3 min"))

# Build database (merge data and add categorical variables)
ovp_dat <- Hmisc::Merge(
  # FIWI data 
  Hmisc::Merge(
    aggregate(BodyTemp ~ Date + ID, ovp_dat, mean, drop = FALSE),
    aggregate(BT_smooth ~ Date + ID, ovp_dat, mean, drop = FALSE),
    aggregate(CollarTemp ~ Date + ID, ovp_dat, mean, drop = FALSE),
    aggregate(ActivityPercent ~ Date + ID, ovp_dat, mean, drop = FALSE),
    aggregate(HeartRate ~ Date + ID, ovp_dat, mean, drop = FALSE),
    aggregate(DurHeadDown ~ Date + ID, ovp_dat, mean, drop = FALSE),
    aggregate(ChangesHeadPos ~ Date + ID, ovp_dat, mean, drop = FALSE),
    all = TRUE, id = ~ Date + ID, verbose = FALSE),
  # Weather data
  Hmisc::Merge(
    aggregate(wind_spd ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(knmi_temp ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(sun_dur ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(prcp_dur ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(prcp_mm ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(cloud_cov ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(prcp ~ Date, ovp_dat, mean, drop = FALSE),
    aggregate(snow ~ Date, ovp_dat, mean, drop = FALSE),
    all = TRUE, id = ~ Date, verbose = FALSE),
  all = TRUE, id = ~ Date, verbose = FALSE) %>% 
  filter(!is.na(BT_smooth), !is.na(HeartRate))

ovp_dat <- ovp_dat %>%
  # Separate date and time 
  mutate(sun_date = as.Date(format(as.POSIXct(Date, "%Y-%m-%d %H:%M:%OS", tz = "CET"), format = "%Y-%m-%d"))) %>%
  # Get sun data from Suncalc 
  group_by(Date) %>%
  mutate(sun = getSunlightTimes(date = sun_date,
                                lat = 52.4573, 
                                lon = 5.4194,
                                tz = "CET",
                                keep = c("sunrise", "sunset", "dusk", "dawn"))) %>%
  unnest_wider(sun) %>% 
  dplyr::select(-sun_date, -date, -lat, -lon) %>% 
  mutate( # Get other variables
    Date = format(as.POSIXct(Date, format='%Y-%m-%d %H:%M:%S')), 
    # Phase
    phase = case_when(
      Date >= dawn & Date < sunrise ~ "dawn",
      Date >= sunrise & Date < sunset ~ "day",
      Date >= sunset & Date < dusk ~ "dusk",
      TRUE ~ "night"),
    # Wind Chill Factor 
    WCF = 13.12+0.621*CollarTemp-11.37*(wind_spd^0.16)+0.3965*CollarTemp*(wind_spd^0.16),
    # Season 
    season = case_when(
      month(Date) %in% 9:11 ~ "Fall",
      month(Date) %in% c(12, 1, 2) ~ "Winter",
      month(Date) %in% 3:5 ~ "Spring",
      TRUE ~ "Summer"),
    # Weight 
    Weight = case_when(
      grepl("OVP_01", ID) ~ 543,
      grepl("OVP_02", ID) ~ 662,
      grepl("OVP_04", ID) ~ 538,
      grepl("OVP_05", ID) ~ 535,
      grepl("OVP_08", ID) ~ 570,
      grepl("OVP_09", ID) ~ 574,
      grepl("OVP_10", ID) ~ 571,
      grepl("OVP_11", ID) ~ 515,
      grepl("OVP_15", ID) ~ 452,
      grepl("OVP_16", ID) ~ 575,
      grepl("OVP_17", ID) ~ 540,
      grepl("OVP_18", ID) ~ 561,
      grepl("OVP_19", ID) ~ 539,
      grepl("OVP_20", ID) ~ 575,
      TRUE ~ NA_real_),
    # Energy expenditure
    EnergyExpHB = 2.907*(HeartRate^0.516)*(Weight^0.777) ,
    # Age
    Age = case_when(
      grepl("OVP_01", ID) ~ 11,
      grepl("OVP_02", ID) ~ 6,
      grepl("OVP_10", ID) ~ 12,
      grepl("OVP_11", ID) ~ 4,
      grepl("OVP_15", ID) ~ 15,
      grepl("OVP_16", ID) ~ 13,
      grepl("OVP_18", ID) ~ 14,
      grepl("OVP_20", ID) ~ 8,
      TRUE ~ NA_integer_),
    # Condition  
    Condition = case_when(
      grepl("OVP_01|OVP_02", ID) ~ "5",
      grepl("OVP_04|OVP_05|OVP_08|OVP_09|OVP_10|OVP_11|OVP_15|OVP_16|OVP_17|OVP_18|OVP_19|OVP_20", ID) ~ "4",
      TRUE ~ NA_character_)) %>% 
  # remove unnecessary data 
  dplyr::select(-dawn, -sunset, -sunrise, -dusk) %>% 
  filter(!is.na(Date))

# >>>>>>>>>> Calculate phase and previous phase temperatures  <<<<<<<<<<

# Here, we calculate the mean/max temperature of the day and night phases per ID and Date, 
# as well as the max and mean of the previous phase per ID and Date. 

temp_dat <- ovp_dat %>%
  # keep relevant columns 
  dplyr::select(Date, ID, phase, CollarTemp) %>% 
  # Get only the date in Y-M-D
  mutate(Date = as.Date(Date)) %>% 
  # Sort by day and ID
  arrange(Date, ID) %>%
  # Keep only day and night phases
  filter(phase %in% c("day", "night")) %>%
  # Calculate mean and max
  group_by(ID, Date, phase) %>%
  summarize(
    phase_mean = mean(CollarTemp, na.rm = TRUE),
    phase_max = max(CollarTemp, na.rm = TRUE)) %>%
  # Calculate prev mean and max 
  ungroup() %>%
  arrange(ID, Date, phase) %>%
  group_by(ID) %>%
  mutate(
    prev_phase_mean = lag(phase_mean, default = NA),
    prev_phase_max = lag(phase_max, default = NA),
    prev_phase_mean = if_else(Date - lag(Date) <= 1 & 
                                phase != lag(phase), 
                              prev_phase_mean, NA_real_),
    prev_phase_max = if_else(Date - lag(Date) <= 1 & 
                               phase != lag(phase), 
                             prev_phase_max, NA_real_)) %>%
  ungroup()

# >>>>>>>>>> Arrange by block variable  <<<<<<<<<<

# Calculate mean per phase per id per day and season 
dat <- ovp_dat %>%
  # Reset all grouping parameters
  ungroup() %>%
  # Get days only and add temp_dat data 
  mutate(Date = as.Date(Date)) %>%
  left_join(temp_dat, by = c("Date", "phase", "ID")) %>%
  # Get cumulative day variable 
  group_by(ID) %>%
  mutate(cum_day = match(Date, seq(min(Date), max(Date), by = "days"))) %>%
  ungroup() %>%
  # group by ID, Date, Phase, Season 
  group_by(ID, Date, phase, season) %>%
  # Summarize numeric values by mean 
  summarize_if(is.numeric, mean, na.rm = TRUE) %>%
  # Get block variable
  dplyr::group_by(ID, phase) %>% 
  dplyr::mutate(block = row_number()) %>% 
  dplyr::select(ID, Date, phase, season, cum_day, block, BodyTemp, BT_smooth, 
                phase_mean, phase_max, prev_phase_mean, prev_phase_max,
                ActivityPercent, HeartRate, DurHeadDown, ChangesHeadPos, 
                EnergyExpHB, Weight, WCF) %>%
  distinct(.keep_all = TRUE) %>%
  # Replace NaN with NA
  mutate_all(~replace(., is.nan(.), NA))

colnames(dat) <- c("ID", "date", "phase", "season", "cum_day", "block", "mean_BT_raw", "mean_BT_smooth", 
                   "phase_mean_CT", "phase_max_CT", "prev_phase_mean_CT", "prev_phase_max_CT", 
                   "mean_activity_percent", "mean_heartrate", "mean_head_down", "mean_head_change", "mean_energy_exp",
                   "weight", "mean_wcf")

# Check if values make sense 
plot_grid(
  
  # Heart Rate 
  ggplot(dat, aes(x=date, y = mean_heartrate, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mean heart rate"),
  
  # Activity
  ggplot(dat, aes(x=date, y = mean_activity_percent, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mean activity percent"), 
  
  # Raw body temperature 
  ggplot(dat, aes(x=date, y = mean_BT_raw, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mean body temperature (RAW)"),
  
  # Clean body temperature 
  ggplot(dat, aes(x=date, y = mean_BT_smooth, group = ID, fill = ID)) + 
    geom_point(size = .2) +
    facet_wrap(~ID) +
    xlab(NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mean body temperature (CLEAN)"),
  
  nrow = 2, ncol = 2)

ggsave(filename = "ovp_values.png",
       path = "/Users/serpent/Desktop",
       plot = last_plot(),
       bg = "white",
       width = 300,
       height = 200,
       units = "mm",
       dpi = 1257)

# Write file 
write_csv(dat, paste0("/Users/serpent/Documents/VHL/RAAK/Data/", "ovp_data_", format(Sys.time(), "%d_%m_%y"), ".csv"))

