library(rabc)
library(plyr)
library(tidyverse)

# prep only DC data (labelled ACC bursts)

mydir <- "/Users/serpent/Documents/VHL/RAAK/Data/dairy campus/DC_001/ACC"
myfiles <- list.files(path=mydir, pattern="*.CSV", full.names=TRUE)
dc_01 <- ldply(myfiles, read_csv)

mydir <- "/Users/serpent/Documents/VHL/RAAK/Data/dairy campus/DC_002/ACC"
myfiles <- list.files(path=mydir, pattern="*.CSV", full.names=TRUE)
dc_02 <- ldply(myfiles, read_csv)

dc <- bind_rows(dc_01, dc_02) %>% select(-prop_Position, -prop_Behaviour, -timestamp, -resolution, -samplerate, -scale)
dc_behaviour <- dc %>% select(-Position) %>% na.omit()

rm(dc, dc_01, dc_02, mydir, myfiles)

dc <- dc_behaviour %>% 
	mutate(Behaviour = ifelse(grepl("f", Behaviour), "Foraging",
														ifelse(grepl("F", Behaviour), "Foraging",
																	 ifelse(grepl("Re", Behaviour), "Resting",
																	 			 	ifelse(grepl("rre", Behaviour), "Resting",
																	 			 				 ifelse(grepl("Ru", Behaviour), "Ruminating", NA)))))) %>%
	na.omit() %>%
	mutate_if(is.numeric, ~ .*9.81) # convert to m/s² (gravitational force)

table(dc_behaviour$Behaviour)
labels <- dc$Behaviour 

# Sort accelerometer data by labels and visualize it 
acc_sorted <- order_acc(df_raw = dc)
plot_acc(df_raw = acc_sorted, axis_num = 3)

# Calculation of time domain features
df_time <- calculate_feature_time(df_raw = acc_sorted, winlen_dba = 10) # 10 second running window 
head(df_time)

# Calculation of frequency domain features 
df_freq <- calculate_feature_freq(df_raw = acc_sorted, samp_freq = 10) # 10 Hz frequency 
head(df_freq)

# Select features 
selection <- select_features(df_feature = cbind(df_time, df_freq), filter = TRUE, vec_label = labels, no_features = 5)

# Train classification model
model_output <- train_model(df = cbind(df_time, df_freq)[ ,selection$features], vec_label = labels) 

# Evaluate model performance using a confusion matrix 
predictions <- plot_confusion_matrix(df_feature = cbind(df_time, df_freq)[ ,selection$features], vec_label = labels)

