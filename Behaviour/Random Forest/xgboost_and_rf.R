##################################################################################################################
#################################### ACC Classification (Random Forest) ####################################
##################################################################################################################

library(randomForest)
library(plyr)
library(tidyverse)
setwd("/Users/serpent/Documents/VHL/RAAK/Stats/Class (MWe)/Random Forest") 

rm(list=ls()) # clean global environment

# Step 1: Get data with labels, merge ACC with GPS, keep only labelled bursts with ACC data, build Random Forest, start wiht DC

# ========== ACC data (labelled) ==========
mydir <- "/Users/serpent/Documents/VHL/RAAK/Data/dairy campus/DC_001/ACC"
myfiles <- list.files(path=mydir, pattern="*.CSV", full.names=TRUE)
dc_01 <- ldply(myfiles, read_csv)
mydir <- "/Users/serpent/Documents/VHL/RAAK/Data/dairy campus/DC_002/ACC"
myfiles <- list.files(path=mydir, pattern="*.CSV", full.names=TRUE)
dc_02 <- ldply(myfiles, read_csv)

dc_01$timestamp <- as.POSIXct(dc_01$timestamp, origin = "1970-01-01")
dc_02$timestamp <- as.POSIXct(dc_02$timestamp, origin = "1970-01-01")
dc_acc <- bind_rows(dc_01, dc_02) %>% select(-prop_Position, -prop_Behaviour, timestamp, -resolution, -samplerate, -scale)
dc_acc <- dc_acc %>% select(-Position) %>% na.omit()
dc_acc <- dc_acc %>% 
	mutate(Behaviour = ifelse(grepl("f", Behaviour), "Foraging",
														ifelse(grepl("F", Behaviour), "Foraging",
																	 ifelse(grepl("Re", Behaviour), "Resting",
																	 			 ifelse(grepl("rre", Behaviour), "Resting",
																	 			 			 ifelse(grepl("Ru", Behaviour), "Ruminating", NA)))))) %>% 
	na.omit() %>%
	mutate_if(is.numeric, ~ .*9.81) # convert to m/s² (gravitational force)


rm(dc_01, dc_02, mydir, myfiles) 

library(randomForest)
library(rabc)

labels <- dc_acc$Behaviour 
dc_acc$timestamp <- NULL
dc_acc <- order_acc(df_raw = dc_acc)
df_time <- calculate_feature_time(df_raw = dc_acc, winlen_dba = 10)
df_freq <- calculate_feature_freq(df_raw = dc_acc, samp_freq = 10) 

acc_desc <- bind_cols(df_time, df_freq)
features <- select_features(df_feature = acc_desc, filter = TRUE, vec_label = labels, no_features = 10)
features <- acc_desc %>% select(any_of(features$features))

# XGBoost
mod1 <- train_model(df = acc_desc, vec_label = labels) 
mod2 <- train_model(df = features, vec_label = labels) 

# Random Forest
dat_rf <- acc_desc
dat_rf$Label <- as.factor(labels)

set.seed(123) 
train_indices <- sample(1:nrow(dat_rf), 0.7 * nrow(dat_rf)) # 70% for training
train_data <- dat_rf[train_indices, ]
test_data <- dat_rf[-train_indices, ]

mod3 <- randomForest(Label ~., data = train_data)
mod3

mean(predict(mod3, newdata = test_data) == test_data$Label) # check accuracy using test dataset 

tuned_mod3 <- tuneRF(train_data[, -29], train_data$Label, ntreeTry = 500, stepFactor = 1.5)
best_mtry <- tuned_mod3$mtry

kmeans_result <- kmeans(acc_desc, centers = 3)
cluster_assignments <- kmeans_result$cluster
print(cluster_assignments)





