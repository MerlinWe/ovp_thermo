### OVP 015 activity ~ ACC behaviour ###

#setwd("C:/Users/weiss_002/OneDrive - Van Hall Larenstein/RAAK/Analysis/OVP015")
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Analysis/OVP015")
options(digits = 3, digits.secs = 1, scipen = 999)

library(tidyverse) # Graphics & Manipulation
library(plyr) # Data manipulation
library(reshape2) # Data manipulation
library(data.table) # Data manipulation
library(gdata) # Data management 
library(Hmisc) # Data management n
library(devtools) # Accessing GitHub packages

#devtools::install_github("YuHuiDeakin/rabc", build_vigmette=TRUE)
library(rabc) # ACC Classification 

list.files() # check file availability

## read GPS & FIWI data of OVP 15 
dat <- subset(read.csv("ovp_gps_movement_activity_01.06.23.csv"), id == "OVP_15")

## ----- Prepare availabel ACC data -----

## Read & prep ACC data (from the solar powered ACC on top!)
mydir <- "ACC" # Import unlabeled data 
myfiles <- list.files(path=mydir, pattern="*.CSV", full.names=TRUE)
acc_15 <- ldply(myfiles, read_csv)[1:150000, ] # keep only a number of rows to prevent R from crashing 

## Check ACC meta data 
meta <- acc_15[ ,c(1:4)] # get acc meta data
meta$Datetime <- as.POSIXct(meta$timestamp, origin = "1970-01-01") # convert time stamp 
meta$timestamp <- NULL # remove times tamp
meta$ID <- "OVP 15" # define animal ID

# Check tracking consistency
ggplot(meta) + # ... in time 
  geom_point(aes(x = Datetime, y = ID))

## Check time lag in minutes
meta %>%
  mutate(time.lag = as.numeric(difftime(Datetime, lag(Datetime), units = "mins"))) %>%
  ggplot() +
  geom_density(aes(x = time.lag, group = ID), fill = "black", alpha = 0.25) +
  scale_x_log10() 

acc_15 <- acc_15[ ,-c(1:4)]/9.81  # transform acceleration to gravitational force 
acc_15$behaviour <- "unknown"     # tell R behaviour is unknown 

## ----- Analyse labeled data -----

## Get labaled data (from all animals, based on the available video material)
acc_labeled <- read.csv("acc_labeled.csv") # Import labeled data to train classification model 
labels <- acc_labeled[ , 301] # get behaviour labels
acc_labeled <- acc_labeled[ ,1:300]/9.81 # Transform acc to gravitational force
acc_labeled$label <- labels # reassign labels 

## Sort accelerometer data by labels and visualize it 
acc_sorted <- order_acc(df_raw = acc_labeled)
plot_acc(df_raw = acc_sorted, axis_num = 3)

## Calculation of time domain features
df_time <- calculate_feature_time(df_raw = acc_sorted, winlen_dba = 10) # 10 second running window 

## Calculation of frequency domain features 
df_freq <- calculate_feature_freq(df_raw = acc_sorted, samp_freq = 10) # 10 Hz frequency 

## Select best fitting descriptive features
selection <- select_features(df_feature = cbind(df_time, df_freq), # without correlation filter
                             filter = FALSE, cutoff = 0.9, vec_label = labels, no_features = 10)

plot_selection_accuracy(results = selection) # Examine accuracy of features 


## ----- Build XGBoost supervised classification model -----

## Train classification model with labeled data
classification_model <- train_model(df = df_time[ ,selection$features], vec_label = labels) 

## Evaluate model performance using a confusion matrix 
predictions <- plot_confusion_matrix(df_feature = df_time[ ,selection$features], vec_label = labels)

## --- Apply model ---

## Calculate features of unlabelled data 
acc_descriptives <- calculate_feature_time(df_raw = acc_15, winlen_dba = 10)

## Apply model to predict labels of unknown behaviours 
predicted_behaviours <- predict(classification_model, acc_descriptives[ , selection$features])

## ----- Add predicted labels to unlabeled data -----

acc_15$behaviour <- predicted_behaviours # add labels
acc_15 <- data.frame(cbind(meta, acc_15, acc_descriptives)) # merge data frames 
acc_15 <- dplyr::select(acc_15, Datetime, behaviour, ODBA) # keep relevant columns

## ----- Calculate time spent per behaviour 

# Set intervals
acc_15$Date_Int <- cut(acc_15$Datetime, breaks = "12 min") # cut in proper intervals
acc_15 <- acc_15[ ,c("behaviour", "Date_Int", "ODBA")]

acc_15_int <- t(apply(xtabs(~ Date_Int + behaviour, data = acc_15), 1, function(z) z/sum(z)))
acc_15_int <- na.omit(cbind(as.data.frame(acc_15_int), Time = rownames(acc_15_int)))
rownames(acc_15_int) <- NULL
colnames(acc_15_int) <- c("%_Foraging", "%_Resting", "%_Walking", "timestamp")

## ----- Merge with fiwi data -----

dat <- Merge(dat, acc_15_int, all = TRUE, id = ~ timestamp)
dat <- na.omit(dat) # clean data 

## ----- Check relationships between activity and behaviours -----

# activity ~ resting
cor(dat$activity, dat$`%_Resting`) 
ggplot(dat, aes(x=activity, y = `%_Resting`)) +
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("activity ~ resting")

# activity ~ walking
cor(dat$activity, dat$`%_Walking`) 
ggplot(dat, aes(x=activity, y = `%_Walking`)) + 
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("activity ~ walking")

# activity ~ foraging
cor(dat$activity, dat$`%_Foraging`) 
ggplot(dat, aes(x=activity, y = `%_Foraging`)) + 
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("activity ~ foraging")

## ----- Check relationships between head positioning and behaviours -----

# head down ~ foraging
cor(dat$head_down, dat$`%_Foraging`) 
ggplot(dat, aes(x=head_down, y = `%_Foraging`)) +
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("activity ~ foraging")



