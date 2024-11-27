######################################################
########## Accelerometer Classification OVP ##########
######################################################

library(tidyverse) # Manipulation & Graphics 
library(devtools) # GitHub Packages
library(lattice) # Supplementary

# Get rabc package
devtools::install_github("YuHuiDeakin/rabc", build_vigmette=TRUE)
library(rabc)

# Prepare work space
setwd("C:/Users/weiss_002/OneDrive - Van Hall Larenstein/Minor Int. Wildlife Management/RAAK Pro/R/Data/Behaviour")
list.files() # check data availability

# Import ACC data 

OVP1_10.3.20 <- read.csv("OVP_001_10.03.20.csv", na.strings = "")[,-c(1:4)]
OVP1_17.3.20 <- read.csv("OVP_001_17.03.20.csv", na.strings = "")[,-c(1:4)]
OVP1_20.3.20 <- read.csv("OVP_001_20.03.20.csv", na.strings = "")[,-c(1:4)]

OVP2_19.3.20 <- read.csv("OVP_002_19.03.20.csv", na.strings = "")[,-c(1:4)]

OVP4_10.3.20 <- read.csv("OVP_004_10.03.20.csv", na.strings = "")[,-c(1:4)]

OVP7_17.3.20 <- read.csv("OVP_007_17.03.20.csv", na.strings = "")[,-c(1:4)]
OVP7_20.3.20 <- read.csv("OVP_007_20.03.20.csv", na.strings = "")[,-c(1:4)]

# Combine data 
acc_all <- data.frame(rbind(OVP1_10.3.20, OVP1_17.3.20, OVP1_20.3.20, OVP2_19.3.20,
                            OVP4_10.3.20,OVP7_17.3.20 ,OVP7_20.3.20))


acc_dat <- na.omit(acc_all) # Get labeled data 
OVP7_17.3.20  <- OVP7_17.3.20[ ,1:300]/9.81 # transform acc to gravitational force 
OVP7_17.3.20$behaviour <- "unknown"

# Concentrate Act_Stand and Pass_Stand to "Resting"
acc_dat <- mutate(acc_dat, behaviour = ifelse(grepl("Act_Stand", behaviour), "Rest",
                                              ifelse(grepl("Pass_Stand", behaviour), "Rest",
                                                     ifelse(grepl("Foraging", behaviour), "Foraging", 
                                                            ifelse(grepl("Walking", behaviour), "Walking", "unknown")))))

# Clean Gl.Env.
rm(OVP1_10.3.20, OVP1_17.3.20, OVP1_20.3.20, OVP2_19.3.20,
   OVP4_10.3.20 ,OVP7_20.3.20, acc_all)

########## Data preparation ##########

# ---------- Analysis of labeled data  ----------

head(acc_dat[ , c(1:6, 301)], n = 5) # Examine data
labels <- acc_dat[ , 301] # get labels
acc_dat <- acc_dat[ ,1:300]/9.81 # Transform acc to gravitational force
acc_dat$label <- labels # reassign labels 

# Sort accelerometer data and visualize it 
acc_sorted <- order_acc(df_raw = acc_dat)
plot_acc(df_raw = acc_sorted, axis_num = 3)

# ----- Calculate descriptive features of labeled data ----- 

# Calculation of time domain features
df_time <- calculate_feature_time(df_raw = acc_sorted, winlen_dba = 10) # 10 second running window 

# Calculation of frequency domain features 
df_freq <- calculate_feature_freq(df_raw = acc_sorted, samp_freq = 10) # 10 Hz frequency 


# ----- Select best fitting descriptive features -----

selection <- select_features(df_feature = cbind(df_time, df_freq), # without correlation filter
                             filter = FALSE, cutoff = 0.9, vec_label = labels, no_features = 10)

plot_selection_accuracy(results = selection) # Examine accuracy of features 

# ----- Plot features -----

plot_feature(df_feature = df_time[, "ODBA", drop = FALSE], vec_label = labels) # dynamic (ODBA)
plot_grouped_feature(df_feature = df_time[, "ODBA", drop = FALSE], vec_label = labels, geom = "boxplot") # static (ODBA)


########## Classification model ##########

# Train classification model with labeled data
model_output <- train_model(df = df_time, vec_label = labels) 

# Calculate features of unlabeled data 
df_time_prediction <- calculate_feature_time(df_raw = OVP7_17.3.20, winlen_dba = 10)

# Apply model to predict labels of unknown behaviors 
predicted_behaviours <- predict(model_output, df_time_prediction)

# Evaluate model performance using a confusion matrix 
predictions <- plot_confusion_matrix(df_feature = df_time, vec_label = labels)

# Plot wrong classifications 
plot_wrong_classifications(df_raw = acc_sorted, df_result = predictions)

########## Finally: View predictions ##########

OVP7_17.3.20$behaviour <- predicted_behaviours
