################################################################################
########## Accelerometer Classification - Supervised machine learning ##########
################################################################################

# load packages (install first)
library(devtools) # GitHub Packages
# devtools::install_github("YuHuiDeakin/rabc", build_vigmette=TRUE)
library(rabc) # ACC Classification

# Prepare work space
setwd("C:/Users/weiss_002/OneDrive - Van Hall Larenstein/Minor Int. Wildlife Management/RAAK Pro/R/ACC Classification")
list.files() # check data availability

########## Prep & Preliminary analysis ##########

# Import & prepare labeled data
acc_dat <- read.csv("acc_labeled.csv")
labels <- acc_dat[ , 301] # get labels
acc_dat <- acc_dat[ ,1:300]/9.81 # Transform acc to gravitational force
acc_dat$label <- labels # reassign labels 

# Import & Prepare data to be labeled (OVP 001 from the 22.04.2020)
ovp_01 <- read.csv("OVP1_22042020.csv")[ ,-c(1:4)]
ovp_01 <- ovp_01/9.81 # transfrom acc to gravitational force
ovp_01$behaviour <- "unknown" 

# ---------- Analysis of labeled data  ----------

head(acc_dat[ , c(1:6, 301)], n = 5) # Examine labeled data

# Sort accelerometer data by labels and visualize it 
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

# --- Build & train model ---

# Train classification model with labeled data
model_output <- train_model(df = df_time[ ,selection$features], vec_label = labels) 

# Evaluate model performance using a confusion matrix 
predictions <- plot_confusion_matrix(df_feature = df_time[ ,selection$features], vec_label = labels)

# --- Apply model ---

# Calculate features of unlabeled data 
df_time_prediction <- calculate_feature_time(df_raw = ovp_01, winlen_dba = 10)

# Apply model to predict labels of unknown behaviours 
predicted_behaviours <- predict(model_output, df_time_prediction[ ,selection$features])

# ---- Add predicted labels to unpredicted data -----

ovp_01$behaviour <- predicted_behaviours

