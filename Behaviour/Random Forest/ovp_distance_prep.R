########## RAAK Proj: OVP data prep GPS ##########

setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS") # set wd to location including raw GPS data

library(tidyverse)      
library(data.table)
library(plyr)
library(reshape2)
library(Hmisc)
library(gdata)
library(lubridate)
library(geosphere)

options(scipem = 999) # disable scientific notation 

# Do GPS prep for OVP 1,2,8,9,10,15,19,20 

## ---------- GPS 01 ----------
mydir <- "GPS_01" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_01 <- ldply(myfiles, read_csv) # read data 
gps_01$id <- "OVP_01" # add ID column 
gps_01 <- dplyr::select(gps_01, id, gpstime, latitude, longitude)
gps_01$gpstime <- as.POSIXct(gps_01$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_01$latitude <- gps_01$latitude / 10^7
gps_01$longitude <- gps_01$longitude / 10^7

# Check time lag 
gps_01 <- gps_01 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_01[, c("longitude", "latitude")][-nrow(gps_01), ],
  p2 = gps_01[, c("longitude", "latitude")][-1, ])

gps_01 <- data.frame(gps_01[-1, ], distance)
gps_01 <- subset(gps_01, distance < 1000)
colnames(gps_01) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_01 <- read.csv("01_FIWI.csv", header = TRUE) # read data
ovp_01$Id <- "OVP_01" # make sure ID is correct
ovp_01$Date <- as.POSIXct(ovp_01$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_01 <- arrange(ovp_01, Date) # arrange chronologically 
colnames(ovp_01) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_01 <- Merge(ovp_01, gps_01, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

summary(ovp_01$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,15,0,0,0), ISOdatetime(2023,3,13,23,57,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_01 <- Merge(time_sq, ovp_01, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_01$timestamp <- cut(ovp_01$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_01, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_01, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_01, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_01, mean, drop = FALSE) # aggregate head changes 

ovp_01 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_01 <- na.omit(ovp_01)
rownames(ovp_01) <- NULL
ovp_01$id <- "OVP_01"
rm(gps_01) # remove redundant data 

## ---------- GPS 02 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_02" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_02 <- ldply(myfiles, read_csv) # read data 
gps_02$id <- "OVP_02" # add ID column 
gps_02 <- dplyr::select(gps_02, id, gpstime, latitude, longitude)
gps_02$gpstime <- as.POSIXct(gps_02$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_02$latitude <- gps_02$latitude / 10^7
gps_02$longitude <- gps_02$longitude / 10^7

# Check time lag 
gps_02 <- gps_02 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_02[, c("longitude", "latitude")][-nrow(gps_02), ],
  p2 = gps_02[, c("longitude", "latitude")][-1, ])

gps_02 <- data.frame(gps_02[-1, ], distance)
gps_02 <- subset(gps_02, distance < 1000)
colnames(gps_02) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_02 <- read.csv("02_FIWI.csv", header = TRUE) # read data
ovp_02$Id <- "OVP_02" # make sure ID is correct
ovp_02$Date <- as.POSIXct(ovp_02$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_02 <- arrange(ovp_02, Date) # arrange chronologically 
colnames(ovp_02) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_02 <- Merge(ovp_02, gps_02, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

summary(ovp_02$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,24,12,0,0), ISOdatetime(2021,5,5,13,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_02 <- Merge(time_sq, ovp_02, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_02$timestamp <- cut(ovp_02$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_02, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_02, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_02, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_02, mean, drop = FALSE) # aggregate head changes 

ovp_02 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_02 <- na.omit(ovp_02)
rownames(ovp_02) <- NULL
ovp_02$id <- "OVP_02"
rm(gps_02) # remove redundant data 


## ---------- GPS 08 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_08" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_08 <- ldply(myfiles, read_csv) # read data 
gps_08$id <- "OVP_08" # add ID column 
gps_08 <- dplyr::select(gps_08, id, gpstime, latitude, longitude)
gps_08$gpstime <- as.POSIXct(gps_08$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_08$latitude <- gps_08$latitude / 10^7
gps_08$longitude <- gps_08$longitude / 10^7

# Check time lag 
gps_08 <- gps_08 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_08[, c("longitude", "latitude")][-nrow(gps_08), ],
  p2 = gps_08[, c("longitude", "latitude")][-1, ])

gps_08 <- data.frame(gps_08[-1, ], distance)
gps_08 <- subset(gps_08, distance < 1000)
colnames(gps_08) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_08 <- read.csv("08_FIWI.csv", header = TRUE) # read data
ovp_08$Id <- "OVP_08" # make sure ID is correct
ovp_08$Date <- as.POSIXct(ovp_08$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_08 <- arrange(ovp_08, Date) # arrange chronologically 
colnames(ovp_08) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_08 <- Merge(ovp_08, gps_08, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

summary(ovp_08$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,15,10,0,0), ISOdatetime(2021,3,25,18,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_08 <- Merge(time_sq, ovp_08, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_08$timestamp <- cut(ovp_08$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_08, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_08, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_08, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_08, mean, drop = FALSE) # aggregate head changes 

ovp_08 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_08 <- na.omit(ovp_08)
rownames(ovp_08) <- NULL
ovp_08$id <- "OVP_08"
rm(gps_08) # remove redundant data 

## ---------- GPS 09 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_09" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_09 <- ldply(myfiles, read_csv) # read data 
gps_09$id <- "OVP_09" # add ID column 
gps_09 <- dplyr::select(gps_09, id, gpstime, latitude, longitude)
gps_09$gpstime <- as.POSIXct(gps_09$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_09$latitude <- gps_09$latitude / 10^7
gps_09$longitude <- gps_09$longitude / 10^7

# Check time lag 
gps_09 <- gps_09 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_09[, c("longitude", "latitude")][-nrow(gps_09), ],
  p2 = gps_09[, c("longitude", "latitude")][-1, ])

gps_09 <- data.frame(gps_09[-1, ], distance)
gps_09 <- subset(gps_09, distance < 1000)
colnames(gps_09) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_09 <- read.csv("09_FIWI.csv", header = TRUE) # read data
ovp_09$Id <- "OVP_09" # make sure ID is correct
ovp_09$Date <- as.POSIXct(ovp_09$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_09 <- arrange(ovp_09, Date) # arrange chronologically 
colnames(ovp_09) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_09 <- Merge(ovp_09, gps_09, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

ovp_09 <- ovp_09[-c(1:179), ] # remove some erroneous timestamp 
summary(ovp_09$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,11,15,10,0,0), ISOdatetime(2021,3,28,11,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_09 <- Merge(time_sq, ovp_09, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_09$timestamp <- cut(ovp_09$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_09, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_09, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_09, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_09, mean, drop = FALSE) # aggregate head changes 

ovp_09 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_09 <- na.omit(ovp_09)
rownames(ovp_09) <- NULL
ovp_09$id <- "OVP_09"
rm(gps_09) # remove redundant data 

## ---------- GPS 10 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_10" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_10 <- ldply(myfiles, read_csv) # read data 
gps_10$id <- "OVP_10" # add ID column 
gps_10 <- dplyr::select(gps_10, id, gpstime, latitude, longitude)
gps_10$gpstime <- as.POSIXct(gps_10$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_10$latitude <- gps_10$latitude / 10^7
gps_10$longitude <- gps_10$longitude / 10^7

# Check time lag 
gps_10 <- gps_10 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_10[, c("longitude", "latitude")][-nrow(gps_10), ],
  p2 = gps_10[, c("longitude", "latitude")][-1, ])

gps_10 <- data.frame(gps_10[-1, ], distance)
gps_10 <- subset(gps_10, distance < 1000)
colnames(gps_10) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_10 <- read.csv("10_FIWI.csv", header = TRUE) # read data
ovp_10$Id <- "OVP_10" # make sure ID is correct
ovp_10$Date <- as.POSIXct(ovp_10$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_10 <- arrange(ovp_10, Date) # arrange chronologically 
colnames(ovp_10) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_10 <- Merge(ovp_10, gps_10, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

summary(ovp_10$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,24,12,0,0), ISOdatetime(2023,3,13,11,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_10 <- Merge(time_sq, ovp_10, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_10$timestamp <- cut(ovp_10$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_10, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_10, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_10, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_10, mean, drop = FALSE) # aggregate head changes 

ovp_10 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_10 <- na.omit(ovp_10)
rownames(ovp_10) <- NULL
ovp_10$id <- "OVP_10"
rm(gps_10) # remove redundant data 

## ---------- GPS 15 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_15" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_15 <- ldply(myfiles, read_csv) # read data 
gps_15$id <- "OVP_15" # add ID column 
gps_15 <- dplyr::select(gps_15, id, gpstime, latitude, longitude)
gps_15$gpstime <- as.POSIXct(gps_15$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_15$latitude <- gps_15$latitude / 10^7
gps_15$longitude <- gps_15$longitude / 10^7

# Check time lag 
gps_15 <- gps_15 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_15[, c("longitude", "latitude")][-nrow(gps_15), ],
  p2 = gps_15[, c("longitude", "latitude")][-1, ])

gps_15 <- data.frame(gps_15[-1, ], distance)
gps_15 <- subset(gps_15, distance < 1000)
colnames(gps_15) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_15 <- read.csv("15_FIWI.csv", header = TRUE) # read data
ovp_15$Id <- "OVP_15" # make sure ID is correct
ovp_15$Date <- as.POSIXct(ovp_15$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_15 <- arrange(ovp_15, Date) # arrange chronologically 
colnames(ovp_15) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_15 <- Merge(ovp_15, gps_15, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

summary(ovp_15$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,24,13,0,0), ISOdatetime(2021,5,05,14,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_15 <- Merge(time_sq, ovp_15, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_15$timestamp <- cut(ovp_15$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_15, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_15, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_15, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_15, mean, drop = FALSE) # aggregate head changes 

ovp_15 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_15 <- na.omit(ovp_15)
rownames(ovp_15) <- NULL
ovp_15$id <- "OVP_15"
rm(gps_15) # remove redundant data 

## ---------- GPS 19 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_19" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_19 <- ldply(myfiles, read_csv) # read data 
gps_19$id <- "OVP_19" # add ID column 
gps_19 <- dplyr::select(gps_19, id, gpstime, latitude, longitude)
gps_19$gpstime <- as.POSIXct(gps_19$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_19$latitude <- gps_19$latitude / 10^7
gps_19$longitude <- gps_19$longitude / 10^7

# Check time lag 
gps_19 <- gps_19 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_19[, c("longitude", "latitude")][-nrow(gps_19), ],
  p2 = gps_19[, c("longitude", "latitude")][-1, ])

gps_19 <- data.frame(gps_19[-1, ], distance)
gps_19 <- subset(gps_19, distance < 1000)
colnames(gps_19) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_19 <- read.csv("19_FIWI.csv", header = TRUE) # read data
ovp_19$Id <- "OVP_19" # make sure ID is correct
ovp_19$Date <- as.POSIXct(ovp_19$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_19 <- arrange(ovp_19, Date) # arrange chronologically 
colnames(ovp_19) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_19 <- Merge(ovp_19, gps_19, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

summary(ovp_19$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,24,12,0,0), ISOdatetime(2023,3,13,11,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_19 <- Merge(time_sq, ovp_19, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_19$timestamp <- cut(ovp_19$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_19, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_19, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_19, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_19, mean, drop = FALSE) # aggregate head changes 

ovp_19 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_19 <- na.omit(ovp_19)
rownames(ovp_19) <- NULL
ovp_19$id <- "OVP_19"
rm(gps_19) # remove redundant data 


## ---------- GPS 20 ----------+
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS")
mydir <- "GPS_20" # extent directory pathway to individual GPS data 
myfiles <- list.files(path = mydir, pattern = "*.CSV", full.names = TRUE) # check file availability
gps_20 <- ldply(myfiles, read_csv) # read data 
gps_20$id <- "OVP_20" # add ID column 
gps_20 <- dplyr::select(gps_20, id, gpstime, latitude, longitude)
gps_20$gpstime <- as.POSIXct(gps_20$gpstime, origin = "1970-01-01")

# Convert latitude and longitude to decimal format
gps_20$latitude <- gps_20$latitude / 10^7
gps_20$longitude <- gps_20$longitude / 10^7

# Check time lag 
gps_20 <- gps_20 %>% mutate(time.lag = as.numeric(difftime(gpstime, lag(gpstime), units = "min")))

# Calculate the distance between consecutive GPS recordings
distance <- distVincentySphere(
  p1 = gps_20[, c("longitude", "latitude")][-nrow(gps_20), ],
  p2 = gps_20[, c("longitude", "latitude")][-1, ])

gps_20 <- data.frame(gps_20[-1, ], distance)
gps_20 <- subset(gps_20, distance < 1000)
colnames(gps_20) <- c("id", "timestamp", "latitude", "logitude", "lag", "distance")

# Access FIWI data (which contains activity and other variables of interest)
setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/data/FIWI") 
ovp_20 <- read.csv("20_FIWI.csv", header = TRUE) # read data
ovp_20$Id <- "OVP_20" # make sure ID is correct
ovp_20$Date <- as.POSIXct(ovp_20$Date, format = "%Y-%m-%d %H:%M:%OS") # set date format
ovp_20 <- arrange(ovp_20, Date) # arrange chronologically 
colnames(ovp_20) <- c("id", "timestamp", "bt", "ct", "activity", "head_down", "head_changes", "hr")

ovp_20 <- Merge(ovp_20, gps_20, id = ~ timestamp + id, all = TRUE) # merge GPS with FIWI data 

ovp_20 <- ovp_20[-c(1:9), ] # remove some erroneous timestamps 
summary(ovp_20$timestamp) # check beginning and starting date
time_sq <- as.data.frame(seq(ISOdatetime(2018,10,24,13,0,0), ISOdatetime(2023,3,13,11,0,0), by=(60*3)))
colnames(time_sq) <- "timestamp"
time_sq$timestamp <- as.POSIXct(time_sq$timestamp, format = "%Y-%m-%d %H:%M:%OS")
ovp_20 <- Merge(time_sq, ovp_20, all = TRUE, id = ~ timestamp) # merge data with time frame
ovp_20$timestamp <- cut(ovp_20$timestamp, breaks = "3 min") # break at 3 min intervals

distance <- aggregate(distance ~ timestamp, ovp_20, mean, drop = FALSE) # aggregate distances
activity <- aggregate(activity ~ timestamp, ovp_20, mean, drop = FALSE) # aggregate activity
head_down <- aggregate(head_down  ~ timestamp, ovp_20, mean, drop = FALSE) # aggregate head down duration
head_changes <- aggregate(head_changes ~ timestamp, ovp_20, mean, drop = FALSE) # aggregate head changes 

ovp_20 <- Merge(distance, activity, head_down, head_changes, all = TRUE, id = ~ timestamp)
ovp_20 <- na.omit(ovp_20)
rownames(ovp_20) <- NULL
ovp_20$id <- "OVP_20"
rm(gps_20) # remove redundant data 

## ---------- Build & clean data base ----------

# Combine data 
ovp_dat <- data.frame(rbind(ovp_01, ovp_02, ovp_08, ovp_09, ovp_10, ovp_15, ovp_19, ovp_20))
gdata::keep(ovp_dat, sure = TRUE) # clean global environment 

# Check for time lag 
ovp_dat %>%
  group_by(id) %>%
  mutate(lag = as.numeric(difftime(timestamp, lag(timestamp), units = "hours"))) %>%
  ggplot() +
  geom_density(aes(x = lag, group = id), fill = "black", alpha = 0.25) +
  scale_x_log10() 

## sort by id and time stamp & write file
ovp_dat <- arrange(gps_dat, id, timestamp)

setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data") # set wd
write.csv(ovp_dat, file = "ovp_gps_movement_activity_01.06.23.csv", row.names = FALSE)
