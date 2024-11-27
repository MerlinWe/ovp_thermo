########## RAAK Proj: OVP data exploration GPS distance ~ Activity ##########

setwd("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/GPS") # set WD to location including raw GPS data

library(tidyverse)      
library(data.table)

options(scipen = 999) # disable scientific notation 

dat <- read.csv("/Users/serpent/OneDrive - Van Hall Larenstein/RAAK/Data/ovp_gps_movement_activity_01.06.23.csv")
table(dat$id) # check available data per ID 

# Make some quick plots

plot(dat$distance)  # distance travelled based on GPS 
plot(dat$activity)  # activity (raw index of locomotion) based on FIWI 
plot(dat$head_down) # duration spend head down based on FIWI
plot(dat$head_changes) # number of head changes based on FIWI 

hist(dat$distance)
hist(log(dat$distance)) # slightly left skewed but better 

hist(dat$activity)
hist(dat$activity) # log gives no improvment 

hist(dat$head_down)
hist(dat$head_changes)

## Write VIF function(s) 

myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

vif <- corvif(dat[ ,c(2:4)]) 

# Correlations of the variables

#             distance    activity   head_down
# distance   1.0000000  0.21246742 -0.01920570
# activity   0.2124674  1.00000000 -0.06693257
# head_down -0.0192057 -0.06693257  1.00000000

# Variance inflation factors

# GVIF
# distance  1.047304
# activity  1.051629
# head_down 1.004526

# Use Spearman's since activity is not normally distributed 
cor.test(dat$distance, dat$activity, method = "spearman")

## There does not seem to be a strong connection between the distance moved and the activity. 

plot(dat$head_down)
hist(dat$head_down)

cor.test(dat$head_down, dat$distance, method = "spearman")
cor.test(dat$head_changes, dat$distance, method = "spearman")
