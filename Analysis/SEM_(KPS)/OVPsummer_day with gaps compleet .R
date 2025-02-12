rm(list=ls())
gc()

library(tidyverse)
dat <- read_csv("/Users/serpent/Documents/VHL/OVP/Data/ovp_data_10_12_24.csv")
Mydata1 <- dat

# Loop through each ID and create a boxplot

df_split1<-split(Mydata1$mean_heartrate,Mydata1$ID)
par(mfrow=c(4,4))
for(i in 1:length(df_split1)) {
 plot(df_split1[[i]], main = names(df_split1[i]))
}

layout(1)

library(feather)
library(data.table)
library(mgcv)
library(car)
library(ggplot2)
library(grid)
library(animation)
library(ggeffects)
library(MASS)
library(glmmTMB)
library(emmeans)
library (gamm4)
library(effects)
library(stats)
library(dplyr)
library(vioplot)
library(lme4)
library(lmerTest)
library(moments)
library(MuMIn)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se= sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])),
      seb=sqrt(mean(x[[col]])/(1-mean(x[[col]]))/length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

nrow(Mydata1) #4714


Mydata<-data.frame(Mydata1)
table(Mydata$season,Mydata$ID)
Mydata$season<-as.factor(Mydata$season)
Mydata$date<-as.factor(Mydata$date)
Mydata$ID<-as.factor(Mydata$ID)
Mydata$phase<-as.factor(Mydata$phase)
Mydata$phase= relevel(Mydata$phase, ref="day")
Mydata$season= relevel(Mydata$season, ref="Fall")
Mydata$BT<-as.numeric(Mydata$mean_BT_raw)
Mydata$TA<-as.numeric(Mydata$phase_mean_CT)
Mydata$HR<-as.numeric(Mydata$mean_heartrate)
Mydata$ACT<-as.numeric(Mydata$mean_activity_percent)
Mydata$WCF<-as.numeric(Mydata$mean_wcf)
Mydata$weight<-as.numeric(Mydata$weight)
Mydata$block<-as.factor(Mydata$block)
Mydata$dayf<-as.factor(Mydata$cum_day)
Mydata$daynr<-as.numeric(Mydata$cum_day)
Mydata$daynr2<-Mydata$daynr^2

#string datum omzetten naar
Mydata$Date1<-as.Date(Mydata$date)
#extract Year form Date1
Mydata$Year <- format(Mydata$Date1, format="%Y")
Mydata$Year <-as.factor(Mydata$Year)
Mydata$IDphasef<-as.factor(interaction(Mydata$phase, Mydata$ID)) 
table(Mydata$IDphasef)

##################################################################
#SELECTION summer day phase

Mydata<-subset(Mydata, season =="Summer" & phase=="day")
nrow(Mydata) #628
n_distinct(Mydata$ID)  #N= 11
n_distinct(Mydata$date)  #N= 177
n_distinct(Mydata$Year)  #N=2
table(Mydata$Year)  #data van 2019 en 2020
table(Mydata$Year,Mydata$daynr) # daynr 193-284 van 2019 en daynr 559-650 van 220
table(Mydata$Year) 
Mydata$Year<-factor(Mydata$Year, levels = c("2019", "2020"))
#create variable serial number per season
Mydata$per_summer<-ifelse(Mydata$daynr<285, "sum19", "sum20")
table(Mydata$per_summer)
table(Mydata$per_summer,Mydata$daynr)
Mydata$day_season<-ifelse(Mydata$daynr<285, Mydata$daynr-192, Mydata$daynr-651)
Mydata$day_season2<-Mydata$day_season^2
table(Mydata$day_season)

Mydata$obs = factor(1:nrow(Mydata))
Mydata$ID.per<-interaction(Mydata$ID,Mydata$per_summer)
####################################################################





########################################################################
#Data exploration
# A Missing values?
# B Outliers in Y / Outliers in X
# C Collinearity X
# D Relationships Y vs X
# E Spatial/temporal aspects of sampling design
# F Interactions (is the quality of the data good enough to include them?)
# G Zero inflation Y
# H Are categorical covariates balanced?
##############################################
# A Missing values?
colSums(is.na(Mydata))
#no missing data

table(Mydata$Year, Mydata$ID) # 
#only 10 observations from OVP01 9 in 2019 and 1 in 2020
#21 observations from OVP08 1 in 2019 and 20 in 2020
#only 9 observations from OVP18 9 in 2019 and 0 in 2019
#zero observations OVP2, ovp16

Mydata$IDn<- gsub('OVP_','',Mydata$ID)
Mydata$IDn<-as.numeric(Mydata$IDn)
table(Mydata$IDn)
#remove ID2 and ID16 
Mydata<-subset(Mydata, IDn %in% c(1,4,5,8,9,10,11,15,18,19,20) )
nrow(Mydata) #628
#remove the 1 observation of OVP08 in 2019
#Mydata<-Mydata[!(Mydata$IDn == 8 & Mydata$Year == 2019),]
#remove the 1 observation of OVP01 in 2020
#Mydata<-Mydata[!(Mydata$IDn == 1 & Mydata$Year == 2020),]
nrow(Mydata) #626

table(Mydata$Year, Mydata$IDn) 
table(Mydata$ID.per,Mydata$day_season)
# 
plot(Mydata$day_season,Mydata$ID.per,col=Mydata$ID,pch=16)
n_distinct(Mydata$ID.per)

#B OUTLIERS in predictors?
#1 cleveland dotplots for detecting outliers (response variable+ continuous covariates)

#response variable
#BT
ggplot(Mydata,aes(x=BT,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$BT,Mydata$obs) #1 outliers case 178
hist(Mydata$BT) #
skewness(Mydata$BT)#0.7431716 ok
kurtosis(Mydata$BT)#5.477701 peaked
#HR
ggplot(Mydata,aes(x=HR,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$HR,Mydata$obs) #no outliers
hist(Mydata$HR) 
skewness(Mydata$HR)#-0.2749446
kurtosis(Mydata$HR)#3.090833 peaked

#ACT
ggplot(Mydata,aes(x=ACT,y=obs,label=obs,col=ID.per))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$ACT,Mydata$obs) #no outliers
hist(Mydata$ACT) #normal distibuted =OK
skewness(Mydata$ACT)# -0.07289144
kurtosis(Mydata$ACT)#3.101462 peaked

#continuous predictors
#TA
ggplot(Mydata,aes(x=TA,y=obs,label=obs,col=ID.per))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$TA,Mydata$obs) #no outliers
hist(Mydata$TA)

#WCF
ggplot(Mydata,aes(x=WCF,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$WCF,Mydata$obs) #0 outliers
hist(Mydata$WCF)

plot(Mydata$day_season,Mydata$obs,col=Mydata$ID.per,pch=16) #no outliers
#wel duidelijk gaps te zien bij de dieren

#C COLLINEARITY between predictors

#1 GVIF per predictor

mod<-lmer(BT~TA+weight+WCF+ACT+day_season+(1|ID.per), data = Mydata, REML=T)
summary(mod)
vif(mod)
#high vif TA and WCF (as expected)
cor(Mydata$WCF,Mydata$TA)
#0.998
cor(Mydata$ACT,Mydata$TA)
#0.126
cor(Mydata$ACT,Mydata$WCF)
#0.1299957

#remove WCF
mod<-lmer(BT~TA+weight+ACT+day_season+(1|ID.per), data = Mydata, REML=T)
summary(mod)

vif(mod) #no problems


###conclusion: strong correlation WCF and TA, I will remove WCF
#############################################################################
#D. Relationships Y vs X
#relation predictors and dependent variable passage

#TA as predictor

#TA-BT
ggplot(Mydata, aes(y=BT,x=TA))+geom_point()+geom_smooth(method=loess)

#TA-HR
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=loess) #dal parabool
#quadratic model TA

#TA-ACT
ggplot(Mydata, aes(y=ACT,x=TA))+geom_point()+geom_smooth(method=loess)
#linear relation

#Conclusion: linear relation TA and ACT; 
#linear relation TA with BT, HR and ACT
#
###################################################
#day_season
#day_season-BT
cor(Mydata$day_season,Mydata$BT)
#[1] -0.1325862
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm1<-gamm(BT~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm1$gam)  #edf=4.575!!
plot(gamm1$gam)
#quadatic relation?

#day_season-HR

ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=loess)
#linear relation

#day_season-ACT
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=loess)

#linear relation

#Conclusion: linear relation day_season and ACT; 
#quadratic relation day_season and BT mountaincparabola;
#linear relation day_season and HR ; 

################################################################################
#predictor ACT for BT and HR response variable

#ACT-BT

ggplot(Mydata, aes(y=BT,x=ACT))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=ACT))+geom_point()+geom_smooth(method=gam)
gamm2<-gamm(BT~s(ACT),random= list(ID.per =~1),data=Mydata)
summary(gamm2$gam)  #edf=1!!
plot(gamm2$gam)


#ACT-HR
ggplot(Mydata, aes(y=HR,x=ACT))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=ACT))+geom_point()+geom_smooth(method=gam)
gamm3<-gamm(HR~s(ACT),random= list(ID.per =~1),data=Mydata)
summary(gamm3$gam)  #edf=1!!
plot(gamm3$gam)
##################################


####################################
#weight
#weight->BT
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=loess)
#no relation 

#weight->HR
ggplot(Mydata, aes(y=HR,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=weight))+geom_point()+geom_smooth(method=gam)
gamm4<-gamm(HR~s(weight),random= list(ID.per =~1),data=Mydata)
summary(gamm4$gam)  #edf=1!!
plot(gamm4$gam)
#linear relation

#weight->ACT
ggplot(Mydata, aes(y=ACT,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=weight))+geom_point()+geom_smooth(method=gam)
gamm5<-gamm(ACT~s(weight),random= list(ID.per =~1),data=Mydata)
summary(gamm5$gam)  #edf=1.9!!
plot(gamm5$gam)
#quadratic??


#################################################################
#E Spatial/temporal/repeated measures aspects of sampling design
table(Mydata$ID)
table(Mydata$IDn)

#controle op autocorrelatie
#ACF en PacF plots voor elke tijdreeks
#BT
getal1=''
getal2=''
for(getal1 in (c(1,4,5,8,9,10,11,15,18,19,20))){
  for(getal2 in c('sum19', 'sum20')){
    test <- subset(Mydata, IDn == getal1 & per_summer==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$BT, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$BT,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}    

#sterke autocorrelatie AR1 

#HR
getal1=''
getal2=''
for(getal1 in (c(1,4,5,8,9,10,11,15,18,19,20))){
  for(getal2 in c('sum19', 'sum20')){
    test <- subset(Mydata, IDn == getal1 & per_summer==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$HR, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$HR,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}   
#ACT
getal1=''
getal2=''
for(getal1 in (c(1,4,5,8,9,10,11,15,18,19,20))){
  for(getal2 in c('sum19', 'sum20')){
    test <- subset(Mydata, IDn == getal1 & per_summer==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$ACT, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$ACT,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}   

#er is sterke autocorrelatie AR1
#Echter er zitten gaps in de waarnemingen:
layout(1)
plot(Mydata$day_season,Mydata$obs,col=Mydata$ID.per,pch=16) #no outliers
#wel duidelijk gaps te zien

#Duidelijke gaps per combinatie van ID en phase
plot(Mydata$day_season,Mydata$obs,col=Mydata$ID,pch=16) #no outliers

plot(Mydata$day_season,Mydata$ID.per,col=Mydata$IDn,pch=16)
table((Mydata$ID.per))
table(as.numeric(Mydata$ID.per))
##################################
#random factor ID.per
#BT
Mydata$ID.pern<-as.numeric(Mydata$ID.per)
boxplot(BT~ID.pern,data=Mydata)
vioplot(BT~ID.pern,data=Mydata) 
#small differences in BT
#HR
boxplot(HR~ID.pern,data=Mydata)
vioplot(HR~ID.pern,data=Mydata)
#large differences in HR between the animals 
#ACT
boxplot(ACT~ID.pern,data=Mydata)
vioplot(ACT~ID.pern,data=Mydata)
#large differences in ACT between the animals in both periods
########################################################################
# F Interactions (is the quality of the data good enough to include them? 

#Interactions
#per_summer*act
#BT
boxplot(ACT~per_summer,data=Mydata) #overlap thus its ok to test for interaction
ggplot(Mydata, aes(x=ACT, y=BT, color=per_summer)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_summer BT



#HR

ggplot(Mydata, aes(x=ACT, y=HR, color=per_summer)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #small interaction (N.S)

mod7 <- lme(HR~per_summer*ACT, random=~1+ACT|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod7) #per_summer:ACT P= 0.0517 
summary(mod7)
AIC(mod7) #REML:3782.437

ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(mod7)), col = "black") +
  geom_smooth(aes(y = predict(mod7), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

#Interactions
#per_summer*TA
boxplot(TA~per_summer,data=Mydata) #overlap thus its ok to test for interaction

#BT

ggplot(Mydata, aes(x=TA, y=BT, color=per_summer)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_summer BT
#no interaction


#HR

ggplot(Mydata, aes(x=TA, y=HR, color=per_summer)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_summer BT
#small interactionut N.S.

mod8 <- lme(HR~per_summer*TA, random=~1+ACT|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod8) #per_summer:TAP= 0.4194
summary(mod8)
AIC(mod8)

#ACT

ggplot(Mydata, aes(x=TA, y=ACT, color=per_summer)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_summer BT
#small interaction N.S.

mod8 <- lme(ACT~per_summer*TA, random=~1|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod8) #per_summer:TA=  0.2426
summary(mod8)
AIC(mod8)


#H. Are categorical covariates balanced?
table(Mydata$IDn)  #no
table(Mydata$per_summer)#no
###################################################################################

#Descriptives;
min(Mydata$TA)  # 12.31125
max(Mydata$TA)  # 39.25
mean(Mydata$TA) #23.37216
sd(Mydata$TA) #4.341544
cvTA <- sd(Mydata$TA) / mean(Mydata$TA) * 100
cvTA #18.57571
quantile(Mydata$TA,c(.025,1-.025))
#2.5%    97.5% 
#15.57464 32.62375 

min(Mydata$ACT) #0.1
max(Mydata$ACT)  #24.1
mean(Mydata$ACT) #11.46532
sd(Mydata$ACT) #3.91261
quantile(Mydata$ACT,c(.025,1-.025))
#2.5%    97.5% 
#3.69375 19.30625 

min(Mydata$BT) #38.57235
max(Mydata$BT)   #40.525
mean(Mydata$BT) #39.21917
sd(Mydata$BT) #0.2163872
quantile(Mydata$BT,c(.025,1-.025))
#2.5%    97.5% 
#38.85875 39.67500 

min(Mydata$HR)  #36.5
max(Mydata$HR)   #95.2
mean(Mydata$HR) #66.9196
sd(Mydata$HR) #10.08715
quantile(Mydata$HR,c(.025,1-.025))
#2.5%    97.5% 
#45.29804 84.87589 
###################################################################################
#modelselection


#HR model

Mydata$per_summer<-as.factor(Mydata$per_summer)
Mydata$per_summer = relevel(Mydata$per_summer, ref="sum20")
#global model
modsel1HR <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)

anova(modsel1HR) 
summary(modsel1HR)
AIC(modsel1HR)#REML3708.028
BIC(modsel1HR)#REML3774.401

ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1HR)), col = "black") +
  geom_smooth(aes(y = predict(modsel1HR), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#duidelijk te zien dat slope ACT verschilt per ID.per
#does global model Modsel1 improve by adding or deleting a random slope variable?


#does global model Modsel1 improve by adding or deleting a random slope variable?
#remove random slope ACT
modsel1HR0 <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR0) #REML 3726.722  higher
BIC(modsel1HR0)  #REML 3784.246 much higher
#global model does not improve if you remove random slope ACT


#adding random slope I(ACT^2) bescause ACT has a quadratic effect
modsel1HR1 <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT+I(ACT^2)|ID.per,method = "REML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)

anova(modsel1HR1) 
summary(modsel1HR1)
AIC(modsel1HR1)#REML3708.822
BIC(modsel1HR1)#REML3788.47
#so I will choose the model without the random slope I(ACT^2)


#adding random slope TA
modsel1HR2 <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT+TA|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR2) #REML 3707.389 delta AIC <2 pts lower
BIC(modsel1HR2)  #REML 3787.036 much higher 
#so I will choose the model without the random slope TA

#does top model improve when adding  random slope day_season
modsel1HR3 <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT+day_season|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR3) #REML 3708.042 higher
BIC(modsel1HR3)#REML 3787.689 much higher
#so I will choose the model without the random slope day_season

#does top model improve when removing the nugget
modsel1HR4 <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
AIC(modsel1HR4) #REML 3839.721 much higher
BIC(modsel1HR4)#REML 3901.669 much higher
#so I will choose the model with nugget


#Modsel1HR  remains global model!!!
#Modsel1HR with ML estimator
modsel1HR <- lme(HR~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "ML",
                 data = Mydata, 
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HR) #3744.992

modsel1HRa <- lme(HR~per_summer+ACT +I(ACT^2)+ TA+I(TA^2)+weight+day_season+I(day_season^2), random=~1+ACT|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HRa) #3744.992

step_modelHR <- stepAIC(modsel1HRa, scope = list(lower = ~ ACT + TA, upper = ~ .),
                      direction = "both")
AIC(step_modelHR) #3742.288
summary(step_modelHR)  #dit is het topmodel

modsel1HRtop1<- lme(HR~ACT+I(ACT^2)+TA+weight+day_season, random=~1+ACT|ID.per,method = "ML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HRtop) #3742.288

modsel1HRtop<- lme(HR~poly(ACT,2)+TA+weight+day_season, random=~1+ACT|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRtop)#3742.288
summary(modsel1HRtop)
modsel1HRtopTST<- lme(HR~ACT+TA+weight+day_season, random=~1+ACT|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRtopTST)#3745.068

anova(modsel1HRtop,modsel1HRtopTST)
AIC(modsel1HRtop,modsel1HRtopTST)
BIC(modsel1HRtop,modsel1HRtopTST)
#########################################
library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1HRtop)
icc(modsel1HRtop)

#grafiek met effect ACT per dier op HR
ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1HRtop)), col = "black") +
  geom_smooth(aes(y = predict(modsel1HRtop), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#my random effects 
summary(modsel1HRtop) #
ranef(modsel1HRtop)
#hier is duidelijk te zien dat het effect van ACT verschilt per dier
plot(allEffects(modsel1HRtop1))


#########################################
modsel1BT <- lme(BT~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                  data = Mydata,
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1BT)
AIC(modsel1BT) #-521.4209
BIC(modsel1BT) #-463.8976

#does global model modsel 1BT improve buy adding or deleting a random slope variable?

#global model modsel1BT with random slope ACT

modsel1BT0<- lme(BT~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
                    data = Mydata,
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT0) #REML -518.2682 higher
BIC(modsel1BT0)  #REML-451.8951higher
#global model does not improve by adding random slope ACT

#does global model improve when adding  random slope TA
modsel1BT1<- lme(BT~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID.per,method = "REML",
                    data = Mydata,
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT1) #REML -518.9754 AIC higher
BIC(modsel1BT1)  #REML -452.6024 BIC higher 
#so I will choose the model without the random slope TA

#does global model improve when adding  random slope day_season
modsel1BT2<- lme(BT~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID.per,method = "REML",
                    data = Mydata,
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT2) #REML -521.5958 0.1 pt lower
BIC(modsel1BT2)#REML -455.2228 much higher
#so I will choose the model without the random slope day_season

#does top model improve when removing the nugget
modsel1BT4<- lme(BT~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                    data = Mydata,
                    correlation = corGaus(form= ~ day_season|ID.per),
                    na.action=na.exclude)
AIC(modsel1BT4) #REML -468.6974  much higher
BIC(modsel1BT4)  #REML-415.599 much higher

#nugget remains in model

#modsel1BT remains gobal model!!!
#ML model
modsel1BT <- lme(BT~per_summer+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1BT) #ML -550.5492
modsel1BTa <- lme(BT~per_summer+ACT+I(ACT^2)+TA+I(TA^2)+weight+day_season+I(day_season^2), random=~1|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1BTa) #-550.5492


step_modelBT <- stepAIC(modsel1BTa, scope = list(lower = ~ ACT + TA, upper = ~ .),
                        direction = "both")

summary(step_modelBT)
modsel1BTtop<- lme(BT~ACT+TA+poly(day_season,2), random=~1|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1BTtop) #ML -557.1377

modsel1BTtop1<- lme(BT~ACT+TA+day_season+I(day_season^2), random=~1|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1BTtop1) #ML -557.1377

library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1BTtop)
icc(modsel1BTtop)

#grafiek met effect ACT per dier op BT
ggplot(Mydata,aes(x=ACT,y=BT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1BTtop)), col = "black") +
  geom_smooth(aes(y = predict(modsel1BTtop), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

plot(allEffects(modsel1BTtop))
########################################################
#ACT model
modsel1ACT <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_summer, random=~1+TA|ID.per,method = "REML",
                  data = Mydata,control =list(msMaxIter = 1000, msMaxEval = 1000),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1ACT) 
summary(modsel1ACT)
AIC(modsel1ACT) #3187.827
BIC(modsel1ACT)  #3245.392

#does global model modsel1ACT improve buy adding or deleting a random slope variable?

#global model modsel1ACT without random slope TA
modsel1ACT0 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_summer, random=~1|ID.per,method = "REML",
                  data = Mydata,control =list(msMaxIter = 1000, msMaxEval = 1000),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1ACT0) 
summary(modsel1ACT0)
AIC(modsel1ACT0) # 3206.4  much higher
BIC(modsel1ACT0) #3255.109 much higher

ggplot(Mydata,aes(x=TA,y=ACT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACT0)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACT0), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

#heeft nugget effect? model without nugget
modsel1ACT1 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_summer, random=~1+TA|ID.per,method = "REML",
                  data = Mydata,
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT1) 
summary(modsel1ACT1)
AIC(modsel1ACT1)#[1] 3194 higher
BIC(modsel1ACT1)#3247.137 higher
#Model with nugget is better

#model met TA2 levert niet een beter topmodel

#ML top global model
modsel1ACT <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_summer, random=~1+TA|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1ACT) 
summary(modsel1ACT)
AIC(modsel1ACT) #3210.42

modsel1ACTa <- lme(ACT~TA+I(TA^2)+weight+day_season+I(day_season^2)+per_summer, random=~1+TA|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1ACTa) 
summary(modsel1ACTa)
AIC(modsel1ACTa) #3210.42


step_modelACT <- stepAIC(modsel1ACTa, scope = list(lower = ~ TA, upper = ~ .),
    direction = "both")
summary(step_modelACT)  
anova(step_modelACT)
AIC(step_modelACT)  #3205.248


modsel1ACTtop <- lme(ACT~TA+per_summer, random=~1+TA|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
AIC(modsel1ACTtop) #3207.082
anova(modsel1ACTtop)
summary(modsel1ACTtop)

plot(allEffects(modsel1ACTtop))

ggplot(Mydata,aes(x=TA,y=ACT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACTtop)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACTtop), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

ranef(modsel1ACTtop)

library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1ACTtop)



ggplot(Mydata,aes(x=TAz,y=ACT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACTtop_1)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACTtop_1), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

plot(allEffects(modsel1ACTtop))

#########################################################################
#plots modellen
#HR


plot(allEffects(modsel1HRtop1, residuals=TRUE), 
     partial.residuals=list(smooth=FALSE,col="black",span=0.75,lty="dashed") )


######################################
#BT

summary(modsel1BTtop)
plot(allEffects(modsel1BTtop, residuals=TRUE), 
     partial.residuals=list(smooth=FALSE,col="black",span=0.75,lty="dashed") )

######################################
#ACT

modsel1ACTtop <- lme(ACT~TA+per_summer, random=~1+TA|ID.per,method = "ML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
summary(modsel1ACTtop)
plot(allEffects(modsel1ACTtop))
plot(allEffects(modsel1ACTtop, residuals=TRUE), 
     partial.residuals=list(smooth=FALSE,col="black",span=0.75,lty="dashed") )

##################################################################################
#controle modeleisen
#normaliteit residuen

#HR model
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1) #-0.5157941
kurtosis(Mydata$E1)  #9.149918
#BTmodel
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)#0.9932993
kurtosis(Mydata$E1) #10.86584
#ACTmodel
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)#0.05364023
kurtosis(Mydata$E1) #4.811322


#autocorrelation residuals model
#HR model
summary(modsel1HRtop)
anova(modsel1HRtop)
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(1,4,5,8,9,10,11,15,18,19,20))){
  for(getal2 in c('sum19', 'sum20')){
    test <- subset(Mydata, IDn == getal1 & per_summer==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}    


#geen autocorrelatie


#BT model
summary(modsel1BTtop)
anova(modsel1BTtop)
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(1,4,5,8,9,10,11,15,18,19,20))){
  for(getal2 in c('sum19', 'sum20')){
    test <- subset(Mydata, IDn == getal1 & per_summer==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}
#geen autocorrelatie (alleen bij dier 5 in 2019 en dier 20 in 2020)

#ACT model
summary(modsel1ACTtop)
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(1,4,5,8,9,10,11,15,18,19,20))){
  for(getal2 in c('sum19', 'sum20')){
    test <- subset(Mydata, IDn == getal1 & per_summer==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}
#geen autocorrelatie
###############################################
#How good is the fit of the model?
#HR model
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
Mydata$F1 <- fitted(modsel1HRtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = Mydata$F1,
     y = Mydata$HR,
     xlab = "Fitted values (with re)",
     ylab = "HR")
abline(a=0,b=1,col=2,lty=2,lwd=4)     

output <- data.frame(observed = Mydata$HR, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)

cor(Mydata$F1, Mydata$HR)
#0.8462317  #prima
########################################################
#BT model
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
Mydata$F1 <- fitted(modsel1BTtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = Mydata$F1,
     y = Mydata$BT,
     xlab = "Fitted values (with re)",
     ylab = "BT")
abline(a=0,b=1,col=2,lty=2,lwd=4)     

output <- data.frame(observed = Mydata$BT, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)

cor(Mydata$F1, Mydata$BT)
#0.5923708

##########################################################
#ACT model
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
Mydata$F1 <- fitted(modsel1ACTtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = Mydata$F1,
     y = Mydata$ACT,
     xlab = "Fitted values (with re)",
     ylab = "ACT")
abline(a=0,b=1,col=2,lty=2,lwd=4)     

output <- data.frame(observed = Mydata$ACT, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)

cor(Mydata$F1, Mydata$ACT)
#0.6577426
###############################################################################
library(piecewiseSEM)
#pathway analysis
modsel1HRtop<- lme(HR~poly(ACT,2)+TA+weight+day_season, random=~1+ACT|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRtop)#3742.288
Mydata$ACT2<-(Mydata$ACT)^2
Mydata$day_season2<-(Mydata$day_season)^2

modsel1HRtop2<- lme(HR~ACT+ACT2+TA+weight+day_season, random=~1+ACT|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRtop2)#3742.288


modsel1BTtop1<- lme(BT~ACT+TA+day_season+I(day_season^2), random=~1|ID.per,method = "ML",
                    data = Mydata,
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BTtop1) #-557.1377
modsel1BTtop2<- lme(BT~ACT+TA+day_season+day_season2, random=~1|ID.per,method = "ML",
                    data = Mydata,
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BTtop2) #-557.1377

path_test<-psem(modsel1HRtop2,modsel1BTtop2,modsel1ACTtop)
summary(path_test)
path_testa<-psem(modsel1HRtop2,modsel1BTtop2,HR%~~%BT,modsel1ACTtop)
summary(path_testa)

plot(path_testa,node_attrs = data.frame(shape = "rectangle", color = "black", 
                                       fillcolor = "orange"),show = "std")
rsquared(path_test, method = NULL)
AIC(path_test) #65.184
############################################################################

#############################
##########################################################################
# plot of pearson residuals against the predicted values
#HR model
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
Mydata$F1 <- fitted(modsel1HRtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = Mydata$F1 , 
     y = Mydata$E1 ,
     xlab = "Fitted values (with re)",
     ylab = "normalized residuals")
abline(h = 0, lty = 2)

output <- data.frame(observed = Mydata$E1, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)
#o.k. 2 outliers

#BT model
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
Mydata$F1 <- fitted(modsel1BTtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = Mydata$F1 , 
     y = Mydata$E1 ,
     xlab = "Fitted values (with re)",
     ylab = "normalized residuals")
abline(h = 0, lty = 2)  #1 outlier

output <- data.frame(observed = Mydata$E1, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)
#o.k. (spread  higher for higher fitted values, 1 outlier!!)

#ACT model
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
Mydata$F1 <- fitted(modsel1ACTtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = Mydata$F1 , 
     y = Mydata$E1 ,
     xlab = "Fitted values (with re)",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

output <- data.frame(observed = Mydata$E1, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)
#o.k. 


