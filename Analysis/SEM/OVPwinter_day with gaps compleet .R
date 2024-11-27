rm(list=ls())
gc()
setwd("C:/Users/kuipe_003/OneDrive - Van Hall Larenstein/project OVP WTN en MeW/analyses per dayphase compleet")
library(readr)
ovp_data_6_2_23 <- read_csv("ovp_data_11.04.23.csv")
nrow(ovp_data_6_2_23) #4708
Mydata1<-ovp_data_6_2_23
df_split1<-split(Mydata1$mean_heartrate,Mydata1$ID)
par(mfrow=c(4,4))
# Loop through each ID and create a boxplot

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

nrow(Mydata1) #4708
#head(Mydata)
Mydata<-data.frame(Mydata1)
Mydata$season<-as.factor(Mydata$season)
table(Mydata$season,Mydata$ID)
Mydata$date<-as.factor(Mydata$date)
Mydata$ID<-as.factor(Mydata$ID)
Mydata$phase<-as.factor(Mydata$phase)
Mydata$phase= relevel(Mydata$phase, ref="day")
Mydata$season= relevel(Mydata$season, ref="Fall")
Mydata$BT<-as.numeric(Mydata$mean_BT_clean)
Mydata$TA<-as.numeric(Mydata$mean_collar_temp)
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
#selection winter season
Mydata<-subset(Mydata, season =="Winter" & phase=="day" )
nrow(Mydata)
#382
n_distinct(Mydata$ID)  #N= 10
n_distinct(Mydata$date)  #N= 166
n_distinct(Mydata$Year)  #N=3
table(Mydata$Year)  #data van winter 2018-2019 en 2019-2020

table(Mydata$Year,Mydata$daynr)

#2 winterperiode periode 1 daynr  11-100 and period 2 daynr 376-466

#create variable serial number per season
Mydata$per_winter<-ifelse(Mydata$daynr<101, "win18-19", "win19-20")
table(Mydata$per_winter,Mydata$daynr)
Mydata$day_season<-ifelse(Mydata$daynr<101, Mydata$daynr-10, Mydata$daynr-375)
Mydata$day_season2<-Mydata$day_season^2

table(Mydata$per_winter,Mydata$day_season)
table(Mydata$per_winter)

Mydata$obs = factor(1:nrow(Mydata))

Mydata$ID.per<-interaction(Mydata$ID,Mydata$per_winter)
table(Mydata$ID.per)
###############################################################
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

Mydata$IDn<- gsub('OVP_','',Mydata$ID)
Mydata$IDn<-as.numeric(Mydata$IDn)
table(Mydata$per_winter, Mydata$IDn) # 
#remove  ID2 (n=3;1)
Mydata<-subset(Mydata, IDn %in% c(4,5,8,9,10,15,18,19,20) )
nrow(Mydata) #378
#remove the 1 observation of OVP08 in 2019-2020
Mydata<-Mydata[!(Mydata$IDn == 8 & Mydata$per_winter == "win19-20" ),]
nrow(Mydata) #377

table(Mydata$per_winter, Mydata$IDn) # 


#B OUTLIERS in predictors?
#1 cleveland dotplots for detecting outliers (response variable+ continuous covariates)



#response variable
#BT
layout(1)
ggplot(Mydata,aes(x=BT,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$BT,Mydata$obs) #4 outliers case 153 t/m 156
hist(Mydata$BT) #
skewness(Mydata$BT)#1.324574  rechtsscheef
kurtosis(Mydata$BT)#7.983484 peaked
#HR
ggplot(Mydata,aes(x=HR,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$HR,Mydata$obs) #no outliers
hist(Mydata$HR) 
skewness(Mydata$HR)#0.8819127
kurtosis(Mydata$HR)#3.56961 peaked

#ACT
ggplot(Mydata,aes(x=ACT,y=obs,label=obs,col=ID.per))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$ACT,Mydata$obs) #no outliers
hist(Mydata$ACT) #normal distibuted =OK
skewness(Mydata$ACT)# 0.4303603
kurtosis(Mydata$ACT)#3.453271 peaked

#continuous predictors
#TA
ggplot(Mydata,aes(x=TA,y=obs,label=obs,col=ID.per))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$TA,Mydata$obs) #4 outliers 342,343,344, 370
hist(Mydata$TA)

#WCF
ggplot(Mydata,aes(x=WCF,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$WCF,Mydata$obs) #4 outliers same as TA
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
#0.9647966
cor(Mydata$ACT,Mydata$TA)
#-0.1584935
cor(Mydata$ACT,Mydata$WCF)
#-0.1456745

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
ggplot(Mydata, aes(y=BT,x=TA))+geom_point()+geom_smooth(method=gam)
gamm1<-gamm(BT~s(TA),random= list(ID.per =~1),data=Mydata)
summary(gamm1$gam)  #edf 3.548 nonlinear relation
plot(gamm1$gam)
#quadratic

#TA-HR
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=loess) #dal parabool
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=gam)
gamm2<-gamm(HR~s(TA),random= list(ID.per =~1),data=Mydata)
summary(gamm2$gam)  #edf=1 linear relation

#TA-ACT
ggplot(Mydata, aes(y=ACT,x=TA))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=TA))+geom_point()+geom_smooth(method=gam)
gamm3<-gamm(ACT~s(TA),random= list(ID.per =~1),data=Mydata)
summary(gamm3$gam) #edf=1
#linear relation

#Conclusion: linear relation TA and ACT; 
# quadratic TA with BT
#
###################################################
#day_season
#day_season-BT
cor(Mydata$day_season,Mydata$BT)
#[1] -0.2863
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm4<-gamm(BT~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm4$gam)  #edf=1!
plot(gamm4$gam)
#linear relation

#day_season-HR

ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm5<-gamm(HR~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm5$gam)  #edf=4.991
plot(gamm5$gam)
#quadratic relation

#day_season-ACT
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm6<-gamm(ACT~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm6$gam)  #edf=1
plot(gamm6$gam)
#linear relation

#Conclusion: linear relation day_season and BT, ACT; 
#quadratic relation day_season and HR valleyparabola;

################################################################################
#predictor ACT for BT and HR response variable

#ACT-BT

ggplot(Mydata, aes(y=BT,x=ACT))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=ACT))+geom_point()+geom_smooth(method=gam)
gamm7<-gamm(BT~s(ACT),random= list(ID.per =~1),data=Mydata)
summary(gamm7$gam)  #edf=1!!
plot(gamm7$gam)
#linear relation

#ACT-HR
ggplot(Mydata, aes(y=HR,x=ACT))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=ACT))+geom_point()+geom_smooth(method=gam)
gamm8<-gamm(HR~s(ACT),random= list(ID.per =~1),data=Mydata)
summary(gamm8$gam)  #edf=1!!
plot(gamm8$gam)

#linear relation
##################################


####################################
#weight
#weight->BT
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=gam)
gamm9<-gamm(BT~s(weight, k=5),random= list(ID.per =~1),data=Mydata)
summary(gamm9$gam)  #edf=1!!
plot(gamm9$gam)
#linear relation 

#weight->HR
ggplot(Mydata, aes(y=HR,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=weight))+geom_point()+geom_smooth(method=gam)
gamm10<-gamm(HR~s(weight,k=5),random= list(ID.per =~1),data=Mydata)
summary(gamm10$gam)  #edf=1!!
plot(gamm10$gam)
#linear relation

#weight->ACT
ggplot(Mydata, aes(y=ACT,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=weight))+geom_point()+geom_smooth(method=gam)
gamm11<-gamm(ACT~s(weight,k=5),random= list(ID.per =~1),data=Mydata)
summary(gamm11$gam)  #edf=1
plot(gamm11$gam)
#linear??


#################################################################
#E Spatial/temporal/repeated measures aspects of sampling design
table(Mydata$ID)
table(Mydata$IDn)

#controle op autocorrelatie
#ACF en PacF plots voor elke tijdreeks
table(Mydata$per_winter)
#BT
getal1=''
getal2=''
for(getal1 in (c(4,5,8,9,10,15,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_winter==getal2) 
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
for(getal1 in (c(4,5,8,9,10,15,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_winter==getal2) 
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
for(getal1 in (c(4,5,8,9,10,15,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_winter==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$ACT, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$ACT,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}   

#er is indere sterke autocorrelatie AR1
#Echter er zitten gaps in de waarnemingen:
layout(1)
plot(Mydata$day_season,Mydata$obs,col=Mydata$ID.per,pch=16) #no outliers
#wel duidelijk gaps te zien

#Duidelijke gaps per combinatie van ID en phase
plot(Mydata$day_season,Mydata$ID.per,col=Mydata$IDn,pch=16)
n_distinct(Mydata$ID.per)
table((Mydata$ID.per))
table(as.numeric(Mydata$ID.per))
##################################
#random factor ID.per
#BT
Mydata$ID.pern<-as.numeric(Mydata$ID.per)
boxplot(BT~ID.pern,data=Mydata)
vioplot(BT~ID.pern,data=Mydata) 
ggplot(Mydata, aes(ID.per, BT, fill=factor(ID))) +
  geom_boxplot()+ coord_flip()

# differences in BT, between years and between antmals within a year
#HR
boxplot(HR~ID.pern,data=Mydata)
vioplot(HR~ID.pern,data=Mydata)
#large differences in HR between the animals 
ggplot(Mydata, aes(ID.per, HR, fill=factor(ID))) +
  geom_boxplot()+ coord_flip()
table(Mydata$IDn)

#ACT
boxplot(ACT~ID.pern,data=Mydata)
vioplot(ACT~ID.pern,data=Mydata)
ggplot(Mydata, aes(ID.per, ACT, fill=factor(ID))) +
  geom_boxplot()+ coord_flip()
table(Mydata$IDn)
#large differences in ACT between the animals in both periods
########################################################################
# F Interactions (is the quality of the data good enough to include them? 

#Interactions
#per_winter*act
#BT
boxplot(ACT~per_winter,data=Mydata) #overlap thus its ok to test for interaction
ggplot(Mydata, aes(x=ACT, y=BT, color=per_winter)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_winter BT



#HR

ggplot(Mydata, aes(x=ACT, y=HR, color=per_winter)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #small interaction (N.S)

mod7 <- lme(HR~per_winter*ACT, random=~1+ACT|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod7) #per_winter:ACT P= 0.8157 
summary(mod7)
AIC(mod7) #REML:2211.077
#no interation

ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(mod7)), col = "black") +
  geom_smooth(aes(y = predict(mod7), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

#Interactions
#per_winter*TA
boxplot(TA~per_winter,data=Mydata) #overlap thus its ok to test for interaction

#BT

ggplot(Mydata, aes(x=TA, y=BT, color=per_winter)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_winter BT
#no interaction


#HR

ggplot(Mydata, aes(x=TA, y=HR, color=per_winter)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_winter BT
#interaction TA*per-winter-> HR ?. NS!!

mod8a <- lme(HR~per_winter*TA, random=~1|ID.per,method = "ML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod8a)  #per_winter:TA P= 0.5643
summary(mod8a)
AIC(mod8a)
#no interaction

#ACT

ggplot(Mydata, aes(x=TA, y=ACT, color=per_winter)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_winter BT
#small interaction N.S.

mod9 <- lme(ACT~per_winter*TA, random=~1|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod9) #per_winter:TA=  0.9640
summary(mod9)
AIC(mod9)


#H. Are categorical covariates balanced?
table(Mydata$IDn)  #no
table(Mydata$per_winter)#no
###################################################################################
###################################################################################
###################################################################################
#Descriptives;
min(Mydata$TA)  #0.9640
max(Mydata$TA)  # 33.38
mean(Mydata$TA) # 9.339567
sd(Mydata$TA) #4.341888
cvTA <- sd(Mydata$TA) / mean(Mydata$TA) * 100
cvTA # 46.48918
quantile(Mydata$TA,c(.025,1-.025))
#2.5%    97.5% 
#1.2194 19.3840 

min(Mydata$ACT) #1.35
max(Mydata$ACT)  #25.5
mean(Mydata$ACT) #10.72463
sd(Mydata$ACT) #4.313494
quantile(Mydata$ACT,c(.025,1-.025))
#2.5%    97.5% 
#3.10 20.81 

min(Mydata$BT) # 38.37111
max(Mydata$BT)   #40.37276
mean(Mydata$BT) #39.03841
sd(Mydata$BT) #0.2440743
quantile(Mydata$BT,c(.025,1-.025))
#2.5%    97.5% 
#38.65096 39.56000

min(Mydata$HR)  #19.8
max(Mydata$HR)   #87.10714
mean(Mydata$HR) #53.68219
sd(Mydata$HR) #11.52309
quantile(Mydata$HR,c(.025,1-.025))
#2.5%    97.5% 
#38.28000 80.63692
###################################################################################
#modelselection


#HR model

Mydata$per_winter<-as.factor(Mydata$per_winter)
Mydata$per_winter = relevel(Mydata$per_winter, ref="win19-20")
#global model
modsel1HR <- lme(HR~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)

anova(modsel1HR) 
summary(modsel1HR)
AIC(modsel1HR)#REML2132.928
BIC(modsel1HR)#2183.733

ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1HR)), col = "black") +
  geom_smooth(aes(y = predict(modsel1HR), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#NIET duidelijk te zien dat slope ACT verschilt per ID.per??

#does global model Modsel1 improve by adding or deleting a random slope variable?
#geeft convergentie problemen als ik poly(ACT,2) in het model heb
#adding random slope ACT
modsel1HR0 <- lme(HR~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
                     data = Mydata,control = lmeControl(opt = list(maxIter = 100)),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
#problemen

#modsel1HR zonder poly(ACT,2)
modsel1HRa <- lme(HR~per_winter+ACT+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)

anova(modsel1HRa) 
summary(modsel1HRa)
AIC(modsel1HRa)#REML2146.072
BIC(modsel1HRa)#2193.002
#Nu met random slope voor ACT
modsel1HRb <- lme(HR~per_winter+ACT+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)

anova(modsel1HRb) 
summary(modsel1HRb)
AIC(modsel1HRb)#REML2147.398 hoger
BIC(modsel1HRb)#REML 2202.15 hoger

#global model does not improve if you add random slope ACT

#adding random slope TA
modsel1HR1 <- lme(HR~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR1) #REML 2133.376 hoger
BIC(modsel1HR1)  #REML 2191.998  hoger 
#so I will choose the model without the random slope TA

#does top model improve when adding  random slope day_season
modsel1HR2 <- lme(HR~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR2) #REML 2163.087 higher
BIC(modsel1HR2)#REML  2221.708 much higher
#so I will choose the model without the random slope day_season

#does model improve when removing the nugget
modsel1HR3 <- lme(HR~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
AIC(modsel1HR3) #REML 2288.849 higher
BIC(modsel1HR3)#REML 2335.746 much higher
#so I will choose the model with the nugget

#Modsel1HR  remains global model!!!

#Modsel1HR with ML estimator
modsel1HR <- lme(HR~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HR) #2165.468
#ACT+ I(ACT^2)+TA+I(TA^2)+weight+day_season+I(day_season^2)
modsel1HRc <- lme(HR~per_winter+ACT+ I(ACT^2)+TA+I(TA^2)+weight+day_season+I(day_season^2), random=~1|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HRc) #2165.468


step_modelHR <- stepAIC(modsel1HRc, scope = list(lower = ~ ACT + TA, upper = ~ .),
                      direction = "both")
summary(step_modelHR)  #dit is het topmodel
AIC(step_modelHR)  #2164.231

#Removing I(TA^2) increases the AIC with 0.4 pts<2 thus I will remove this variable 
#(I prefer the simpler model).      


modsel1HRtop<- lme(HR~per_winter+ACT+TA+weight+poly(day_season,2), random=~1|ID.per,method = "ML",
                     data = Mydata,
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HRtop)
r.squaredGLMM(modsel1HRtop)

#########################################
library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1HRtop)
icc(modsel1HRtop)

library(rptR)
#install.packages(("rptR"))
library(devtools)

r.squaredGLMM(modsel1HRtop)
#R2m       R2c
#[1,] 0.6219341 0.6662584
plot(allEffects(modsel1HRtop))


###########################################################

modsel1BT <- lme(BT~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2),
                  random=~1|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
AIC(modsel1BT) #-152.6428
BIC(modsel1BT) #-95.57595
summary(modsel1BT)

#Does global model modsel 1BT improve buy adding or deleting a random slope variable?
#global model modsel1BT with random slope ACT

modsel1BT0<- lme(BT~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT0) #REML-149.7055 higher
BIC(modsel1BT0)  #REML-91.08421 higher
#global model does not improve by adding random slope ACT

#does global model improve when adding  random slope TA
modsel1BT1<- lme(BT~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT1) #REML -148.6466 AIC higher
BIC(modsel1BT1)  #REML -90.02533 BIC higher 
#so I will choose the model without the random slope TA

#does global model improve when adding  random slope day_season
modsel1BT2<- lme(BT~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT2) #REML -149.0293 higher
BIC(modsel1BT2)#REML -90.40801 higher
#so I will choose the model without the random slope day_season

#does top model improve when removing the nugget
modsel1BT4<- lme(BT~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per),
                    na.action=na.exclude)
AIC(modsel1BT4) #REML -126.5783 much higher
BIC(modsel1BT4)  #REML-79.68128 much higher

#modsel1BT remains gobal model!!!
#ML model
modsel1BT <- lme(BT~per_winter+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "ML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1BT) #-180.5832

modsel1BTc <- lme(BT~per_winter+ACT+ I(ACT^2)+TA+I(TA^2)+weight+day_season+I(day_season^2), random=~1|ID.per,method = "ML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1BTc) #-180.5832

step_modelBTc <- stepAIC(modsel1BTc, scope = list(lower = ~ ACT + TA+day_season, upper = ~ .),
                        direction = "both")

summary(step_modelBT)
AIC(step_modelBT) #-182.3503
#Removing I(ACT^2) increase the AIC with 0.1 pts , thus this variable will be removed 
#from the model.
modsel1BTtopA<- lme(BT~ACT+TA+I(day_season^2), random=~1|ID.per,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1BTtopA) #-184.2604
summary(modsel1BTtopA)
coef(summary(modsel1BTtopA))
#Is a model with day_season^2 better than a model with day_season?
modsel1BTtopB<- lme(BT~ACT+TA+day_season, random=~1|ID.per,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1BTtopB) #-183.9626
anova(modsel1BTtopA,modsel1BTtopB)

modsel1BTtop<- lme(BT~ACT+TA+day_season, random=~1|ID.per,method = "ML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BTtop)
library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1BTtop,digits = 4)
icc(modsel1BTtop)

plot(allEffects(modsel1BTtop))
########################################################
#ACT model
modsel1ACT <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_winter, random=~1|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT) 
summary(modsel1ACT)
AIC(modsel1ACT) #2065.374
BIC(modsel1ACT)  #2104.51

#does global model modsel1ACT improve buy adding or deleting a random slope variable?

#global model modsel1ACT with random slope TA
modsel1ACT0 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_winter, random=~1+TA|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT0) 
summary(modsel1ACT0)
AIC(modsel1ACT0) #2067.763  higher
BIC(modsel1ACT0) #2112.725  higher

#global model modsel1ACT with random slope day_season
modsel1ACT1 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_winter, random=~1+day_season|ID.per,method = "REML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                   na.action=na.exclude)
anova(modsel1ACT1) 
summary(modsel1ACT1)
AIC(modsel1ACT1) #2062.359 lower 3.1 pts (was 2065.374 )
BIC(modsel1ACT1) #2109.321  higher (was 2104)
#no random slope for day_season

ggplot(Mydata,aes(x=day_season,y=ACT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACT1)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACT1), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#no random slope for day_season

#heeft nugget effect? model with nugget
modsel1ACT2 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_winter, random=~1|ID.per,method = "REML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
anova(modsel1ACT2) 
summary(modsel1ACT2)
AIC(modsel1ACT2)#2067.374 lower
BIC(modsel1ACT2)#2110.423 higher
#I choose a model without nugget

#ML top global model
modsel1ACT <- lme(ACT~poly(TA,2)+weight+poly(day_season,2)+per_winter, random=~1|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT) 
summary(modsel1ACT)
AIC(modsel1ACT) #2080.63

modsel1ACTc <- lme(ACT~TA+I(TA^2)+weight+day_season+I(day_season^2)+per_winter, random=~1|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACTc) 
summary(modsel1ACTc)
AIC(modsel1ACTc)

step_modelACT <- stepAIC(modsel1ACTc, scope = list(lower = ~ TA, upper = ~ .),
    direction = "both")
summary(step_modelACT)
AIC(step_modelACT) #2077.392

modsel1ACTtop <- lme(ACT~TA, random=~1|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
AIC(modsel1ACTtop) #2077.392
anova(modsel1ACTtop)
summary(modsel1ACTtop)


library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1ACTtop)
r.squaredGLMM(modsel1ACTtop)


ggplot(Mydata,aes(x=TAz,y=ACT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACTtop_1)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACTtop_1), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

plot(allEffects(modsel1ACTtop))

#########################################################################
#plots modellen
#HR

plot(allEffects(modsel1HRtop, residuals=TRUE), 
     partial.residuals=list(smooth=FALSE,col="black",span=0.75,lty="dashed") )


######################################
#BT


plot(allEffects(modsel1BTtop, residuals=TRUE), 
     partial.residuals=list(smooth=FALSE,col="black",span=0.75,lty="dashed") )

######################################
#ACT


plot(allEffects(modsel1ACTtop))
plot(allEffects(modsel1ACTtop, residuals=TRUE), 
     partial.residuals=list(smooth=FALSE,col="black",span=0.75,lty="dashed") )

##################################################################################
#controle modeleisen
#normaliteit residuen

#HR model
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)
kurtosis(Mydata$E1)  #7.568504
#BTmodel
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)#1.185421
kurtosis(Mydata$E1) #13.03594
#ACTmodel
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)#0.558323
kurtosis(Mydata$E1) #3.94697


#autocorrelation residuals model
#HR model
summary(modsel1HRtop)
anova(modsel1HRtop)
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(4,5,8,9,10,15,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_winter==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}    
#no autocorrelation

Mydata$per_winter
#BT model
summary(modsel1BTtop)
anova(modsel1BTtop)
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(4,5,8,9,10,15,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_winter==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}
#no autocorrelation

#ACT model
summary(modsel1ACTtop)
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(4,5,8,9,10,15,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_winter==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}
#no autocorrelation
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
#high  fitted values have pos residuals ; I misssed a non-linear correlation? !
output <- data.frame(observed = Mydata$HR, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)

cor(Mydata$F1, Mydata$HR)
#0.8644974  #prima
########################################################
#BT model
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
Mydata$F1 <- fitted(modsel1BTtop,type="normalized")

par(mfrow = c(1,1), cex.lab = 1.5, mar = c(5,5,2,2))
plot(x = Mydata$F1,
     y = Mydata$BT,
     xlab = "Fitted values (with re)",
     ylab = "BT")
abline(a=0,b=1,col=2,lty=2,lwd=4)  #4 outliers   

output <- data.frame(observed = Mydata$BT, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)

cor(Mydata$F1, Mydata$BT)
#0.5617718

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
#0.5799971
###############################################################################
library(piecewiseSEM)
#pathway analysis
modsel1HRtop1<- lme(HR~per_winter+ACT+TA+weight+day_season+day_season2, random=~1|ID.per,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)

AIC(modsel1HRtop1,modsel1HRtop)
path_test<-psem(modsel1HRtop1,modsel1BTtop,modsel1ACTtop)
summary(path_test)
path_testa<-psem(modsel1HRtop1,modsel1BTtop,HR%~~%BT,modsel1ACTtop)
summary(path_testa)

plot(path_testa,node_attrs = data.frame(shape = "rectangle", color = "black", 
                                       fillcolor = "orange"),show = "std")
rsquared(path_test, method = NULL)
AIC(path_test) #65.184
nrow(Mydata)
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


