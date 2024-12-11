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
#selection spring season
Mydata<-subset(Mydata, season =="Spring" & phase=="day" )
nrow(Mydata)
#646
n_distinct(Mydata$ID)  #N= 12
n_distinct(Mydata$date)  #N= 182
n_distinct(Mydata$Year)  #N=2
table(Mydata$Year)  #data of spring 2019 en 2020

table(Mydata$Year,Mydata$daynr)

#2 spring periode 2019 daynr=101 t/m 192 and period 2 daynr = 467 t/ 558

#create variable serial number per season
Mydata$per_spring<-ifelse(Mydata$daynr<195, "spring2019", "spring2020")
table(Mydata$per_spring,Mydata$daynr)

table(Mydata$per_spring,Mydata$Year)
Mydata$day_season<-ifelse(Mydata$daynr<195, Mydata$daynr-101, Mydata$daynr-467)
Mydata$day_season2<-Mydata$day_season^2

table(Mydata$per_spring,Mydata$day_season)
table(Mydata$per_spring)

Mydata$obs = factor(1:nrow(Mydata))

Mydata$ID.per<-interaction(Mydata$ID,Mydata$per_spring)
table(Mydata$ID.per)
n_distinct(Mydata$ID.per)
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
table(Mydata$per_spring, Mydata$IDn) # 
#remove  nothing
nrow(Mydata) #646

#B OUTLIERS in predictors?
#1 cleveland dotplots for detecting outliers (response variable+ continuous covariates)

#response variable
#BT
layout(1)
ggplot(Mydata,aes(x=BT,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$BT,Mydata$obs) #4 outliers case 153 t/m 156
hist(Mydata$BT) #
skewness(Mydata$BT)#0.3991502  
kurtosis(Mydata$BT)#3.105796 peaked
#HR
ggplot(Mydata,aes(x=HR,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$HR,Mydata$obs) #1 outlier
hist(Mydata$HR) 
skewness(Mydata$HR)#-0.009156316
kurtosis(Mydata$HR)#2.45279 peaked

#ACT
ggplot(Mydata,aes(x=ACT,y=obs,label=obs,col=ID.per))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$ACT,Mydata$obs) #no outliers
hist(Mydata$ACT) #normal distibuted =OK
skewness(Mydata$ACT)# 0.3207268
kurtosis(Mydata$ACT)#3.509585 peaked

#continuous predictors
#TA
ggplot(Mydata,aes(x=TA,y=obs,label=obs,col=ID.per))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$TA,Mydata$obs) #0 outliers 
hist(Mydata$TA)

#WCF
ggplot(Mydata,aes(x=WCF,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$WCF,Mydata$obs) #4 outliers same as TA
hist(Mydata$WCF)
table(Mydata$ID.per,Mydata$ID)
plot(Mydata$day_season,Mydata$ID.per,col=Mydata$IDn,pch=16) #no outliers
#wel duidelijk gaps te zien bij de dieren

#C COLLINEARITY between predictors

#1 GVIF per predictor

mod<-lmer(BT~TA+weight+WCF+ACT+day_season+(1|ID.per), data = Mydata, REML=T)
summary(mod)
vif(mod)
#high vif TA and WCF (as expected)
cor(Mydata$WCF,Mydata$TA)
# 0.9941949
cor(Mydata$ACT,Mydata$TA)
#-0.03917299
cor(Mydata$ACT,Mydata$WCF)
#-0.02871649

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
summary(gamm1$gam)  #edf 1.908 nonlinear relation
plot(gamm1$gam)
#quadratic

#TA-HR
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=loess) #dal parabool
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=gam)
gamm2<-gamm(HR~s(TA),random= list(ID.per =~1),data=Mydata)
summary(gamm2$gam)  #edf=3.847  quadatic relation

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
#[1] 0.3251193
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm4<-gamm(BT~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm4$gam)  #edf=7.1412
plot(gamm4$gam)
#paar onverklaarbare piekjes linear relation

#day_season-HR

ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm5<-gamm(HR~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm5$gam)  #edf=4.392
plot(gamm5$gam)
#linear relation

#day_season-ACT
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm6<-gamm(ACT~s(day_season),random= list(ID.per =~1),data=Mydata)
summary(gamm6$gam)  #edf=1
plot(gamm6$gam)
#linear relation

#Conclusion: linear relation day_season and BT, ACT, HR; 


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
summary(gamm8$gam)  #edf=3.116!!
plot(gamm8$gam)

#quadratic relation
##################################


####################################
#weight
#weight->BT
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=gam)
gamm9<-gamm(BT~s(weight, k=5),random= list(ID.per =~1),data=Mydata)
summary(gamm9$gam)  #edf=1.457!!
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
#linear


#################################################################
#E Spatial/temporal/repeated measures aspects of sampling design
table(Mydata$ID)
table(Mydata$IDn)

#controle op autocorrelatie
#ACF en PacF plots voor elke tijdreeks
table(Mydata$per_spring)
#BT
getal1=''
getal2=''
for(getal1 in (c(1,2,4,5,8,9,10,15,16,18,19,20))){
  for(getal2 in c('spring2019', 'spring2020')){
    test <- subset(Mydata, IDn == getal1 & per_spring==getal2) 
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
for(getal1 in (c(1,2,4,5,8,9,10,15,16,18,19,20))){
  for(getal2 in c('spring2019', 'spring2020')){
    test <- subset(Mydata, IDn == getal1 & per_spring==getal2) 
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
for(getal1 in (c(1,2,4,5,8,9,10,15,16,18,19,20))){
  for(getal2 in c('win18-19', 'win19-20')){
    test <- subset(Mydata, IDn == getal1 & per_spring==getal2) 
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
#per_spring*act
#BT
boxplot(ACT~per_spring,data=Mydata) #overlap thus its ok to test for interaction
ggplot(Mydata, aes(x=ACT, y=BT, color=per_spring)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_spring BT



#HR

ggplot(Mydata, aes(x=ACT, y=HR, color=per_spring)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #small interaction (N.S)

mod7 <- lme(HR~per_spring*ACT, random=~1+ACT|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod7) #per_spring:ACT P= 0.3924 
summary(mod7)
AIC(mod7) #ML:3842.866
#no interation

ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(mod7)), col = "black") +
  geom_smooth(aes(y = predict(mod7), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)

#Interactions
#per_spring*TA
boxplot(TA~per_spring,data=Mydata) #overlap thus its ok to test for interaction

#BT

ggplot(Mydata, aes(x=TA, y=BT, color=per_spring)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_spring BT
#no interaction
mod8 <- lme(BT~per_spring*ACT, random=~1+ACT|ID.per,method = "ML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod8) #per_spring:ACT P= 0.5422
summary(mod8)
AIC(mod8) #

#HR

ggplot(Mydata, aes(x=TA, y=HR, color=per_spring)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_spring BT
#interaction TA*per-spring-> HR ?. NS!!

mod8a <- lme(HR~per_spring*TA, random=~1+TA|ID.per, method = "ML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod8a)  #per_spring:TA P= 0.0441
summary(mod8a)
AIC(mod8a)
ranef(mod8a)
#no interaction

ggplot(Mydata,aes(x=TA,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(mod8a)), col = "black") +
  geom_smooth(aes(y = predict(mod8a), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)



#ACT

ggplot(Mydata, aes(x=TA, y=ACT, color=per_spring)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) #no interaction  per_spring BT
#small interaction N.S.

mod9 <- lme(ACT~per_spring*TA, random=~1|ID.per,method = "ML",
            data = Mydata,
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)
anova(mod9) #per_spring:TA=  0.1745
summary(mod9)
AIC(mod9)


#H. Are categorical covariates balanced?
table(Mydata$IDn)  #no
table(Mydata$per_spring)#no
###################################################################################
###################################################################################
###################################################################################
#Descriptives;
min(Mydata$TA)  #2.81
max(Mydata$TA)  #29.62
mean(Mydata$TA) #16.42255
sd(Mydata$TA) #5.171882
cvTA <- sd(Mydata$TA) / mean(Mydata$TA) * 100
cvTA # 31.49256
quantile(Mydata$TA,c(.025,1-.025))
#     2.5%     97.5% 
# 7.386667 26.562333 
#############################
min(Mydata$ACT) #0.3
max(Mydata$ACT)  #23.9
mean(Mydata$ACT) #10.67907
sd(Mydata$ACT) #3.92034
quantile(Mydata$ACT,c(.025,1-.025))
#2.5%    97.5% 
# 3.333333 19.252083 

min(Mydata$BT) # 38.52
max(Mydata$BT)   #39.93333
mean(Mydata$BT) #39.0946
sd(Mydata$BT) #0.2334456
quantile(Mydata$BT,c(.025,1-.025))
#2.5%    97.5% 
#38.68912 39.59875

min(Mydata$HR)  #29.03846
max(Mydata$HR)   #87.55714
mean(Mydata$HR) #56.86895
sd(Mydata$HR) #11.16285
quantile(Mydata$HR,c(.025,1-.025))
#2.5%    97.5% 
#36.29375 77.59097 
###################################################################################
#modelselection


#HR model

Mydata$per_spring<-as.factor(Mydata$per_spring)
Mydata$per_spring = relevel(Mydata$per_spring, ref="spring2020")
#global model
modsel1HR <- lme(HR~per_spring+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
            na.action=na.exclude)

anova(modsel1HR) 
summary(modsel1HR)
AIC(modsel1HR)#REML3737.435
BIC(modsel1HR)#3795.373

ggplot(Mydata,aes(x=ACT,y=HR,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1HR)), col = "black") +
  geom_smooth(aes(y = predict(modsel1HR), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#NIET duidelijk te zien dat slope ACT verschlt per ID.per??

#does global model ModselHR1 improve by adding or deleting a random slope variable?
#adding random slope ACT
modsel1HR0 <- lme(HR~per_spring+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR0) #REML  3754.41 higher
BIC(modsel1HR0)  #REML 3821.262 higher
#global model does not improve if you add random slope ACT

#adding random slope TA
modsel1HR1 <- lme(HR~per_spring+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR1) #REML  3741.4351 higher
BIC(modsel1HR1)  #REML 3808.286 BIC higher 
#so I will choose the model without the random slope TA

#does top model improve when adding  random slope day_season
modsel1HR2 <- lme(HR~per_spring+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID.per,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR2) #REML 3741.429 higher
BIC(modsel1HR2)#REML  3808.281  higher
#so I will choose the model without the random slope day_season

#does model improve when removing the nugget
modsel1HR3 <- lme(HR~per_spring+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
AIC(modsel1HR3) #REML 3918.367 much higher
BIC(modsel1HR3)#REML 3971.848 much higher
#so I will choose the model with the nugget

#Modsel1HR  remains global model!!!

#Modsel1HR with ML estimator
modsel1HR <- lme(HR~per_spring+poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HR) #3772.603
summary(modsel1HR)

modsel1HRa <- lme(HR~per_spring+ACT+I(ACT^2)+TA +I(TA^2)+weight+day_season+I(day_season^2), random=~1|ID.per,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HRa) #3772.603
summary(modsel1HRa)


step_modelHR <- stepAIC(modsel1HRa, scope = list(lower = ~ ACT + TA, upper = ~ .),
                      direction = "both")
summary(step_modelHR)  #dit is het topmodel
AIC(step_modelHR)  #3770.124

#per_spring and I(day_season^2) results in a decrease of the AIC with less than 2 pts. 
#So I will remove these variables.

modsel1HRtop<- lme(HR~ACT+TA+day_season, random=~1|ID.per,method = "ML",
                     data = Mydata,
                     correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HRtop) #3770.903
summary(modsel1HRtop)
r.squaredGLMM(modsel1HRtop)

#########################################
library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1HRtop)


library(rptR)
#install.packages(("rptR"))
library(devtools)

r.squaredGLMM(modsel1HRtop)
#R2m       R2c
#[1,] 0.4081319 0.6154324
plot(allEffects(modsel1HRtop))


###########################################################

modsel1BT <- lme(BT~per_spring+ACT+TA+weight+poly(day_season,2),
                  random=~1|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
AIC(modsel1BT) #-737.4079
BIC(modsel1BT) #-688.3489
summary(modsel1BT)


#for random model with day_season+day_season2 does not give same results when using REML 
#It does with ML!!

#Does global model modsel 1BT improve buy adding or deleting a random slope variable?
#global model modsel1BT with random slope ACT

modsel1BT0<- lme(BT~per_spring+ACT+TA+weight+poly(day_season,2), random=~1+day_season+ACT|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT0) #REML-734.7412 higher
BIC(modsel1BT0)  #REML-676.7625 higher
#global model does not improve by adding random slope ACT

#does global model improve when adding  random slope TA
modsel1BT1<- lme(BT~per_spring+ACT+TA+weight+poly(day_season,2), random=~1+TA|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT1) #REML -733.4079AIC higher
BIC(modsel1BT1)  #REML -675.4291 BIC higher 
#so I will choose the model without the random slope TA

#does global model improve when adding  random slope day_season
modsel1BT2<- lme(BT~per_spring+ACT+TA+weight+poly(day_season,2), random=~1+day_season|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT2) #REML -742.0891  4.6 pts lower!
BIC(modsel1BT2)#REML -684.1103 higher
#so I will choose the model without the random slope day_season (because BIC is hiher)

ggplot(Mydata,aes(x=day_season,y=BT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1BT2)), col = "black") +
  geom_smooth(aes(y = predict(modsel1BT2), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#no big differences


#does top model improve when removing the nugget
modsel1BT4<- lme(BT~per_spring+ACT+TA+weight+poly(day_season,2), random=~1|ID.per,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID.per),
                    na.action=na.exclude)
AIC(modsel1BT4) #REML -635.15969 much higher
BIC(modsel1BT4)  #REML-590.5606 much higher

#modsel1BT remains gobal model!!!
#ML model
modsel1BT <- lme(BT~per_spring+ACT+TA+weight+poly(day_season,2), random=~1|ID.per,method = "ML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1BT) # -783.0847
summary(modsel1BT)

#does model improve if you remove the quadratic term of day_season?
modsel1BTa <- lme(BT~per_spring+ACT+TA+weight+day_season, random=~1|ID.per,method = "ML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1BTa) #-783.7516
summary(modsel1BT)
#YES!!

step_modelBT <- stepAIC(modsel1BTa, scope = list(lower = ~ ACT + TA, upper = ~ .),
                        direction = "both")

summary(step_modelBT)
AIC(step_modelBT) #-785.6631
summary(step_modelBT)
#does model improve if you only use the linear term for day_season(quadratic term is N.S.) 
#(yes!!, because)

modsel1BTtop<- lme(BT~ACT+TA+day_season, random=~1|ID.per,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1BTtop) #-785.6631
summary(modsel1BTtop)
coef(summary(modsel1BTtop))


library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1BTtop,digits = 4)
icc(modsel1BTtop)

plot(allEffects(modsel1BTtop))
########################################################
#ACT model (without nugget)
modsel1ACTa <- lme(ACT~TA+weight+poly(day_season,2)+per_spring, random=~1|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACTa) 
summary(modsel1ACTa)
AIC(modsel1ACTa) # 3526.168
BIC(modsel1ACTa)  #3566.321



#does global model modsel1ACT improve buy adding or deleting a random slope variable?

#global model modsel1ACT with random slope TA
modsel1ACT0 <- lme(ACT~TA+weight+poly(day_season,2)+per_spring, random=~1++TA|ID.per,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT0) 
summary(modsel1ACT0)
AIC(modsel1ACT0) #3528.116  higher
BIC(modsel1ACT0) #3577.192 higher



#global model modsel1ACT with random slope day_season
modsel1ACT1 <- lme(ACT~TA+weight+poly(day_season,2)+per_spring, random=~1+day_season|ID.per,method = "REML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                   na.action=na.exclude)
anova(modsel1ACT1) 
summary(modsel1ACT1)
AIC(modsel1ACT1) #3525.657  lower 0.5  pts (was 3526.168 )
BIC(modsel1ACT1) #3574.734 higher
#no random slope for day_season

ggplot(Mydata,aes(x=day_season,y=ACT,colour=ID.per))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACT1)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACT1), colour = Mydata$ID.per), method = "lm") + facet_wrap(~Mydata$ID.per)
#no random slope for day_season

#Does nugget have effect? model with nugget
modsel1ACT2 <- lme(ACT~TA+weight+poly(day_season,2)+per_spring, random=~1|ID.per,method = "REML",
                  data = Mydata,
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1ACT2) 
summary(modsel1ACT2)
AIC(modsel1ACT2)#[1] 3524.569 higher
BIC(modsel1ACT2)#3569.184 higher
#I choose a model with nugget

#ML top global model
modsel1ACT <- lme(ACT~TA+weight+poly(day_season,2)+per_spring, random=~1|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT) 
summary(modsel1ACT)
AIC(modsel1ACT) #3524.276

#removing quadratic effect day_season improves the model?
modsel1ACT2 <- lme(ACT~TA+weight+day_season+per_spring, random=~1|ID.per,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                   na.action=na.exclude)

anova(modsel1ACT2) 
summary(modsel1ACT2)
AIC(modsel1ACT2) #3522.389  lower
BIC(modsel1ACT2)  #3558.155  lower
#YES

step_modelACT <- stepAIC(modsel1ACT2, scope = list(lower = ~ TA, upper = ~ .),
    direction = "both")
summary(step_modelACT)
AIC(step_modelACT) #2077.392
#removing day_season (N.S.) decreases the AIC with <2pts, thus it is removed from the topmodel
modsel1ACTtop <- lme(ACT~TA, random=~1|ID.per,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID.per, nugget=FALSE),
                  na.action=na.exclude)
AIC(modsel1ACTtop) #3521.786
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
skewness(Mydata$E1) #-0.4743562
kurtosis(Mydata$E1)  #7.144179
#BTmodel
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)#-0.01384719
kurtosis(Mydata$E1) #4.126723
#ACTmodel
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
hist(Mydata$E1)
skewness(Mydata$E1)# 0.471474
kurtosis(Mydata$E1) # 3.963528


#autocorrelation residuals model
#HR model
summary(modsel1HRtop)
anova(modsel1HRtop)
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(1,2,4,5,8,9,10,15,16,18,19,20))){
  for(getal2 in c('spring2019', 'spring2020')){
    test <- subset(Mydata, IDn == getal1 & per_spring==getal2) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1),.(getal2))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1),.(getal2))))}
    layout(1)}
}    
#no autocorrelation

Mydata$per_spring
#BT model
summary(modsel1BTtop)
anova(modsel1BTtop)
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(1,2,4,5,8,9,10,15,16,18,19,20))){
  for(getal2 in c('spring2019', 'spring2020')){
    test <- subset(Mydata, IDn == getal1 & per_spring==getal2) 
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
for(getal1 in (c(1,2,4,5,8,9,10,15,16,18,19,20))){
  for(getal2 in c('spring2019', 'spring2020')){
    test <- subset(Mydata, IDn == getal1 & per_spring==getal2) 
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
#0.8444339  #prima
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
#0.7836466

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
#0.436975
###############################################################################
library(piecewiseSEM)
#pathway analysis

path_test<-psem(modsel1HRtop,modsel1BTtop,modsel1ACTtop)
summary(path_test)
path_testa<-psem(modsel1HRtop,modsel1BTtop,HR%~~%BT,modsel1ACTtop)
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


