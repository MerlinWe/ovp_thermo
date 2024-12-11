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
Mydata$IDn<- gsub('OVP_','',Mydata$ID)
Mydata$IDn<-as.numeric(Mydata$IDn)
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
table(Mydata$season)
##################################################################
#selection spring season
Mydata<-subset(Mydata, season =="Fall" & phase=="day" )
nrow(Mydata)
#356
n_distinct(Mydata$ID)  #N= 12
n_distinct(Mydata$date)  #N= 117
n_distinct(Mydata$Year)  #N=3
table(Mydata$date)
table(Mydata$Year,Mydata$daynr)
table(Mydata$Year) 
table(Mydata$Year,Mydata$IDn) 
#data of Fall 2018 van 21/11-30/11 (n=41)  
#2019 van 1/9-30/11 (n=299) en 2020 1/9-6/10 (n=16)
#Fall 2018 daynr 1-10; 2019 daynr 285-375;  2020 daynr 651-686
#Ik neem voor de anlayse allen 2019 mee (andere jaren niet volledig gesampled)
Mydata<-subset(Mydata, Year == "2019")
#Dier 8 wordt niet meegnomen (maar 1 waarneming)
Mydata<-subset(Mydata, IDn != 8)
nrow(Mydata)
Mydata$day_season<- Mydata$daynr-284
Mydata$obs = factor(1:nrow(Mydata))

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

table(Mydata$IDn) # 
#remove  nothing
nrow(Mydata) #298

#B OUTLIERS in predictors?
#1 cleveland dotplots for detecting outliers (response variable+ continuous covariates)

#response variable
#BT
layout(1)
ggplot(Mydata,aes(x=BT,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$BT,Mydata$obs) #1 outliers case 146
hist(Mydata$BT) #
skewness(Mydata$BT)#-0.2061893
kurtosis(Mydata$BT)#4.100146 peaked
#HR
ggplot(Mydata,aes(x=HR,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$HR,Mydata$obs) #0 outliers

hist(Mydata$HR) 
skewness(Mydata$HR)#0.0241975
kurtosis(Mydata$HR)#2.684683 peaked

#ACT
ggplot(Mydata,aes(x=ACT,y=obs,label=obs,col=ID))+geom_point()+geom_text(hjust=0, vjust=0)
plot(Mydata$ACT,Mydata$obs) #no outliers
hist(Mydata$ACT) #normal distibuted =OK
skewness(Mydata$ACT)#0.4971346
kurtosis(Mydata$ACT)#4.426463 peaked

#continuous predictors
#TA
ggplot(Mydata,aes(x=TA,y=obs,label=obs,col=ID))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$TA,Mydata$obs) #0 outliers 
hist(Mydata$TA)

#WCF
ggplot(Mydata,aes(x=WCF,y=obs,label=obs))+geom_point()+geom_text(hjust=0, vjust=0) #no outliers
plot(Mydata$WCF,Mydata$obs) #0outliers
hist(Mydata$WCF)

table(Mydata$IDn)
plot(Mydata$day_season,Mydata$IDn,col=Mydata$IDn,pch=16) #no outliers
#wel duidelijk gaps te zien bij de dieren

#C COLLINEARITY between predictors

#1 GVIF per predictor

mod<-lmer(BT~TA+weight+WCF+ACT+day_season+(1|ID), data = Mydata, REML=T)
summary(mod)
vif(mod)
#high vif TA and WCF (as expected)
cor(Mydata$WCF,Mydata$TA)
# 0.9927069
cor(Mydata$ACT,Mydata$TA)
#0.2207239
cor(Mydata$ACT,Mydata$WCF)
#0.2262178

#remove WCF
mod<-lmer(BT~TA+weight+ACT+day_season+(1|ID), data = Mydata, REML=T)
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
gamm1<-gamm(BT~s(TA),random= list(ID =~1),data=Mydata)
summary(gamm1$gam)  #edf 1 linear relation
plot(gamm1$gam)


#TA-HR
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=loess) #dal parabool
ggplot(Mydata, aes(y=HR,x=TA))+geom_point()+geom_smooth(method=gam)
gamm2<-gamm(HR~s(TA),random= list(ID =~1),data=Mydata)
summary(gamm2$gam)  #edf=5.076  non linear relation
plot(gamm2$gam)
gam.check(gamm2$gam)
#TA-ACT
ggplot(Mydata, aes(y=ACT,x=TA))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=TA))+geom_point()+geom_smooth(method=gam)
gamm3<-gamm(ACT~s(TA),random= list(ID =~1),data=Mydata)
summary(gamm3$gam) #non linear relation edf=1.471
plot(gamm3$gam)
#linear relation

#Conclusion: linear relation TA and BT; 
# quadratic TA with HR and ACT?
#
###################################################
#day_season
#day_season-BT
cor(Mydata$day_season,Mydata$BT)
#[1] -0.01674475
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm4<-gamm(BT~s(day_season),random= list(ID =~1),data=Mydata)
summary(gamm4$gam)  #edf=1
plot(gamm4$gam)
# linear relation

#day_season-HR

ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm5<-gamm(HR~s(day_season),random= list(ID =~1),data=Mydata)
summary(gamm5$gam)  #edf=1
plot(gamm5$gam)
#linear relation

#day_season-ACT
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=day_season))+geom_point()+geom_smooth(method=gam)
gamm6<-gamm(ACT~s(day_season),random= list(ID =~1),data=Mydata)
summary(gamm6$gam)  #edf=1
plot(gamm6$gam)
#linear relation

#Conclusion: linear relation day_season and BT, ACT, HR; 


################################################################################
#predictor ACT for BT and HR response variable

#ACT-BT

ggplot(Mydata, aes(y=BT,x=ACT))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=ACT))+geom_point()+geom_smooth(method=gam)
gamm7<-gamm(BT~s(ACT),random= list(ID =~1),data=Mydata)
summary(gamm7$gam)  #edf=1.973!!
plot(gamm7$gam)
#non linear relation

#ACT-HR
ggplot(Mydata, aes(y=HR,x=ACT))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=ACT))+geom_point()+geom_smooth(method=gam)
gamm8<-gamm(HR~s(ACT),random= list(ID =~1),data=Mydata)
summary(gamm8$gam)  #edf=4.525!!
plot(gamm8$gam)
gam.check(gamm8$gam)
#non linear relation
##################################


####################################
#weight
#weight->BT
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=BT,x=weight))+geom_point()+geom_smooth(method=gam)
gamm9<-gamm(BT~s(weight, k=5),random= list(ID =~1),data=Mydata)
summary(gamm9$gam)  #edf=1!!
plot(gamm9$gam)
#linear relation 

#weight->HR
ggplot(Mydata, aes(y=HR,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=HR,x=weight))+geom_point()+geom_smooth(method=gam)
gamm10<-gamm(HR~s(weight,k=5),random= list(ID =~1),data=Mydata)
summary(gamm10$gam)  #edf=1!!
plot(gamm10$gam)
#linear relation

#weight->ACT
ggplot(Mydata, aes(y=ACT,x=weight))+geom_point()+geom_smooth(method=loess)
ggplot(Mydata, aes(y=ACT,x=weight))+geom_point()+geom_smooth(method=gam)
gamm11<-gamm(ACT~s(weight,k=5),random= list(ID =~1),data=Mydata)
summary(gamm11$gam)  #edf=1
plot(gamm11$gam)
#linear


#################################################################
#E Spatial/temporal/repeated measures aspects of sampling design
table(Mydata$ID)
table(Mydata$IDn)

#controle op autocorrelatie
#ACF en PacF plots voor elke tijdreeks
#ID 8 weggelaten (was maar 1 waarneming)
#BT
getal1=''
for(getal1 in (c(1,2,4,5,9,10,15,18,19))){
      test <- subset(Mydata, IDn == getal1)  
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$BT, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1))))
      pacf(test$BT,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1))))}
    layout(1)
}    

#sterke autocorrelatie AR1 

#HR
getal1=''
for(getal1 in (c(1,2,4,5,9,10,15,18,19))){
    test <- subset(Mydata, IDn == getal1) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$HR, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1))))
      pacf(test$HR,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1))))}
    layout(1)
}   
#ACT
getal1=''
for(getal1 in (c(1,2,4,5,9,10,15,18,19))){
    test <- subset(Mydata, IDn == getal1) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$ACT, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1))))
      pacf(test$ACT,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1))))}
    layout(1)
}   

#er isautocorrelatie AR1
#Echter er zitten gaps in de waarnemingen en daar houden deze plots geen rekening mee:
layout(1)
plot(Mydata$day_season,Mydata$obs,col=Mydata$ID,pch=16) #no outliers
#wel duidelijk gaps te zien

#Duidelijke gaps per combinatie van ID en phase
plot(Mydata$day_season,Mydata$ID,col=Mydata$IDn,pch=16)
n_distinct(Mydata$ID)
table((Mydata$ID))
table(as.numeric(Mydata$ID))
##################################
#random factor ID
#BT
Mydata$IDn<-as.numeric(Mydata$ID)
boxplot(BT~IDn,data=Mydata)
vioplot(BT~IDn,data=Mydata) 
ggplot(Mydata, aes(ID, BT, fill=factor(ID))) +
  geom_boxplot()+ coord_flip()

# differences in BT, between years and between antmals within a year
#HR
boxplot(HR~IDn,data=Mydata)
vioplot(HR~IDn,data=Mydata)
#large differences in HR between the animals 
ggplot(Mydata, aes(ID, HR, fill=factor(ID))) +
  geom_boxplot()+ coord_flip()
table(Mydata$IDn)

#ACT
boxplot(ACT~IDn,data=Mydata)
vioplot(ACT~IDn,data=Mydata)
ggplot(Mydata, aes(ID, ACT, fill=factor(ID))) +
  geom_boxplot()+ coord_flip()
table(Mydata$IDn)
#large differences in ACT between the animals in both periods
########################################################################
# F Interactions (is the quality of the data good enough to include them? 

#Interactions

#H. Are categorical covariates balanced?
table(Mydata$IDn)  #no
###################################################################################
###################################################################################
###################################################################################
#Descriptives;
min(Mydata$TA)  #3.81
max(Mydata$TA)  #30
mean(Mydata$TA) #16.58972
sd(Mydata$TA) #5.050667
cvTA <- sd(Mydata$TA) / mean(Mydata$TA) * 100
cvTA # 30.44456
quantile(Mydata$TA,c(.025,1-.025))
#     2.5%     97.5% 
# 6.610625 26.874917 
#############################
min(Mydata$ACT) #0
max(Mydata$ACT)  #21.8
mean(Mydata$ACT) #10.12333
sd(Mydata$ACT) #3.477113
quantile(Mydata$ACT,c(.025,1-.025))
#2.5%    97.5% 
#3.54250 18.60875 

min(Mydata$BT) # 38.39083
max(Mydata$BT)   #39.47
mean(Mydata$BT) #39.02995
sd(Mydata$BT) #0.1450928
quantile(Mydata$BT,c(.025,1-.025))
#2.5%    97.5% 
#38.75 39.31 

min(Mydata$HR)  #40.13509
max(Mydata$HR)   #85.01429
mean(Mydata$HR) #61.42318
sd(Mydata$HR) #9.89919
quantile(Mydata$HR,c(.025,1-.025))
#2.5%    97.5% 
#43.08011 81.37667  
###################################################################################
#modelselection


#HR model

#global model
modsel1HR <- lme(HR~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "REML",
            data = Mydata,control=lmeControl(opt='optim'),
            correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
            na.action=na.exclude)

anova(modsel1HR) 
summary(modsel1HR)
AIC(modsel1HR)#REML1653.119
BIC(modsel1HR)#1697.157


ggplot(Mydata,aes(x=ACT,y=HR,colour=ID))+  geom_point() +
  geom_point(aes(y = predict(modsel1HR)), col = "black") +
  geom_smooth(aes(y = predict(modsel1HR), colour = Mydata$ID), method = "lm") + facet_wrap(~Mydata$ID)
#NIET duidelijk te zien dat slope ACT verschlt per ID??

#does global model ModselHR1 improve by adding or deleting a random slope variable?
#adding random slope ACT
modsel1HR0 <- lme(HR~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+ACT|ID,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR0) #REML  1657.12 higher
BIC(modsel1HR0)  #REML 1708.498 higher
#global model does not improve if you add random slope ACT

#adding random slope TA
modsel1HR1 <- lme(HR~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR1) #REML  1657.101 higher
BIC(modsel1HR1)  #REML 1708.479 BIC higher 
#so I will choose the model without the random slope TA

#does top model improve when adding  random slope day_season
modsel1HR2 <- lme(HR~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID,method = "REML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HR2) #REML 
BIC(modsel1HR2)#REML  
#modle convergeert niet dus heb laat ik het weg

#does model improve when removing the nugget
modsel1HR3 <- lme(HR~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                  na.action=na.exclude)
AIC(modsel1HR3) #REML 1752.928 much higher
BIC(modsel1HR3)#REML 1793.297 much higher
#so I will choose the model with the nugget

#Modsel1HR  remains global model!!!

#Modsel1HR with ML estimator
modsel1HR <- lme(HR~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "ML",
                 data = Mydata,
                 correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                 na.action=na.exclude)
AIC(modsel1HR) #1684.719
BIC(modsel1HR) #1729.084
summary(modsel1HR)



modsel1HRc <- lme(HR~ACT+ I(ACT^2)+TA+I(TA^2)+weight+day_season+I(day_season^2), random=~1|ID,method = "ML",
                  data = Mydata,
                  correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                  na.action=na.exclude)
AIC(modsel1HRc) #1684.719
step_modelHR <- stepAIC(modsel1HRc, scope = list(lower = ~  TA+ACT, upper = ~ .),
                      direction = "both")
summary(step_modelHR)  #dit is het topmodel
AIC(step_modelHR)  #3770.124


modsel1HRtop<- lme(HR~ACT+I(ACT^2)+ TA, random=~1|ID,method = "ML",
                     data = Mydata,
                     correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                     na.action=na.exclude)
AIC(modsel1HRtop) #1680.337
summary(modsel1HRtop)
r.squaredGLMM(modsel1HRtop)

#or
modsel1HRtop1<- lme(HR~poly(ACT,2)+ TA, random=~1|ID,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRtop1) #1680.337
summary(modsel1HRtop1)
r.squaredGLMM(modsel1HRtop1)

gammHRtop<-gamm(HR~s(TA)+s(ACT),random= list(ID =~1),data=Mydata, method="ML")
summary(gammHRtop$gam)  #edf=5.076  non linear relation
plot(gammHRtop$gam)

gam.check(gammHRtop$gam)
AIC(gammHRtop)#2176.219

gammHRtop1<-gamm(HR~s(TA)+ACT,random= list(ID =~1),data=Mydata, method="ML")
summary(gammHRtop1$gam)  #edf=5.076  non linear relation
plot(gammHRtop1$gam)
gam.check(gammHRtop1$gam)
AIC(gammHRtop1)#2175.823

gammHRtop2<-gamm(HR~TA+ACT,random= list(ID =~1),data=Mydata, method="ML")
summary(gammHRtop2$gam)  #edf=5.076  non linear relation

AIC(gammHRtop2) #2188.507

#########################################
library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1HRtop1,digits = 2)


library(rptR)
#install.packages(("rptR"))
library(devtools)

r.squaredGLMM(modsel1HRtop)
#R2m       R2c
#[1,] 0.4081319 0.6154324
plot(allEffects(modsel1HRtop1))


###########################################################

modsel1BT <- lme(BT~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2),
                  random=~1|ID,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID),
                  na.action=na.exclude)
AIC(modsel1BT) # -334.5266
BIC(modsel1BT) # -294.1579
summary(modsel1BT)


#Does global model modsel 1BT improve buy adding or deleting a random slope variable?
#global model modsel1BT with random slope ACT

modsel1BT0<- lme(BT~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season+ACT|ID,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID),
                    na.action=na.exclude)
AIC(modsel1BT0) #REML -329.2435higher
BIC(modsel1BT0)  #REML-270.5254 higher
#global model does not improve by adding random slope ACT

#does global model improve when adding  random slope TA
modsel1BT1<- lme(BT~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID),
                    na.action=na.exclude)
AIC(modsel1BT1) #REML -334.9418   lower (<2 pts)
BIC(modsel1BT1)  #REML -287.2333   much higher
#global model does not improve by adding random slope TA

#does global model improve when adding  random slope day_season
modsel1BT2<- lme(BT~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID),
                    na.action=na.exclude)
AIC(modsel1BT2) #REML -331.9985 higher!
BIC(modsel1BT2)#REML -284.29 higher
#global model does not improve by adding random slope day_season

ggplot(Mydata,aes(x=day_season,y=BT,colour=ID))+  geom_point() +
  geom_point(aes(y = predict(modsel1BT2)), col = "black") +
  geom_smooth(aes(y = predict(modsel1BT2), colour = Mydata$ID), method = "lm") + facet_wrap(~Mydata$ID)
#no big differences


#does top model improve when adding the nugget
modsel1BT4<- lme(BT~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "REML",
                    data = Mydata,control=lmeControl(opt='optim'),
                    correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                    na.action=na.exclude)
AIC(modsel1BT4) #REML -333.4057 higher
BIC(modsel1BT4)  #REML-289.3671 higher

#modsel1BT remains gobal model!!!
#ML model
modsel1BT <- lme(BT~poly(ACT,2)+poly(TA,2)+weight+poly(day_season,2),
                 random=~1|ID,method = "ML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID),
                 na.action=na.exclude)
AIC(modsel1BT)  #ML -363.8026

modsel1BTa <- lme(BT~ACT+ I(ACT^2)+TA+I(TA^2)+weight+day_season+ I(day_season^2), random=~1|ID,method = "ML",
                 data = Mydata,control=lmeControl(opt='optim'),
                 correlation = corGaus(form= ~ day_season|ID),
                 na.action=na.exclude)
AIC(modsel1BTa) # -363.8026
summary(modsel1BTa)

step_modelBT <- stepAIC(modsel1BTa, scope = list(lower = ~ ACT + TA, upper = ~ .),
                        direction = "both")

summary(step_modelBT)
AIC(step_modelBT) #-785.6631
summary(step_modelBT)
#I will remove the quadratic term of ACT (delta AIC=0.2 <2 pts)

modsel1BTtop<- lme(BT~ACT+TA, random=~1|ID,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1BTtop) #-369.0502
summary(modsel1BTtop)
coef(summary(modsel1BTtop))


library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1BTtop,digits = 4)
icc(modsel1BTtop)
r.squaredGLMM(modsel1BTtop)

plot(allEffects(modsel1BTtop))
########################################################
#ACT model (without nugget)
modsel1ACTa <- lme(ACT~poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACTa) 
summary(modsel1ACTa)
AIC(modsel1ACTa) #  1526.676
BIC(modsel1ACTa)  #1559.767



#does global model modsel1 ACT improve buy adding or deleting a random slope variable?

#global model modsel1ACT with random slope TA
modsel1ACT0 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2), random=~1+TA|ID,method = "REML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT0) 
summary(modsel1ACT0)
AIC(modsel1ACT0) #1530.668  higher
BIC(modsel1ACT0) # 1571.113higher
#global model does not improve by adding random slope TA

#global model modsel1ACT with random slope day_season
modsel1ACT1 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2), random=~1+day_season|ID,method = "REML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                   na.action=na.exclude)
anova(modsel1ACT1) 
summary(modsel1ACT1)
AIC(modsel1ACT1) #1530.044  higher
BIC(modsel1ACT1) #1570.488 higher
##global model does not improve by adding random slope day_season

ggplot(Mydata,aes(x=day_season,y=ACT,colour=ID))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACT1)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACT1), colour = Mydata$ID), method = "lm") + facet_wrap(~Mydata$ID)
#no random slope for day_season

#Does nugget have effect? model with nugget
modsel1ACT2 <- lme(ACT~poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "REML",
                  data = Mydata,
                  correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                  na.action=na.exclude)
anova(modsel1ACT2) 
summary(modsel1ACT2)
AIC(modsel1ACT2)#1526.466 0.2 lowe
BIC(modsel1ACT2)#1563.234 higher
#I choose a model without nugget

#ML top global model
modsel1ACT <- lme(ACT~poly(TA,2)+weight+poly(day_season,2), random=~1|ID,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT) 
summary(modsel1ACT)
AIC(modsel1ACT) #1539.359

modsel1ACT2 <- lme(ACT~TA+I(TA^2)+weight+day_season+I(day_season^2), random=~1|ID,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                  na.action=na.exclude)
anova(modsel1ACT2) 
summary(modsel1ACT2)
AIC(modsel1ACT2) #1539.359

step_modelACT <- stepAIC(modsel1ACT2, scope = list(lower = ~ TA, upper = ~ .),
    direction = "both")
summary(step_modelACT)
anova(step_modelACT)
AIC(step_modelACT) #1538.235

modsel1ACTtop <- lme(ACT~poly(TA,2)+poly(day_season,2), random=~1|ID,method = "ML",
                  data = Mydata,control=lmeControl(opt='optim'),
                  correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                  na.action=na.exclude)
AIC(modsel1ACTtop) #1538.235
anova(modsel1ACTtop)
summary(modsel1ACTtop)


library(sjPlot) # table functions
library(sjmisc) # sample data
tab_model(modsel1ACTtop)
r.squaredGLMM(modsel1ACTtop)


ggplot(Mydata,aes(x=TA,y=ACT,colour=ID))+  geom_point() +
  geom_point(aes(y = predict(modsel1ACTtop)), col = "black") +
  geom_smooth(aes(y = predict(modsel1ACTtop), colour = Mydata$ID), method = "lm") + facet_wrap(~Mydata$ID)

plot(allEffects(modsel1ACTtop))

#########################################################################
#plots modellen
#HR

plot(allEffects(modsel1HRtop1, residuals=TRUE), 
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
skewness(Mydata$E1) #0.16353162
kurtosis(Mydata$E1)   #6.093241
#BTmodel
Mydata$E2 <- resid(modsel1BTtop,type="normalized")
hist(Mydata$E2)
skewness(Mydata$E2)# -0.2710051
kurtosis(Mydata$E2) #4.07476
#ACTmodel
Mydata$E3 <- resid(modsel1ACTtop,type="normalized")
hist(Mydata$E3)
skewness(Mydata$E3)# 0.211321
kurtosis(Mydata$E3) # 5.024607


#autocorrelation residuals model
#HR model
summary(modsel1HRtop)
anova(modsel1HRtop)
Mydata$E1 <- resid(modsel1HRtop,type="normalized")
getal1=''
Mydata$IDn
for(getal1 in (c(1,2,4,5,9,10,15,18,19))){
    test <- subset(Mydata, IDn == getal1) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=10)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1))))}
    layout(1)}
#no autocorrelation


#BT model
summary(modsel1BTtop)
anova(modsel1BTtop)
Mydata$E1 <- resid(modsel1BTtop,type="normalized")
getal1=''

for(getal1 in (c(1,2,4,5,9,10,15,18,19))){
    test <- subset(Mydata, IDn == getal1) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1))))}
    layout(1)}
#no autocorrelation

#ACT model
summary(modsel1ACTtop)
Mydata$E1 <- resid(modsel1ACTtop,type="normalized")
getal1=''
getal2=''
for(getal1 in (c(1,2,4,5,9,10,15,18,19))){
    test <- subset(Mydata, IDn == getal1) 
    if(nrow(test) > 0){
      layout(matrix(1:2, ncol = 2))
      acf(test$E1, main = "",lag=40)
      title(cex = 1.5, bquote(paste('ACF = ', .(getal1))))
      pacf(test$E1,main = "")
      title(cex = 1.5, bquote(paste('PacF = ', .(getal1))))}
    layout(1)}

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
#!
output <- data.frame(observed = Mydata$HR, fitted =  Mydata$F1)
ggplot(output, aes(fitted, observed)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_smooth(method="loess")+geom_abline(intercept=0, slope=1,col=2,lty=2)

cor(Mydata$F1, Mydata$HR)
#-0.07522144  #heel slecht
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
#0.5375262

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
#0.5143059
###############################################################################
library(piecewiseSEM)
#pathway analysis

modsel1HRtop<- lme(HR~ACT+I(ACT^2)+ TA, random=~1|ID,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRtop)#1680.337
Mydata$ACT2<-(Mydata$ACT)^2
modsel1HRsem<- lme(HR~ACT+ACT2+ TA, random=~1|ID,method = "ML",
                   data = Mydata,
                   correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                   na.action=na.exclude)
AIC(modsel1HRsem)#1680.337


modsel1BTtop<- lme(BT~ACT+TA, random=~1|ID,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID, nugget=TRUE),
                   na.action=na.exclude)

modsel1ACTtop <- lme(ACT~TA+I(TA^2)+day_season+I(day_season^2), random=~1|ID,method = "ML",
                   data = Mydata,control=lmeControl(opt='optim'),
                   correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                   na.action=na.exclude)
AIC(modsel1ACTtop) #1538.235
Mydata$TA2<-(Mydata$TA)^2
Mydata$day_season2<-(Mydata$day_season)^2
modsel1ACTsem<- lme(ACT~TA+TA2+day_season+day_season2, random=~1|ID,method = "ML",
                     data = Mydata,control=lmeControl(opt='optim'),
                     correlation = corGaus(form= ~ day_season|ID, nugget=FALSE),
                     na.action=na.exclude)
summary(modsel1ACTsem)
AIC(modsel1ACTsem)
path_test<-psem(modsel1HRsem,modsel1BTtop,modsel1ACTsem)
summary(path_test)
path_testa<-psem(modsel1HRsem,modsel1BTtop,HR%~~%BT,modsel1ACTsem)
summary(path_testa)

table(Mydata$IDn)
nrow(Mydata)
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
#o.k. 

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
#

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


