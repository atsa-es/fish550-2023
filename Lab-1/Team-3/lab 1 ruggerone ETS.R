ruggerone_data <- readRDS(here::here("Lab-1", "Data_Images", "ruggerone_data.rds"))

library(tidyverse)
library(forecast)


#create forecast model for sockeye salmon across all regions

#Is the data stationary?
sockeye.dat<-subset(ruggerone_data, species=='sockeye')

sockeye.dat %>%
  group_by(year) %>%
  summarize(total = sum(returns, na.rm=T)) %>%
ggplot(aes(x=year, y=log(total))) +
  geom_line() +
  ylab('Log (Returns)') +
  xlab('Year')

total.sockeye<-sockeye.dat %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

sockeye.ts<-ts(total.sockeye$lntotal, start=total.sockeye$year[1])
forecast::ndiffs(sockeye.ts, test='adf')
#1
forecast::ndiffs(sockeye.ts, test='kpss')
#1
#requires d=1
plot(diff(sockeye.ts))

train.sockeye<-window(sockeye.ts, start=1952, end=2010)
test.sockeye<-window(sockeye.ts, start=2011, end=2015)
fit <- forecast::auto.arima(train.sockeye, trace=T)
# Best model: ARIMA(0,1,2)  
#AIC = 1.111738
fit.final<-forecast::auto.arima(train.sockeye, approximation = F, stepwise = F)
ARIMA(0,1,2) 

Coefficients:
  ma1      ma2
-0.3538  -0.2708
s.e.   0.1248   0.1244

sigma^2 = 0.05491:  log likelihood = 2.67
AIC=0.67   AICc=1.11   BIC=6.85

#check ACF and PACF
acf(train.sockeye)
#significant lags through 6, then 9-11
pacf(train.sockeye)
#significant lag at 9

residuals(fit.final)
plot(residuals(fit.final))
#these look like white noise
acf(residuals(fit.final))
#sig correlation at 7
forecast::checkresiduals(fit.final)
forecast::checkresiduals(fit.final, plot=F)
#p = 0.01 so null is rejected and model is not a good fit
#see Eli's notes/code on adding seasonal component since there is still strong correlation at lag of 7 years

forecast::accuracy(fit.final)
accuracy(forecast(fit.final, h=5), test.sockeye)

                     ME      RMSE       MAE        MPE     MAPE      MASE
Training set  0.031180976 0.2283027 0.1682948  0.5082170 4.058891 0.9116935
Test set     -0.009636128 0.1742628 0.1711543 -0.3561161 3.828493 0.9271840
ACF1 Theils U
Training set -0.00179341        NA
Test set      0.30234708 0.9239656

#plot forecast all fish
fit.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.sockeye))

#Model forecast for sockeye returns by region
unique(sockeye.dat$region)
#"ci"    "e_kam" "japan" "kod"   "korea" "m_i"   "nbc"   "pws"   "sbc"  "seak"  "s_pen" "wak"   "wc"    "w_kam"
sockeye.dat %>%
  group_by(year, region) %>%
  summarize(total = sum(returns, na.rm=T)) %>%
  ggplot(aes(x=year, y=log(total))) +
  geom_line() +
  ylab('Log (Returns)') +
  xlab('Year') +
  facet_wrap(~region)

#Is the data stationary?
#These regions look stationary: m_i, nbc, pws, sbc, seak, wc
#ci
sock.ci<-subset(sockeye.dat, region=='ci')
total.ci<-sock.ci %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

ci.ts<-ts(total.ci$lntotal, start=total.ci$year[1])
forecast::ndiffs(ci.ts, test='adf')
#1
forecast::ndiffs(ci.ts, test='kpss')
#1
#requires d=1
plot(diff(ci.ts))

train.ci<-window(ci.ts, start=1952, end=2010)
test.ci<-window(ci.ts, start=2011, end=2015)
fit.ci <- forecast::auto.arima(train.ci, trace=T)
# Best model: ARIMA(1,1,1) 
#AIC = 55.19379
ci.final<-forecast::auto.arima(train.ci, approximation = F, stepwise = F)

ARIMA(4,1,1) 

Coefficients:
  ar1      ar2      ar3      ar4     ma1
-1.0117  -0.4855  -0.4679  -0.4269  0.7036
s.e.   0.1820   0.1697   0.1631   0.1164  0.1797

sigma^2 = 0.1264:  log likelihood = -20.22
AIC=52.44   AICc=54.09   BIC=64.81

#check ACF and PACF
acf(train.ci)
#significant positive lags through 11
pacf(train.ci)
#significant negative lag at 6

residuals(ci.final)
forecast::checkresiduals(ci.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(ci.final, plot=F)
#p = 0.4443 so null is not rejected and model is a good fit

forecast::accuracy(ci.final)
accuracy(forecast(ci.final, h=5), test.ci)

                   ME      RMSE       MAE       MPE     MAPE      MASE
Training set 0.02471079 0.3370305 0.2702999 -38.06705 62.87923 0.8967555
Test set     0.34919358 0.4082214 0.3491936  17.84353 17.84353 1.1584957
ACF1 Theils U
Training set  0.008040561        NA
Test set     -0.080588177  1.829245

#plot forecast
ci.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci))

##ekam
sock.ekam<-subset(sockeye.dat, region=='e_kam')
total.ekam<-sock.ekam %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

ekam.ts<-ts(total.ekam$lntotal, start=total.ekam$year[1])
forecast::ndiffs(ekam.ts, test='adf')
#0
forecast::ndiffs(ekam.ts, test='kpss')
#0

train.ekam<-window(ekam.ts, start=1952, end=2010)
test.ekam<-window(ekam.ts, start=2011, end=2015)
fit.ekam <- forecast::auto.arima(train.ekam, trace=T)
# Best model: ARIMA(1,0,0) with non-zero mean
#AIC = 91.28832
ekam.final<-forecast::auto.arima(train.ekam, approximation = F, stepwise = F)

ARIMA(1,0,0) with non-zero mean 

Coefficients:
  ar1    mean
0.7744  1.0348
s.e.  0.0948  0.2721

sigma^2 = 0.2514:  log likelihood = -42.43
AIC=90.85   AICc=91.29   BIC=97.08

#check ACF and PACF
acf(train.ekam)
#significant positive lags through 2
pacf(train.ekam)
#no sig lags

residuals(ekam.final)
forecast::checkresiduals(ekam.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(ekam.final, plot=F)
#p = 0.5239 so null is not rejected and model is a good fit

forecast::accuracy(ekam.final)
accuracy(forecast(ekam.final, h=5), test.ekam)

                 ME      RMSE       MAE        MPE      MAPE     MASE
Training set 0.04091151 0.4928187 0.3428350 -238.67853 281.31610 1.020966
Test set     0.46145661 0.5037707 0.4614566   24.49621  24.49621 1.374223
ACF1 Theils U
Training set 0.0188609        NA
Test set     0.1718058  3.439219

#plot forecast
ekam.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ekam))

#m_i
sock.mi<-subset(sockeye.dat, region=='m_i')
total.mi<-sock.mi %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

mi.ts<-ts(total.mi$lntotal, start=total.mi$year[1])
forecast::ndiffs(mi.ts, test='adf')
#0
forecast::ndiffs(mi.ts, test='kpss')
#0

train.mi<-window(mi.ts, start=1952, end=2010)
test.mi<-window(mi.ts, start=2011, end=2015)
fit.mi <- forecast::auto.arima(train.mi, trace=T)
#  Best model: ARIMA(1,0,0) with non-zero mean 
#AIC = 74.21702
mi.final<-forecast::auto.arima(train.mi, approximation = F, stepwise = F)

ARIMA(1,0,0) with non-zero mean 

Coefficients:
  ar1     mean
0.4790  -1.2783
s.e.  0.1145   0.1055

sigma^2 = 0.1903:  log likelihood = -33.89
AIC=73.78   AICc=74.22   BIC=80.01

#check ACF and PACF
acf(train.mi)
#significant positive lag at 1
pacf(train.mi)
#no sig lags

residuals(mi.final)
forecast::checkresiduals(mi.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(mi.final, plot=F)
#p = 0.8835 so null is not rejected and model is a good fit

forecast::accuracy(mi.final)
accuracy(forecast(mi.final, h=5), test.mi)

                  ME      RMSE       MAE        MPE     MAPE      MASE
Training set 0.004688289 0.4288126 0.3196422  -0.438757 44.06885 0.9489701
Test set     0.529338744 0.5445771 0.5293387 -89.086879 89.08688 1.5715281
ACF1 Theils U
Training set  0.003297423        NA
Test set     -0.518165481  1.711458

#plot forecast
mi.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.mi))

#wkam
sock.wkam<-subset(sockeye.dat, region=='w_kam')
total.wkam<-sock.wkam %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

wkam.ts<-ts(total.wkam$lntotal, start=total.wkam$year[1])
forecast::ndiffs(wkam.ts, test='adf')
#1
forecast::ndiffs(wkam.ts, test='kpss')
#1
#requires d=1
plot(diff(wkam.ts))

train.wkam<-window(wkam.ts, start=1952, end=2010)
test.wkam<-window(wkam.ts, start=2011, end=2015)
fit.wkam <- forecast::auto.arima(train.wkam, trace=T)
# Best model: ARIMA(1,1,0) 
#AIC = 75.77518
wkam.final<-forecast::auto.arima(train.wkam, approximation = F, stepwise = F)

ARIMA(3,1,2) 

Coefficients:
  ar1      ar2      ar3     ma1     ma2
-0.5191  -0.7677  -0.4684  0.2400  0.8618
s.e.   0.1359   0.1434   0.1170  0.0977  0.1427

sigma^2 = 0.1787:  log likelihood = -30.46
AIC=72.92   AICc=74.57   BIC=85.28

#check ACF and PACF
acf(train.wkam)
#significant positive lags through 8
pacf(train.wkam)
#significant positive lag at 2 and negative lag at 6

residuals(wkam.final)
forecast::checkresiduals(wkam.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(wkam.final, plot=F)
#p = 0.6775 so null is not rejected and model is a good fit

forecast::accuracy(wkam.final)
accuracy(forecast(wkam.final, h=5), test.wkam)

                    ME      RMSE       MAE      MPE     MAPE      MASE
Training set 0.01512549 0.4006757 0.3122417 27.42798 76.59597 0.8458469
Test set     0.33419927 0.3528808 0.3341993 13.05811 13.05811 0.9053289
ACF1 Theils U
Training set -0.04056413        NA
Test set     -0.39439701  1.888705

#plot forecast
wkam.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wkam))

  #geom_line(aes(x=total.wkam$year[60:64], y=total.wkam$lntotal[60:64]), colour='red')

#wak
sock.wak<-subset(sockeye.dat, region=='wak')
total.wak<-sock.wak %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

wak.ts<-ts(total.wak$lntotal, start=total.wak$year[1])
forecast::ndiffs(wak.ts, test='adf')
#0
forecast::ndiffs(wak.ts, test='kpss')
#1
plot(wak.ts)
plot(diff(wak.ts))
#visually differencing seems like a good idea (at least once)
#1950-1980 seems to have large seasonal effect that then disappears

train.wak<-window(wak.ts, start=1952, end=2010)
test.wak<-window(wak.ts, start=2011, end=2015)
fit.wak <- forecast::auto.arima(train.wak, trace=T)
# Best model: ARIMA(2,1,2)
#AIC = 65.61292
wak.final<-forecast::auto.arima(train.wak, approximation = F, stepwise = F)

ARIMA(2,1,2) 

Coefficients:
  ar1      ar2      ma1     ma2
0.5083  -0.9240  -0.6982  0.6790
s.e.  0.0576   0.0554   0.1363  0.1231

sigma^2 = 0.1554:  log likelihood = -27.23
AIC=64.46   AICc=65.61   BIC=74.76

#check ACF and PACF
acf(train.wak)
#significant positive lags: 1, 4, 5, 9, 10, 14, 15
pacf(train.wak)
#significant lags: 1(pos), 2( neg), 3 (pos), 4(pos), 6 (neg)

residuals(wak.final)
forecast::checkresiduals(wak.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(wak.final, plot=F)
#p = 0.2666 so null is not rejected and model is a good fit

forecast::accuracy(wak.final)
accuracy(forecast(wak.final, h=5), test.wak)

accuracy(forecast(wak.final, h=5), test.wak)
                    ME      RMSE       MAE        MPE      MAPE      MASE
Training set  0.012449325 0.3771240 0.3029073 -0.8791146 10.093141 0.6867955
Test set     -0.005148165 0.2581967 0.2177408 -0.7186798  5.944272 0.4936936
ACF1 Theils U
Training set -0.05606583        NA
Test set      0.15124484  0.807023

#plot forecast
wak.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wak))

#what if training set only contained periodic data?
train.wak<-window(wak.ts, start=1952, end=1979)
test.wak<-window(wak.ts, start=1980, end=2015)
fit.wak <- forecast::auto.arima(train.wak, trace=T)
# Best model: ARIMA(0,0,1) with non-zero mean 
#AIC = 50.83221
wak.final<-forecast::auto.arima(train.wak, approximation = F, stepwise = F)

ARIMA(2,0,0) with non-zero mean 

Coefficients:
  ar1      ar2    mean
0.3887  -0.6754  2.8226
s.e.  0.1460   0.1356  0.0653

sigma^2 = 0.2081:  log likelihood = -16.81
AIC=41.61   AICc=43.35   BIC=46.94

#check ACF and PACF
acf(train.wak)
#significant positive lags: 5, 9
#neg: 2, 3, 7
pacf(train.wak)
#significant lags: 2 (neg)

residuals(wak.final)
forecast::checkresiduals(wak.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(wak.final, plot=F)
#p = 0.2629 so null is not rejected and model is a good fit

forecast::accuracy(wak.final)
accuracy(forecast(wak.final, h=36), test.wak)

                     ME      RMSE       MAE       MPE     MAPE      MASE
Training set -0.001068959 0.4310956 0.3592882 -2.754584 13.42142 0.5449189
Test set      0.880912148 0.9401431 0.8809121 23.327872 23.32787 1.3360463
ACF1 Theils U
Training set -0.04541685        NA
Test set      0.47679133  2.657052

#plot forecast
wak.final %>%
  forecast(h=36) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wak))

#s_pen
sock.spen<-subset(sockeye.dat, region=='s_pen')
total.spen<-sock.spen %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

spen.ts<-ts(total.spen$lntotal, start=total.spen$year[1])
forecast::ndiffs(spen.ts, test='adf')
#1
forecast::ndiffs(spen.ts, test='kpss')
#1
#requires d=1
plot(diff(spen.ts))

train.spen<-window(spen.ts, start=1952, end=2010)
test.spen<-window(spen.ts, start=2011, end=2015)
fit.spen <- forecast::auto.arima(train.spen, trace=T)
# Best model: ARIMA(2,1,2)  
#AIC = 44.73109
spen.final<-forecast::auto.arima(train.spen, approximation = F, stepwise = F)

ARIMA(0,1,3) 

Coefficients:
  ma1      ma2     ma3
-0.5261  -0.4831  0.3698
s.e.   0.1376   0.1176  0.1227

sigma^2 = 0.1089:  log likelihood = -17.01
AIC=42.02   AICc=42.77   BIC=50.26

#check ACF and PACF
acf(train.spen)
#significant positive lags through 11
pacf(train.spen)
#significant positive lags at 1, 3, 5

residuals(spen.final)
forecast::checkresiduals(spen.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(spen.final, plot=F)
#p = 0.6163 so null is not rejected and model is a good fit

forecast::accuracy(spen.final)
accuracy(forecast(spen.final, h=5), test.ci)

                   ME      RMSE       MAE       MPE      MAPE      MASE
Training set 0.06673541 0.3185482 0.2470833 208.81892 278.10005 0.7374556
Test set     1.21459658 1.2289063 1.2145966  64.01319  64.01319 3.6251372
ACF1 Theils U
Training set -0.02292016        NA
Test set      0.19715418  8.227514

#plot forecast
spen.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.spen))

#kod
sock.kod<-subset(sockeye.dat, region=='kod')
total.kod<-sock.kod %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

kod.ts<-ts(total.kod$lntotal, start=total.kod$year[1])
forecast::ndiffs(kod.ts, test='adf')
#1
forecast::ndiffs(kod.ts, test='kpss')
#1
#requires d=1
plot(diff(kod.ts))
#maybe needs second differencing?

train.kod<-window(kod.ts, start=1952, end=2010)
test.kod<-window(kod.ts, start=2011, end=2015)
fit.kod <- forecast::auto.arima(train.kod, trace=T)
# Best model: ARIMA(0,1,1)  
#AIC =56.34381
kod.final<-forecast::auto.arima(train.kod, approximation = F, stepwise = F)

ARIMA(3,1,0) 

Coefficients:
  ar1      ar2      ar3
-0.4678  -0.3059  -0.4194
s.e.   0.1183   0.1258   0.1160

sigma^2 = 0.1335:  log likelihood = -22.73
AIC=53.46   AICc=54.22   BIC=61.7

#check ACF and PACF
acf(train.kod)
#significant positive lags through 12 or 13
pacf(train.kod)
#significant pos lag at 1 and 4

residuals(kod.final)
forecast::checkresiduals(kod.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(kod.final, plot=F)
#p = 0.7926 so null is not rejected and model is a good fit

forecast::accuracy(kod.final)
accuracy(forecast(kod.final, h=5), test.ci)

                   ME      RMSE       MAE      MPE     MAPE      MASE
Training set 0.01967286 0.3527077 0.2796002 26.77484 82.84994 0.8560982
Test set     0.97584135 0.9857948 0.9758413 51.47752 51.47752 2.9878950
ACF1 Theils U
Training set -0.02639163        NA
Test set      0.20214107  6.667924

#plot forecast
kod.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.kod))

#pws
sock.pws<-subset(sockeye.dat, region=='pws')
total.pws<-sock.pws %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

pws.ts<-ts(total.pws$lntotal, start=total.pws$year[1])
forecast::ndiffs(pws.ts, test='adf')
#0
forecast::ndiffs(pws.ts, test='kpss')
#1
#maybe requires d=1
plot(pws.ts)
plot(diff(pws.ts))
#I like the differencing

train.pws<-window(pws.ts, start=1952, end=2010)
test.pws<-window(pws.ts, start=2011, end=2015)
fit.pws <- forecast::auto.arima(train.pws, trace=T)
# Best model: ARIMA(0,0,3) with non-zero mean 
#AIC = 67.8246
pws.final<-forecast::auto.arima(train.pws, approximation = F, stepwise = F)

ARIMA(0,0,3) with non-zero mean 

Coefficients:
  ma1     ma2      ma3    mean
0.4085  0.4070  -0.2924  0.3883
s.e.  0.1240  0.1193   0.1444  0.0765

sigma^2 = 0.1608:  log likelihood = -28.35
AIC=66.69   AICc=67.82   BIC=77.08

#check ACF and PACF
acf(train.pws)
#significant positive lag 1
pacf(train.pws)
#significant positive lag at 1
#This would make me choose a (1,1,1)

residuals(pws.final)
forecast::checkresiduals(pws.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(pws.final, plot=F)
#p = 0.3683 so null is not rejected and model is a good fit

forecast::accuracy(pws.final)
accuracy(forecast(pws.final, h=5), test.pws)

                    ME      RMSE       MAE      MPE      MAPE      MASE
Training set -1.232416e-05 0.3872030 0.2951372 66.22481 157.84719 0.7871924
Test set      4.170751e-01 0.4753133 0.4311048 45.18425  46.70421 1.1498463
ACF1 Theils U
Training set  0.003127207        NA
Test set     -0.042746070  2.288412

#plot forecast
pws.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.pws))

#seak
sock.seak<-subset(sockeye.dat, region=='seak')
total.seak<-sock.seak %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

seak.ts<-ts(total.seak$lntotal, start=total.seak$year[1])
forecast::ndiffs(seak.ts, test='adf')
#1
forecast::ndiffs(seak.ts, test='kpss')
#1
#requires d=1
plot(diff(seak.ts))

train.seak<-window(seak.ts, start=1952, end=2010)
test.seak<-window(seak.ts, start=2011, end=2015)
fit.seak <- forecast::auto.arima(train.seak, trace=T)
# Best model: ARIMA(0,1,1) 
#AIC = 8.581275
seak.final<-forecast::auto.arima(train.seak, approximation = F, stepwise = F)

ARIMA(0,1,1) 

Coefficients:
  ma1
-0.4705
s.e.   0.1299

sigma^2 = 0.06396:  log likelihood = -2.18
AIC=8.36   AICc=8.58   BIC=12.48

#check ACF and PACF
acf(train.seak)
#significant positive lags through 11
pacf(train.seak)
#maybe significant positive lag at 8
#suggests ARMA(p,q)?

residuals(seak.final)
forecast::checkresiduals(seak.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(seak.final, plot=F)
#p = 0.4304 so null is not rejected and model is a good fit

forecast::accuracy(seak.final)
accuracy(forecast(seak.final, h=5), test.seak)

                   ME      RMSE       MAE        MPE       MAPE     MASE
Training set -0.003354505 0.2485703 0.1909454 9229.38968 9287.87728 0.860063
Test set      0.243291865 0.2805902 0.2432919   55.55413   55.55413 1.095844
ACF1 Theils U
Training set 0.03286126        NA
Test set     0.21205641  1.359235

#plot forecast
seak.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak))

#nbc
sock.nbc<-subset(sockeye.dat, region=='nbc')
total.nbc<-sock.nbc %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

nbc.ts<-ts(total.nbc$lntotal, start=total.nbc$year[1])
forecast::ndiffs(nbc.ts, test='adf')
#0
forecast::ndiffs(nbc.ts, test='kpss')
#0

train.nbc<-window(nbc.ts, start=1952, end=2010)
test.nbc<-window(nbc.ts, start=2011, end=2015)
fit.nbc <- forecast::auto.arima(train.nbc, trace=T)
# Best model: ARIMA(0,0,1) with non-zero mean
#AIC = 65.47524
nbc.final<-forecast::auto.arima(train.nbc, approximation = F, stepwise = F)

ARIMA(0,0,1) with non-zero mean 

Coefficients:
  ma1    mean
0.4764  1.3359
s.e.  0.1231  0.0761

sigma^2 = 0.1641:  log likelihood = -29.52
AIC=65.04   AICc=65.48   BIC=71.27

#check ACF and PACF
acf(train.nbc)
#significant positive lags 1,5, negative lag 16
pacf(train.nbc)
#significant negative lag at 1 and 16
#maybe p=1

residuals(nbc.final)
forecast::checkresiduals(nbc.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(nbc.final, plot=F)
#p = 0.5427 so null is not rejected and model is a good fit

forecast::accuracy(nbc.final)
accuracy(forecast(nbc.final, h=5), test.nbc)

                     ME      RMSE       MAE        MPE      MAPE      MASE
Training set  0.004804253 0.3982042 0.3286109  -10.86895  29.39848 0.8642957
Test set     -0.371518749 0.6260534 0.4524202 -131.22412 137.24046 1.1899327
ACF1 Theils U
Training set -0.01478231        NA
Test set     -0.34558261  0.560803

#plot forecast
nbc.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.nbc))

#sbc
sock.sbc<-subset(sockeye.dat, region=='sbc')
total.sbc<-sock.sbc %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

sbc.ts<-ts(total.sbc$lntotal, start=total.sbc$year[1])
forecast::ndiffs(sbc.ts, test='adf')
#0
forecast::ndiffs(sbc.ts, test='kpss')
#0

train.sbc<-window(sbc.ts, start=1952, end=2010)
test.sbc<-window(sbc.ts, start=2011, end=2015)
fit.sbc <- forecast::auto.arima(train.sbc, trace=T)
# Best model: ARIMA(3,0,2) with non-zero mean 
#AIC = 108.4809
sbc.final<-forecast::auto.arima(train.sbc, approximation = F, stepwise = F)

ARIMA(3,0,2) with non-zero mean 

Coefficients:
  ar1      ar2     ar3      ma1     ma2    mean
0.5211  -0.9869  0.4196  -0.3068  0.8199  2.1514
s.e.  0.1577   0.0570  0.1374   0.1330  0.1250  0.0965

sigma^2 = 0.2993:  log likelihood = -46.14
AIC=106.28   AICc=108.48   BIC=120.83

#check ACF and PACF
acf(train.sbc)
#significant positive lags 4,8, negative 6
pacf(train.sbc)
#significant positive lag at 4
#suggests (4,0,0)

residuals(sbc.final)
forecast::checkresiduals(sbc.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(sbc.final, plot=F)
#p = 0.4291 so null is not rejected and model is a good fit

forecast::accuracy(sbc.final)
accuracy(forecast(sbc.final, h=5), test.sbc)

                   ME      RMSE       MAE        MPE     MAPE      MASE
Training set -0.002416181 0.5185349 0.4041240  -9.928464 24.80507 0.7166195
Test set     -0.470945358 0.6434687 0.4709454 -41.815610 41.81561 0.8351116
ACF1 Theils U
Training set -0.01417896        NA
Test set     -0.49636332 0.6833409

#plot forecast
sbc.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.sbc))

#wc
sock.wc<-subset(sockeye.dat, region=='wc')
total.wc<-sock.wc %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

wc.ts<-ts(total.wc$lntotal, start=total.wc$year[1])
forecast::ndiffs(wc.ts, test='adf')
#0
forecast::ndiffs(wc.ts, test='kpss')
#0

train.wc<-window(wc.ts, start=1952, end=2010)
test.wc<-window(wc.ts, start=2011, end=2015)
fit.wc <- forecast::auto.arima(train.wc, trace=T)
# Best model: ARIMA(1,0,0) with non-zero mean 
#AIC = 103.0188
wc.final<-forecast::auto.arima(train.wc, approximation = F, stepwise = F)

ARIMA(1,0,2) with non-zero mean 

Coefficients:
  ar1     ma1     ma2     mean
-0.9328  1.2620  0.4296  -1.1193
s.e.   0.0740  0.1397  0.1227   0.0922

sigma^2 = 0.28:  log likelihood = -44.45
AIC=98.9   AICc=100.04   BIC=109.29

#check ACF and PACF
acf(train.wc)
#significant positive lag 4
pacf(train.wc)
#significant positive lag at 4
#not sure how to interpret

residuals(wc.final)
forecast::checkresiduals(wc.final)
#residuals resemble white noise and demonstrate no sig correlation
forecast::checkresiduals(ci.final, plot=F)
#p = 0.4443 so null is not rejected and model is a good fit

forecast::accuracy(wc.final)
accuracy(forecast(wc.final, h=5), test.wc)

                    ME      RMSE       MAE       MPE     MAPE      MASE
Training set 0.002607753 0.5109361 0.4067843 -44.99390 68.70901 0.7026772
Test set     0.319794343 0.4046153 0.3549579 -54.83908 57.32173 0.6131525
ACF1 Theils U
Training set -0.01807646        NA
Test set      0.12345269  1.153447

#plot forecast
wc.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wc))
