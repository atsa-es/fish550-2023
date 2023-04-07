ruggerone_data <- readRDS(here::here("Lab-1", "Data_Images", "ruggerone_data.rds"))

library(tidyverse)



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