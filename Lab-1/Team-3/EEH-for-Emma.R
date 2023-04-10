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

sockeye.ts<-ts(total.sockeye$lntotal, 
               start=total.sockeye$year[1])
plot(diff(sockeye.ts))
# there is a negative autocorrelation at 7 years in the data
acf(diff(sockeye.ts))

# We could try including that 7 year "seasonality" in our model
y <- msts(sockeye.ts, seasonal.periods=7)
# Decomposition
y %>% mstl() %>% autoplot()
# Horribly we don't keep the year numbering when we add the season
train.sockeye<-window(y, start=c(1,1), end=c(7,7))
test.sockeye<-window(y, start=c(8,1), end=c(10,1))
fit.final<-forecast::auto.arima(train.sockeye)
forecast::checkresiduals(fit.final)
fit.final %>%
  forecast(h=15) %>%
  autoplot() + 
    geom_point(aes(x=x, y=y), data=fortify(test.sockeye)) +
    # this stuff is to fix the x-axis back to years
    scale_x_continuous(breaks = seq(1-2/7,10,5/7),
                     labels = seq(1950,2015,5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("")

# Another approach is to use a cyclic covariate. 
# We will cover this in 2 weeks
# https://otexts.com/fpp2/complexseasonality.html
fit <- auto.arima(train.sockeye, seasonal=FALSE, 
                  xreg=fourier(train.sockeye, K=c(1)))
fit %>%
  forecast(xreg=fourier(test.sockeye, K=c(1), h=15)) %>%
  autoplot() + 
  geom_point(aes(x=x, y=y), data=fortify(test.sockeye)) +
  # this stuff is to fix the x-axis back to years
  scale_x_continuous(breaks = seq(1-2/7,10,5/7),
                     labels = seq(1950,2015,5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("")

# No seasonality
train.sockeye<-window(sockeye.ts, start=1952, end=2000)
test.sockeye<-window(sockeye.ts, start=2001, end=2015)
fit.final<-forecast::auto.arima(train.sockeye)
forecast::checkresiduals(fit.final)
fit.final %>%
  forecast(h=15) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.sockeye))
