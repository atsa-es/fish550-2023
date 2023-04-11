#Madi's Exploratory Script 

library(tidyverse)
library(forecast)

#Code for testing accuracy with training and test data 
#Read data in
ruggerone_data <- readRDS("C:/GitHub/fish550-2023/Lab-1/Data_Images/ruggerone_data.rds")
#Filter by species (Pink)
dat <- ruggerone_data %>%  filter(species=="pink" & region=="ci") %>% mutate(log.returns = log(returns)) %>% select(year, log.returns)
#Create a time series 
dat <- ts(dat$log.returns, start=dat$year[1])
#Training window
train <- window(dat, start=1961, end=1980)
#Testing window
test <- window(dat, start=1981, end=1985)
#Fit Forecast
fit <- forecast::auto.arima(train)
#Assess acuracy
accuracy(fit) # fit within the train data
accuracy(forecast(fit, h=5), test) # comparison of a forecast (5 year) to the test data (5 year)

#The salmon residuals are going to look patterned because they return on a cycle

#=======================
#Code that Eli provided for Emma to address cyclic nature of Salmon

#create forecast model for Pink salmon across all regions

#Is the data stationary?
pink.dat<-subset(ruggerone_data, species=='pink')

pink.dat %>%
  group_by(year) %>%
  summarize(total = sum(returns, na.rm=T)) %>%
  ggplot(aes(x=year, y=log(total))) +
  geom_line() +
  ylab('Log (Returns)') +
  xlab('Year')

total.pink<-pink.dat %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

pink.ts<-ts(total.pink$lntotal, 
               start=total.pink$year[1])
plot(diff(pink.ts))
# there is a negative autocorrelation at all years in the data
    #Pinks have a two year life cycle, and the ACF shows significance the entire time series. 
    #This makes sense.
acf(diff(pink.ts)) 
# We could try including "seasonality" in our model, I started with 10.....no great reason
y <- msts(pink.ts, seasonal.periods=10)
# Decomposition
y %>% mstl() %>% autoplot()
# Horribly we don't keep the year numbering when we add the season
train.pink<-window(y, start=c(1,1), end=c(6,6))
test.pink<-window(y, start=c(4,1), end=c(6,1))
fit.final<-forecast::auto.arima(train.pink)
forecast::checkresiduals(fit.final)
fit.final %>%
  forecast(h=15) %>%
  autoplot() + 
  geom_point(aes(x=x, y=y), data=fortify(test.pink)) +
  # this stuff is to fix the x-axis back to years
  scale_x_continuous(breaks = seq(1-2/7,10,5/7),
                     labels = seq(1950,2015,5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("")

# Another approach is to use a cyclic covariate. 
# We will cover this in 2 weeks
# https://otexts.com/fpp2/complexseasonality.html
fit <- auto.arima(train.pink, seasonal=FALSE, 
                  xreg=fourier(train.pink, K=c(1)))
fit %>%
  forecast(xreg=fourier(test.pink, K=c(1), h=15)) %>%
  autoplot() + 
  geom_point(aes(x=x, y=y), data=fortify(test.pink)) +
  # this stuff is to fix the x-axis back to years
  scale_x_continuous(breaks = seq(1-2/7,10,5/7),
                     labels = seq(1950,2015,5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("")

# No seasonality
train.pink<-window(pink.ts, start=1952, end=2000)
test.pink<-window(pink.ts, start=2001, end=2015)
fit.final<-forecast::auto.arima(train.pink)
forecast::checkresiduals(fit.final)
fit.final %>%
  forecast(h=15) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.pink))

