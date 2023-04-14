#Madi's Exploratory Script 

library(tidyverse)
library(forecast)

#Autoregressive integrated moving average (ARIMA) model process 

# A. ARIMA(p,d,q) Model Selection
#   1) Evaluate Stationarity
#   2) Fix any issues and select a differencing level (d)
#   3) Selection of the AR level (p)
#   4) Selection of the MA level (q)

# B. Parameter Estimation
# C. Model Checking 
#   1)Test model residuals for normality (or other distribution assumptions)
#   2) Test model residuals for temporal correlation

#Code for testing accuracy with training and test data 
#Read data in
ruggerone_data <- readRDS("C:/GitHub/fish550-2023/Lab-1/Data_Images/ruggerone_data.rds")

#Filter by species (Pink)
dat <- ruggerone_data %>%  
  filter(species=="pink" & region=="ci") %>% 
  mutate(log.returns = log(returns)) %>% 
  select(year, log.returns)

#ID years 
unique(dat$year) #1952 start year, 2015 end year

#Plot by region
ruggerone_data %>% 
  filter(species=="pink") %>% 
  ggplot(aes(x=year, y=log(returns))) + 
  geom_line() + 
  ggtitle("pink salmon log abundance by region") +
  facet_wrap(~region)

#Note, no data in Korea 

PinkByRegion<-ruggerone_data %>%
  filter(region != "korea") %>% #Remove WA too, it's trouble 
  group_by(species, region, year) %>%
  summarize(total = sum(returns, na.rm=TRUE)) %>% 
  mutate(smoothedreturns = total +1)%>% ### TOO JANKY
  mutate(lnreturns = log(total)) %>%
  mutate(lnsmooth = log(smoothedreturns))%>%
  filter(species == "pink")%>% 
  print(n=10)

#Make Inf 0 and move the decimal place in returns over? 
#Or replace with NAs
#Post a question on the discussion board 


#check to see if start and end years are all the same 
PinkByRegion %>% group_by(region) %>% summarise(startyear = min(year), endyear = max(year))

#ID Stationarity 

#All data 

PinkByRegion %>%
  group_by(year) %>%
  summarize(total = sum(total, na.rm=T)) %>%
  ggplot(aes(x=year, y=log(total))) +
  geom_line() +
  ylab('Log (Returns)') +
  xlab('Year')

total.pink<-PinkByRegion %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(total, na.rm=T)))

pink.ts<-ts(total.pink$lntotal, 
            start=total.pink$year[1])
plot(diff(pink.ts)) #something odd happened between 1990 and 2005
acf(diff(pink.ts)) #ruh roh, ACF correlation for entire series 


#Differenced plots for all Regions 
PinkByRegion %>%
  group_by(region) %>%
  mutate(diff_total = c(NA, diff(total))) %>%
  ggplot(aes(x = year, y = diff_total)) +
  geom_line() +
  facet_wrap(~region, scales = "free_y") +
  ylab("Difference in Total Returns") +
  xlab("Year") +
  ggtitle("Diff by Region") 
  #ggfortify::ggstat_acf(method = "ma", na.action = na.pass)


#Compare ADF and KPSS
  #Note for the ADF null hypothesis is that the system is non-stationary (we want to reject)
  #The KPSS test null hypothesis is that there is stationarity

#Augmented Dicky Fuller 

tseries::adf.test(pink.ts, k=0)

tseries::kpss.test(pink.ts, null = c("Level", "Trend"))

# WAT DO HERE 

#By region
#regions vector
regions<-unique(PinkByRegion$region)
#regions key
regionskey<-c("Cook Inlet", "E. Kamchatka", "Japan", "Kodiak", "Russia", "N.British Columbia",
              "Prince William Sound", "S. Alaska Pen.", "S. British Columbia", "SE Alaska", "W. Kamchatka", "Washington", "W. Alaska")
names(regionskey)<-regions #for plotting

#Functions from Zoe! 

#======================================================

ACFandPACF<-function(reg){
  Pinkdat<-PinkByRegion %>% filter(region == reg)
  #create time series
  #datts <- ts(Pinkdat$lnreturns, start=Pinkdat$year[1]) 
  datts <- ts(Pinkdat$lnsmooth, start=Pinkdat$year[1])  #Changed this to lnsmooth because of issues
  return(list(a = acf(datts), p = pacf(datts)))
}

FitModFunction<-function(reg, forelevel){
  #filter region
  Pinkdat<-PinkByRegion %>% filter(region == reg)
  #create time series
  #datts <- ts(Pinkdat$lnreturns, start=Pinkdat$year[1])} 
  datts <- ts(Pinkdat$lnsmooth, start=Pinkdat$year[1])  #this assumes the first year in data is the start of the time series (they are in order) 
  cutoff<-2015-forelevel
  train <- window(datts, 1952, cutoff)
  test <- window(datts, cutoff+1, 2015)
  
  mod <- auto.arima(train)
  mod
  
  res<-accuracy(forecast(mod, h=forelevel), test)[2,"MASE"] #test set MASE
  
  return(list(Fit = mod, MASE = res))
}

#================================================================

#loop through regions/levels
DiagPlots<-lapply(regions, ACFandPACF) #HAD TO CHANGE lnRETURNS with Smoothing
   names(DiagPlots)<-regions

#ACF and PACF

#ACF plots for each region
par(mfrow=c(3,5))
for(r in 1:length(regions)){
  plot(DiagPlots[[r]][[1]], main = paste0("Region: ", regionskey[r]))
}

#PACF plots for each region
par(mfrow=c(3,5))
for(r in 1:length(regions)){
  plot(DiagPlots[[r]][[2]], main = paste0("Region: ", regionskey[r]))
}

#ARIMA model and Forcasting

#forecast levels
forecastlevels<-c(5, 10, 20)
#all combinations
Allcombs<-expand_grid(regions, forecastlevels)

#Added a column above as a smoother otherwise No suitable ARIMA model found error. 
RegionMods<-mapply(FitModFunction, Allcombs$regions, Allcombs$forecastlevels, SIMPLIFY = FALSE)
head(RegionMods)
names(RegionMods) #should be three for each region

#getting MASE
RegionMASE<-sapply(RegionMods, function(x){y<-x$MASE})
RegionBestMod<-sapply(RegionMods, function(x){y<-as.character(x$Fit)})
#combine into tables
ResultsTable<-Allcombs %>% add_column(Model = RegionBestMod, MASE = RegionMASE)
ResultsTable

#plot results
ggplot(ResultsTable) + 
  geom_bar(aes(x = regions, y = MASE, fill = as.factor(forecastlevels)), stat = "identity", position = "dodge") + 
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  scale_x_discrete(labels = as_labeller(regionskey)) +
  labs(fill = "Forecast Levels", x = "Region") + 
  ggtitle("Pink") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))

#Plot Forecasts

#One good, one med, one bad 

#Good: Washington
pink.wa<-subset(PinkByRegion, region=='wa')
wa.ts<-ts(pink.wa$lnsmooth, start=pink.wa$year[1])
#create training and test datasets for the 5, 10, and 20 year forecasts
train.wa5<-window(wa.ts, start=1952, end=2010)
test.wa5<-window(wa.ts, start=2011, end=2015)
wa.final5<-forecast::auto.arima(train.wa5, approximation = F, stepwise = F)

train.wa10<-window(wa.ts, start=1952, end=2005)
test.wa10<-window(wa.ts, start=2006, end=2015)
wa.final10<-forecast::auto.arima(train.wa10, approximation = F, stepwise = F)

train.wa20<-window(wa.ts, start=1952, end=1995)
test.wa20<-window(wa.ts, start=1996, end=2015)
wa.final20<-forecast::auto.arima(train.wa20, approximation = F, stepwise = F)

#Plots
par(mfrow=c(3,1))

wa.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wa5))

wa.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wa10))

wa.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.wa20))



############

RegionMods[1] %>%
  forecast(xreg=fourier(PinkByRegion, K=c(1), h=15)) %>%
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




#
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

#plot by region 

pink.dat %>%
  filter(region != "korea") %>%
  group_by(region, year) %>%
  summarize(total = sum(returns, na.rm=T)) %>%
  ggplot(aes(x = year, y = log(total))) %>%
  geom_line() %>%
  facet_wrap(~ region) %>%
  ylab('Log (Returns)') %>%
  xlab('Year') -> NewPink

#Look at species in total
total.pink<-pink.dat %>%
  group_by(year, region) %>%
  summarize(lntotal=log(sum(returns, na.rm=T)))

#Create time series object for total pinks
pink.ts<-ts(total.pink$lntotal, 
               start=total.pink$year[1])
#Plot the difference and ACD
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

