#Lab 1: Ruggerone Sockeye Returns
#code adapted from Zoe

ruggerone_data <- readRDS(here::here("Lab-1", "Data_Images", "ruggerone_data.rds"))

library(tidyverse)
library(forecast)

sockeye.dat<-subset(ruggerone_data, species=='sockeye')

#Forecast returns for full dataset
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

train.sockeye<-window(sockeye.ts, start=1952, end=2010)
test.sockeye<-window(sockeye.ts, start=2011, end=2015)

fit <- forecast::auto.arima(train.sockeye, trace=T)
# Best model: ARIMA(0,1,2)  
#AIC = 1.111738
fit.final<-forecast::auto.arima(train.sockeye, approximation = F, stepwise = F)
#ARIMA(0,1,2) 
#check ACF and PACF
acf(train.sockeye)
#significant lags through 6, then 9-11
pacf(train.sockeye)
#significant lag at 9

forecast::checkresiduals(fit.final)
forecast::checkresiduals(fit.final, plot=F)
#p = 0.01 so null is rejected and model is not a good fit
accuracy(forecast(fit.final, h=5), test.sockeye)
#MASE training = 0.912, test = 0.927

#plot forecast
fit.final %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.sockeye))

#Below is adapted from Zoe's work
#removing Korea Japan because there's no data
SockByRegion<-ruggerone_data %>%
  filter(region != "japan") %>%
  filter(region != "korea") %>%
  group_by(species, region, year) %>%
  summarize(total = sum(returns, na.rm=TRUE)) %>% 
  mutate(lnreturns = log(total)) %>%
  filter(species == "sockeye")
head(SockByRegion)

#regions vector
regions<-unique(SockByRegion$region)
#regions key
regionskey<-c("Cook Inlet", "E. Kamchatka", "Kodiak", "Russia", "N.British Columbia",
              "Prince William Sound", "S. Alaska Pen.", "S. British Columbia", "SE Alaska", "W. Kamchatka", "Washington", "W. Alaska")
names(regionskey)<-regions #for plotting
#forecast levels
forecastlevels<-c(5, 10, 20)
#all combinations
Allcombs<-expand_grid(regions, forecastlevels)

#ACF and PACF
ACFandPACF<-function(reg){
  Sockdat<-SockByRegion %>% filter(region == reg)
  #create time series
  datts <- ts(Sockdat$lnreturns, start=Sockdat$year[1])
  return(list(a = acf(datts), p = pacf(datts)))
}

#function for ARIMA models
FitModFunction<-function(reg, forelevel){
  #filter region
  Sockdat<-SockByRegion %>% filter(region == reg)
  #create time series
  datts <- ts(Sockdat$lnreturns, start=Sockdat$year[1]) #this assumes the first year in data is the start of the time series (they are in order) 
  cutoff<-2015-forelevel
  train <- window(datts, 1952, cutoff)
  test <- window(datts, cutoff+1, 2015)
  
  mod <- auto.arima(train)
  mod
  
  res<-accuracy(forecast(mod, h=forelevel), test)[2,"MASE"] #test set MASE
  
  return(list(Fit = mod, MASE = res))
}


#loop through regions/levels
DiagPlots<-lapply(regions, ACFandPACF)
names(DiagPlots)<-regions

#ACF plots for each region
par(mfrow=c(3,4))
for(r in 1:length(regions)){
  plot(DiagPlots[[r]][[1]], main = paste0("Region: ", regionskey[r]))
}

#PACF plots for each region
par(mfrow=c(3,4))
for(r in 1:length(regions)){
  plot(DiagPlots[[r]][[2]], main = paste0("Region: ", regionskey[r]))
}


RegionMods<-mapply(FitModFunction, Allcombs$regions, Allcombs$forecastlevels, SIMPLIFY = FALSE)

head(RegionMods)
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
  ggtitle("Sockeye") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))
                                       