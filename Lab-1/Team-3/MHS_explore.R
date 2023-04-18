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
  filter(region != "korea") %>% 
  filter(region != "wa") %>%#Remove WA too, it's trouble 
  group_by(species, region, year) %>%
  summarize(total = sum(returns, na.rm=TRUE)) %>% 
  mutate(lnreturns = log(total)) %>%
  mutate(lnreturns = ifelse(lnreturns == -Inf, NA, lnreturns)) %>%
  filter(species == "pink")%>% 
  print(n=10)


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

#The ACF looks like there is a lot of corelation, that's probably because 
#Pinks have a very consistent two year cycle

#Let's try a forecast model to see what happens

#Let's train and test with a 10 year period
train.pink<-window(pink.ts, start=1952, end=2005)
test.pink<-window(pink.ts, start=2006, end=2015)

fit <- forecast::auto.arima(train.pink, trace=T)
fit.final.pink<-forecast::auto.arima(train.pink, approximation = F, stepwise = F)

#30 year forcast, not so believeable
fit.final.pink %>%
  forecast(h=15) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.pink))

#Lets parse this into two pieces, Even and Odd years 
#======================================================

#Even Years: Totals 

PinkByRegion_even<-PinkByRegion %>% 
  filter(year %% 2 == 0)

#Odd Years: Totals 
PinkByRegion_odd<-PinkByRegion %>% 
  filter(year %% 2 == 1)

#Trends
PinkByRegion_even %>%
  group_by(year) %>%
  summarize(total = sum(total, na.rm=T)) %>%
  ggplot(aes(x=year, y=log(total))) +
  geom_line() +
  ylab('Log (Returns)') +
  xlab('Year') +
  ggtitle('Total Pinks (Even Years)')

PinkByRegion_odd %>%
  group_by(year) %>%
  summarize(total = sum(total, na.rm=T)) %>%
  ggplot(aes(x=year, y=log(total))) +
  geom_line() +
  ylab('Log (Returns)') +
  xlab('Year') +
  ggtitle('Total Pinks (Odd Years)')


#Even Years -- Total

total.pink_even<-PinkByRegion_even %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(total, na.rm=T)))

pink.ts_even<-ts(total.pink_even$lntotal, 
            start=total.pink$year[1], frequency = 0.5)
plot(diff(pink.ts_even)) #Looks pretty stationary
acf(diff(pink.ts_even)) #This looks much better

#Train and test for a 10 year period
train.pink_even<-window(pink.ts_even, start=1952, end=2004)
test.pink_even<-window(pink.ts_even, start=2006, end=2014)

fit <- forecast::auto.arima(train.pink_even, trace=T)
fit.final.pink_even<-forecast::auto.arima(train.pink_even, approximation = F, stepwise = F)

fit.final.pink_even %>%
  forecast(h=15) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.pink_even))

#That is a straight ass line..... not very good.....

#Odd Years -- Total

total.pink_odd<-PinkByRegion_odd %>%
  group_by(year) %>%
  summarize(lntotal=log(sum(total, na.rm=T)))

pink.ts_odd<-ts(total.pink_odd$lntotal, 
                 start=total.pink_odd$year[1], frequency = 0.5)
plot(diff(pink.ts_odd)) #Looks pretty stationary
acf(diff(pink.ts_odd)) #This also looks better

#Train and test for a 10 year period
train.pink_odd<-window(pink.ts_odd, start=1953, end=2005)
test.pink_odd<-window(pink.ts_odd, start=2007, end=2015)

fit <- forecast::auto.arima(train.pink_odd, trace=T)
fit.final.pink_odd<-forecast::auto.arima(train.pink_odd, approximation = F, stepwise = F)

fit.final.pink_odd %>%
  forecast(h=15) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.pink_odd))

#This is also a straight ass line..... not very good.....


#======= Regional Parts =================#

#Even Years First 
#Differenced plots for all Regions 
PinkByRegion_even %>%
  group_by(region) %>%
  mutate(diff_total = c(NA, diff(total))) %>%
  ggplot(aes(x = year, y = diff_total)) +
  geom_line() +
  facet_wrap(~region, scales = "free_y") +
  ylab("Difference in Total Returns") +
  xlab("Year") +
  ggtitle("Diff by Region (Even Years)") 
  #ggfortify::ggstat_acf(method = "ma", na.action = na.pass)

# Odd Years First 
#Differenced plots for all Regions 
PinkByRegion_odd %>%
  group_by(region) %>%
  mutate(diff_total = c(NA, diff(total))) %>%
  ggplot(aes(x = year, y = diff_total)) +
  geom_line() +
  facet_wrap(~region, scales = "free_y") +
  ylab("Difference in Total Returns") +
  xlab("Year") +
  ggtitle("Diff by Region (Odd Years)") 

#Compare ADF and KPSS
  #Note for the ADF null hypothesis is that the system is non-stationary (we want to reject)
  #The KPSS test null hypothesis is that there is stationarity

#Augmented Dicky Fuller 

tseries::adf.test(pink.ts_even, k=0)
tseries::adf.test(pink.ts_odd, k=0)

tseries::kpss.test(pink.ts_even, null = c("Level", "Trend"))
tseries::kpss.test(pink.ts_odd, null = c("Level", "Trend"))

#Great, both these tests are now showing stationarity 


#Functions from Zoe! 

#======================================================

ACFandPACF_even<-function(reg){
  Pinkdat<-PinkByRegion_even %>% filter(region == reg)
  #create time series
  datts <- ts(Pinkdat$lnreturns, start=Pinkdat$year[1]) 
  return(list(a = acf(datts), p = pacf(datts)))
}

ACFandPACF_odd<-function(reg){
  Pinkdat<-PinkByRegion_odd %>% filter(region == reg)
  #create time series
  datts <- ts(Pinkdat$lnreturns, start=Pinkdat$year[1]) 
  return(list(a = acf(datts), p = pacf(datts)))
}

FitModFunction_even<-function(reg, forelevel){
  #filter region
  Pinkdat<-PinkByRegion_even %>% filter(region == reg)
  #create time series
  datts <- ts(Pinkdat$lnreturns, start=Pinkdat$year[1], frequency = 0.5) #set Frequency for 0.5
  cutoff<-2014-forelevel
  train <- window(datts, Pinkdat$year[1], cutoff)
  test <- window(datts, cutoff+1, 2014)
  
  mod <- auto.arima(train)
  #testing to be sure that this is the best model (is the best mode the simplest if it is within 2 AIC values?)
  trace <- capture.output({
    # assign so it doesn't pollute the output
    model <- auto.arima(datts, trace = TRUE)
  })
  con    <- textConnection(trace)
  models <- read.table(con, sep=":")
  close(con)
  
  #getting the "best models" that are within 2 AIC units
  BestMods<-models%>% filter(row_number() != nrow(models)) %>% mutate(AIC = replace(V2, V2 == "Inf", 99999), AIC = as.numeric(AIC), DeltaAIC = AIC-min(AIC)) %>% filter(DeltaAIC <= 2.0)
  for(i in 1:nrow(BestMods)){
    BestMods$Mod[i]<-strsplit(strsplit(strsplit(BestMods$V1[i], "[(]")[[1]][2], "[)]")[[1]][1],"[,]")
    BestMods$npar[i]<-sum(as.numeric(BestMods$Mod[i][[1]][c(1,3)]))
    if(strsplit(strsplit(BestMods$V1[i], "[(]")[[1]][2], "[)]")[[1]][2] == " with drift         "){
      BestMods$npar[i] = BestMods$npar[i] + 1
    }
  }
  
  New<-BestMods %>% filter(npar == min(npar))
  if(0 %in% New$DeltaAIC){
    #auto arima picked the best model
    res<-accuracy(forecast(mod, h=forelevel), test)[2,"MASE"] #test set MASE
  }else{
    #of the models with the fewest parameters, pick the lowest AIC
    newmod<-New %>% filter(AIC == min(AIC)) %>% select(Mod)
    mod<-Arima(train, order = as.numeric(strsplit(newmod$Mod[[1]], "[,]")), include.constant = TRUE)
    res<-accuracy(forecast(mod, h=forelevel), test)[2,"MASE"] #test set MASE
  }
  
  return(list(Fit = mod, MASE = res, Bm = BestMods)) #include best mods for testing to see that it's doing what I want
  
}

FitModFunction_odd<-function(reg, forelevel){
  #filter region
  Pinkdat<-PinkByRegion_odd %>% filter(region == reg)
  #create time series
  datts <- ts(Pinkdat$lnreturns, start=Pinkdat$year[1], frequency = 0.5) 
  cutoff<-2015-forelevel
  train <- window(datts, Pinkdat$year[1], cutoff)
  test <- window(datts, cutoff+1, 2015)
  
  mod <- auto.arima(train)
  #testing to be sure that this is the best model (is the best mode the simplest if it is within 2 AIC values?)
  trace <- capture.output({
    # assign so it doesn't pollute the output
    model <- auto.arima(datts, trace = TRUE)
  })
  con    <- textConnection(trace)
  models <- read.table(con, sep=":")
  close(con)
  
  #getting the "best models" that are within 2 AIC units
  BestMods<-models%>% filter(row_number() != nrow(models)) %>% mutate(AIC = replace(V2, V2 == "Inf", 99999), AIC = as.numeric(AIC), DeltaAIC = AIC-min(AIC)) %>% filter(DeltaAIC <= 2.0)
  for(i in 1:nrow(BestMods)){
    BestMods$Mod[i]<-strsplit(strsplit(strsplit(BestMods$V1[i], "[(]")[[1]][2], "[)]")[[1]][1],"[,]")
    BestMods$npar[i]<-sum(as.numeric(BestMods$Mod[i][[1]][c(1,3)]))
    if(strsplit(strsplit(BestMods$V1[i], "[(]")[[1]][2], "[)]")[[1]][2] == " with drift         "){
      BestMods$npar[i] = BestMods$npar[i] + 1
    }
  }
  
  New<-BestMods %>% filter(npar == min(npar))
  if(0 %in% New$DeltaAIC){
    #auto arima picked the best model
    res<-accuracy(forecast(mod, h=forelevel), test)[2,"MASE"] #test set MASE
  }else{
    #of the models with the fewest parameters, pick the lowest AIC
    newmod<-New %>% filter(AIC == min(AIC)) %>% select(Mod)
    mod<-Arima(train, order = as.numeric(strsplit(newmod$Mod[[1]], "[,]")), include.constant = TRUE)
    res<-accuracy(forecast(mod, h=forelevel), test)[2,"MASE"] #test set MASE
  }
  
  return(list(Fit = mod, MASE = res, Bm = BestMods)) #include best mods for testing to see that it's doing what I want
  
}

#================================================================

#Regional Considerations 

#By region
#regions vector
regions<-unique(PinkByRegion$region)
#regions key
regionskey<-c("Cook Inlet", "E. Kamchatka", "Japan", "Kodiak", "Russia", "N.British Columbia",
              "Prince William Sound", "S. Alaska Pen.", "S. British Columbia", "SE Alaska", "W. Kamchatka", "W. Alaska")
names(regionskey)<-regions #for plotting

#loop through regions/levels
#Even
#par(mfrow=c(3,4))
DiagPlots_even<-lapply(regions, ACFandPACF_even)
   names(DiagPlots_even)<-regions
   
#Odd
DiagPlots_odd<-lapply(regions, ACFandPACF_odd)
   names(DiagPlots_odd)<-regions   

#ACF and PACF
#=====================
#ACF plots for each region

#Even Years 
####################
par(mfrow=c(3,4))
for(r in 1:length(regions)){
  plot(DiagPlots_even[[r]][[1]], main = paste0("Region: ", regionskey[r]))
}

#PACF plots for each region
par(mfrow=c(3,4))
for(r in 1:length(regions)){
  plot(DiagPlots_even[[r]][[2]], main = paste0("Region: ", regionskey[r]))
}

#Odd Years 
##################
par(mfrow=c(3,4))
for(r in 1:length(regions)){
  plot(DiagPlots_odd[[r]][[1]], main = paste0("Region: ", regionskey[r]))
}

#PACF plots for each region
par(mfrow=c(3,4))
for(r in 1:length(regions)){
  plot(DiagPlots_odd[[r]][[2]], main = paste0("Region: ", regionskey[r]))
}

#================================================
#ARIMA model and Forcasting
#================================================

#forecast levels
forecastlevels<-c(5, 10, 20)
#all combinations
Allcombs<-expand_grid(regions, forecastlevels)

#================================================

#Split into even and odd years  
#Even
RegionMods_even<-mapply(FitModFunction_even, Allcombs$regions, Allcombs$forecastlevels, SIMPLIFY = FALSE)
head(RegionMods_even)
names(RegionMods_even) #should be three for each region

#Odd 
RegionMods_odd<-mapply(FitModFunction_odd, Allcombs$regions, Allcombs$forecastlevels, SIMPLIFY = FALSE)
head(RegionMods_odd)
names(RegionMods_odd) #should be three for each region

#getting MASE
#Even
RegionMASE_even<-sapply(RegionMods_even, function(x){y<-x$MASE})
RegionBestMod_even<-sapply(RegionMods_even, function(x){y<-as.character(x$Fit)})
#Odd
RegionMASE_odd<-sapply(RegionMods_odd, function(x){y<-x$MASE})
RegionBestMod_odd<-sapply(RegionMods_odd, function(x){y<-as.character(x$Fit)})
#combine into tables
#even
ResultsTable_even<-Allcombs %>% add_column(Model = RegionBestMod_even, MASE = RegionMASE_even)
ResultsTable_even

#Odd
ResultsTable_odd<-Allcombs %>% add_column(Model = RegionBestMod_odd, MASE = RegionMASE_odd)
ResultsTable_odd

#plot results


ggplot(ResultsTable_even) + 
  geom_bar(aes(x = regions, y = MASE, fill = as.factor(forecastlevels)), stat = "identity", position = "dodge") + 
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  scale_x_discrete(labels = as_labeller(regionskey)) +
  labs(fill = "Forecast Levels", x = "Region") + 
  ggtitle("Pinks Even Years") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))

ggplot(ResultsTable_odd) + 
  geom_bar(aes(x = regions, y = MASE, fill = as.factor(forecastlevels)), stat = "identity", position = "dodge") + 
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  scale_x_discrete(labels = as_labeller(regionskey)) +
  labs(fill = "Forecast Levels", x = "Region") + 
  ggtitle("Pinks Odd Years") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))


#=================================================
#Plot Forecasts

#==============================================
#Cook Inlet (bad in both even and odd years)
#Japan, even (medium) and odd (pretty good) 
#SE AK, (Good for both even and odd years)
#---------------------------------------------
#One good, one med, one bad 

#BAD: Cook Inlet
#----------------------------------------------------
#Even Years
pink.ci_even<-subset(PinkByRegion_even, region=='ci')
ci.ts_even<-ts(pink.ci_even$lnreturns, start=pink.ci_even$year[1], frequency = 0.5)

#create training and test datasets for the 5 year forecast 
train.ci_even5<-window(ci.ts_even, start=1952, end=2010)
test.ci_even5<-window(ci.ts_even, start=2011, end=2014)
ci_even.final5<-forecast::auto.arima(train.ci_even5, approximation = F, stepwise = F)

#create training and test datasets for the 10 year forecast 
train.ci_even10<-window(ci.ts_even, start=1952, end=2005)
test.ci_even10<-window(ci.ts_even, start=2006, end=2014)
ci_even.final10<-forecast::auto.arima(train.ci_even10, approximation = F, stepwise = F)

#create training and test datasets for the 20 year forecast 
train.ci_even20<-window(ci.ts_even, start=1952, end=1995)
test.ci_even20<-window(ci.ts_even, start=1996, end=2014)
ci_even.final20<-forecast::auto.arima(train.ci_even20, approximation = F, stepwise = F)

#------------------------------
#Odd Years
pink.ci_odd<-subset(PinkByRegion_odd, region=='ci')
ci.ts_odd<-ts(pink.ci_odd$lnreturns, start=pink.ci_odd$year[1], frequency = 0.5)

#create training and test datasets for the 5 year forecast 
train.ci_odd5<-window(ci.ts_odd, start=1952, end=2010)
test.ci_odd5<-window(ci.ts_odd, start=2011, end=2014)
ci_odd.final5<-forecast::auto.arima(train.ci_odd5, approximation = F, stepwise = F)

#create training and test datasets for the 10 year forecast 
train.ci_odd10<-window(ci.ts_odd, start=1952, end=2005)
test.ci_odd10<-window(ci.ts_odd, start=2006, end=2014)
ci_odd.final10<-forecast::auto.arima(train.ci_odd10, approximation = F, stepwise = F)

#create training and test datasets for the 20 year forecast 
train.ci_odd20<-window(ci.ts_odd, start=1952, end=1995)
test.ci_odd20<-window(ci.ts_odd, start=1996, end=2014)
ci_odd.final20<-forecast::auto.arima(train.ci_odd20, approximation = F, stepwise = F)

#Plots

library(cowplot)

plot_even_5 <- ci_even.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci_even5))

plot_even_10 <- ci_even.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci_even10))

plot_even_20 <- ci_even.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci_even20))

plot_odd_5 <- ci_odd.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci_odd5))

plot_odd_10 <- ci_odd.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci_odd10))

plot_odd_20 <- ci_odd.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.ci_odd20))

plot_grid(plot_even_5, plot_even_10, plot_even_20, plot_odd_5, plot_odd_10, plot_odd_20, ncol = 3, nrow = 2)

#---------------------------------------------------------------------

#--Medium Performance, Japan 

#----------------------------------------------------
#Even Years
pink.j_even<-subset(PinkByRegion_even, region=='japan')
j.ts_even<-ts(pink.j_even$lnreturns, start=pink.j_even$year[1], frequency = 0.5)

#create training and test datasets for the 5 year forecast 
train.j_even5<-window(j.ts_even, start=1952, end=2010)
test.j_even5<-window(j.ts_even, start=2011, end=2014)
j_even.final5<-forecast::auto.arima(train.j_even5, approximation = F, stepwise = F)

#create training and test datasets for the 10 year forecast 
train.j_even10<-window(j.ts_even, start=1952, end=2005)
test.j_even10<-window(j.ts_even, start=2006, end=2014)
j_even.final10<-forecast::auto.arima(train.j_even10, approximation = F, stepwise = F)

#create training and test datasets for the 20 year forecast 
train.j_even20<-window(j.ts_even, start=1952, end=1995)
test.j_even20<-window(j.ts_even, start=1996, end=2014)
j_even.final20<-forecast::auto.arima(train.j_even20, approximation = F, stepwise = F)

#------------------------------
#Odd Years
pink.j_odd<-subset(PinkByRegion_odd, region=='japan')
j.ts_odd<-ts(pink.j_odd$lnreturns, start=pink.j_odd$year[1], frequency = 0.5)

#create training and test datasets for the 5 year forecast 
train.j_odd5<-window(j.ts_odd, start=1952, end=2010)
test.j_odd5<-window(j.ts_odd, start=2011, end=2014)
j_odd.final5<-forecast::auto.arima(train.j_odd5, approximation = F, stepwise = F)

#create training and test datasets for the 10 year forecast 
train.j_odd10<-window(j.ts_odd, start=1952, end=2005)
test.j_odd10<-window(j.ts_odd, start=2006, end=2014)
j_odd.final10<-forecast::auto.arima(train.j_odd10, approximation = F, stepwise = F)

#create training and test datasets for the 20 year forecast 
train.j_odd20<-window(j.ts_odd, start=1952, end=1995)
test.j_odd20<-window(j.ts_odd, start=1996, end=2014)
j_odd.final20<-forecast::auto.arima(train.j_odd20, approximation = F, stepwise = F)

#Plots 

plot_even_5 <- j_even.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.j_even5))

plot_even_10 <- j_even.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.j_even10))

plot_even_20 <- j_even.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.j_even20))

plot_odd_5 <- j_odd.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.j_odd5))

plot_odd_10 <- j_odd.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.j_odd10))

plot_odd_20 <- j_odd.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.j_odd20))

plot_grid(plot_even_5, plot_even_10, plot_even_20, plot_odd_5, plot_odd_10, plot_odd_20, ncol = 3, nrow = 2)

#---------------------------------------------------------------------

#--Good Performance, SE AK

#----------------------------------------------------
#Even Years
pink.seak_even<-subset(PinkByRegion_even, region=='seak')
seak.ts_even<-ts(pink.seak_even$lnreturns, start=pink.seak_even$year[1], frequency = 0.5)

#create training and test datasets for the 5 year forecast 
train.seak_even5<-window(seak.ts_even, start=1952, end=2010)
test.seak_even5<-window(seak.ts_even, start=2011, end=2014)
seak_even.final5<-forecast::auto.arima(train.seak_even5, approximation = F, stepwise = F)

#create training and test datasets for the 10 year forecast 
train.seak_even10<-window(seak.ts_even, start=1952, end=2005)
test.seak_even10<-window(seak.ts_even, start=2006, end=2014)
seak_even.final10<-forecast::auto.arima(train.seak_even10, approximation = F, stepwise = F)

#create training and test datasets for the 20 year forecast 
train.seak_even20<-window(seak.ts_even, start=1952, end=1995)
test.seak_even20<-window(seak.ts_even, start=1996, end=2014)
seak_even.final20<-forecast::auto.arima(train.seak_even20, approximation = F, stepwise = F)

#------------------------------
#Odd Years
pink.seak_odd<-subset(PinkByRegion_odd, region=='seak')
seak.ts_odd<-ts(pink.seak_odd$lnreturns, start=pink.seak_odd$year[1], frequency = 0.5)

#create training and test datasets for the 5 year forecast 
train.seak_odd5<-window(seak.ts_odd, start=1952, end=2010)
test.seak_odd5<-window(seak.ts_odd, start=2011, end=2014)
seak_odd.final5<-forecast::auto.arima(train.seak_odd5, approximation = F, stepwise = F)

#create training and test datasets for the 10 year forecast 
train.seak_odd10<-window(seak.ts_odd, start=1952, end=2005)
test.seak_odd10<-window(seak.ts_odd, start=2006, end=2014)
seak_odd.final10<-forecast::auto.arima(train.seak_odd10, approximation = F, stepwise = F)

#create training and test datasets for the 20 year forecast 
train.seak_odd20<-window(seak.ts_odd, start=1952, end=1995)
test.seak_odd20<-window(seak.ts_odd, start=1996, end=2014)
seak_odd.final20<-forecast::auto.arima(train.seak_odd20, approximation = F, stepwise = F)

#Plots 

plot_even_5 <- seak_even.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak_even5))

plot_even_10 <- seak_even.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak_even10))

plot_even_20 <- seak_even.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak_even20))

plot_odd_5 <- seak_odd.final5 %>%
  forecast(h=5) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak_odd5))

plot_odd_10 <- seak_odd.final10 %>%
  forecast(h=10) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak_odd10))

plot_odd_20 <- seak_odd.final20 %>%
  forecast(h=20) %>%
  autoplot() + geom_point(aes(x=x, y=y), data=fortify(test.seak_odd20))

plot_grid(plot_even_5, plot_even_10, plot_even_20, plot_odd_5, plot_odd_10, plot_odd_20, ncol = 3, nrow = 2)

# Looking at Stationarity 

#Even Years 
Ndiff_even<-sapply(RegionBestMod_even, function(x){
  a<-strsplit(strsplit(strsplit(x, "[(]")[[1]][2], "[)]")[[1]][1],"[,]")
  return(a[[1]][2])}
)

tibble(Ndiff = Ndiff_even, region = Allcombs$regions, level = Allcombs$forecastlevels) %>%
  ggplot() + geom_bar(aes(x = region, y = Ndiff, fill = as.factor(level)), stat = "identity", position = "dodge") +
  scale_x_discrete(labels = as_labeller(regionskey)) +
  labs(fill = "Forecast Levels", x = "Region", y = "Number of Differences") + 
  ggtitle("Number of differences to achieve stationarity (Pink-Even Years)") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))

#Odd Years 
Ndiff_odd<-sapply(RegionBestMod_odd, function(x){
  a<-strsplit(strsplit(strsplit(x, "[(]")[[1]][2], "[)]")[[1]][1],"[,]")
  return(a[[1]][2])}
)

tibble(Ndiff = Ndiff_odd, region = Allcombs$regions, level = Allcombs$forecastlevels) %>%
  ggplot() + geom_bar(aes(x = region, y = Ndiff, fill = as.factor(level)), stat = "identity", position = "dodge") +
  scale_x_discrete(labels = as_labeller(regionskey)) +
  labs(fill = "Forecast Levels", x = "Region", y = "Number of Differences") + 
  ggtitle("Number of differences to achieve stationarity (Pink-Odd Years)") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))

#Residuals 

#Even
ac_mods_Pink_even<-c(19, 20) 

for(i in 1:length(ac_mods_Pink_even)){
print(paste(regionskey[ResultsTable_even$regions[ac_mods_Pink_even[i]]], ResultsTable_even$forecastlevels[ac_mods_Pink_even[i]]))
checkresiduals(RegionMods_even[[ac_mods_Pink_even[i]]]$Fit)
}

#Odd
ac_mods_Pink_odd<-c(8) 

for(i in 1:length(ac_mods_Pink_odd)){
  print(paste(regionskey[ResultsTable_odd$regions[ac_mods_Pink_odd[i]]], ResultsTable_odd$forecastlevels[ac_mods_Pink_odd[i]]))
  checkresiduals(RegionMods_odd[[ac_mods_Pink_odd[i]]]$Fit)
}

##Extra Unused Code ##################################

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

