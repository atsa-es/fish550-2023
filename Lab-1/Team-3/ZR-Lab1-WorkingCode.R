#Zoe's Work For Lab 1
#4/6/22
library(tidyverse)
library(forecast)
library(here)

# Reading In Data ---------------------------------------------------------
ruggerone_data <- readRDS(here::here("Lab-1", "Data_Images", "ruggerone_data.rds"))

#Exploratory plots
ruggerone_data %>% 
  filter(species=="chum") %>% 
  ggplot(aes(x=year, y=log(returns))) + 
  geom_line() + 
  ggtitle("chum salmon log abundance by region") +
  facet_wrap(~region)

ruggerone_data %>% 
  group_by(species, year) %>%
  summarize(total = sum(returns, na.rm=TRUE)) %>%
  ggplot(aes(x=year, y=log(total))) + 
  geom_line() + 
  ggtitle("log abundance by species") +
  facet_wrap(~species)



# Analysis ----------------------------------------------------------------


#total chum returns over year
ChumByYear<-ruggerone_data %>% 
  group_by(species, year) %>%
  summarize(total = sum(returns, na.rm=TRUE)) %>% 
  mutate(lnreturns = log(total)) %>%
  filter(species == "chum")

ChumByYear$year[1] #1952
ChumByYear$year[nrow(ChumByYear)] #2015

#create timeseries object
datts <- ts(ChumByYear$lnreturns, start=ChumByYear$year[1]) #this assumes the first year in data is the start of the time series (they are in order) 
train <- window(datts, 1952, 2010)
test <- window(datts, 2011, 2015)


#acf
acf(datts) #significant autocorrelation up to lag 4ish
pacf(datts) #significant up to lag 2



#ARIMA model
mod <- auto.arima(train, trace = TRUE)
mod #picked a best model but a lot are within 2 AIC units of each other
plot(forecast(mod)) 
points(test) #just visually seems like there's an increase in the last few years that are not caught by the model




# Update 4/12 Region and Forecast Level Differences -----------------------
ruggerone_data %>% 
  filter(species=="chum") %>% 
  ggplot(aes(x=year, y=log(returns))) + 
  geom_line() + 
  ggtitle("chum salmon log abundance by region") +
  facet_wrap(~region)


#removing Korea Japan because there's no data
ChumByRegion<-ruggerone_data %>%
  filter(region != "japan") %>%
  filter(region != "korea") %>%
  group_by(species, region, year) %>%
  summarize(total = sum(returns, na.rm=TRUE)) %>% 
  mutate(lnreturns = log(total)) %>%
  filter(species == "chum")
head(ChumByRegion)

#making sure all the regions cover all the years (or at least start and end)
ChumByRegion %>% group_by(region) %>% summarise(startyear = min(year), endyear = max(year))
#all start in 1952 and end in 2015

#regions vector
regions<-unique(ChumByRegion$region)
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
  Chumdat<-ChumByRegion %>% filter(region == reg)
  #create time series
  datts <- ts(Chumdat$lnreturns, start=Chumdat$year[1])
  return(list(a = acf(datts), p = pacf(datts)))
}


Chumdat<-ChumByRegion %>% filter(region == reg)
#create time series
datts <- ts(Chumdat$lnreturns, start=Chumdat$year[1])
trace <- capture.output({
  # assign so it doesn't pollute the output
  model <- auto.arima(datts, trace = TRUE)
})
con    <- textConnection(trace)
models <- read.table(con, sep=":")
close(con)

BestMods<-models%>% filter(row_number() != nrow(models)) %>% mutate(AIC = replace(V2, V2 == "Inf", 99999), AIC = as.numeric(AIC), DeltaAIC = AIC-min(AIC)) %>% filter(DeltaAIC <= 2.0)
for(i in 1:nrow(BestMods)){
  BestMods$Mod[i]<-strsplit(strsplit(strsplit(BestMods$V1[i], "[(]")[[1]][2], "[)]")[[1]][1],"[,]")
  BestMods$npar[i]<-sum(as.numeric(BestMods$Mod[i][[1]][c(1,3)]))
  if(strsplit(strsplit(BestMods$V1[i], "[(]")[[1]][2], "[)]")[[1]][2] == " with drift         "){
    BestMods$npar[i] = BestMods$npar[i] + 1
  }
}

BestMods$Mod
BestMods$npar

#function for ARIMA models
FitModFunction<-function(reg, forelevel){
  #filter region
  Chumdat<-ChumByRegion %>% filter(region == reg)
  #create time series
  datts <- ts(Chumdat$lnreturns, start=Chumdat$year[1]) #this assumes the first year in data is the start of the time series (they are in order) 
  cutoff<-2015-forelevel
  train <- window(datts, 1952, cutoff)
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
  ggtitle("Chum") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))



#get stationarity results

Ndiff<-sapply(RegionBestMod, function(x){
  a<-strsplit(strsplit(strsplit(x, "[(]")[[1]][2], "[)]")[[1]][1],"[,]")
  return(a[[1]][2])}
)

tibble(Ndiff = Ndiff, region = Allcombs$regions, level = Allcombs$forecastlevels) %>%
  ggplot() + geom_bar(aes(x = region, y = Ndiff, fill = as.factor(level)), stat = "identity", position = "dodge") +
  scale_x_discrete(labels = as_labeller(regionskey)) +
  labs(fill = "Forecast Levels", x = "Region", y = "Number of Differences") + 
  ggtitle("Number of differences to achieve stationarity (Chum)") + theme_bw() + theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.5 ))


#checking residuals
for(i in 1:36){
  print(i)
  print(checkresiduals(RegionModsChum[[i]]$Fit))
}
