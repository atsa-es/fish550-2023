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
