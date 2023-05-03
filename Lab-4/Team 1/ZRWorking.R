stoplight <- read.csv(here::here("Lab-4", "stoplight.csv"))
df <- rbind.data.frame(stoplight[stoplight$Ecosystem.Indicators == "Copepod_richness",],stoplight[stoplight$Ecosystem.Indicators == "N_copepod",],stoplight[stoplight$Ecosystem.Indicators == "S_copepod",])
head(df)
library(tidyverse)
newdf<-df%>%pivot_longer(-c(Ecosystem.Indicators, Type), names_to = "Year") %>%
  mutate(Year = as.numeric(gsub("X","",Year)))
head(newdf)

ggplot(newdf) + geom_line(aes(x = Year, y = value, color = Ecosystem.Indicators)) + theme_classic()

df2<-df %>% select(-Type) %>%
  gather(Year, value, -Ecosystem.Indicators) %>% 
  spread(Ecosystem.Indicators, value) 

df2$Year<-as.numeric(gsub("X","",df2$Year))
df2

set.seed(123)
library(depmixS4)

mod<-depmix(list(Copepod_richness ~1, N_copepod ~1, S_copepod ~1), 
           nstates = 2, 
           family = list(gaussian(), gaussian(), gaussian()), 
           data = df2)
fitmod<-fit(mod)
summary(fitmod)

#making sure it gets the best value
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best <- 1e10
best_model <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod <- fit(mod)
  if(AIC(fitmod)< best){
    best_model <- fitmod
    best <- AIC(fitmod)
  }
}

summary(best_model)

#plotting
prstates<-apply(posterior(fitmod)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model)[,1]
pred <- tibble("Year" = df2$Year, "Fit" = mu[prstates])
ggplot() + geom_point(data = newdf, aes(x = Year, y = value, color = Ecosystem.Indicators)) + 
  geom_line(data = pred, aes(x = Year, y = Fit)) + theme_classic()
