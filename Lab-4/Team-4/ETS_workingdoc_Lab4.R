library(depmixS4)
library(ggplot2)

stoplight <- read.csv(here::here("Lab-4", "stoplight.csv"))

#we will focus on copepod richness, N cop. biomass, S cop. biomass
#rows 9-11
cop.dat<-stoplight[9:11,]
cop.row<-cop.dat[,1]
cop.dat<-cop.dat[,3:26]
rownames(cop.dat)<-cop.row

#look at the data over time
cop.t<-t(cop.dat)
rows <- gsub("X","", rownames(cop.t))
row.names(cop.t)<-rows
cop.df<-data.frame(cop.t)
cop.df$Year<-row.names(cop.df)

plot(x=cop.df$Year, y=cop.df$Copepod_richness, type='l')
lines(x=cop.df$Year, y=cop.df$N_copepod, col='blue')
lines(x=cop.df$Year, y=cop.df$S_copepod, col='magenta')

#building the 2 state model
set.seed(1)
mod.2st<-depmix(list(Copepod_richness ~ 1, N_copepod ~ 1, S_copepod ~ 1),
                nstates = 2,
                family=list(gaussian(), gaussian(),gaussian()),
                data=cop.df)
fitmod.2st <- fit(mod.2st)
#converged at iteration 8 with logLik: -36.64126 

#most probable states
prstates2 <- apply(posterior(fitmod.2st) [,c("S1", "S2")], 1, which.max)
plot(prstates2, type='b', xlab='Time', ylab='State')

#estimated data
mu2 <- summary(fitmod.2st)[,1]
pred2 <- data.frame("year"=seq(min(cop.df$Year), max(cop.df$Year)), "fit" = mu2[prstates2])

cop.df$Year<-as.numeric(cop.df$Year)

ggplot(cop.df) + 
  geom_point(aes(x=Year, y=Copepod_richness), col="black", size=1.5) + 
  geom_point(aes(x=Year, y=N_copepod), col="blue", size=1.5) +
  geom_point(aes(x=Year, y=S_copepod), col="magenta", size=1.5) +
  geom_line(data=pred2, aes(year, fit)) + 
  ggtitle("Copepods - raw data and predictions") + 
  ylab("Richness or Biomass") + xlab("Year") + theme_bw()

#3 state model
set.seed(314)
mod.3st<-depmix(list(Copepod_richness ~ 1, N_copepod ~ 1, S_copepod ~ 1),
                nstates = 3,
                family=list(gaussian(), gaussian(),gaussian()),
                data=cop.df)
fitmod.3st <- fit(mod.3st)
#converged at iteration 14 with logLik: -14.21542

#most probable states
prstates3 <- apply(posterior(fitmod.3st) [,c("S1", "S2")], 1, which.max)
plot(prstates3, type='b', xlab='Time', ylab='State')

#estimated data
mu3 <- summary(fitmod.3st)[,1]
pred3 <- data.frame("year"=seq(min(cop.df$Year), max(cop.df$Year)), "fit" = mu3[prstates3])

ggplot(cop.df) + 
  geom_point(aes(x=Year, y=Copepod_richness), col="black", size=1.5) + 
  geom_point(aes(x=Year, y=N_copepod), col="blue", size=1.5) +
  geom_point(aes(x=Year, y=S_copepod), col="magenta", size=1.5) +
  geom_line(data=pred3, aes(year, fit)) + 
  ggtitle("Copepods - raw data and predictions") + 
  ylab("Richness or Biomass") + xlab("Year") + theme_bw()
