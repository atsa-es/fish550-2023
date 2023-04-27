#Zoe Rand Working Code
#4/26/23
library(MARSS)
library(tidyverse)
library(zoo)

## load the raw data (there are 3 datasets contained here)
data(lakeWAplankton, package = "MARSS")

## we want `lakeWAplanktonTrans`, which has been transformed
## so the 0's are replaced with NA's and the data z-scored
all_dat <- lakeWAplanktonTrans

head(all_dat)
cn<-colnames(all_dat)

#PLOTTING 

all_dat %>% as_tibble() %>%
  mutate(Date = as.yearmon(paste(Year, Month), "%Y %m")) %>%
  select(c("Date", "Cyclops", "Diaptomus", "Epischura", "Neomysis")) %>% 
  pivot_longer(-c(Date), names_to = "Species", values_to = "Vals") %>%
  ggplot(aes(x = Date, y = Vals, color = Species)) + geom_point() + 
  geom_line() + facet_wrap(~Species) + scale_x_yearmon(format = "%Y") + 
  theme_minimal() + labs(y = "Abundance Index", x = "Year") + theme(legend.position =  "none")

#Neomysis is very patchy so potentially could cause problems? 
#using data up to 1985
spp<-c("Cyclops", "Diaptomus", "Epischura", "Neomysis")
dat_1985<-all_dat[all_dat[, "Year"] <= 1985,]
dat_zoo_1985<-dat_1985[, spp]
head(dat_zoo_1985)
dat_zoo<-t(dat_zoo_1985)

## get number of time series
N_ts <- dim(dat_zoo)[1]

## get length of time series
TT <- dim(dat_zoo)[2] 

## mean of each taxon
y_bar <- apply(dat_zoo, 1, mean, na.rm = TRUE)

## subtract the means
dat <- dat_zoo - y_bar

## assign new column names
rownames(dat) <- spp

head(dat)


#for covariates
#temp, pH, TP, cryptomonas, diatoms, greens, unicells, other.algae 
covs<-colnames(dat_1985)[c(3:8,10:11)]
cov_dat<-t(dat_1985[,covs])
head(cov_dat)


pairs(dat_1985[,covs])
library(corrplot)
corrplot(cor(dat_1985[,covs], use = "pairwise.complete.obs"), addCoef.col = "blue")
