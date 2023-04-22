#Working doc for Emma Lab 2
#Lower Columbia River Chinook
#I am in charge of question 2 and finalizing the lab report
library(tidyverse)
library(MARSS)
library(reshape)

load(here::here("Lab-2", "Data_Images", "columbia-river.rda"))

esu <- unique(columbia.river$esu_dps)

plotesu <- function(esuname){
  df <- columbia.river %>% subset(esu_dps %in% esuname)
  ggplot(df, aes(x=spawningyear, y=log(value), color=majorpopgroup)) + 
    geom_point(size=0.2, na.rm = TRUE) + 
    theme(strip.text.x = element_text(size = 3)) +
    theme(axis.text.x = element_text(size = 5, angle = 90)) +
    facet_wrap(~esapopname) +
    ggtitle(paste0(esuname, collapse="\n"))
}
#Chinook
plotesu(esu[5])



df <- columbia.river %>% subset(species == "Chinook salmon")
ggplot(df, aes(x=spawningyear, y=log(value), color=run)) + 
  geom_point(size=0.2, na.rm = TRUE) +
  theme(strip.text.x = element_text(size = 3)) +
  theme(axis.text.x = element_text(size = 5, angle = 90)) + 
  facet_wrap(~esapopname)

#check start and end date for each esapopname for Chinook
min_vals<-aggregate(spawningyear~esapopname, data=df, FUN=min)
max_vals<-aggregate(spawningyear~esapopname, data=df, FUN=max)
#data go from 1964-2021 but there are a lot of missing values, designated NA

#count NA values for each esapopname
na_counts<-aggregate(value~esapopname, data=df, FUN=function(x) sum(is.na(x)))
#this does not work. I just get a table with a 0 count for each esapopname.

#adjust esapopname so that it does not include Salmon, Chinook (Lower Columbia River ESU)
#df$esapopname<-gsub('Salmon, Chinook (Lower Columbia River ESU)', "", df$esapopname)
#didn't work

#Create estimates of spawner abundance for all missing years and provide estimates of the decline from the historical abundance.
#From the user manual: Section 3.4.7, discusses at length that the standard missing values correction leads to an inexact likelihood when there are missing values.
#see page 21: The {MARSS} package uses instead the exact likelihood correction for missing values

#SEE CHAPTER 11 IN ATSA FOR PREDICTING MISSING VALUES

#variables: major pop group; esa pop name; common pop name; run
#we could assume correlation at any of these levels
unique(df$majorpopgroup) #6 levels
unique(df$esapopname) #17 levels
unique(df$commonpopname) #17 levels
unique(df$run) #3 levels

#hypothesis 1a: there are 3 independent populations (fall, late fall, spring)
ggplot(df, aes(x=spawningyear, y=log(value), color=run)) + 
  geom_point(size=0.2, na.rm = TRUE) +
  theme(strip.text.x = element_text(size = 3)) +
  theme(axis.text.x = element_text(size = 5, angle = 90)) + 
  facet_wrap(~run)

#format data for MARSS
#time needs to be column headers
run.dat<-aggregate(log(value)~run + spawningyear, data=df, FUN=sum) #na.rm = TRUE
#replace the 0s with NAs
colnames(run.dat)<-c('run', 'spawningyear', 'logvalue')
#replace -inf with 0
run.dat$logvalue[!is.finite(run.dat$logvalue)]<-0
run.wide<-pivot_wider(run.dat, names_from=run, values_from=logvalue)

dat<-t(run.wide)
years<-dat[1,]
dat<-dat[2:nrow(dat),]
colnames(dat)<-years

U.model <- 'unequal' #diff pop growth rates
Q.model<- 'diagonal and unequal' #diff process errors
R.model <- 'diagonal and equal' #same observation error variance
Z.model <- matrix(1,3,1)

mod.list.1a <- list( Q = Q.model, Z=Z.model, R=R.model, U=U.model, 
                    A = "scaling",  
                   x0 = matrix("mu"), tinitx = 0)
fit.1a <- MARSS(dat, model = mod.list.1a)

autoplot(fit.1a)

#hypothesis 1b: process errors are correlated
U.model <- 'unequal'
Q.model<- 'unequal'
R.model <- 'diagonal and equal'
Z.model
#hypothesis 1c: same process underlying everything
U.model <- 'unequal'
Q.model<- 'diagonal and equal'
R.model <- 'diagonal and equal'
Z.model
#hypothesis 2: there are 6 independent populations (major pop group)
#hypothesis 3: there are 17 independent populations (esa group)

#Evaluate support for the major population groups. Are the populations in the groups more correlated than outside the groups?

#Evaluate the evidence of cycling in the data. We will talk about how to do this on the Tuesday after lab.
#in lab, include covariate of seasonality/cycling - looking for evidence that there is some kind of cycling (for Chinook should be ~ 5 year cycle based on generation time); ACF on raw data will pick up natural cycle to be included as covariate

#You can assume that R="diagonal and equal" and A="scaling". Assume that “historical” means the earliest years available for your group.
