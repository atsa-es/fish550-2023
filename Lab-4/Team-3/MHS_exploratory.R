stoplight <- read.csv(here::here("Lab-4", "stoplight.csv"))

head(stoplight)
names(stoplight)

## add some code here
class(stoplight)
names(stoplight)

#ecosystem indicators
stoplight[,1]  
years<-seq(1998,2021,1)

#ecosystem indicators to pull
PDO_DecMarch<-unlist(as.vector(stoplight[1,3:26])) 
SST<-unlist(as.vector(stoplight[4,3:26])) 
Copepod_richness<-unlist(as.vector(stoplight[9,3:26])) 
Ichthy_community_index<-unlist(as.vector(stoplight[14,3:26]))

dat<-matrix(nrow=length(years), ncol=5)
dat[,1]<-years
dat[,2]<-PDO_DecMarch
dat[,3]<-SST
dat[,4]<-Copepod_richness
dat[,5]<-Ichthy_community_index
colnames(dat)<-c("years","PDO", "SST", "CopeRich", "FishInx")
rownames(dat)<-seq(1,24,1)

dat<-as.data.frame(dat)
print(dat)
#five times series in response and keep things to your block 

#EX code 

library(hmmTMB)

fitmod1 <- depmix(dat=surv_dat, 
                  logit.s~1,
                  nstates = 2,
                  transition = ~ 1,
                  family = gaussian())


#Playing around 
set.seed(123)

mod = depmix(list(CopeRich ~1, PDO~1, SST~1),
             nstates = 2, 
             family = list(gaussian(),gaussian(),gaussian()),
             data=dat)
             
fitmod = fit(mod)

summary(fitmod)

#Copepod richness as related to PDO 
mod_PDO = depmix(CopeRich ~ 1 + PDO, 
                  nstates = 2, 
                  transition = ~ 1, # the =~ means no covariates in transition
                  family = , #defaults to gaussian (response variable)
                  data=dat)

fitmod_PDO = fit(mod_PDO)

summary(fitmod_PDO)
test_PDO<-getpars(fitmod_PDO)

init_states1<-test_PDO[1:2]
mat1<-matrix(test_PDO[3:6],2,2,byrow=TRUE)

plot(ts(posterior(fitmod_PDO, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


#More likely to stay in state one by a little bit


#Copepod richness as realted to SST 
mod_SST = depmix(CopeRich ~ 1 + SST, 
                  nstates = 2, 
                  transition = ~ 1, # the =~ means no covariates in transition
                  family = , #defaults to gaussian (response variable)
                  data=dat)

fitmod_SST = fit(mod_SST)

summary(fitmod_SST)
test_SST<-getpars(fitmod_SST)

init_states2<-test_SST[1:2]
mat2<-matrix(test_SST[3:6],2,2,byrow=TRUE)

plot(ts(posterior(fitmod_SST, type="smoothing")[,1], start=c(2019,2), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

#State 2 is more stable than state 1 as it pertains to copepod richness 
#states are being informed by SST? 


#Play with Transition assumptions 

#Copepod richness as related to PDO 
mod_PDO_test = depmix(CopeRich ~ 1, 
                 nstates = 2, 
                 transition = ~PDO, # the =~ means no covariates in transition
                 family = , #defaults to gaussian (response variable)
                 data=dat)

fitmod_PDO_test = fit(mod_PDO_test)

summary(fitmod_PDO_test)
test_PDO2<-getpars(fitmod_PDO_test)

init_states1<-test_PDO2[1:2]
mat_PDO_test<-matrix(test_PDO2[3:6],2,2,byrow=TRUE)

plot(ts(posterior(fitmod_PDO_test, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


# Looking at inputs as a list 

set.seed(123)

multimod1 = depmix(list(CopeRich ~1, PDO~1, SST~1),
             nstates = 2, 
             family = list(gaussian(),gaussian(),gaussian()),
             data=dat)

fitmultimod = fit(multimod1)

summary(fitmultimod)











mod_test = depmix(CopeRich ~ 1
                 nstates = 2, 
                 transition = SST, # the =~ means no covariates in transition
                 family = , #defaults to gaussian (response variable)
                 data=dat)

fitmod_SST = fit(mod_SST)

summary(fitmod_SST)
test_SST<-getpars(fitmod_SST)

init_states2<-test_SST[1:2]
mat2<-matrix(test_SST[3:6],2,2,byrow=TRUE)



#Extreme example?
mod = depmix(list(SST ~1, PDO ~ 1), 
             nstates = 2, 
             family = list(gaussian(),gaussian()),
             data=calcofi)
fitmod = fit(mod)


#pick some plankton, pick some drivers, 




