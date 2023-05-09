library(depmixS4)
library(tidyverse)
library(kableExtra)
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


#exploring the data more 
#rule 1 plot your data 

ggplot(dat, aes(x = years)) +
  geom_line(aes(y = PDO, color = "PDO"), size = 1) +
  geom_line(aes(y = SST, color = "SST"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green", FishInx = "purple")) +
  labs(x = "Year", y = "Value", color = "Variable") +
  theme_minimal()


# Plotting all
ggplot(dat, aes(x = years)) +
  geom_line(aes(y = PDO, color = "PDO"), size = 1) +
  geom_line(aes(y = SST, color = "SST"), size = 1) +
  geom_line(aes(y = CopeRich, color = "CopeRich"), size = 1) +
  geom_line(aes(y = FishInx, color = "FishInx"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green", FishInx = "purple")) +
  labs(x = "Year", y = "Value", color = "Variable") +
  theme_minimal()

# Plotting Coperich and env cov
ggplot(dat, aes(x = years)) +
  geom_line(aes(y = PDO, color = "PDO"), size = 1) +
  geom_line(aes(y = SST, color = "SST"), size = 1) +
  geom_line(aes(y = CopeRich, color = "CopeRich"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green", FishInx = "purple")) +
  labs(x = "Year", y = "Value", color = "Variable") +
  theme_minimal()

# Plotting Coperich and env cov
ggplot(dat, aes(x = years)) +
  geom_line(aes(y = PDO, color = "PDO"), size = 1) +
  geom_line(aes(y = CopeRich, color = "CopeRich"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green", FishInx = "purple")) +
  labs(x = "Year", y = "Value", color = "Variable") +
  theme_minimal()

# PDO states 
ggplot(dat, aes(x = years)) +
  geom_line(aes(y = PDO, color = "PDO"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green", FishInx = "purple")) +
  labs(x = "Year", y = "Value", color = "Variable") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 2.5, linetype = "dotted") +
  geom_hline(yintercept = -2.5, linetype = "dotted") +
  theme_minimal()

#Think of weak vs strong PDO?

#CopeRichness tracks very closely with PDO
#SST seems to have aa negligable corelation 

ggplot(dat, aes(x = PDO, y = CopeRich, color = CopeRich < 0)) +
  geom_point() +
  geom_line(aes(x = PDO, y = CopeRich)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("red", "blue"), labels = c("Above 0", "Below 0")) +
  labs(x = "PDO", y = "CopeRich") +
  theme_minimal()
  
#fish index all over the place, let's not use it. 

#Playing around 
set.seed(123)

#Copepods with PDO

mod = depmix(list(CopeRich ~1, PDO~1),
             nstates = 2, 
             transition = ~PDO, # the =~ means no covariates in transition?
             family = list(gaussian(),gaussian()),
             data=dat)
             
#Nstates=how many possible states does the system have 
#Transition= how quickly you jump between states
  # 1= instantanious (one time step)
  #0.5 = slower -- slope of the transition
  #Put in PDO for the transition, the data from PDO is mediating the transition
  #between copepod richness state 

#AIC more states = better 
  #need a biological rational for the number of states you have 


fitmod = fit(mod)

summary(fitmod)

fit_pars<-getpars(fitmod) #understand these parameters more.

init_states1<-fit_pars[1:2]
mat1<-matrix(fit_pars[3:6],2,2,byrow=TRUE)

plot(ts(posterior(fitmod, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


#=== model options
#=======
#   1
#=======
#Just copepods
mod1 = depmix(CopeRich ~1,
                   nstates = 2, 
                   transition = ~1, 
                   family = gaussian(),
                   data=dat)

fitmod1 = fit(mod1)
summary(fitmod1)

#Find best fit 
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best1 <- 1e10
best_model1 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod1 <- fit(mod1)
  if(AIC(fitmod1)< best1){
    best_model1 <- fitmod1
    best1 <- AIC(fitmod1)
  }
}

plot(ts(posterior(best_model1, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

prstates<-apply(posterior(best_model1)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model1)[,1]
#mu2<-summary(fitmod1)[,3]
pred <- tibble("Year" = dat$years, "Copepod Richness" = mu[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Copepod Richness", values_to = "Fit")
ggplot() + geom_point(data = dat, aes(x = years, y = CopeRich, color = "Copepod Richness")) + 
  geom_line(data = pred2, aes(x = years, y = Fit, group = "Copepod Richness", color = "Copepod Richness")) + theme_classic()

#=======
#   2
#=======
#Copepods where PDO informs the transition 
mod2 = depmix(CopeRich ~1,
              nstates = 2, 
              transition = ~PDO, 
              family = gaussian(),
              data=dat)

fitmod2 = fit(mod2)
summary(fitmod2)

#Find best fit 
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best2 <- 1e10
best_model2 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod2 <- fit(mod2)
  if(AIC(fitmod2)< best2){
    best_model2 <- fitmod2
    best2 <- AIC(fitmod2)
  }
}

plot(ts(posterior(best_model2, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

prstates<-apply(posterior(best_model2)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model2)[,1]
#mu2<-summary(fitmod2)[,3]
pred <- tibble("Year" = dat$years, "Copepod Richness" = mu[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Copepod Richness", values_to = "Fit")
ggplot() + geom_point(data = dat, aes(x = years, y = CopeRich, color = "Copepod Richness")) + 
  geom_line(data = pred2, aes(x = years, y = Fit, group = "Copepod Richness", color = "Copepod Richness")) + theme_classic()

#=======
#   3
#=======
#Copepods and PDO as intercepts 
mod3 = depmix(list(CopeRich ~1, PDO~1),
             nstates = 2, 
             transition = ~1, 
             family = list(gaussian(),gaussian()),
             data=dat)

fitmod3 = fit(mod3)
summary(fitmod3)

#Find best fit 
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best3 <- 1e10
best_model3 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod3 <- fit(mod3)
  if(AIC(fitmod3)< best3){
    best_model3 <- fitmod3
    best3 <- AIC(fitmod3)
  }
}

plot(ts(posterior(best_model3, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

prstates<-apply(posterior(best_model3)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model3)[,1]
mu2<-summary(best_model3)[,3]
pred <- tibble("Year" = dat$years, "Copepod Richness" = mu[prstates],"PDO" = mu2[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Covariates", values_to = "Fit")
ggplot() +
  geom_point(data = dat, aes(x = years, y = CopeRich, color = "Copepod Richness")) +
  geom_point(data = dat, aes(x = years, y = PDO, color = "PDO")) +
  geom_line(data = pred2, aes(x = Year, y = Fit, group = Covariates, color = Covariates)) +
  theme_classic()

#=======
#   4
#=======
#Copepods and PDO as intercepts where the transition is informed by PDO
mod4 = depmix(list(CopeRich ~1, PDO~1),
              nstates = 2, 
              transition = ~PDO, 
              family = list(gaussian(),gaussian()),
              data=dat)
fitmod4 = fit(mod4)
summary(fitmod4)

#Find best fit 
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best4 <- 1e10
best_model4 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod4 <- fit(mod4)
  if(AIC(fitmod4)< best4){
    best_model4 <- fitmod4
    best4 <- AIC(fitmod4)
  }
}

plot(ts(posterior(best_model4, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

prstates<-apply(posterior(best_model4)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model4)[,1]
#mu2<-summary(fitmod2)[,3]
pred <- tibble("Year" = dat$years, "Copepod Richness" = mu[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Copepod Richness", values_to = "Fit")
ggplot() + geom_point(data = dat, aes(x = years, y = CopeRich, color = "Copepod Richness")) + 
  geom_line(data = pred2, aes(x = years, y = Fit, group = "Copepod Richness", color = "Copepod Richness")) + theme_classic()

#=======
#   5
#=======
#Copepods as they relate to PDO
mod5 = depmix(list(CopeRich ~PDO),
              nstates = 2, 
              transition = ~1, 
              family = list(gaussian()),
              data=dat)
fitmod5 = fit(mod5)
summary(fitmod5)

#Find best fit 
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best5 <- 1e10
best_model5 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod5 <- fit(mod5)
  if(AIC(fitmod5)< best5){
    best_model5 <- fitmod5
    best5 <- AIC(fitmod5)
  }
}

plot(ts(posterior(best_model5, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

prstates<-apply(posterior(best_model5)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model5)[,1]
#mu2<-summary(fitmod2)[,3]
pred <- tibble("Year" = dat$years, "Copepod Richness" = mu[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Copepod Richness", values_to = "Fit")
ggplot() + geom_point(data = dat, aes(x = years, y = CopeRich, color = "Copepod Richness")) + 
  geom_line(data = pred2, aes(x = years, y = Fit, group = "Copepod Richness", color = "Copepod Richness")) + theme_classic()

#=======
#   6
#=======
#Copepods as they relate to PDO where the transition is informed by PDO
mod6 = depmix(list(CopeRich ~PDO),
              nstates = 2, 
              transition = ~PDO, 
              family = list(gaussian()),
              data=dat)
fitmod6 = fit(mod6)
summary(fitmod6)

#Find best fit 
iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best6 <- 1e10
best_model6 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod6 <- fit(mod6)
  if(AIC(fitmod6)< best6){
    best_model6 <- fitmod6
    best6 <- AIC(fitmod6)
  }
}

plot(ts(posterior(best_model6, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

prstates<-apply(posterior(best_model6)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model6)[,1]
#mu2<-summary(fitmod2)[,3]
pred <- tibble("Year" = dat$years, "Copepod Richness" = mu[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Copepod Richness", values_to = "Fit")
ggplot() + geom_point(data = dat, aes(x = years, y = CopeRich, color = "Copepod Richness")) + 
  geom_line(data = pred2, aes(x = years, y = Fit, group = "Copepod Richness", color = "Copepod Richness")) + theme_classic()

# Looking for best model
######################################################
#AIC

mods <- c("1","2","3","4","5","6")
aic <- c(AIC(best_model1), AIC(best_model2), AIC(best_model3), AIC(best_model4), AIC(best_model5),AIC(best_model6))
daic <- aic-min(aic)
tab <- cbind.data.frame(mods, aic, daic)





#======OLD CODE ==============
#Two state Assumptions

multimod2 = depmix(list(CopeRich ~1, PDO~1),
             nstates = 2, 
             transition = ~PDO, 
             family = list(gaussian(),gaussian()),
             data=dat)

fitmod = fit(multimod2)

iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best <- 1e10
best_model <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod <- fit(multimod2)
  if(AIC(fitmod)< best){
    best_model <- fitmod
    best <- AIC(fitmod)
  }
}

summary(best_model)

fit@response
fit@transition

test<-getpars(best_model)

init_states1<-test[1:2]
mat_test<-matrix(test[3:6],2,2,byrow=TRUE)

plot(ts(posterior(best_model, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)

#=========================
#Two states modded
#========================
#Two state Assumptions

multimod2 = depmix(CopeRich~PDO,
                   nstates = 2, 
                   transition = ~1, 
                   family = gaussian(),
                   data=dat)

fitmod = fit(multimod2)

iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best <- 1e10
best_model <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod <- fit(multimod2)
  if(AIC(fitmod)< best){
    best_model <- fitmod
    best <- AIC(fitmod)
  }
}

summary(best_model)

test<-getpars(best_model)

fit@response
fit@transition

init_states1<-test[1:2]
mat_test<-matrix(test[3:6],2,2,byrow=TRUE)

plot(ts(posterior(best_model, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


#===================================
#Alt 2  state Assumptions
#===================================

multimod3 = depmix(CopeRich ~1,
                   nstates = 2, 
                   transition = ~PDO, 
                   family = gaussian(),
                   data=dat)

fitmod3 = fit(multimod3)

fitmod3@response
fitmod3@transition

iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best3 <- 1e10
best_model <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod <- fit(multimod3)
  if(AIC(fitmod3)< best3){
    best_model <- fitmod3
    best <- AIC(fitmod3)
  }
}

summary(best_model)

test<-getpars(best_model)

init_states1<-test[1:2]
mat_test<-matrix(test[3:6],2,2,byrow=TRUE)

plot(ts(posterior(best_model, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


#Plot fits to data? 
#plotting

#FROM ZOE RAND
prstates<-apply(posterior(fitmod)[,c("S1", "S2")],1, which.max)
mu<-summary(best_model)[,1]
mu2<-summary(best_model)[,3]
mu3<-summary(best_model)[,5]
pred <- tibble("Year" = df2$Year, "Copepod_richness" = mu[prstates], "N_copepod" = mu2[prstates], "S_copepod" = mu3[prstates])
pred2<-pred %>% pivot_longer(-Year, names_to = "Ecosystem.Indicators", values_to = "Fit")
ggplot() + geom_point(data = newdf, aes(x = Year, y = value, color = Ecosystem.Indicators)) + 
  geom_line(data = pred2, aes(x = Year, y = Fit, group = Ecosystem.Indicators, color = Ecosystem.Indicators)) + theme_classic()


#States overlayed on data (plots from above)
#Diagnostics 
#comparative tables for 2 v 3 states 








