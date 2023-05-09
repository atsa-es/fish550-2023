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


# Looking for best model
######################################################

#Two state Assumptions

multimod2 <- depmix(list(CopeRich ~1, PDO~1),
             nstates = 2, 
             transition = ~PDO, 
             family = list(gaussian(),gaussian()),
             data=dat)

# for 2 states
fitmod2 <- fit(multimod2)

iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best2 <- 1e10
best_model2 <- NA

for(i in 1:iter){
  set.seed(seeds[i])
  fitmod2 <- fit(multimod2)
  if(AIC(fitmod2)< best2){
    best_model2 <- fitmod2
    best2 <- AIC(fitmod2)
  }
}

summary(best_model2)

test<-getpars(best_model2)

init_states2<-test[1:2]
mat_test2<-matrix(test[3:6],2,2,byrow=TRUE)

plot(ts(posterior(best_model2, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


#Plot model output
#slide 50

# get the states
prstates2 = apply(posterior(fitmod2)[,c("S1","S2")], 
                 1, which.max)
#plot the states
plot(prstates2, type="b", xlab="Time", ylab="State", main = "Two States")

#plot the states against data
# get the means
mu2 <- summary(fitmod2)[,1]
# then predict estimated data
pred2 <- data.frame("year"=seq(min(years), 
                              max(years)),
                   "fit" = mu2[prstates2])


ggplot(pred2, aes(x = years))+
  geom_line(aes( y=fit, color = "predicted_states"), size = 1)+
  geom_line(aes(y = PDO_DecMarch, color = "PDO" ), size = 1)+
  geom_line(aes(y = Copepod_richness, color = "CopeRich"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green",
                                FishInx = "purple", predicted_states = "orange")) +
  labs(title = "Two States", x = "Year", y = "Value", color = "Variable") +
  theme_minimal()

#===================================
#Three state Assumptions
#===================================

####3####

multimod3 = depmix(list(CopeRich ~1, PDO~1),
                   nstates = 3, 
                   transition = ~PDO, 
                   family = list(gaussian(),gaussian()),
                   data=dat)

fitmod3 = fit(multimod3)

iter <-100
seeds<-sample(100:1000, size = iter, replace = F)
best3 <- 1e10
best_model3 <- NA
for(i in 1:iter){
  set.seed(seeds[i])
  fitmod3 <- fit(multimod3)
  if(AIC(fitmod3)< best3){
    best_model3 <- fitmod3
    best3 <- AIC(fitmod3)
  }
}

summary(best_model3)

test<-getpars(best_model3)

init_states3<-test[1:2]
mat_test3<-matrix(test[3:6],2,2,byrow=TRUE)

plot(ts(posterior(best_model3, type="smoothing")[,1], start=c(2003,5), deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (copepods good?).",
     frame=FALSE)


#Plot model output
#slide 50

# get the states
prstates3 <- apply(posterior(fitmod)[,c("S1","S2")], 
                 1, which.max)
#plot the states
plot(prstates3, type="b", xlab="Time", ylab="State", main = "Three States")

#create estimated data from the model

# get the means
mu3 <- summary(fitmod3)[,1]
# then predict estimated data
pred3 <- data.frame("year"=seq(min(years), 
                             max(years)),
                  "fit" = mu3[prstates3])

# Three states does not look good, it only shows two and they 
# only show up at the upper range of the data ?
ggplot(pred, aes(x = years))+
  geom_line(data=pred, aes( y=fit, color = "predicted_states"), size = 1)+
  geom_line(aes(y = PDO_DecMarch, color = "PDO" ), size = 1)+
  geom_line(aes(y = Copepod_richness, color = "CopeRich"), size = 1) +
  scale_color_manual(values = c(PDO = "blue", SST = "red", CopeRich = "green",
                                FishInx = "purple", predicted_states = "orange")) +
  labs(title = "Three States", x = "Year", y = "Value", color = "Variable") +
  theme_minimal()

####3####

#Diagnostics using two state model from above

summary(fitmod2)

#these give same values as summary, but split up?
fitmod2@response #gives the 4 intercepts with sd
fitmod2@transition # gives coefficients and transition probabilities?

#comparative tables for 2 v 3 states 


convert_probs <- function(fit, coef_value = 0) {
  # extract coefficients
  p1 <- fit@transition[[1]]@parameters$coefficients
  p2 <- fit@transition[[2]]@parameters$coefficients
  
  # assemble transition matrix in mlogit space
  m <- matrix(0, 2, 2)
  m[1,] <- p1[1,] + p1[2,] * coef_value #intercept + slope * coefficient
  m[2,] <- p2[1,] + p2[2,] * coef_value
  
  # exponentiate
  m <- exp(m)
  
  # normalize
  m[1,] <- m[1,] / sum(m[1,])
  m[2,] <- m[2,] / sum(m[2,])
  return(m)
}

#returns the transition matrix
convert_probs(fit = fitmod2)



