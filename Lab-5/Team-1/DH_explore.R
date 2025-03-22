library(MARSS)
data(KvichakSockeye, package="atsalibrary")
SR_data <- KvichakSockeye

#create log ratio for yt
SR_data$y <- log(SR_data$recruits/SR_data$spawners)

#date to columns
dat <- t(SR_data)
colnames(dat) <- SR_data$brood_year
dat <- dat[-1,]
dat <- dat[c(5,1,3,4,2),]

mod1 <- list(
  Z = "identity",
  U="zero",
  R=matrix("r",nrow=1),
  B="identity",
  A="zero",
  Q=matrix("q",nrow=1)
)
m1 <- MARSS(dat[1,],model = mod1)
#plot alpha values
alpha1 <- as.numeric(m1$states)
alpha.se1 <- as.numeric(m1$states.se)
plot(alpha1~SR_data$brood_year,type='l', ylim=c(-3,3))
lines(alpha1+2*alpha.se1~SR_data$brood_year, lty="dashed")
lines(alpha1-2*alpha.se1~SR_data$brood_year, lty="dashed")

#AIC value
m1$AICc

#plot model residuals diagnostics. this simplified model should show signs of 
#autocorrelation in residuals
res <- residuals(m1)[,7]
tmp <- which(!is.na(res))
res <- res[tmp]
acf(res)
#or 
autoplot.marssMLE(m1)

######################
#                    #
#      PART 2        #
#                    #
######################

#first model we assume there is no density dependence
#we model the underlying states of alpha (brood-year productivity) and beta
#beta is prevented from changing.
#spawn_z <- matrix(scale(dat[2,]),nrow=1)  #z score spawners
#spawn_z <- spawn_z[,-c(1:4)]              #remove the NAs
spawn <- matrix(dat[2,-c(1:4)],nrow=1)    #not sure if we should scale this or not
ratio <- SR_data$y[-c(1:4)]               #remove responses that corresponded to Spawner = NA
# z-score the predictor variable
m <- 2                     #number of parameters in process 
TT <- length(spawn)      #number of data points

B <-  diag(m)                         #"identity"
U <-  matrix(0, nrow = m, ncol = 1)   #"zero"
Q <- matrix(list("q.alpha", 0, 0, 0), nrow = 2)  # to have characters and numbers in same matrix use a list
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2,] <- spawn            ## Nx1; predictor variable 
A = matrix(0)            #"zero"
R <-  matrix("r")

inits_list <- list(x0 = matrix(c(0, 0), nrow = m))     
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

#this model uses beta/spawners as a covariate rather than a time stable variable in the z matrix
mod3 <- list(
  Z=matrix(1,1),
  Q=matrix("q",1),
  R=matrix("r",1),
  U="zero",
  A="zero",
  B=matrix(1,1),
  D=matrix("d",1),
  d=spawn
)
init=list(x0=matrix(0,1))
m3 <- MARSS(ratio,mod3,init)
beta3 <- as.numeric(m3$coef)  

plot(dat[5,]~dat[2,])
## fit the model with MARSS this model uses beta/spawner as time invarying in the z matrix
m2 <- MARSS(ratio, mod_list, inits = inits_list)

hist(spawn)
alpha <- as.numeric(m2$states[1,])
beta <- as.numeric(m2$states[2,])
alpha.se <- as.numeric(m2$states.se[1,])
beta.se <- as.numeric(m2$states.se[2,])
plot(alpha~SR_data$brood_year[-c(1:4)],type='l', ylim=c(-3,3))
lines(alpha+2*alpha.se~SR_data$brood_year[-c(1:4)], lty="dashed")
lines(alpha-2*alpha.se~SR_data$brood_year[-c(1:4)], lty="dashed")


#from first model
plot(alpha1~SR_data$brood_year,type='l', ylim=c(-3,3), xlim= c(1950, 2005))
lines(alpha1+2*alpha.se1~SR_data$brood_year, lty="dashed")
lines(alpha1-2*alpha.se1~SR_data$brood_year, lty="dashed")
#second model
lines(alpha~SR_data$brood_year[-c(1:4)],type='l',col="red")
lines(alpha+2*alpha.se~SR_data$brood_year[-c(1:4)], lty="dashed", col="red")
lines(alpha-2*alpha.se~SR_data$brood_year[-c(1:4)], lty="dashed", col="red")

m1$states[-c(1:4)] #alpha from model 1 trimmed to same length as model 2 for comparison
m2$states[1,]      #alpha from model 2
m2$states[2,]     # beta from model 2

#the estimated values of alpha are very similar in each model. I was expecting the alphas
#from model 2 to be smaller than model 1.

#AIC value
m1$AICc
m2$AICc  #AICc is within a single point.


#plot model diagnostics. ACF still showing a lag at 5 and 10
autoplot.marssMLE(m2)

#############
#           #
#  Part 3   #
#############


#we model the underlying states of alpha (brood-year productivity) and beta

#spawn_z <- matrix(scale(dat[2,]),nrow=1)  #z score spawners
#spawn_z <- spawn_z[,-c(1:4)]              #remove the NAs
pdo_s <- matrix(dat[3,-c(1:4)], nrow=1)
spawn <- matrix(dat[2,-c(1:4)],nrow=1)    #not sure if we should scale this or not
            
# z-score the predictor variable
m <- 3                     #number of parameters in process 
TT <- length(spawn)      #number of data points

B <-  diag(m)                         #"identity"
U <-  matrix(0, nrow = m, ncol = 1)   #"zero"
Q <- matrix(list(0, 0,0,              #we no longer want alpha to vary
                 0, 0,0,              #hold beta steady as well
                 0,0,"q.pdo_s"), nrow = 3, byrow = TRUE)  # to have characters and numbers in same matrix use a list
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2,] <- spawn           ## Nx1; predictor variable spawners
Z[1,3,] <-  pdo_s            ## Nx1; predictor PDO summer
A = matrix(0)                #"zero"
R <-  matrix("r")

inits_list2 <- list(x0 = matrix(c(0,0,0), nrow = m))     
mod_list2 <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
cont_list <- list(maxit=10000)
# run model
m4 <- MARSS(ratio, mod_list2, inits = inits_list2)

alpha3 <- as.numeric(m4$states[1,])  #static number
beta3 <- as.numeric(m4$states[2,])   #should be static as well
pdo_s <- as.numeric(m4$states[3,])

beta.se3 <- as.numeric(m4$states.se[2,])
pdo_s.se3 <- as.numeric(m4$states.se[3,])


#from first model
plot(alpha1~SR_data$brood_year,type='l', ylim=c(-3,3), xlim= c(1950, 2005))
lines(alpha1+2*alpha.se1~SR_data$brood_year, lty="dashed")
lines(alpha1-2*alpha.se1~SR_data$brood_year, lty="dashed")
#second model
lines(alpha~SR_data$brood_year[-c(1:4)],type='l',col="red")
lines(alpha+2*alpha.se~SR_data$brood_year[-c(1:4)], lty="dashed", col="red")
lines(alpha-2*alpha.se~SR_data$brood_year[-c(1:4)], lty="dashed", col="red")
#third model
lines(alpha3~SR_data$brood_year[-c(1:4)],type='l', col="blue")


m1$states[-c(1:4)] #alpha from model 1 trimmed to same length as model 2 for comparison
m2$states[1,]      #alpha from model 2
m4$states[1,]     #alpha from model 3
m2$states[2,]     # beta from model 2
m4$states[2,]     #beta from model 3
m4$states[3,]    #estimate of pdo_s coeficient was allowed to vary in Q but was estimated as static

#AIC value
m1$AICc
m2$AICc  #AICc is within a single point.
m4$AICc  # getting worse

autoplot.marssMLE(m4)
#more autocorrelations in acf plot

#########model with PDO winter
pdo_w <- matrix(dat[4,-c(1:4)], nrow=1)
spawn <- matrix(dat[2,-c(1:4)],nrow=1)    #not sure if we should scale this or not

# z-score the predictor variable
m <- 3                     #number of parameters in process 
TT <- length(spawn)      #number of data points

B <-  diag(m)                         #"identity"
U <-  matrix(0, nrow = m, ncol = 1)   #"zero"
Q <- matrix(list(0, 0,0,              #we no longer want alpha to vary
                 0, 0,0,              #hold beta steady as well
                 0,0,"q.pdo_w"), nrow = 3, byrow = TRUE)  # to have characters and numbers in same matrix use a list
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2,] <- spawn           ## Nx1; predictor variable spawners
Z[1,3,] <-  pdo_w            ## Nx1; predictor PDO winter
A = matrix(0)                #"zero"
R <-  matrix("r")

inits_list2 <- list(x0 = matrix(c(0,0,0), nrow = m))     
mod_list2 <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

# run model
m5 <- MARSS(ratio, mod_list2, inits = inits_list2)

m1$AICc  #model 1: time varying alpha only
m2$AICc  #model 2: time varying alpha, time stable beta
m4$AICc  #model 3: time stable alpha and beta, time varying PDO_summer
m5$AICc  #model 4: time stable alpha and beta, time varying PDO_winter
