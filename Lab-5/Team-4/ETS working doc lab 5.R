library(MARSS)
library(ggplot2)
## get data
data(KvichakSockeye, package="atsalibrary")
SR_data <- KvichakSockeye
plot.ts(x=SR_data$brood_year, SR_data$spawners)
lines(x=SR_data$brood_year,SR_data$recruits, col='blue')

#plot.ts(SR_data$recruits)

plot(x=SR_data$spawners, y=SR_data$recruits)

#look at 1956-1998
SR_56to98<-subset(SR_data, brood_year > 1955 & brood_year < 1999)

plot.ts(x=SR_56to98$brood_year, SR_56to98$spawners, ylim=c(min(SR_56to98$recruits), max(SR_56to98$recruits)))
lines(x=SR_56to98$brood_year,SR_56to98$recruits, col='blue')

#alpha is intercept and beta is slope
#log(Rt/St) = at - bSt + vt
#where at is a RW
#need to estimate a and b (as hidden states)
TT<-length(SR_56to98$brood_year)
SR_56to98$RtovSt<-log(SR_56to98$recruits/SR_56to98$spawners)
dat <- matrix(SR_56to98$RtovSt, nrow = 1)
dat.z<-zscore(dat)

plot(x=SR_56to98$brood_year, y=SR_56to98$RtovSt, type='l')

#covariate/predictor variable
#start with pdo summer
pdosum <- SR_56to98$pdo_summer_t2
## z-score the CUI
pdosum_z <- matrix(zscore(pdosum), nrow = 1)
## number of regr params (slope + intercept)
m <- dim(pdosum_z)[1] + 1

##### Estimate 1 regression parameter
B <- 'identity'  
U <- matrix(0, nrow = 1, ncol = 1)  ## 2x1; both elements = 0
Q <- matrix(list(0), 1, 1)  ## 2x2; all 0 for now
diag(Q) <- c("q.alpha")

#Miranda wants no covariates
MirandaZ<-"identity"
A <- matrix(0)  ## 1x1; scalar = 0
R <- matrix("r")  ## 1x1; scalar = r

## only need starting values for regr parameters
inits_list <- list(x0 = 0)

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = MirandaZ, A = A, R = R)

#DLM without covariates
dlm_1 <- MARSS(dat.z, inits = inits_list, model = mod_list)

plot(x=SR_56to98$brood_year, y=dlm_1$states, type='l')
lines(x=SR_56to98$brood_year, y=zscore(SR_56to98$RtovSt), col='blue')

dlm_1$AICc
#121.0344

autoplot(dlm_1)
#bit of autocorrelation at 5 and 10 years

#### Estimate 2 regression parameters
B <- diag(2)  ## 2x2; Identity
U <- matrix(0, nrow = 2, ncol = 1)  ## 2x1; both elements = 0
#Q <- matrix(list(0), 1, 1)
Q <- matrix(list(0), 2, 2)  ## 2x2; all 0 for now
diag(Q) <- c("q.alpha", "q.beta")

Z<-matrix(1, nrow = 1, ncol = 2)
A <- matrix(0)  ## 1x1; scalar = 0
R <- matrix("r")  ## 1x1; scalar = r

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = 2))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

#DLM without covariates
dlm_2 <- MARSS(dat.z, inits = inits_list, model = mod_list)
autoplot(dlm_2)
#AICc: 126.0491

######time-varying alpha; static beta
B <- diag(2)  ## 2x2; Identity
U <- matrix(0, nrow = 2, ncol = 1)  ## 2x1; both elements = 0

alpha <- log(SR_56to98$recruits)
## z-score
alpha_z <- matrix(zscore(alpha), nrow = 1)
Q <- array(NA, c(2, 2, TT))  ## NxMxT; empty for now
Q[1, 1, ] <- alpha_z
Q[1, 2, ] <- 0
Q[2, 1, ] <- 0
Q[2, 2, ] <- 1

Z<-matrix(1, nrow = 1, ncol = 2)
A <- matrix(0)  ## 1x1; scalar = 0
R <- matrix("r")  ## 1x1; scalar = r

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = 2))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

#DLM without covariates
dlm_3 <- MARSS(dat.z, inits = inits_list, model = mod_list, control=list(maxit=1000))
  