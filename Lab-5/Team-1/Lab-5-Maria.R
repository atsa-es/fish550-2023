library(here)
library(devtools)
library(MARSS)


devtools::install_github("nwfsc-timeseries/atsalibrary")

## get data
data(KvichakSockeye, package="atsalibrary")
SR_data <- KvichakSockeye


#create log ratio for yt
SR_data$y <- log(SR_data$recruits/SR_data$spawners)

#date to columns
dat <- t(SR_data)
colnames(dat) <- SR_data$brood_year
dat <- dat[-1,]
dat <- dat[c(5,1,3,4,2),]


spawn_z <- matrix(scale(dat[2,]),nrow=1)  #z score spawners
spawn_z <- spawn_z[,-c(1:4)]              #remove the NAs
pdo_s <- matrix(dat[3,-c(1:4)], nrow=1)
pdo_s <- scale(pdo_s[1,])
ratio <- scale(SR_data$y)[-c(1:4)]


# z-score the predictor variable
m <- 3                     #number of parameters in process 
TT <- length(ratio)      #number of data points

B <-  diag(m)                         #"identity"
U <-  matrix(0, nrow = m, ncol = 1)   #"zero"
Q <- matrix(list(0, 0,0,              #we no longer want alpha to vary
                 0, "q.beta",0,
                 0,0,"q.pdo_s"), nrow = 3, byrow = TRUE)  # to have characters and numbers in same matrix use a list
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2,] <- spawn_z           ## Nx1; predictor variable spawners
Z[1,3,] <-  pdo_s            ## Nx1; predictor PDO summer
A = matrix(0)                #"zero"
R <-  matrix("r")



inits_list2 <- list(x0 = matrix(c(0,0,0), nrow = m))     
mod_list2 <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
cont_list <- list(maxit=10000)
# run model
m3 <- MARSS(ratio, mod_list2, inits = inits_list2, control = cont_list, method = "BFGS")




######################
#                    #
#      PART 2        #
#                    #
######################

#first model we assume there is no density dependence
#we model the underlying states of alpha (brood-year productivity) and beta
#beta is prevented from changing.
spawn_z <- matrix(scale(dat[2,]),nrow=1)  #z score spawners
spawn_z <- spawn_z[,-c(1:4)]              #remove the NAs
spawn <- matrix(dat[2,-c(1:4)],nrow=1)    #not sure if we should scale this or not
ratio <- SR_data$y[-c(1:4)]               #remove responses that corresponded to Spawner = NA
# z-score the predictor variable
m <- 2                     #number of parameters in process 
TT <- length(spawn_z)      #number of data points

B <-  matrix(1)                         #"identity"
U <-  matrix(0, nrow = 1, ncol = 1)   #"zero"
Q <- matrix("q")  # to have characters and numbers in same matrix use a list
Z <- matrix(1) # array(NA, c(1, m, TT))  ## NxMxT; empty for now
# Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
# Z[1,2,] <- spawn_z             ## Nx1; predictor variable 
A = matrix(0)            #"zero"
R <-  matrix("r")
d <- spawn
D <- matrix("beta")

inits_list <- list(x0 = matrix(c(0), nrow = 1))     
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, d = d, D = D)

## fit the model with MARSS
m2 <- MARSS(ratio, mod_list, inits = inits_list)

