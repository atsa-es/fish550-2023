
# load MARSS for data and analyses
library(MARSS)

## load the raw data (there are 3 datasets contained here)
data(lakeWAplankton, package = "MARSS")

## we want `lakeWAplanktonTrans`, which has been transformed
## so the 0's are replaced with NA's and the data z-scored
all_dat <- lakeWAplanktonTrans

head(all_dat)

## use only the 10 years from 1980-1989
yr_frst <- 1962
yr_last <- 1994
plank_dat <- all_dat[all_dat[, "Year"] >= yr_frst & 
                       all_dat[, "Year"] <= yr_last,]

## create vector of phytoplankton group names
phytoplankton1 <- c("Cryptomonas", "Diatoms", "Greens",
                   "Bluegreens", "Unicells")
phytoplankton2 <- c("Other.algae", "Conochilus", "Cyclops",
                   "Daphnia", "Diaptomus")
phytoplankton3 <- c("Epischura", "Leptodora", "Neomysis",
                   "Non.daphnid.cladocerans", "Non.colonial.rotifers")

## get only the phytoplankton
dat_init <- plank_dat[, phytoplankton3]

## transpose data so time goes across columns
dat_init <- t(dat_init)

## get number of time series
N_ts <- dim(dat_init)[1]

## get length of time series
TT <- dim(dat_init)[2] 

## mean of each taxon
y_bar <- apply(dat_init, 1, mean, na.rm = TRUE)

## subtract the means
dat <- dat_init - y_bar

## assign new column names
spp <- rownames(dat_init)
rownames(dat) <- spp


## set plot colors
clr <- c("brown", "blue", "darkgreen", "darkred", "purple")

## initialize a counter
cnt <- 1

## set up plotting space & make plots
par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))

for(i in spp){
  plot(dat[i,],  bty = "L", xaxt = "n", pch = 16,
       xlab = "",
       ylab = "Abundance index", col=clr[cnt], type="b")
  axis(1, 12 * (0:dim(dat_init)[2]) + 1, yr_frst + 0:dim(dat_init)[2])
  title(i)
  cnt <- cnt + 1
}


#================================
#Official RMarkdown start 
#=================================
# set the window we wish to examine
yr_frst <- 1978
yr_last <- 1985
init_dat <- all_dat[all_dat[, "Year"] >= yr_frst & 
                      all_dat[, "Year"] <= yr_last,]

# put the 5 response variables in the yt matrix, can be replace whenever
# This represents multiple trophic levels
# Primary prod: Diatoms & Other Algae
# Grazers: Daphnia & Cyclops
# Predators: Epischura
plank_taxa <- c("Diatoms","Other.algae","Daphnia","Cyclops","Epischura")
dat <- init_dat[,plank_taxa] %>%
  t()
## mean of each taxon
y_bar <- apply(dat, 1, mean, na.rm = TRUE)
## subtract the means
plank_dat <- dat - y_bar
## assign new column names
spp <- rownames(plank_dat)
rownames(plank_dat) <- spp

# create a covariates matrix
covnames <- c("Temp","TP") #pH has missing data points so it will only work in certain date ranges
covar <- init_dat[,covnames] %>%
  t()

#==========
#Plot Taxa
#==========

## set plot colors
clr <- c("brown", "blue", "darkgreen", "darkred", "purple")

## initialize a counter
cnt <- 1

## set up plotting space & make plots
par(mfrow = c(nrow(plank_dat), 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))
mm = ncol(plank_dat)

for(i in spp){
  plot(plank_dat[i,],  bty = "L", xaxt = "n", pch = 16,
       xlab = "",
       ylab = "Abundance index", col=clr[cnt], type="b")
  axis(1, 12 * (0:mm) + 1, yr_frst + 0:mm)
  title(i)
  cnt <- cnt + 1
}

#================
#Plot covariates 
#===============

## set plot colors
clr <- c("brown", "blue", "darkgreen","purple")

## initialize a counter
cnt <- 1

## set up plotting space & make plots
par(mfrow = c(nrow(covar), 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))

for(i in row.names(covar)){
  plot(covar[i,],  bty = "L", xaxt = "n", pch = 16,
       xlab = "",
       ylab = "Value", col=clr[cnt], type="b")
  axis(1, 12 * (0:mm) + 1, yr_frst + 0:mm)
  title(i)
  cnt <- cnt + 1
}

#================
#Model Set up 
#================
model.states <- 2:4
names(model.states) <- c("two_states","three_states","four_states")

# Also testing different observation error matrices
R.models <- c("diagonal and unequal","equalvarcov", "diagonal and equal")

# Ignoring pH for now because it does not have a complete dataset. May add seasonal variables later. From early results it looks like best fits all have covariation. I have taken out a few of the options to reduce the number of models
# Also the covariates must be z-scored if the response variables are demeaned
d.models <- list( 
  # d1 = "zero", # No covariates
  d2 = zscore(t(init_dat[,"Temp"])), # Temperature
  d3 = zscore(t(init_dat[,"TP"])), # Total Phosphorus
  d4 = zscore(t(init_dat[,"pH"])), # ph
  d5 = zscore(t(init_dat[,c("Temp","TP")])),
  d6 = zscore(t(init_dat[,c("TP", "pH")])),
  d7 = zscore(t(init_dat[,c("Temp", "pH")])),
  d8 = zscore(t(init_dat[,c("Temp","TP", "pH")])) # Temperature and Total Phosphorus
)
names(d.models) <- c("Temp", "TP", "pH", "Temp and TP", "TP and pH", "Temp and pH", "Temp,TP, and pH")


# Setting fixed portion of mod list, some of these we can mess around with but we will already be testing 48 models and I didn't want to eat up a bunch of processing time
mod.list = list(
  A = "zero", # this is set to zero because the data has been demeaned
  Q = "identity", # Q is set to the default value, this can also be "diagonal and equal" or "diagonal and unequal"
  x0 = "zero" # x0 can be also be set to "unconstrained" or "unequal" it might be useful to try these later
)

#DFA form for mars model Q--reduces option to streamline runtime

# all other params are left to their default settings

start.time <- Sys.time()


out.tab <- NULL
fits <- list()
#for(i in model.states){
  for (j in 1:length(d.models)){
    for(R.model in R.models){
      fit.model = c(list(m=i, R=R.model), mod.list)
      
      #MARSS will throw an error if covariates are not a matrix, This will not pass the covariate arg if it does not exist
      #      if(j==1) {
      #      fit <- MARSS(plank_dat, model=fit.model, 
      #                form = "dfa", z.score = FALSE, silent = TRUE,
      #                control=list(maxit=2000, allow.degen=TRUE))
      #      }
      #      else {
      
      fit <- MARSS(plank_dat, model=fit.model,
                   form = "dfa", z.score = FALSE, silent = TRUE,
                   control=list(maxit=800, allow.degen = TRUE),
                   covariates = d.models[[j]])

      
      #      }
      
      out=data.frame(States=names(model.states)[i-1],d=names(d.models)[j],R=R.model, # 
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab=rbind(out.tab,out)
      fits=c(fits,list(fit))
      
    }
  }
}
end.time <- Sys.time()

elapsed.time <- round((end.time - start.time), 3)
print(elapsed.time)
min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]
out.tab.1

#===============
#Diagnostics 
#===============

# Rotating the model, plotting the states and loadings
best_fit <- fits[[19]]

# Get Z Estimations and rotate the data
## get the estimated ZZ
Z_est <- coef(best_fit, type = "matrix")$Z

## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv   

## rotate processes
proc_rot = solve(H_inv) %*% best_fit$states

N_ts <- nrow(plank_dat)
m = ncol(Z_est)
## plot labels
ylbl <- plank_taxa
w_ts <- seq(mm)

## set up plot area
layout(matrix(1:4, m, 2), widths = c(3,2))
par(mai = c(0.3, 0.3, 0.3, 0.1), omi = c(0, 0, 0, 0))

## plot the processes
for(i in 1:m) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i,], lwd = 2)
  lines(w_ts, proc_rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i), side = 3, line = 0.5)
  axis(1, 12 * (0:mm) + 1, yr_frst + 0:mm)
}
## plot the loadings
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for(i in 1:m) {
  plot(x = c(1:N_ts)[abs(Z_rot[,i])>minZ],
       y = as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]),
       type = "h",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylm,
       xlim = c(0.5, N_ts + 0.5), col = clr)
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, col = clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, col = clr[j])}
    abline(h = 0, lwd = 1.5, col = "gray")
  } 
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}

## plot CCF's
ccf(proc_rot[1,],proc_rot[2,], lag.max = 15, main="")
title(main = "State 1 vs 2")

####
#More plots 

###

get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
  ## empty list for results
  fits <- list()
  ## extra stuff for var() calcs
  Ey <- MARSS:::MARSShatyt(MLEobj)
  ## model params
  ZZ <- coef(MLEobj, type = "matrix")$Z
  ## number of obs ts
  nn <- dim(Ey$ytT)[1]
  ## number of time steps
  TT <- dim(Ey$ytT)[2]
  ## get the inverse of the rotation matrix
  H_inv <- varimax(ZZ)$rotmat
  ## check for covars
  if (!is.null(dd)) {
    DD <- coef(MLEobj, type = "matrix")$D
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
  } else {
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states
  }
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for (tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                         , tt] %*% t(ZZ)
    SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
      t(MLEobj$states[, tt, drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
  fits$lo <- qnorm(alpha/2) * SE + fits$ex
  return(fits)
}


clr <- c("brown", "blue", "darkgreen", "darkred", "purple")
N_ts <- dim(plank_dat)[1]
w_ts <- seq(dim(plank_dat)[2])
## get model fits & CI's
mod_fit <- get_DFA_fits(best_fit)
## plot the fits
ylbl <- plank_taxa
par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                             0, 0, 0))
for (i in 1:N_ts) {
  up <- mod_fit$up[i, ]
  mn <- mod_fit$ex[i, ]
  lo <- mod_fit$lo[i, ]
  plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
       cex.lab = 1.2, ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(plank_dat)[2]) + 1, yr_frst + 0:dim(plank_dat)[2])
  points(w_ts, plank_dat[i, ], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}

#=====================
#Cycling Assumptions 
#=====================

TT <- dim(plank_dat)[2]
period <- 6
cos.t <- cos(2 * pi * seq(TT)/period)
sin.t <- sin(2 * pi * seq(TT)/period)
c.Four <- rbind(cos.t, sin.t)

model.states <- 2
names(model.states) <- c("two_states")

# Also testing different observation error matrices
R.models <- c("diagonal and unequal")

# Ignoring pH for now because it does not have a complete dataset. May add seasonal variables later. From early results it looks like best fits all have covariation. I have taken out a few of the options to reduce the number of models
# Also the covariates must be z-scored if the response variables are demeaned
d10 = zscore(t(init_dat[,c("Temp","TP", "pH")]))
d9 = rbind(c.Four,d10)

d.models <- list(d10, d9)
  

names(d.models) <- c("Temp,TP, and pH", "with seasonality")


# Setting fixed portion of mod list, some of these we can mess around with but we will already be testing 48 models and I didn't want to eat up a bunch of processing time
mod.list = list(
  A = "zero", # this is set to zero because the data has been demeaned
  Q = "identity", # Q is set to the default value, this can also be "diagonal and equal" or "diagonal and unequal"
  x0 = "zero" # x0 can be also be set to "unconstrained" or "unequal" it might be useful to try these later
)


# Added in another for loop to run through cycling options, be warned, this means there are 80 total models. It takes a long time to run
start.time <- Sys.time()


out.tab <- NULL
fits <- list()
for(i in model.states){
  for (j in 1:length(d.models)){
    for(R.model in R.models){
      fit.model = c(list(m=i, R=R.model), mod.list)
      
      #MARSS will throw an error if covariates are not a matrix, This will not pass the covariate arg if it does not exist
      #      if(j==1) {
      #      fit <- MARSS(plank_dat, model=fit.model, 
      #                form = "dfa", z.score = FALSE, silent = TRUE,
      #                control=list(maxit=2000, allow.degen=TRUE))
      #      }
      #      else {
      fit <- MARSS(plank_dat, model=fit.model, 
                   form = "dfa", z.score = FALSE, silent = TRUE,
                   control=list(maxit=2000, allow.degen=TRUE), 
                   covariates = d.models[[j]])
      #      }
      
      out=data.frame(States=names(model.states)[i-1],d=names(d.models)[j],R=R.model, # 
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab=rbind(out.tab,out)
      fits=c(fits,list(fit))
      
    }
  }
}
end.time <- Sys.time()

elapsed.time <- round((end.time - start.time), 3)
print(elapsed.time)
min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]
out.tab.1


######

#Notes on CCF 

#three states, out to ~40 steps were correlation 

#Hole where values are 

## 50 day lag between season change and temperature change 

