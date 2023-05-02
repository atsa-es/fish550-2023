## load MARSS for data and analyses
library(MARSS)

## load the raw data (there are 3 datasets contained here)
data(lakeWAplankton, package = "MARSS")

## we want `lakeWAplanktonTrans`, which has been transformed
## so the 0's are replaced with NA's and the data z-scored
all_dat <- lakeWAplanktonTrans

all_dat


# set the window we wish to examine
yr_frst <- 1974
yr_last <- 1994
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

model.states <- 2:4
names(model.states) <- c("two_states","three_states","four_states")

# Also testing different observation error matrices
R.models <- c("diagonal and unequal","equalvarcov","unconstrained")

# Ignoring pH for now because it does not have a complete dataset. May add seasonal variables later. From early results it looks like best fits all have covariation. I have taken out a few of the options to reduce the number of models
# Also the covariates must be z-scored if the response variables are demeaned
d.models <- list( 
  #  d1 = "zero", # No covariates
  #  d2 = t(init_dat[,"Temp"]), # Temperature
  #  d3 = t(init_dat[,"TP"]), # Total Phosphorus
  d4 = zscore(t(init_dat[,c("Temp","TP")])) # Temperature and Total Phosphorus
)
names(d.models) <- c("Temp and TP")


# Setting fixed portion of mod list, some of these we can mess around with but we will already be testing 48 models and I didn't want to eat up a bunch of processing time
mod.list = list(
  A = "zero", # this is set to zero because the data has been demeaned
  Q = "identity", # Q is set to the default value, this can also be "diagonal and equal" or "diagonal and unequal"
  x0 = "zero" # x0 can be also be set to "unconstrained" or "unequal" it might be useful to try these later
)

# Added in another for loop to run through cycling options, be warned, this means there are 80 total models. It takes a long time to run
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

min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]
out.tab.1

# all other params are left to their default settings


# Rotating the model, plotting the states and loadings
best_fit <- fits[[3]]

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
layout(matrix(1:6, m, 2), widths = c(3,2))
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


## set up plotting area
par(mai = c(0.9,0.9,0.3,0.1))

## plot CCF's
ccf(proc_rot[1,],proc_rot[2,], lag.max = 30, main="")
title(main = "State 1 vs 2")
