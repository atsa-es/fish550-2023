library(MARSS)
library(ggplot2)
library(tidyr)
library(forecast)
data(lakeWAplankton, package = "MARSS")
all_dat <- lakeWAplanktonTrans

colnames(all_dat)[6:20]
for(i in 6:20){
  name <- colnames(all_dat)[i]
  plot(all_dat[,i] ~ all_dat[,1], main = paste0(name))
}
plot(all_dat[,5]~all_dat[,1])

#we only want 5 of them
crit <- all_dat[c(61:300),c(1:5,7,8,10,13,20)]

need <- t(crit[,6:10])
colnames(need) <- crit[,1]
rownames(need) <- c("Diatoms", "Greens", "Unicells", "Cyclops", "Rotifers")


#get the number of time series (y)
n_ts <- nrow(need)
#get length of time series
TT <- ncol(need)

#find the mean of each time series and remove it
y_bar <- apply(need, 1, mean, na.rm = TRUE)
need.av <- need - y_bar

#plot each group
yr_start <- 1967
yr_end <- 1986
spp <- rownames(need.av)
clr <- c("brown", "blue", "darkgreen", "darkred", "purple")
cnt <- 1
par(mfrow = c(n_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                             0, 0, 0))
for (i in spp) {
  plot(need.av[i, ], xlab = "", ylab = "Abundance index", bty = "L", 
       xaxt = "n", pch = 16, col = clr[cnt], type = "b")
  axis(side = 1, at =  12 * (0:dim(need.av)[2]) + 1,
       labels = yr_start + 0:dim(need.av)[2])
  title(i)
  cnt <- cnt + 1
}
#set up observation model parameters
#create z loading matrix. assume 3 underlying processes
Z_mat <- matrix(c("z11",  0,    0,
                "z21","z22",  0,
                "z31","z32","z33",
                "Z41","z42","z43",
                "z51","z52","z53"),5,3, byrow = TRUE)

#A scaling 
a <- matrix(c(0,0,0,"a1","a2"),5,1)

#covariates blank for now
D <- "zero"
d <- "zero"

#R matrix
r <- "diagonal and equal"

#Set up the process models
#number of underlying states
mm <- 3

#B matrix will be a diagonal (random walks)
B <- diag(3)

#U matrix for drift
u <- "zero"

#Covariates for underlying process
C <- "zero"
c <- "zero"

#Q matrix will be fixed as diagonal and equal
Q <- "diagonal and equal"

#create matrices list for MARSS model
mod_list <- list(Z = Z_mat, A = a, D = D, d = d, R = r, B = B, 
                 U = u, C = C, c = c, Q = Q)
#tell MARSS what values to use to start with
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))

#run the model
dfa_1 <- MARSS(y = need.av, model = mod_list, inits = init_list)

## get the estimated ZZ
Z_est <- coef(dfa_1, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat
# rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% dfa_1$states

######################################
# plot the states and their loading  #
######################################

ylbl <- rownames(need.av)
w_ts <- seq(dim(need.av)[2])
layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi =
## c(0,0,0,0))
par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
## plot the processes
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(need.av)[2]) + 1, yr_start + 0:dim(need.av)[2])
}
## plot the loadings
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:n_ts)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, 
                                                                    i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, n_ts + 0.5), col = clr)
  for (j in 1:n_ts) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
           col = clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
           col = clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}
dev.off()
#look at cross correlation between states
ccf(proc_rot[1, ], proc_rot[2, ], lag.max = 12, main = "")
ccf(proc_rot[1, ], proc_rot[3, ], lag.max = 12, main = "")
ccf(proc_rot[2, ], proc_rot[3, ], lag.max = 12, main = "")

#get model fits this custom function does the rotation process automatically
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

## get model fits & CI's
mod_fit <- get_DFA_fits(dfa_1)
## plot the fits
ylbl <- rownames(need.av)
par(mfrow = c(n_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                             0, 0, 0))
for (i in 1:n_ts) {
  up <- mod_fit$up[i, ]
  mn <- mod_fit$ex[i, ]
  lo <- mod_fit$lo[i, ]
  plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
       cex.lab = 1.2, ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(need.av)[2]) + 1, yr_start + 0:dim(need.av)[2])
  points(w_ts, need.av[i, ], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}

###now see if the covariates help
#create a matrix for each of the covariates
temp <- t(crit[,3])
Phos <- t(crit[,4])
ph <- t(crit[,5])
#make a season matrix
cos_t <- cos(2 * pi * seq(TT)/12)
sin_t <- sin(2 * pi * seq(TT)/12)
season <- rbind(cos_t,sin_t)
dim(season)

#try it with all covariates, arg. form="dfa" creates model list for you besides m and R

mod_list2 <- list(m=3, R = "diagonal and equal")
cont_list <- list(maxit = 3000, allow.degen = TRUE)
dfa_global <- MARSS(need.av, model = mod_list2,control = cont_list, inits = init_list, form = "dfa",
                    z.score = FALSE, covariates = rbind(temp,Phos,ph,season))

#get model fits for global model and plot them
mod_fit2 <- get_DFA_fits(dfa_global, dd = rbind(temp,Phos,ph,season))
par(mfrow = c(n_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1),
        omi = c(0,0, 0, 0))
for (i in 1:n_ts) {
  up <- mod_fit2$up[i, ]
  mn <- mod_fit2$ex[i, ]
  lo <- mod_fit2$lo[i, ]
  plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
       cex.lab = 1.2, ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(need.av)[2]) + 1, yr_start + 0:dim(need.av)[2])
  points(w_ts, need.av[i, ], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}

#########look at AIC###########
print(cbind(model = c("no covars", "global"), 
            AICc = round(c(dfa_1$AICc, dfa_global$AICc))), 
      quote = FALSE)


################look at Z and D matrix########################
coef(dfa_global, type = "matrix")$D
coef(dfa_global, type = "matrix")$Z

########################
#   Plot the loadings  #
########################

## get the estimated ZZ
Z_est <- coef(dfa_global, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat
# rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% dfa_global$states

ylbl <- rownames(need.av)
w_ts <- seq(dim(need.av)[2])
layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai = c(0.5,0.5,0.5,0.1), omi =
## c(0,0,0,0))
par(mai = c(0.25, 0.5, 0.25, 0.1), omi = c(0, 0, 0, 0))

## plot the processes
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(need.av)[2]) + 1, yr_start + 0:dim(need.av)[2])
}
## plot the loadings
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:n_ts)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, 
                                                                    i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, n_ts + 0.5), col = clr)
  for (j in 1:n_ts) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
           col = clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
           col = clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}

dev.off()
#look at cross correlation between states
ccf(proc_rot[1, ], proc_rot[2, ], lag.max = 12, main = "")
ccf(proc_rot[1, ], proc_rot[3, ], lag.max = 12, main = "")
ccf(proc_rot[2, ], proc_rot[3, ], lag.max = 12, main = "")




#make unique combinations of covariate lists
covar <- rbind(temp,Phos,ph,season)
#rownames(covar) <- c("temp", "Phos", "ph", "season_cos", "season_sin") # do NOT name the rows it adds parameters?

library(rje)
combo <- powerSet(1:4)
covar[combo[[6]],]
#make sure cos and sin are always together as season
for (i in 9:16) {
  combo[[i]] <- c(combo[[i]],5)
}
combo  # [[1]] is empty so dont loop it

#test this works

dfa_global_test0 <- MARSS(need.av, model = mod_list2,control = cont_list, inits = init_list, form = "dfa",
                         z.score = FALSE, covariates = rbind(temp,Phos,ph,season))
dfa_global_test <- MARSS(need.av, model = mod_list2,control = cont_list, inits = init_list, form = "dfa",
                    z.score = FALSE, covariates = covar[combos[[16]],])
dfa_global_test2 <- MARSS(need.av, model = mod_list2,control = cont_list, inits = init_list, form = "dfa",
                         z.score = FALSE, covariates = covar)


