## load MARSS for data and analyses
library(MARSS)
library(ggplot2)

## load the raw data (there are 3 datasets contained here)
data(lakeWAplankton, package = "MARSS")

## we want `lakeWAplanktonTrans`, which has been transformed
## so the 0's are replaced with NA's and the data z-scored
all_dat <- lakeWAplanktonTrans
colnames(all_dat)

#look at the different organisms
plot.ts(all_dat[,6]) #Cryptomonas: missing a lot of first half of ts
plot.ts(all_dat[,7]) #Diatoms: ts mostly complete
plot.ts(all_dat[,8]) #Greens: missing a big of beginning of ts
plot.ts(all_dat[,9]) #Bluegreens: most of last half of ts missing
plot.ts(all_dat[,10]) #unicells: ts mostly complete
plot.ts(all_dat[,11]) #other.algae: ts mostly complete
plot.ts(all_dat[,12]) #conochilus: some time points in first half missing
plot.ts(all_dat[,13]) #cyclops: ts mostly complete
plot.ts(all_dat[,14]) #daphnia: first half of ts mostly missing
plot.ts(all_dat[,15]) #diaptomus: ts mostly complete
plot.ts(all_dat[,16]) #epischura: ts mostly complete
plot.ts(all_dat[,17]) #leptodora: very gappy ts
plot.ts(all_dat[,18]) #neomysis: almost no data
plot.ts(all_dat[,19]) #non daphnic cladocera: mostly completely ts
plot.ts(all_dat[,20]) #rotifers: mostly completely ts
plot.ts(all_dat[,3]) #temp: complete
plot.ts(all_dat[,4]) #TP: complete
plot.ts(all_dat[,5]) #pH : mostly complete

#good data: diatoms, unicells, other.algae, cyclops, diaptomus, epischura, non daphnic cladocera, rotifers
#what eats algae: cyclops, diaptomus, rotifers
clr <- c("red", "orange",  "green", "blue", "purple")
#subset data to taxa of interest and pH, temp, and TP

plank.for.DFA<-subset(all_dat, select=c("Year", 'Unicells', 'Other.algae', 'Cyclops', 'Diaptomus', 'Non.colonial.rotifers'))

#transpose data
plank.t<-t(plank.for.DFA)
year.cols<-plank.t[1,]
plank.t<-plank.t[2:6,]
colnames(plank.t)<-year.cols
N.ts<-nrow(plank.t)
TT<-ncol(plank.t)
# mean of each taxon
y_bar <- apply(plank.t, 1, mean)
# subtract the means
dat <- plank.t - y_bar
yr_frst <- 1962
yr_last <- 1994

#since we did not subset the data we do not need to take a new z-score
#5 time series and assume 3 hidden trends
Z.vals<-list(
  "z11", 0, 0,
  "z21", "z22", 0,
  "z31", "z32", "z33",
  "z41", "z42", "z43",
  "z51", "z52", "z53"
)

Z<-matrix(Z.vals, nrow=N.ts, ncol=3, byrow=T)

#Q and B are equal to identity matrix
Q<-B<-diag(1,3) # "identity"

#assume different observation variance for R matrix
R<-"diagonal and unequal"

x0<-U<-A<-"zero" # also init_list <- list(x0 = matrix(rep(0,m),m,1)) where m = # processes
V0<-diag(5,3)

dfa.model1<-list(Z=Z, A="zero", R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)
cntl.list<-list(maxit=100) #may want to set this higher since error that convergence was not reached

#fit the model
dfa.fit1<-MARSS(plank.t, model=dfa.model1, control=cntl.list)

autoplot(dfa.fit1, plot.type='acf.std.model.resids.ytt1')
#ACF shows quite a bit of autocorrelation - this needs to be addressed

#compare AICc for 2, 3, or 4 underlying trends
model.list<-list(m=2, R="diagonal and unequal")
dfa.fit2<-MARSS(plank.t, model=model.list, z.score=T, form='dfa', control=cntl.list)

model.list<-list(m=4, R="diagonal and unequal")
dfa.fit3<-MARSS(plank.t, model=model.list, z.score=T, form='dfa', control=cntl.list)


print(cbind(model=c("3 trends", "2 trends", "4 trends"), AICc=round(c(dfa.fit1$AICc, dfa.fit2$AICc, dfa.fit3$AICc))), quote=F)

#4 trends has lowest AIC
#model    AICc
#[1,] 2 trends 4818
#[2,] 3 trends 4930
#[3,] 4 trends 4629

# It looks like a 4 state model is best. 



# get estimated ZZ
Z_est <- coef(dfa.fit3, type = "matrix")$Z
H_inv <- varimax(Z_est)$rotmat

# rotate factor loadings
Z.rot = Z_est %*% H_inv
trends.rot = solve(H_inv) %*% dfa.fit3$states

#####################################################
#####################################################
w_ts <- seq(dim(plank.t)[2])
mm <- 4
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(trends.rot[i,]))
  ## set up plot area
  plot(w_ts,trends.rot[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, trends.rot[i,], lwd = 2)
  lines(w_ts, trends.rot[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(plank.t)[2]) + 1, yr_frst + 0:dim(plank.t)[2])
}

# plot the factor loadings
spp <- rownames(plank.t)
minZ <- 0.05
m <- dim(trends.rot)[1]
ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
par(mfrow = c(ceiling(m / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:m) {
  plot(c(1:N.ts)[abs(Z.rot[, i]) > minZ], as.vector(Z.rot[abs(Z.rot[, i]) > minZ, i]),
       type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1)
  )
  for (j in 1:N.ts) {
    if (Z.rot[j, i] > minZ) {
      text(j, -0.05, spp[j], srt = 90, adj = 1, cex = 0.9)
    }
    if (Z.rot[j, i] < -minZ) {
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9)
    }
    abline(h = 0, lwd = 1, col = clr[i])
  } # end j loop
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5)
} # end i loop
#####################################################
#####################################################

# If there were no missing values, this function will return the fits and CIs
getDFAfits <- function(MLEobj, alpha = 0.05, dd = NULL) { #dd = covariates
  fits <- list()
  Ey <- MARSShatyt(MLEobj) # for var() calcs
  ZZ <- coef(MLEobj, type = "matrix")$Z # estimated Z
  nn <- dim(Ey$ytT) # number of obs ts
 # mm <- ncol(ZZ) # number of factors/states
  TT <- ncol(Ey$ytT) # number of time steps
  ## check for covars
  if (!is.null(dd)) {
    DD <- coef(MLEobj, type = "matrix")$D
    cov_eff <- DD %*% dd
  } else {
    cov_eff <- matrix(0, nn, TT)
  }
  ## model expectation
  fits$ex <- ZZ %*% H_inv %*% MLEobj$states + cov_eff
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for (tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, , tt] %*% t(ZZ)
    SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% t(MLEobj$states[, tt, drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1 - alpha / 2) * SE + fits$ex
  fits$lo <- qnorm(alpha / 2) * SE + fits$ex
  return(fits)
}

# examine fits
mod_fit <- getDFAfits(dfa.fit3)

############################################
# plot the fits

ylbl <- rownames(plank.t)
#w_ts <- seq(dim(plank.t)[2])
par(mfcol = c(5, 1), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:N.ts) {
  up <- mod_fit$up[i, ]
  mn <- mod_fit$ex[i, ]
  lo <- mod_fit$lo[i, ]
  plot(w_ts, mn,
       xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", cex.lab = 1.2,
       ylim = c(min(lo), max(up))
  )
  axis(1, 12 * (0:dim(plank.t)[2]) + 1, 1980 + 0:dim(plank.t)[2])
  points(w_ts, plank.t[i, ], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}

#########################
# This was in manual but I don't think we need it
require(ggplot2)
alpha <- 0.05
theme_set(theme_bw())
d <- residuals(dfa.fit3, type = "tT")
d$up <- qnorm(1 - alpha / 2) * d$.sigma + d$.fitted
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted
ggplot(data = subset(d, name=="model")) +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) +
  xlab("Time Step") +
  ylab("Count")
#########################

# Covariates
temp <- t(all_dat[,"Temp", drop = FALSE])
TP <- t(all_dat[,"TP", drop = FALSE])

mod_list = list(m = 4, R = "diagonal and unequal")
dfa_temp <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                  control=cntl.list, covariates = temp)
dfa_TP <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                control=cntl.list, covariates = TP)
dfa_both <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                  control=cntl.list, covariates = rbind(temp, TP))

print(cbind(model = c("no covars", "Temp", "TP", "Temp & TP"),
            AICc = round(c(dfa.fit3$AICc,
                           dfa_temp$AICc,
                           dfa_TP$AICc,
                           dfa_both$AICc))),
      quote = FALSE)

# model     AICc
# [1,] no covars 4629
# [2,] Temp      4508
# [3,] TP        4607
# [4,] Temp & TP 4477


# It looks like temp and TP together are the best fit model.



## create dummy sine and cosine waves
cos_t <- cos(2 * pi * seq(TT) / 12)
sin_t <- sin(2 * pi * seq(TT) / 12)
dd <- rbind(cos_t, sin_t)

## fit model
dfa_seas <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                  control = cntl.list, covariates = dd)

dfa_TPseas <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                  control=cntl.list, covariates = rbind(TP, dd))

dfa_TPtempseas <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                    control=cntl.list, covariates = rbind(temp,TP, dd))

dfa_tempseas <- MARSS(plank.t, model = mod_list, form = "dfa", z.score = FALSE,
                        control=cntl.list, covariates = rbind(temp, dd))

#AICc for TP + seasons and TP + temp + seasons are similar and lowest

## get model fits & CI's
mod_fit <- getDFAfits(dfa_TPseas)
mod_fit2 <- getDFAfits(dfa_TPtempseas)

## plot the fits
par(mfrow = c(N.ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))
for(i in 1:N.ts) {
  up <- mod_fit$up[i,]
  mn <- mod_fit$ex[i,]
  lo <- mod_fit$lo[i,]
  plot(w_ts, mn, type = "n",
       xlab = "", ylab = ylbl[i],
       xaxt = "n", cex.lab = 1.2,
       ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(plank.t)[2]) + 1, yr_frst + 0:dim(plank.t)[2])
  points(w_ts, plank.t[i,], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}

for(i in 1:N.ts) {
  up <- mod_fit2$up[i,]
  mn <- mod_fit2$ex[i,]
  lo <- mod_fit2$lo[i,]
  plot(w_ts, mn, type = "n",
       xlab = "", ylab = ylbl[i],
       xaxt = "n", cex.lab = 1.2,
       ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(plank.t)[2]) + 1, yr_frst + 0:dim(plank.t)[2])
  points(w_ts, plank.t[i,], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}


#plots of both fits appear visually the same

print(cbind(model = c("no covars", "Temp", "TP", "Temp & TP", "Seasonality", "TP Seas", "TP & temp & Seas", "temp & Seas"),
            AICc = round(c(dfa.fit3$AICc,
                           dfa_temp$AICc,
                           dfa_TP$AICc,
                           dfa_both$AICc,
                           dfa_seas$AICc,
                           dfa_TPseas$AICc,
                           dfa_TPtempseas$AICc,
                           dfa_tempseas$AICc))),
      quote = FALSE)

# It looks like TP and seasonality together are the best fit model.

#####################
#The first 10 years of the data seem to have a stronger seasonality. Does the model fit better to just the first 10 years?
all_dat2 <- lakeWAplanktonRaw
plank.10yr<-subset(all_dat2, select=c("Year", 'Unicells', 'Other.algae', 'Cyclops', 'Diaptomus', 'Non.colonial.rotifers'), all_dat[,1]<1973)

#transpose data
plank.t10<-t(plank.10yr)
year.cols<-plank.t10[1,]
plank.t10<-plank.t10[2:6,]
colnames(plank.t10)<-year.cols
N.ts<-nrow(plank.t10)
TT<-ncol(plank.t10)
# re-do z-score
dat.z<-zscore(plank.t10)

temp.10yr <- t(all_dat2[1:132,"Temp", drop = FALSE])
TP.10yr <- t(all_dat2[1:132,"TP", drop = FALSE])

temp.z<-zscore(temp.10yr)
TP.z<-zscore(TP.10yr)

mod_list = list(m = 4, R = "diagonal and unequal")

cos_t <- cos(2 * pi * seq(TT) / 12)
sin_t <- sin(2 * pi * seq(TT) / 12)
dd.10 <- rbind(cos_t, sin_t)

#dfa_TPseas.10 <- MARSS(plank.t10, model = mod_list, form = "dfa", z.score = TRUE,
 #                control=cntl.tmp, covariates = rbind(TP.z,dd.10))
#Error: Stopped in MARSS.dfa() due to specification problem(s).
#I cannot figure out why the TP does not work here

dfa_tempseas.10 <- MARSS(plank.t10, model = mod_list, form = "dfa", z.score = TRUE,
                        control=cntl.list, covariates = rbind(temp.z, dd.10))

mod_fit.tempseas.10 <- getDFAfits(dfa_tempseas.10)

autoplot(dfa_tempseas.10)

w_ts_10 <- seq(dim(plank.t10)[2])

par(mfrow = c(N.ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))
for(i in 1:N.ts) {
  up <- mod_fit.tempseas.10$up[i,]
  mn <- mod_fit.tempseas.10$ex[i,]
  lo <- mod_fit.tempseas.10$lo[i,]
  plot(w_ts_10, mn, type = "n",
       xlab = "", ylab = ylbl[i],
       xaxt = "n", cex.lab = 1.2,
       ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(plank.t10)[2]) + 1, yr_frst + 0:dim(plank.t10)[2])
  points(w_ts_10, dat.z[i,], pch = 16, col = clr[i])
  lines(w_ts_10, up, col = "darkgray")
  lines(w_ts_10, mn, col = "black", lwd = 2)
  lines(w_ts_10, lo, col = "darkgray")
}

Z_est.10 <- coef(dfa_tempseas.10, type = "matrix")$Z
H_inv.10 <- varimax(Z_est.10)$rotmat

# rotate factor loadings
Z.rot.10 = Z_est.10 %*% H_inv.10
trends.rot.10 = solve(H_inv.10) %*% dfa_tempseas.10$states

m10 <- dim(trends.rot.10)[1]
ylims <- c(-1.1 * max(abs(Z.rot.10)), 1.1 * max(abs(Z.rot.10)))
par(mfrow = c(ceiling(m10 / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:m10) {
  plot(c(1:N.ts)[abs(Z.rot.10[, i]) > minZ], as.vector(Z.rot.10[abs(Z.rot.10[, i]) > minZ, i]),
       type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1)
  )
  for (j in 1:N.ts) {
    if (Z.rot.10[j, i] > minZ) {
      text(j, -0.05, spp[j], srt = 90, adj = 1, cex = 0.9)
    }
    if (Z.rot.10[j, i] < -minZ) {
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9)
    }
    abline(h = 0, lwd = 1, col = clr[i])
  } # end j loop
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5)
}

for(i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(trends.rot.10[i,]))
  ## set up plot area
  plot(w_ts_10,trends.rot.10[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts_10, trends.rot.10[i,], lwd = 2)
  lines(w_ts_10, trends.rot.10[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(plank.t10)[2]) + 1, yr_frst + 0:dim(plank.t10)[2])
}

#OR!!! Try this after the big TP dump into the lake - does this change effect of TP in the model?
all_dat2 <- lakeWAplanktonRaw
plot.ts(all_dat2[,4])
TP.sub<-subset(all_dat2, select=c("Year", "TP"), all_dat[,1]<1991)
TP.sub<-subset(all_dat2, select=c("Year", "TP"), all_dat[,1]>1979)
plot(x=TP.sub[,1], y=TP.sub[,2])
plank.10yr<-subset(all_dat2, select=c("Year", 'Unicells', 'Other.algae', 'Cyclops', 'Diaptomus', 'Non.colonial.rotifers'), all_dat[,1]<1991)
plank.10yr<-subset(all_dat2, select=c("Year", 'Unicells', 'Other.algae', 'Cyclops', 'Diaptomus', 'Non.colonial.rotifers'), all_dat[,1]>1979)


#transpose data
plank.t10<-t(plank.10yr)
year.cols<-plank.t10[1,]
plank.t10<-plank.t10[2:6,]
colnames(plank.t10)<-year.cols
N.ts<-nrow(plank.t10)
TT<-ncol(plank.t10)
# re-do z-score
dat.z<-zscore(plank.t10)

temp.10yr <- t(all_dat2[1:180,"Temp", drop = FALSE])
TP.10yr <- t(all_dat2[1:180,"TP", drop = T])

temp.z<-zscore(temp.10yr)
TP.z<-zscore(TP.10yr)
#replace NA with 0 in TP.z
TP.z[is.na(TP.z)] <- 0

mod_list = list(m = 3, R = "diagonal and unequal")

dfa_TPseas.10 <- MARSS(plank.t10, model = mod_list, form = "dfa", z.score = TRUE,
              control=cntl.tmp, covariates = rbind(TP.z,dd.10))


dfa_temp.10 <- MARSS(plank.t10, model = mod_list, form = "dfa", z.score = TRUE,
                         control=cntl.tmp, covariates = rbind(temp.z))

mod_fit.temp.10 <- getDFAfits(dfa_temp.10)
mod_fit.TPseas.10 <- getDFAfits(dfa_TPseas.10)
#Error in ZZ %*% H_inv : non-conformable arguments
#this may be because of missing values?

autoplot(dfa_temp.10)
autoplot(dfa_TPseas.10)

w_ts_10 <- seq(dim(plank.t10)[2])

par(mfrow = c(N.ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 0, 0, 0))
for(i in 1:N.ts) {
  up <- mod_fit.TPseas.10$up[i,]
  mn <- mod_fit.TPseas.10$ex[i,]
  lo <- mod_fit.TPseas.10$lo[i,]
  plot(w_ts_10, mn, type = "n",
       xlab = "", ylab = ylbl[i],
       xaxt = "n", cex.lab = 1.2,
       ylim = c(min(lo), max(up)))
  axis(1, 12 * (0:dim(plank.t10)[2]) + 1, yr_frst + 0:dim(plank.t10)[2])
  points(w_ts_10, dat.z[i,], pch = 16, col = clr[i])
  lines(w_ts_10, up, col = "darkgray")
  lines(w_ts_10, mn, col = "black", lwd = 2)
  lines(w_ts_10, lo, col = "darkgray")
}

#doesn't look like a great fit

Z_est.10 <- coef(dfa_TPseas.10, type = "matrix")$Z
H_inv.10 <- varimax(Z_est.10)$rotmat

# rotate factor loadings
Z.rot.10 = Z_est.10 %*% H_inv.10
trends.rot.10 = solve(H_inv.10) %*% dfa_temp.10$states

m10 <- dim(trends.rot.10)[1]
ylims <- c(-1.1 * max(abs(Z.rot.10)), 1.1 * max(abs(Z.rot.10)))
par(mfrow = c(ceiling(m10 / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:m10) {
  plot(c(1:N.ts)[abs(Z.rot.10[, i]) > minZ], as.vector(Z.rot.10[abs(Z.rot.10[, i]) > minZ, i]),
       type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1)
  )
  for (j in 1:N.ts) {
    if (Z.rot.10[j, i] > minZ) {
      text(j, -0.05, spp[j], srt = 90, adj = 1, cex = 0.9)
    }
    if (Z.rot.10[j, i] < -minZ) {
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9)
    }
    abline(h = 0, lwd = 1, col = clr[i])
  } # end j loop
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5)
}

mm<-3
for(i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(trends.rot.10[i,]))
  ## set up plot area
  plot(w_ts_10,trends.rot.10[i,], type = "n", bty = "L",
       ylim = ylm, xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts_10, trends.rot.10[i,], lwd = 2)
  lines(w_ts_10, trends.rot.10[i,], lwd = 2)
  ## add panel labels
  mtext(paste("State",i), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(plank.t10)[2]) + 1, yr_frst + 0:dim(plank.t10)[2])
}

# Results and Discussion

# We first determined the optimal number of states by fitting multiple models, 
# each with different numbers of states. The models were then compared via AICc to determine 
# what number of states best fit the time series. 
# 
# The best fitting model had 4 separate states over the full period of time, there were no other competing
# models based on delta AICc.
# 
# 
# We then chose covariates and compared models with different combinations off
# temperature, total phosphorus, a combination of temperature and total phosphorus, 
# dummy sine and cosine waves to simulate seasonality, and total phosphorus combined with seasonality.
# We did not fit a model with temperature and seasonality, as it was assumed that these two 
# variables would be highly correlated, and thus shouldn't be modeled together as their effects 
# wouldn't be able to be distinguished.
# 
# A model with total phosphorus combined with seasonality was the best fit according to delta AICc. 
# 

# However, the model fits do not appear to fit the data well when plotted. This perhaps indicates that 
# Across the entire time series, there is no model that performs well over all areas. Breaking the time series up 
# into smaller lengths and modeling those individually may help alleviate this problem, and highlight 
# how the influence of different covariates changes over time. The lake has shifted away from a state off
# eutrophication after raw sewage stopped bying pummped into the lake. The covariates that had dominate roles 
# during eutrophication may not be the same afterwards. 
# # 





      