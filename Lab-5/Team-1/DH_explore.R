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

#first model we assume there is no density dependence
#we model the underlying state of alpha (brood-year productivity)
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
alpha <- as.numeric(m1$states)
alpha.se <- as.numeric(m1$states.se)
plot(alpha~SR_data$brood_year,type='l', ylim=c(-3,3))
lines(alpha+2*alpha.se~SR_data$brood_year, lty="dashed")
lines(alpha-2*alpha.se~SR_data$brood_year, lty="dashed")

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


