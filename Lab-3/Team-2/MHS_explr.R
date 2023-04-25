
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
