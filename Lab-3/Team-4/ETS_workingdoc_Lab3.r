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

#subset data to taxa of interest and pH, temp, and TP

plank.for.DFA<-subset(all_dat, select=c("Year", 'Unicells', 'Other.algae', 'Cyclops', 'Diaptomus', 'Non.colonial.rotifers'))

#transpose data
plank.t<-t(plank.for.DFA)
year.cols<-plank.t[1,]
plank.t<-plank.t[2:6,]
colnames(plank.t)<-year.cols
N.ts<-nrow(plank.t)
TT<-ncol(plank.t)

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
Q<-B<-diag(1,3)

#assume same observation variance for R matrix
R<-"diagonal and equal"

x0<-U<-A<-"zero"
V0<-diag(5,3)

dfa.model1<-list(Z=Z, A="zero", R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)
cntl.list<-list(maxit=100) #increase this later to 1000

#fit the model
dfa.fit1<-MARSS(plank.t, model=dfa.model1, control=cntl.list)

autoplot(dfa.fit1)
#ACF shows quite a bit of autocorrelation - this needs to be addressed

#compare AICc for 2, 3, or 4 underlying trends
model.list<-list(m=2, R="diagonal and equal")
dfa.fit2<-MARSS(plank.t, model=model.list, z.score=T, form='dfa', control=cntl.list)

#this was in the MARSS user guide but the object saved.res was never defined so I am skipping it
if(!saved.res) {
  model.list<-list(m=2, R="diagonal and equal")
  dfa.fit2<-MARSS(plank.t, model=model.list, z.score=T, form="dfa", control=big.maxit.cntl.list)
}

model.list<-list(m=4, R="diagonal and equal")
dfa.fit3<-MARSS(plank.t, model=model.list, z.score=T, form='dfa', control=cntl.list)


print(cbind(model=c("2 trends", "3 trends", "4 trends"), AICc=round(c(dfa.fit1$AICc, dfa.fit2$AICc, dfa.fit3$AICc))), quote=F)

#4 trends has lowest AIC
model    AICc
[1,] 2 trends 4882
[2,] 3 trends 5016
[3,] 4 trends 4735
