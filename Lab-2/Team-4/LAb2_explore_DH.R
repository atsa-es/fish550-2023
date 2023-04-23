#dependencies
library(tidyverse)
library(dplyr)
library(MARSS)
library(corrplot)

#load data
load(here::here("Lab-2", "Data_Images", "columbia-river.rda"))
dat <- columbia.river

#we are only interested in the Middle Columbia
#plot the unique Rivers in Middle Columbia
esuname <- unique(dat$esu_dps)
plotesu <- function(esuname){
  df <- dat %>% subset(esu_dps %in% esuname)
  ggplot(df, aes(x=spawningyear, y=log(value), color=majorpopgroup)) + 
    geom_point(size=1, na.rm = TRUE) + 
    theme(strip.text.x = element_text(size = 8)) +
    theme(axis.text.x = element_text(size = 8, angle = 90)) +
    facet_wrap(~esapopname) +
    ggtitle(paste0(esuname, collapse="\n"))
}
plotesu(esuname[1])

#Prepare our data so the columns are the years and rows are unique rivers
esuname <- esuname[1]
dat <- columbia.river %>% 
  subset(esu_dps == esuname) %>% # get only this ESU
  mutate(log.spawner = log(value)) %>% # create a column called log.spawner
  dplyr::select(esapopname, spawningyear, log.spawner) %>% # get just the columns that I need
  pivot_wider(names_from = "esapopname", values_from = "log.spawner") %>% 
  column_to_rownames(var = "spawningyear") %>% # make the years rownames
  as.matrix() %>% # turn into a matrix with year down the rows
  t() # make time across the columns
# MARSS complains if I don't do this
dat[is.na(dat)] <- NA
dat
#clean the rownames I am going to leave "Summer/Winter" on them
tmp <- rownames(dat)
tmp <- stringr::str_replace(tmp, "Steelhead [(]Middle Columbia River DPS[)]", "")
tmp <- stringr::str_replace(tmp, " - summer", "-S")
tmp <- stringr::str_replace(tmp, " - winter", "-W")

tmp <- stringr::str_trim(tmp)
rownames(dat) <- tmp

#fill any blanks with NA
any(is.null(dat))
any(is.infinite(dat))
dat[is.infinite(dat)] <- NA
#need to make your plot window as big as you can to print this
#nothing too wild going on in the data
layout_mat <- matrix(c(1:15,0),4,4)
layout(layout_mat)
for (i in 1:15) {
  hist(dat[i,])
}
dev.off()
dat
############################################################3

#Hypothesis 2.1: hypothesis that the four main groups form subpopulations. Random walk is allowed
#to drift uniquely in each of the 4 hidden states based on U. The Q matrix for
#variance of process errors is diagonal and equal meaning all have the same variance
#but they are not correlated to each other

#give U values names to make it easier to read results
#this hypothesis has 4 hidden states based on major groups
U_mat <- matrix(c("Cascades","JohnDay","Walla","Yakima"),4,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0,0,0),3),
                  rep(c(0,1,0,0),5),
                  rep(c(0,0,1,0),3),
                  rep(c(0,0,0,1),4)),15,4, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "diagonal and equal",
  Z = Z_mat
)
m2.1 <- MARSS(dat, model = mod.list1)
#notes: again yakima is the only positive drift value on the hidden state's random walk
#       Only one Q value output as all hidden states have the same one and there is no
#       covariance/correlation between states
autoplot(m2.1)
#huge ballooned CIs on the missing data in hidden states
#fitted CI have the same balloon shaped CIs on missing data, I would predict this
#model does poorly when we compare AICs
#maybe some structure in residuals as well


#look at corrplot
Q2.1 <- coef(m2.1, type = "matrix")$Q
corrmat2.1 <- diag(1/sqrt(diag(Q2.1))) %*% Q2.1 %*% diag(1/sqrt(diag(Q2.1)))
corrplot(corrmat2.1)
#As expected output displays only diagonal as we told it diagonal and equal


############################################################
#Hypothesis 2.2: hypothesis that the four main groups form subpopulations. Random walk is allowed
#to drift uniquely in each of the 4 hidden states based on U. The Q matrix for
#variance of process errors is diagonal and unequal meaning variance can be different
#but they are not correlated to each other

#give U values names to make it easier to read results
#this hypothesis has 4 hidden states based on major groups
U_mat <- matrix(c("Cascades","JohnDay","Walla","Yakima"),4,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0,0,0),3),
                  rep(c(0,1,0,0),5),
                  rep(c(0,0,1,0),3),
                  rep(c(0,0,0,1),4)),15,4, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  Z = Z_mat
)
m2.2 <- MARSS(dat, model = mod.list1)
#notes: again yakima is the only positive drift value on the hidden state's random walk
#       4 varied Q values output as all hidden states have a different one and there is no
#       covariance/correlation between states

autoplot(m2.2)
#huge ballooned CIs on the missing data in hidden states
#fitted CI have the same balloon shaped CIs on missing data and is quite large 
#on the Walla Walla group

#look at corrplot
Q2.2 <- coef(m2.2, type = "matrix")$Q
corrmat2.2 <- diag(1/sqrt(diag(Q2.2))) %*% Q2.2 %*% diag(1/sqrt(diag(Q2.2)))
corrplot(corrmat2.2)

###############################################################################################

#Hypothesis 2.3: hypothesis that the four main groups form subpopulations. Random walk is allowed
#to drift uniquely in each of the 4 hidden states based on U. The Q matrix for
#variance of process errors have equal variance and covariance so they fluctuate
#together

#give U values names to make it easier to read results
#this hypothesis has 4 hidden states based on major groups
U_mat <- matrix(c("Cascades","JohnDay","Walla","Yakima"),4,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0,0,0),3),
         rep(c(0,1,0,0),5),
         rep(c(0,0,1,0),3),
         rep(c(0,0,0,1),4)),15,4, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "equalvarcov",
  Z = Z_mat
)
m2.3 <- MARSS(dat, model = mod.list1)
#Notes: Yakima group's hidden state is the only one with a positive U, indicating this random 
#       walk has a positive drift while the others have a negative

autoplot(m2.3)

#look at corrplot
Q2.3 <- coef(m2.3, type = "matrix")$Q
corrmat2.3 <- diag(1/sqrt(diag(Q2.3))) %*% Q2.3 %*% diag(1/sqrt(diag(Q2.3)))
corrplot(corrmat2.3)
#corrplot mirrors what we told MARSS to use as a Q matrix (equal variance and covariance)

############################################################################
#Hypothesis 2.4: hypothesis that the four main groups form subpopulations. Random walk is allowed
#to drift uniquely in each of the 4 hidden states based on U. The Q matrix for
#variance of process errors is unconstrained and each hidden state can vary separately.
#this will report the covariances as well so we can see how the major groups relate to each other

#give U values names to make it easier to read results
#this hypothesis has 4 hidden states based on major groups
U_mat <- matrix(c("Cascades","JohnDay","Walla","Yakima"),4,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0,0,0),3),
                  rep(c(0,1,0,0),5),
                  rep(c(0,0,1,0),3),
                  rep(c(0,0,0,1),4)),15,4, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "unconstrained",
  Z = Z_mat
)
m2.4 <- MARSS(dat, model = mod.list1)
#notes: again yakima is the only positive drift value on the hidden state's random walk
#       John Day and Yakima have strongest correlation and then John Day and Cascades
#       I expect this to be seen in the corrplot
autoplot(m2.4)
#fitted CI look decent, Deschutes and Fifteenmile are the worst otherwise
#the CI doesn't blow up on others where data was missing. Apparently not a good
#for the Cascades group
#Standardized residuals may show a bit of structure in a few graphs but otherwise 
#white noise

#look at corrplot
Q2.4 <- coef(m2.4, type = "matrix")$Q
corrmat1 <- diag(1/sqrt(diag(Q2.4))) %*% Q2.4 %*% diag(1/sqrt(diag(Q2.4)))
corrplot(corrmat1)
#not at all what I was expecting from the model output, 
#All groups are highly correlated with each other

#################################################################################
#Hypothesis 4: Salmon of the Yakima group have to swim the furthest to reach their 
#spawning ground, including a large bend in the river that heads back west. They
#are the most isolated group and thus may have their own hidden state while the 
#other 3 major population groups maybe more more closely linked to each other.
#################################################################################

#Hypothesis 4.1: Yakima salmon have separate hidden state and Q matrix is 
#diagonal and equal, this will prevent the correlation of the two hidden states
#but they will have the same variance

#give U values names to make it easier to read results
#this hypothesis has 2 hidden states based spatial distribution of major groups
U_mat <- matrix(c("West_group","Yakima"),2,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0),3),
                  rep(c(1,0),5),
                  rep(c(1,0),3),
                  rep(c(0,1),4)),15,2, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "diagonal and equal",
  Z = Z_mat
)
m4.1 <- MARSS(dat, model = mod.list1)
#notes: yakima continues its positive drift, collective group shows negative drift
# single Q value reported as "diagonal and equal"
autoplot(m4.1)
#HUGE CIs on hidden state of Yakima group where data is missing. Less missing data 
#for the other hidden state but CI do get big at the last time steps where data 
#CIs for fitted values show the same pattern
#qqplots show more variation in the Yakima group then I have noticed before.
#North Fork John Daly is continually the worst performer

#look at corrplot
Q4.1 <- coef(m4.1, type = "matrix")$Q
corrmat4.1 <- diag(1/sqrt(diag(Q4.1))) %*% Q4.1 %*% diag(1/sqrt(diag(Q4.1)))
corrplot(corrmat4.1)
#As expected for this Q call

##############################################################################

#Hypothesis 4.2: Yakima salmon have separate hidden state and Q matrix is 
#diagonal and unequal, this will prevent the correlation of the two hidden states
#and the same variance, they were still very similar in the unconstrained hypothesis
#so maybe see the same thing

#give U values names to make it easier to read results
#this hypothesis has 2 hidden states based spatial distribution of major groups
U_mat <- matrix(c("West_group","Yakima"),2,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0),3),
                  rep(c(1,0),5),
                  rep(c(1,0),3),
                  rep(c(0,1),4)),15,2, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  Z = Z_mat
)
m4.2 <- MARSS(dat, model = mod.list1)
#notes: yakima continues its positive drift, collective group shows negative drift
# Q values are similar but not as close as they were in the unconstrained model
autoplot(m4.2)
#HUGE, HUGE! CIs on hidden state of Yakima group where data is missing. Less missing data 
#for the other hidden state but CI do get big at the last time steps where data 
#CIs for fitted values show the same pattern
#qqplots show more variation in the Yakima group then I have noticed before.
#North Fork John Daly is continually the worst performer

#look at corrplot
Q4.2 <- coef(m4.2, type = "matrix")$Q
corrmat4.2 <- diag(1/sqrt(diag(Q4.2))) %*% Q4.2 %*% diag(1/sqrt(diag(Q4.2)))
corrplot(corrmat4.2)
#As expected for this Q call

####################################################################################

#Hypothesis 4.3: Yakima salmon have separate hidden state and Q matrix is 
#Equal variance and covariance,

#maybe be very similar to outcome of Hyp4.1


#give U values names to make it easier to read results
#this hypothesis has 2 hidden states based spatial distribution of major groups
U_mat <- matrix(c("West_group","Yakima"),2,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0),3),
                  rep(c(1,0),5),
                  rep(c(1,0),3),
                  rep(c(0,1),4)),15,2, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "equalvarcov",
  Z = Z_mat
)
m4.3 <- MARSS(dat, model = mod.list1)
#notes: yakima continues its positive drift, collective group shows negative drift
# Both the variance and covariance almost have the same value. likely will have same
#results as hypothesis 4.1
autoplot(m4.3)
#positive and negative drifts apparent in States graphs. Decent looking CIs too
#Pretty good looking CIs on the fitted values
#maybe some structure in residuals, qqplots ok: North Fork John Daly has heavy tails
#results look very similar to 4.1

#look at corrplot
Q4.3 <- coef(m4.3, type = "matrix")$Q
corrmat4.3 <- diag(1/sqrt(diag(Q4.3))) %*% Q4.3 %*% diag(1/sqrt(diag(Q4.3)))
corrplot(corrmat4.3)
#As expected almost the exact results of the 4.1 hypothesis

#################################################################################
#Hypothesis 4.4: Yakima salmon have separate hidden state and Q matrix is unconstrained
#We will be estimating a 2X2 Q matrix in which variance and covariance is independent


#give U values names to make it easier to read results
#this hypothesis has 2 hidden states based spatial distribution of major groups
U_mat <- matrix(c("West_group","Yakima"),2,1)
#make Z matrix correspond to 4 hidden states
Z_mat <- matrix(c(rep(c(1,0),3),
                  rep(c(1,0),5),
                  rep(c(1,0),3),
                  rep(c(0,1),4)),15,2, byrow=TRUE)

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "unconstrained",
  Z = Z_mat
)
m4.4 <- MARSS(dat, model = mod.list1)
#notes: yakima continues its positive drift, collective group shows negative drift
# Q values are very similar across variance and covariance

autoplot(m4.4)
#positive and negative drifts apparent in States graphs. Decent looking CIs too
#Pretty good looking CIs on the fitted values
#maybe some structure in residuals, qqplots ok: North Fork John Daly has heavy tails

#look at corrplot
Q4.4 <- coef(m4.4, type = "matrix")$Q
corrmat4.4 <- diag(1/sqrt(diag(Q4.4))) %*% Q4.4 %*% diag(1/sqrt(diag(Q4.4)))
corrplot(corrmat4.4)
#As hinted at by the model out put very high correlation even though this was an
#unconstrained model

##################################################################################

#compare the AIC values
mods <- c("2.1","2.2","2.3","2.4","4.1","4.2","4.3","4.4")
aic <- c(m2.1$AICc, m2.2$AICc, m2.3$AICc, m2.4$AICc, m4.1$AICc, m4.2$AICc, m4.3$AICc, m4.4$AICc)
daic <- aic-min(aic)
tab <- cbind.data.frame(mods, aic, daic)
kable(tab, col.names = c("Hypothesis", "AICc", "delta AICc"))
