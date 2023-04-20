
library(tidyverse)
options(dplyr.summarise.inform = FALSE)


load(here::here("Lab-2", "Data_Images", "columbia-river.rda"))

esu <- unique(columbia.river$esu_dps)
esu


plotesu <- function(esuname){
  df <- columbia.river %>% subset(esu_dps %in% esuname)
  ggplot(df, aes(x=spawningyear, y=log(value), color=majorpopgroup)) + 
    geom_point(size=0.2, na.rm = TRUE) + 
    theme(strip.text.x = element_text(size = 3)) +
    theme(axis.text.x = element_text(size = 5, angle = 90)) +
    facet_wrap(~esapopname) +
    ggtitle(paste0(esuname, collapse="\n"))
}


#John Day Lower Mainstream Tributaries --test on this data and see how well you can track the other John Day pops

#Cascades 1980's....what do we connect to? Disjunt from John Day 

#Start with using John Day data to inform Cascade interpolation

plotesu(esu[1])


esuname <- esu[1]

dat <- columbia.river %>% 
  subset(esu_dps == esuname) %>% # get only this ESU
  mutate(log.spawner = log(value)) %>% # create a column called log.spawner
  select(esapopname, spawningyear, log.spawner) %>% # get just the columns that I need
  pivot_wider(names_from = "esapopname", values_from = "log.spawner") %>% 
  column_to_rownames(var = "spawningyear") %>% # make the years rownames
  as.matrix() %>% # turn into a matrix with year down the rows
  t() # make time across the columns
# MARSS complains if I don't do this
dat[is.na(dat)] <- NA
dat[is.nan(dat)] <- NA

any(is.nan(dat))
dat[is.infinite(dat)]<- NA #Make any -inf 
any(is.infinite(dat))


tmp <- rownames(dat)
tmp <- stringr::str_replace(tmp, "Steelhead [(]Middle Columbia River DPS[)]", "")
tmp <- stringr::str_replace(tmp, "River - summer", "")
tmp <- stringr::str_trim(tmp)
rownames(dat) <- tmp


# Hypthesis 1 


mod.list1 <- list(
  U = "unequal", #each of the rivers are estimated separately (different U)
  R = "diagonal and equal", #Process errors are all assumed to be the same 
  Q = "diagonal and equal" #Observation error 
)

#We will eventually create matrices that link the population groups that we think
#are related for model fitting 


library(MARSS)
fit1 <- MARSS(dat, model=mod.list1, method="BFGS")
autoplot(fit1)

fit1$states

fit1_smooth<-tsSmooth(fit1)

#look at corrplot
Q1 <- coef(fit1, type = "matrix")$Q
corrmat1 <- diag(1/sqrt(diag(Q1))) %*% Q1 %*% diag(1/sqrt(diag(Q1)))
corrplot(corrmat1)

##################################
#Hypthesis 2
#hypothesis that the four main groups form subpopulations. Random walk is allowed
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

mod.list2 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "equalvarcov",
  Z = Z_mat
)
m2.1 <- MARSS(dat, model = mod.list1)
autoplot(m1)

#look at corrplot
Q2.1 <- coef(m2.1, type = "matrix")$Q
corrmat2.1 <- diag(1/sqrt(diag(Q2.1))) %*% Q2.1 %*% diag(1/sqrt(diag(Q2.1)))
corrplot(corrmat2.1)
#corrplot mirrors what we told MARSS to use as a Q matrix (equal variance and covariance)


#Hypthesis 2.2

U_mat2.2 <- matrix(c("Cascades","JohnDay","Walla","Yakima"),4,1)
#make Z matrix correspond to 4 hidden states
Z_mat2.2 <- matrix(c(rep(c(1,0,0,0),3),
                  rep(c(0,1,0,0),5),
                  rep(c(0,0,1,0),3),
                  rep(c(0,0,0,1),4)),15,4, byrow=TRUE)

mod.list2.2 <- list(
  U = U_mat2,
  R = "diagonal and equal",
  Q = "unconstrained",
  Z = Z_mat2
)
m2.2 <- MARSS(dat, model = mod.list1)
autoplot(m2.2)

#Notes. Four states. Q values on top of line. 
#Confindence Intervals are not emcompassing data, 
#but from the description doesn't seem like an issue

m2.2$states

m2.2_smooth<-tsSmooth(fit)

#look at corrplot
Q2.2 <- coef(m2.2, type = "matrix")$Q
corrmat2 <- diag(1/sqrt(diag(Q2.2))) %*% Q2.2 %*% diag(1/sqrt(diag(Q2.2)))
corrplot(corrmat2)

#######################################
# Hypothesis 3 

#East vs West underlying states. John Day and Cascades vs Walla Walla and Yakima 
# 3.1: Q matrix = equal variance covariance 
# 3.2: Q matrix = unconstrained 
# 3.3: Q matrix = diagonal and equal 


# 3.1: Q matrix = equal variance covariance 
U_mat3 <- matrix(c("North","South"),2,1)
#make Z matrix correspond to 2 hidden states
Z_mat3 <- matrix(c(rep(c(0,1),8),
                     rep(c(1,0),7)),15,2, byrow=TRUE)
                     

mod.list3.1 <- list(
  U = U_mat3,
  R = "diagonal and equal",
  Q = "equalvarcov",
  Z = Z_mat3
)
m3.1 <- MARSS(dat, model = mod.list3.1)
autoplot(m3.1)

# 3.2: Q matrix = unconstrained 

mod.list3.2 <- list(
  U = U_mat3,
  R = "diagonal and equal",
  Q = "unconstrained",
  Z = Z_mat3
)
m3.2 <- MARSS(dat, model = mod.list3.2)
autoplot(m3.2)

# 3.3: Q matrix = diagonal and equal 
mod.list3.3 <- list(
  U = U_mat3,
  R = "diagonal and equal",
  Q = "diagonal and equal",
  Z = Z_mat3
)
m3.3 <- MARSS(dat, model = mod.list3.3)
autoplot(m3.3)

#1. Create estimates of spawner abundance for all missing years and provide estimates of the decline from the historical abundance.
