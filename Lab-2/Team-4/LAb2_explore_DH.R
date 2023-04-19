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

mod.list1 <- list(
  U = U_mat,
  R = "diagonal and equal",
  Q = "equalvarcov",
  Z = Z_mat
)
m1 <- MARSS(dat, model = mod.list1)
autoplot(m1)

#look at corrplot
Q1 <- coef(m1, type = "matrix")$Q
corrmat1 <- diag(1/sqrt(diag(Q1))) %*% Q1 %*% diag(1/sqrt(diag(Q1)))
corrplot(corrmat1)
#corrplot mirrors what we told MARSS to use as a Q matrix (equal variance and covariance)