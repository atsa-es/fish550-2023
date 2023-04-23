#loading data

load(here::here("Lab-2", "Data_Images", "columbia-river.rda"))

#choosing only lower columbia chinook
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
#Chinook
plotesu(esu[5])


#looking at only cascade populations
chin_c_r<-columbia.river %>% subset(esu_dps == esu[5])
chin_c_r %>% 
  filter(majorpopgroup %in% c("Cascade fall", "Cascade late fall", "Cascade spring")) %>% ggplot(aes(x=spawningyear, y=log(value), color=majorpopgroup)) + 
  geom_point(size=2, na.rm = TRUE) + 
  theme(strip.text.x = element_text(size = 3)) +
  theme(axis.text.x = element_text(size = 5, angle = 90)) +
  facet_wrap(~esapopname) +
  ggtitle(paste0(esu[5], collapse="\n"))


chin_newdat <- chin_c_r %>% 
  filter(majorpopgroup %in% c("Cascade fall", "Cascade late fall", "Cascade spring")) %>% #just looking at cascade populations
  mutate(log.spawner = log(value)) %>% # create a column called log.spawner
  select(esapopname, spawningyear, log.spawner) %>% # get just the columns that I need
  pivot_wider(names_from = esapopname, values_from = log.spawner) %>% 
  column_to_rownames(var = "spawningyear") %>% # make the years rownames
  as.matrix() %>% # turn into a matrix with year down the rows
  t() # make time across the columns
# MARSS complains if I don't do this
chin_newdat[is.na(chin_newdat)] <- NA

#clean up row names
tmp <- rownames(chin_newdat)
tmp <- stringr::str_replace(tmp, "Salmon, Chinook [(]Lower Columbia River ESU[)]", "")
tmp <- stringr::str_trim(tmp)
rownames(chin_newdat) <- tmp

#look at data
print(chin_newdat[,1:5])



#hypothesis 1 - 11 independednt populations

mod.list1 <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)
fit1<-MARSS(chin_newdat, model=mod.list1, method = "BFGS")

plot(fit1, plot.type = "fitted.ytT")

#residuals

plot(fit1, plot.type = "model.resids.ytt1")
plot(fit1, plot.type = "qqplot.std.model.resids.ytt1")
plot(fit1, plot.type = "acf.std.model.resids.ytt1")

# trying a corr plot

mod.list4.unconst <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "unconstrained")
  
fit4.unconst<-MARSS(chin_newdat, model=mod.list4.unconst, method = "BFGS")

library(corrplot)
Q4.unconst <- coef(fit4.unconst, type="matrix")$Q
corrmat4.unconst <- diag(1/sqrt(diag(Q4.unconst))) %*% Q4.unconst %*% diag(1/sqrt(diag(Q4.unconst)))
corrplot(corrmat4.unconst)


#Let's look at ACF of the data for each river

#make ts object

for(i in 1:length(chin_newdat[,1])){
  ts <- ts(chin_newdat[i,], start = 1964, end = 2021)
  acf(ts, na.action = na.pass)
}
#except for the last river, all rivers have some significant acf at lag>0
#most of them are 4

#let's do this for #Lower Cowlitz river - fall

# setting up the data

yt <- chin_newdat[3,]
TT <- length(yt)
p <- 4

#specifying Z

Z <- array(1, dim = c(1, 3, TT))
Z[1, 2, ] <- sin(2 * pi * (1:TT)/p)
Z[1, 3, ] <- cos(2 * pi * (1:TT)/p)

#specifying model list
mod.list <- list(U = "zero", Q = "diagonal and unequal", Z = Z, 
                 A = "zero")

m <- dim(Z)[2]
fit <- MARSS(yt, model = mod.list, inits = list(x0 = matrix(0, 
                                                            m, 1)))
plot(fit, plot.type = "xtT")
#no cycling?

plot(chin_newdat[11,],type = 'l')

fitriver <- function(i, p = 5) {
  yt <- chin_newdat[i,]
  TT <- length(yt)
  Z <- array(1, dim = c(1, 3, TT))
  Z[1, 2, ] <- sin(2 * pi * (1:TT)/p)
  Z[1, 3, ] <- cos(2 * pi * (1:TT)/p)
  mod.list <- list(U = "zero", Q = "diagonal and unequal", 
                   Z = Z, A = "zero")
  fit <- MARSS(yt, model = mod.list, inits = list(x0 = matrix(0, 
                                                              3, 1)), silent = TRUE)
  return(fit)
}

fits <- list()
for (i in 1:11){
  fits[[i]] <- fitriver(i)
}


dfz <- data.frame()
for (i in 1:11) {
  fit <- fits[[i]]
  tmp <- data.frame(amplitude = sqrt(fit$states[2, ]^2 + fit$states[3, 
  ]^2), trend = fit$states[1, ], river = rownames(chin_newdat)[i], year = as.numeric(colnames(chin_newdat)))
  dfz <- rbind(dfz, tmp)
}

ggplot(dfz, aes(x = year, y = amplitude)) + geom_line() + 
  facet_wrap(~river, scales = "free_y") + ggtitle("Cycle Amplitude")


ggplot(dfz, aes(x = year, y = trend)) + 
  geom_line() + facet_wrap(~river, scales = "free_y") + ggtitle("Stochastic Level")
