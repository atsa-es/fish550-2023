#Task 1: Estimate the missing values and come up with estimates of "historical" population size 
  #Can compare with different MARSS model formulations
  #will need to make some assumption for how these are related

#Task 2: Decompose the VCV to look at correlation in best model

library(tidyverse)
load(here::here("Lab-2", "Data_Images", "columbia-river.rda"))
esu <- unique(columbia.river$esu_dps)
esu #Group 1: Lower Columbia River Chinook which is index 5

#just chinook
chin_c_r<-columbia.river %>% subset(esu_dps == esu[5])
summary(chin_c_r)


chin_newdat <- chin_c_r %>% 
  mutate(log.spawner = log(value)) %>% # create a column called log.spawner
  select(esapopname, spawningyear, value, log.spawner) #%>% # get just the columns that I need
  #group_by(spawningyear) %>%
  #pivot_wider(names_from = esapopname, values_from = log.spawner) #%>% 
  #column_to_rownames(var = "spawningyear") %>% # make the years rownames
  #as.matrix() %>% # turn into a matrix with year down the rows
  #t() # make time across the columns
# MARSS complains if I don't do this

#some esapopnames and years have two values...but all of these are either both NA, or only have one count, so can probably reduce down

duplicates<-
  chin_c_r %>% 
  mutate(log.spawner = log(value)) %>% # create a column called log.spawner
  select(esapopname, spawningyear, value, log.spawner) %>% group_by(esapopname, spawningyear) %>%
  summarise(n = n()) %>% filter(n > 1)

for(i in 1:nrow(duplicates)){
  print(chin_newdat %>% filter(esapopname == duplicates$esapopname[i]) %>% 
    filter(spawningyear == duplicates$spawningyear[i]))
}

ifelse(sum(is.na(.x)) == 2, NA, sum(.x, na.rm = TRUE))

test<-chin_newdat %>%  select(esapopname, spawningyear, log.spawner) %>%
  pivot_wider(names_from = esapopname, values_from = log.spawner, values_fn = ~ifelse(sum(is.na(.x)) == 2, "NA", sum(.x, na.rm = TRUE))) #%>% 
#column_to_rownames(var = "spawningyear") %>% # make the years rownames
#as.matrix() %>% # turn into a matrix with year down the rows
#t() # make time across the columns

#removing this population
"Salmon, Chinook (Lower Columbia River ESU) Lower Gorge Tributaries - fall"
#replace the one 0 with 1
#then na.rm = TRUE
dat[is.na(dat)] <- NA



#We assume that the time series from each population has identical and independent observation error variance
#this means that diagonal R matrix with one variance term on the diagonal


#plots: From Emma's code: 


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


#subset to just look at Cascade major population group

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

#hypotheses: 
#1) All independent populations:
#1x for 1y and then 
#diagonal and unequal process errors - environment impacts all populations differently
#daigonal and equal process errors  -environment impacts all populations the same
#2) 3 populations according to major popualation group
#1x for each run time 
#3) All one population
#just one process error
#each will allow for independent observation errors


mod.list3 <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "equalvarcov"
)
fit3<-MARSS(chin_newdat, model=mod.list3, control = list(maxit = 2000))

autoplot(fit3, plot.type="fitted.ytT")
autoplot(fit3, plot.type = "residuals") #some of these residual plots don't look great

aic <- c(fit1$AICc, fit2$AICc, fit3$AICc)
#delta aic
aic-min(aic) #model 1 has the best AICc of these


#3 populations
mod.list4<-mod.list1
rownames(chin_newdat)
mod.list4$Z<-factor(c(rep("fa", 5), rep("l_fa", 2), rep("Sp", 4)))
mod.list4

fit4<-MARSS(chin_newdat, model = mod.list4, method = "BFGS")
autoplot(fit4, plot.type="fitted.ytT")
autoplot(fit4, plot.type = "residuals") #some of these residual plots don't look great but most don't seem too bad?

mod.list5<-mod.list4
mod.list5$Q<-"diagonal and equal"
fit5<-MARSS(chin_newdat, model = mod.list5, control = list(maxit = 2000))
autoplot(fit5, plot.type="fitted.ytT")
autoplot(fit5, plot.type = "residuals") #some of these residual plots don't look great but most don't seem too bad?


mod.list6<-mod.list4
mod.list6$Q<-"equalvarcov"
fit6<-MARSS(chin_newdat, model = mod.list6, control = list(maxit = 2000))
autoplot(fit6, plot.type="fitted.ytT")
autoplot(fit6, plot.type = "residuals") #some of these residual plots don't look great but most don't seem too bad?

aic <- c(fit1$AICc, fit2$AICc, fit3$AICc, fit4$AICc, fit5$AICc, fit6$AICc)
mods<-seq(1,6)
aic.names<-paste("Model", mods, sep = " ")
names(aic)<-aic.names
aic-min(aic) #model 1 has the best AICc of these
a<-matrix(c(aic, aic-min(aic)), nrow = 6, byrow = F)
rownames(a)<-aic.names
a<-a[order(a[,2],decreasing=FALSE),]
rownames(a)<-aic.names[order(a[,2],decreasing=FALSE),]
knitr::kable(a, col.names = c("AICc", "Delta AIC"))

fit1$states
str(fit1$states)
#%change from historical abundance 
#using log scale-change in log abundance 
hist.abund<-tibble(Pop = rownames(chin_newdat), 
                   Year1 = fit1$states[,1],
                   Year58 = fit1$states[,58])
                   
hist.abund<-hist.abund %>% mutate(PChange = ((Year58 - Year1)/abs(Year1)) * 100)

ggplot(hist.abund) + geom_histogram(aes(x = Pop, y = PChange), stat = "identity")

hist.abund %>% filter(Pop != "Upper Cowlitz River - fall") %>%
  ggplot() + geom_histogram(aes(x = Pop, y = PChange, fill = factor(c(rep("fa", 5), rep("l_fa", 2), rep("Sp", 3)))), stat = "identity") + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  scale_fill_discrete(labels = c("Fall", "Late Fall", "Spring"), name = "Run") + 
  labs(x = "Population", y="% Change From Hist. Abund.") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 270, hjust = -0.1))

        