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
