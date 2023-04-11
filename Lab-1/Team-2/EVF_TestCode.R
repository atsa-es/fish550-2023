library(tidyselect)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
# Load in and rename the Bristol Bay Data
bb_data <- readRDS(here::here("Lab-1", "Data_Images", "bristol_bay_data_plus_covariates.rds"))

#Plot of just the 2.2 age group, separated by river system
bb_data %>% 
  filter(age_group == "2.2") %>% 
  ggplot(aes(x=ret_yr, y=log(ret))) + 
  geom_line() + 
  ggtitle("2.2 Age Group Returns") +
  facet_wrap(~system)

#Plot of just the 1.3 age group, separated by river system
bb_data %>% 
  filter(age_group == "1.3") %>% 
  ggplot(aes(x=ret_yr, y=log(ret))) + 
  geom_line() + 
  ggtitle("1.3 Age Group Returns") +
  facet_wrap(~system)

#Merge the 1.3 and 2.2 age groups together into a 4-year subgroup
subdata <- bb_data %>% 
  mutate(
    fish_age = case_match(
      age_group, 
      c("2.2", "1.3") ~ "4-yr-salmon",
      .default = age_group
    )) %>%
  group_by(system, fish_age, ret_yr) %>%
  summarize(lntotal = log(sum(ret, na.rm=TRUE))) #Unsure about the point of this line
head(subdata) 

#Plot the subgroup of data
subdata %>% 
  filter(fish_age == "4-yr-salmon") %>%
  ggplot(aes(x=ret_yr, y=lntotal)) +
  geom_line() +
  ggtitle("4 Year Old Fish Returns") +
  facet_wrap(~system)

#Create a time series with training and testing groups
dat <- subdata %>%
  filter(system == "Kvichak", fish_age == "4-yr-salmon") %>%
  mutate(year = ret_yr) %>%
  select(year, lntotal)
datts <- ts(dat$lntotal, start=dat$year[1])
train <- window(datts, dat$year[1], dat$year[1]+46)
test <- window(datts, dat$year[1]+47, dat$year[1]+47+10)

plot(train)
plot(test)