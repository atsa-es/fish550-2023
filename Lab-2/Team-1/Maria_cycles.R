
library(tidyverse)
library(MARSS)

load(here::here("Lab-2", "Data_Images", "columbia-river.rda"))

#team 2 needs to work on lower columbia chinook
esu <- unique(columbia.river$esu_dps)
esu

#esu[5] is lower columbia chinook

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

chin_c_r<-columbia.river %>% subset(esu_dps == esu[5])

chin_c_r %>% 
  filter(majorpopgroup %in% c("Cascade fall", "Cascade late fall", "Cascade spring")) %>% ggplot(aes(x=spawningyear, y=log(value), color=majorpopgroup)) + 
  geom_line(size=2, na.rm = TRUE) + 
  theme(strip.text.x = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 5, angle = 90), axis.title = element_text(size = 14),
        strip.text = element_text(size = 20)) +
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

# let's look at only lower cowlitz, lewis fall and lewis late fall up untill 2000

chin_subset <- chin_newdat[c(2,3,6),1:45]
print(chin_subset)


par(mfrow=c(2,2))
for(i in 1:3){
  acf(chin_subset[i,], na.action=na.pass, main=rownames(chin_subset)[i])
}


TT <- dim(chin_subset)[2] #number of time steps
covariates <- rbind(
  forecast::fourier(ts(1:TT, freq=3), K=1) |> t(),
  forecast::fourier(ts(1:TT, freq=4), K=1) |> t(),
  forecast::fourier(ts(1:TT, freq=5), K=1) |> t()
)


mod.list_cycle <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  D = "unconstrained",
  d = covariates
)

fit_cycle <- MARSS(chin_subset, model=mod.list_cycle, method = "BFGS")


library(broom)
df <- tidy(fit_cycle) %>%
  subset(stringr::str_sub(term,1,1)=="D")
df$lag <- as.factor(rep(3:5, each=6))
df$river <- as.factor(rep(rownames(chin_subset),3))
df$sc <- rep(rep(c("S","C"), each=3), 3)
df$type <- paste0(df$sc,df$lag)

ggplot(df, aes(x=type, y=estimate, col=lag)) + 
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=.2, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~river) +
  ggtitle("The cycle estimates with CIs")
