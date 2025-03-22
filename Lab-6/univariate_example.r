library(TMB)
library(ggplot2)
compile("Lab-6/univariate.cpp")
dyn.load(dynlib("Lab-6/univariate"))

x = cumsum(rnorm(30))
y = x + rnorm(length(x), 0, 0.01)
estimate_drift = TRUE # U in MARSS
estimate_rho = FALSE # AR(1) parameter, b in MARSS

parameters <- list(
  log_obs_sd = 0,
  log_pro_sd = 0,
  x = rep(0, length(y)),
  u = 0,
  #x0 = 0,
  logit_rho = 0
)

# Map off parameters not being estimted
tmb_map <- list(x = as.factor(c(NA,1:(length(y)-1))))
if(estimate_drift == FALSE) tmb_map <- c(tmb_map, list(u = factor(NA)))
if(estimate_rho == FALSE) tmb_map <- c(tmb_map, list(logit_rho = factor(NA)))

# Create TMB data
data_list <- list(Y = y, n = length(y),
                  est_drift = as.numeric(estimate_drift),
                  est_rho = as.numeric(estimate_rho),
                  keep = ifelse(!is.na(y),1,0))

# Create object for fitting
obj <- TMB::MakeADFun(
  data = data_list,
  map = tmb_map,
  random = "x",
  parameters = parameters,
  DLL = "univariate",
  silent = TRUE
)

# Do the fitting with stats::nlminb, sometimes need to change default control args if not converging
pars <- stats::nlminb(
  start = obj$par, objective = obj$fn,
  gradient = obj$gr
)

par_summary <- summary(sdreport(obj))

indx <- grep("pred", rownames(par_summary))
df <- data.frame(
  pred = as.numeric(par_summary[indx,"Estimate"]),
  se = as.numeric(par_summary[indx,"Std. Error"]),
  y = y,
  t = 1:length(y)
)

ggplot(df, aes(t, pred)) + 
  geom_ribbon(aes(ymin=pred-2*se, ymax = pred+2*se),alpha=0.5) + 
  geom_line() + 
  geom_point(aes(t,y),col="red",alpha=0.5) + 
  xlab("Time") + ylab("Data") + 
  theme_bw()
