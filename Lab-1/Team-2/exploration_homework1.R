library(forecast)
bb_data <- readRDS(here::here("Lab-1", "Data_Images", "bristol_bay_data_plus_covariates.rds"))
head(bb_data)
unique(bb_data$system)
#creating list of rivers
rivers <- unique(bb_data$system)
#adding a log abundance and full age column to dataframe
bb_data$log_ret <- log(bb_data$ret)
bb_data$full_age <- bb_data$fw_age + bb_data$o_age

#just data for 4year old fish from a particular river
wood_dt <- bb_data[(bb_data$system == rivers[2]) & (bb_data$full_age == 4),]
unique(wood_dt$o_age)
#split the dataframe by years in the ocean
wood3 <- wood_dt[wood_dt$o_age==3,]
wood2 <- wood_dt[wood_dt$o_age==2,]

#make them time series
wood3.ts <- ts(wood3$log_ret, start = 1963, frequency = 1)
wood2.ts <- ts(wood2$log_ret, start = 1963, frequency = 1)

#check the plot
plot.ts(wood3.ts, ylab = "log(fish counted x1000)", main = "4 yr olds that spent 3 years in Ocean" )
#looks like there is a positive trend to the data and the variance is increasing: not stationary
plot.ts(wood2.ts, ylab = "log(counted fish x1000)", main = "4 yr olds that spent 2 years in Ocean" )
#looks like it could potentially be stationary, not as clear to me

#we need to make our data stationary before we can fit models to it. it can be stationary
#with mean of zero or non-zero, or stationary around a trend

#wood3.ts testing
#run tests for stationarity
tseries::adf.test(wood3.ts)                  # fails to reject: Non-stationary
tseries::adf.test(wood3.ts,k = 0)            #forced it to test AR(1) it does reject: Stationary
tseries::kpss.test(wood3.ts, null = "Level") #not sure this test makes sense, visually there is a trend in the data opposite null hypothesis for kpss it rejects: Non-stationary
tseries::kpss.test(wood3.ts, null = "Trend") #this tested to see if data is stationary around a trend
                                             # just barelt rejects stationarity: meaning it is non-stationary
#wood3 has disagreement on stationairity, so we need to fix that by differencing the data
#figure out how many differences you need to reach stationairity
forecast::ndiffs(wood3.ts, test = "kpss")    # 1 difference needed
forecast::ndiffs(wood3.ts, test = "adf")     # 0 differences needed

#wood2.ts testing
tseries::adf.test(wood2.ts)                  # Rejects: stationary
tseries::adf.test(wood2.ts,k = 0)            #forced it to test AR(1) it does reject: Stationary
tseries::kpss.test(wood2.ts, null = "Level") #opposite null hypothesis for kpss it fails to reject: Stationary
tseries::kpss.test(wood2.ts, null = "Trend") #not sure that this test makes sense, visually data does not appear to have a trend.it fails to reject: stationary
#wood2 appears like its going to be an AR1 stationary
#see if ndiffs says 0, i expect it to
forecast::ndiffs(wood2.ts, test = "kpss")    # 0 needed
forecast::ndiffs(wood2.ts, test = "adf")     # 0 needed

#look at acf and pacf plots for each time series using differencing if needed
#if it is not a combined ARMA model we should see one of these plots decay over a longer time period
#and the other will have just a few significant lags
#telling us the order of AR() or MA()

#acf and pacf for differenced wood3
layout_matrix_1 <- matrix(c(1,2),ncol=2)
layout(layout_matrix_1)
acf(diff(wood3.ts))
pacf(diff(wood3.ts))
#pacf looks like a slow decay indicating it will be and MA model. signif lags at 1, 5, 10 on ACF

#acf and pacf for wood2
layout_matrix_1 <- matrix(c(1,2),ncol=2)
layout(layout_matrix_1)
acf(wood2.ts)
pacf(wood2.ts)
dev.off()
#no significant lags in either plot. This is white noise, no autocorrelation. expecting ARIMA(0,0,0)

#auto fit the models to data
forecast::auto.arima(diff(wood3.ts)) #MA(1) only not sure what to make of the lags at 5 and 10
forecast::auto.arima(wood3.ts)       # is in agreement with need for differencing and MA(1) with drift may be in agreement with non-zero mean
forecast::auto.arima(wood2.ts)       #my expectation panned out that this is just white noise but has non-zero mean which makes sense as this is count data

#look at residuals
#fit the models
fit.wood3 <- forecast::auto.arima(wood3.ts)
fit.wood2 <- forecast::auto.arima(wood2.ts)

#plot of our models
layout_mat_2 <- matrix(c(1,1,2,2), nrow =2, byrow = TRUE)
layout(layout_mat_2)
plot.ts(wood3.ts)
lines(fitted(fit.wood3), col ="red")
plot.ts(wood2.ts)
lines(fitted(fit.wood2), col ="red")

#box test
forecast::checkresiduals(fit.wood3)  #All good news. fail to reject null hypothesis of Ljung-Box, which means not autocorrelated. the ACF plot agrees with no significant lags and the histogram looks like a normal distribution
forecast::checkresiduals(fit.wood2)  #All good news again. fail to reject null of Ljung-Box. no significant lags on ACF. Histogram looks pretty good, two large negative values but otherwise good.

#in the github discussion Eli says we will see cyclic behavior that these models
#will not be able to handle that will present itself as autocorrelation
#She advises that we say in the report that we see that, but we don't see it in these smaller subsets of data.


#creating the training and testing periods for models
train.wood2 <- window(wood2.ts, start = 1963, end = 2010)
test.wood2 <- window(wood2.ts, start = 2011, end = 2020)
fit.train <- auto.arima(train.wood2)
accuracy(fit.train)
accuracy(forecast(fit.train, h = 10), test.wood2)

train.wood3 <- window(wood3.ts, start = 1963, end = 2010)
test.wood3 <- window(wood3.ts, start = 2011, end = 2020)
fit.train3 <- auto.arima(train.wood3)
for3x <- predict(fit.train3, n.ahead = 10)
plot.ts(test.wood3)
points(for3x$pred)
accuracy(fit.train3)
accuracy(forecast(fit.train3, h = 10), test.wood3)
#the training set on wood3 resulted in a lower Root Mean Squared Error than the training model for 
#wood2. I think this makes sense as wood2 was just white noise so there is no structure for the model to use to make an accurate prediction. 