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
plot.ts(wood2.ts, ylab = "log(counted fish x1000)", main ="4 year old sockeye returning to Wood River\nthat spent 2 years in Ocean" )
#looks like it could potentially be stationary
plot.ts(wood3.ts, ylab = "log(fish counted x1000)", main ="4 year old sockeye returning to Wood River\nthat spent 3 years in Ocean" )
#looks like there is a positive trend to the data: not stationary, or could be stationary around a trend


#we need to make our data stationary before we can fit models to it. it can be stationary
#with mean of zero or non-zero, or stationary around a trend

#wood3.ts testing
#run tests for stationarity
tseries::adf.test(wood3.ts)                  # fails to reject: Non-stationary
tseries::adf.test(wood3.ts,k = 0)            #forced it to test AR(1) it does reject: Stationary. So conflicts with the test above
tseries::kpss.test(wood3.ts, null = "Level") #not sure this test makes sense, visually there is a trend in the data. opposite null hypothesis for kpss it rejects: Non-stationary
tseries::kpss.test(wood3.ts, null = "Trend") #this tested to see if data is stationary around a trend
# just barely rejects stationarity: meaning it is non-stationary
#wood3 has disagreement on stationairity, so we need to fix that by differencing the data
#figure out how many differences you need to reach stationairity
forecast::ndiffs(wood3.ts, test = "kpss")    # 1 difference needed
forecast::ndiffs(wood3.ts, test = "adf")     # 0 differences needed, again disagreement. Likely will need a difference but we can compare to the ARIMA model that gets automatically fit

#wood2.ts testing
tseries::adf.test(wood2.ts)                  # Rejects: stationary
tseries::adf.test(wood2.ts,k = 0)            #forced it to test AR(1) it does reject: Stationary
tseries::kpss.test(wood2.ts, null = "Level") #opposite null hypothesis for kpss it fails to reject: Stationary
tseries::kpss.test(wood2.ts, null = "Trend") #not sure that this test makes sense (talked to Mark about this he said it can still be good to check as the trend could be very subtle), visually data does not appear to have a trend.it fails to reject: stationary
#wood2 has full agreement on stationairity from all tests used. Very likely is stationary and will not require differencing
#see if ndiffs says 0, i expect it to
forecast::ndiffs(wood2.ts, test = "kpss")    # 0 needed
forecast::ndiffs(wood2.ts, test = "adf")     # 0 needed. Agrees with the unit root tests above that found stationairity of the data

#look at acf and pacf plots for each time series using differencing if needed
#if it is not a combined ARMA model we should see one of these plots decay over a longer time period
#and the other will have just a few significant lags
#telling us the order of AR() or MA()

#acf and pacf for differenced wood3
layout_matrix_1 <- matrix(c(1,2),ncol=2)  #setting up how the plots will display
layout(layout_matrix_1)                   #setting up how the plots will display
acf(diff(wood3.ts))
pacf(diff(wood3.ts))
#pacf looks like a slow decay indicating it will be and MA model. signif lags at 1, 5, 10 on ACF

#acf and pacf for wood2
#carries the plot display set up from the above over to these plots
acf(wood2.ts)
pacf(wood2.ts)
dev.off()        #resets the plot window so it should only show a single plot at a time
#no significant lags in either plot. This is white noise, no autocorrelation. expecting ARIMA(0,0,0)

#auto fit the models to data
forecast::auto.arima(diff(wood3.ts)) #MA(1) only not sure what to make of the lags at 5 and 10. could be the cyclic nature of the fish Eli has mentioned (see note further down)
forecast::auto.arima(wood3.ts)       # is in agreement with need for differencing and MA(1). I think "with drift" in the non-differenced data is the equivalent of "with non-zero mean" in the pre-differenced model above.
forecast::auto.arima(wood2.ts)       #my expectation panned out that this is just white noise but has non-zero mean which makes sense as this is count data

#look at residuals
#fit and save the models
fit.wood3 <- forecast::auto.arima(wood3.ts)
fit.wood2 <- forecast::auto.arima(wood2.ts)

#plot of our models
layout_mat_2 <- matrix(c(1,1,2,2), nrow =2, byrow = TRUE)
layout(layout_mat_2)
plot.ts(wood2.ts, main ="4 year old sockeye returning to Wood River\nthat spent 2 years in Ocean", ylab = "log(counted fish x1000")
lines(fitted(fit.wood2), col ="red")
plot.ts(wood3.ts, main ="4 year old sockeye returning to Wood River\nthat spent 3 years in Ocean", ylab = "log(counted fish x1000")
lines(fitted(fit.wood3), col ="red")

#our fitted models are able to put a line through the data. The wood3 model was
#the undifferenced data and so our model has drift. The wood2 model is just a straight line that I think
#we will see down below that the predicted values is just the mean around which the says the data is white noise (4.6514)

#box test
forecast::checkresiduals(fit.wood3)  #All good news. fail to reject null hypothesis of Ljung-Box, which means the residuals are not autocorrelated. the ACF plot agrees with no significant lags and the histogram looks like a normal distribution
forecast::checkresiduals(fit.wood2)  #All good news again. fail to reject null of Ljung-Box. no significant lags on ACF. Histogram looks pretty good, two large negative values but otherwise good.

#in the github discussion Eli says we will see cyclic behavior that these models
#will not be able to handle that will present itself as autocorrelation
#She advises that we say in the report that we see that. That may have been the significant lags seen in the
#wood3 pacf() plot. but it doesn't seem to be showing itself in the residuals of our models


#creating the training and testing periods for models
train.wood2 <- window(wood2.ts, start = 1963, end = 2010)
test.wood2 <- window(wood2.ts, start = 2011, end = 2020)
fit.train <- auto.arima(train.wood2)
pred.wood2 <- forecast(fit.train, h = 10)
plot(pred.wood2,ylab = "log(fish counted x1000)", xlab = "Time")
points(test.wood2, pch=19, col="red")
#the actual observations do fall within the 80% prediction interval
#This data is treated as white noise with the ARIMA(0,0,0) model so 
#I imagine the prediction intervals are tied very closely to the variance 
#of the residuals themselves.
residuals <- (test.wood2[1:10] - pred.wood2$mean[1:10])
RMSE.wood2 <- sqrt(mean(residuals^2))
accuracy(forecast(fit.train, h = 10), test.wood2)
#the RMSE done manually matches with the RMSE value that is found by the 
#model in row 2, column 2. This must be how the model does when compared to the 
#test data. So we have RMSE = 0.731

train.wood3 <- window(wood3.ts, start = 1963, end = 2010)
test.wood3 <- window(wood3.ts, start = 2011, end = 2020)
fit.train3 <- auto.arima(train.wood3)

pred.wood3 <- forecast(fit.train3, h = 10)
plot(pred.wood3,ylab = "log(fish counted x1000)", xlab = "Time")
points(test.wood3, pch=19, col="red")
#we see that the prediction intervals widen as the forecast moves forward
#likely due to the autocorrelation of the errors, because this is an MA(1) model,
#which the model has to estimate with each successive prediction. The fitted points
#are a flat line again because this is an MA(1) model, we would write the model as
# Xt = Xt-1 + error. When we fitted a model to the entire dataset, auto.arima() 
#selected an ARIMA(0,1,1) with drift, however with just the training data auto.arima() 
#selected an ARIMA(0,1,1) because the model does not have drift, the predictions stay as a flat line
# rather than having a trend.
residuals.3 <- test.wood3[1:10] - pred.wood3$mean[1:10]
sqrt(mean(residuals.3^2))
accuracy(forecast(fit.train3, h = 10), test.wood3)
#the Root Mean Squared Error for this model is 0.532. This value is less
#than the value for the wood2 model. This is a bit surprising as they both
#put a straight line on the graph for predictions. These values indicate that 
#the model for wood.3 was more accurate. The returns of fish that spent less time in fresh
#water and longer in the ocean were more accurately predicted than of the fish who spent
#two years in both freshwater and the ocean.
dev.off()



#look at one more river system because I have time.
#subset data from Kvichak region
Kvichak_dt <- bb_data[(bb_data$system == rivers[4]) & (bb_data$full_age == 4),]
unique(Kvichak_dt$o_age)
kvi.3 <- Kvichak_dt[Kvichak_dt$o_age == 3,]
kvi.2 <- Kvichak_dt[Kvichak_dt$o_age == 2,]

#make the ts
kvi.3.ts <- ts(kvi.3$log_ret, start = 1963, end = 2020)
kvi.2.ts <- ts(kvi.2$log_ret, start = 1963, end = 2020)

#plots
layout(layout_mat_2)
plot(kvi.2.ts, ylab = "log(fish counted x1000)", main = "4 year old sockey returning to Kvichak River\nthat spent 2 years in Ocean" )
#looks fairly stationary to me
plot(kvi.3.ts, ylab = "log(fish counted x1000)", main = "4 year old sockeye returning to Kvichak River\nthat spent 3 years in Ocean" )
#one huge spike downward at start otherwise variance looks similar througout. maybe stationary around a trend if the initial spike doesn't ruin it

dev.off()

#test stationairity kvi.2
tseries::adf.test(kvi.2.ts)                  # fails to Reject: non-stationary
tseries::adf.test(kvi.2.ts,k = 0)            #forced it to test AR(1) it does reject: Stationary. disagreement between them
tseries::kpss.test(kvi.2.ts, null = "Level") #opposite null hypothesis for kpss it fails to reject: Stationary
tseries::kpss.test(kvi.2.ts, null = "Trend") #not sure that this test makes sense (talked to Mark about this he said it can still be good to check as the trend could be very subtle), visually data does not appear to have a trend.it fails to reject: stationary
#there is disagreement from adf and kpss, may require differencing

forecast::ndiffs(kvi.2.ts, test = "kpss")    # 0 needed
forecast::ndiffs(kvi.2.ts, test = "adf")     # 0 needed. So I am expecting this to be an AR(1) as the full adf test was the only one to say non-stationary

#test stationairity kvi.3
tseries::adf.test(kvi.3.ts)                  #  Rejects just barely: stationary
tseries::adf.test(kvi.3.ts,k = 0)            #forced it to test AR(1) it does reject: Stationary.
tseries::kpss.test(kvi.3.ts, null = "Level") #opposite null hypothesis for kpss it  rejects: non-Stationary
tseries::kpss.test(kvi.3.ts, null = "Trend") #.it rejects:non- stationary
#there is disagreement from adf and kpss, may require differencing

forecast::ndiffs(kvi.3.ts, test = "kpss")    # 1 needed
forecast::ndiffs(kvi.3.ts, test = "adf")     #0 needed. again disagreement. likely will be differenced by auto.arima()

#plot acf and pacf
layout(layout_matrix_1)
acf(kvi.2.ts)
pacf(kvi.2.ts)
#neither shows any clear lag or decay signs. maybe an combined ARMA model

acf(kvi.3.ts)
pacf(kvi.3.ts)
#both have a significant lag at 1. maybe just a differencing will help
acf(diff(kvi.3.ts))
pacf(diff(kvi.3.ts))
#looks like slow decay on pacf and a significant lag at 1 on acf. expecting an ARIMA(0,1,1)

#fit models
fit.kvi.2 <- forecast::auto.arima(kvi.2.ts)
#ARIMA(0,0,1) with non-zero mean. significant pacf lag at 7 could be cyclic nature of date obscuring the MA() pattern
fit.kvi.3 <- forecast::auto.arima(kvi.3.ts)
#ARIMA(0,1,1) as expected from acf and pacf plots
#plot
dev.off()
layout(layout_mat_2)
plot(kvi.2.ts, main= "4 year old sockeye returning to Kvichak River\nthat spent 2 years in ocean", ylab = "log(count of fish x1000)")
lines(fitted(fit.kvi.2), col = "red")
plot(kvi.3.ts, main= "4 year old sockeye returning to Kvichak River\nthat spent 3 years in ocean", ylab = "log(count of fish x1000)")
lines(fitted(fit.kvi.3), col= "red")

#check residuals
forecast::checkresiduals(fit.kvi.2)
#here we see the that the residuals are autocorrelated as Eli had mentioned
#we might find. we fail to reject Ljung-Box telling us the same thing. Histogram looks
#good. Just mention we see the autocorrelation
forecast::checkresiduals(fit.kvi.3)
#this one does much better. fails to reject Ljung-Box and acf plot both telling us no autocorrelation. 
#histogram looks ok. we see the giant negative spike in the tail. lots of variance in the first few years

#train and test predictions
train.kvi.2 <- window(kvi.2.ts, start= 1963, end=2010)
test.kvi.2 <- window(kvi.2.ts, start= 2011, end= 2020)
fit.train.k2 <- forecast::auto.arima(train.kvi.2)
#chooses same model as the full dataset model
pred.kvi2 <- forecast(fit.train.k2, h=10)
plot(pred.kvi2,ylab= "log(counted fish x1000", xlab ="Time")
points(test.kvi.2, pch=19, col ="red")
resd.k2 <- test.kvi.2[1:10] - pred.kvi2$mean[1:10]
RMSE.k2 <- sqrt(mean(resd.k2^2))
accuracy(pred.kvi2, test.kvi.2)
#manual calc matches accuracy(). RMSE = 2.088. points fell outside of the
#prediction intervals with this model. Same straight line provided by the MA() 
#model but there was more variance in this data which corresponds to the poor 
#prediction performance

train.kvi.3 <- window(kvi.3.ts, start= 1963, end=2010)
test.kvi.3 <- window(kvi.3.ts, start= 2011, end= 2020) 
fit.train.k3 <- forecast::auto.arima(train.kvi.3)
#chooses same exact model as the full dataset model
pred.kvi3 <- forecast(fit.train.k3, h=10)
plot(pred.kvi3, ylab= "log(counted fish x1000", xlab ="Time")
points(test.kvi.3, pch=19, col ="red" ) 
resd.k3 <- test.kvi.3[1:10] - pred.kvi3$mean[1:10]
RMSE.k3 <- sqrt(mean(resd.k3^2))
accuracy(pred.kvi3, test.kvi.3)
#RMSE = 0.65 this time we see the prediction line with a positive trend
#this is because the model is ARIMA(0,1,1) with drift. less variance in
#this data set so the points fall much closer to the straight line prediction
#again we see that returns for fish that spent longer in the ocean were more accurately 
#predicted.