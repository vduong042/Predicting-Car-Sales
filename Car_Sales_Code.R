# Load necessary libraries
library(TSA)
library(forecast)
library(lmtest)
library(fUnitRoots)
library(tseries)
library(knitr)
library(dLagM)
library(lattice)
library(bestglm)
library(leaps)
library(ltsa)
library(FitAR)
library(CombMSC)
library(fGarch)
library(zoo)
library(astsa)

# Load data
cars <- read.csv("cars.csv", header = TRUE)
head(cars)

# Create time series object
cars.ts <- ts(as.vector(t(matrix(cars$Sales, nrow = 108, ncol = 1))), start = c(1960, 1), end = c(1968, 12), frequency = 12)
plot(cars.ts, type = 'o', ylab = 'Car Sales')

# ACF and PACF of the time series
acf(cars.ts, lag.max = 36)
pacf(cars.ts, lag.max = 36)

# Fit initial ARIMA model
m1.cars <- arima(cars.ts, order = c(0, 0, 0), seasonal = list(order = c(0, 1, 0), period = 12))
res.m1 <- residuals(m1.cars)
par(mfrow = c(1, 1))
plot(res.m1, xlab = 'Time', ylab = 'Residuals')

# ACF and PACF of residuals for m1.cars
par(mfrow = c(1, 2))
acf(res.m1, lag.max = 36)
pacf(res.m1, lag.max = 36)

# Fit ARIMA(0,0,0)x(1,1,1) model
m2.cars <- arima(cars.ts, order = c(0, 0, 0), seasonal = list(order = c(1, 1, 1), period = 12))
res.m2 <- residuals(m2.cars)
par(mfrow = c(1, 1))
plot(res.m2, xlab = 'Time', ylab = 'Residuals', main = "Time series plot of the residuals")

# ACF and PACF of residuals for m2.cars
par(mfrow = c(1, 2))
acf(res.m2, lag.max = 36)
pacf(res.m2, lag.max = 36)

# Fit ARIMA(0,0,0)x(1,1,2) model
m3.cars <- arima(cars.ts, order = c(0, 0, 0), seasonal = list(order = c(1, 1, 2), period = 12))
res.m3 <- residuals(m3.cars)
par(mfrow = c(1, 2))
acf(res.m3, lag.max = 36)
pacf(res.m3, lag.max = 36)

# Log transformation of the time series
log.cars.ts <- log(cars.ts)
par(mfrow = c(1, 1))
plot(log.cars.ts, ylab = 'log of sales count', xlab = 'Year', type = 'o')

# Fit ARIMA model to log-transformed data
m4.cars <- arima(log.cars.ts, order = c(0, 0, 0), seasonal = list(order = c(1, 1, 2), period = 12))
res.m4 <- residuals(m4.cars)
plot(res.m4, xlab = 'Time', ylab = 'Residuals', main = "Time series plot of the residuals")

# ACF of residuals for m4.cars
acf(res.m4, lag.max = 36)

# PACF of residuals for m4.cars
pacf(res.m4, lag.max = 36)

# Fit SARIMA(0,1,0)x(1,1,2) model
m5.cars <- arima(log.cars.ts, order = c(0, 1, 0), seasonal = list(order = c(1, 1, 2), period = 12))
res.m5 <- residuals(m5.cars)
par(mfrow = c(1, 2))
acf(res.m5, lag.max = 36)
pacf(res.m5, lag.max = 36)

# ADF test for residuals of m5.cars
adf.test(res.m5)

# EACF analysis for residuals of m5.cars
eacf(res.m5)

# Fit various ARIMA models and evaluate
model2.cars <- arima(log.cars.ts, order = c(0, 1, 3), seasonal = list(order = c(1, 1, 2), period = 12), method = "ML")
coeftest(model2.cars)

res.model2 <- residuals(model2.cars)
par(mfrow = c(1, 2))
acf(res.model2, lag.max = 36)
pacf(res.model2, lag.max = 36)

model1.cars <- arima(log.cars.ts, order = c(0, 1, 4), seasonal = list(order = c(1, 1, 2), period = 12), method = "ML")
coeftest(model1.cars)

res.model1 <- residuals(model1.cars)
par(mfrow = c(1, 2))
acf(res.model1, lag.max = 36)
pacf(res.model1, lag.max = 36)

model3.cars <- arima(log.cars.ts, order = c(3, 1, 4), seasonal = list(order = c(1, 1, 2), period = 12), method = "ML")
coeftest(model3.cars)

res.model3 <- residuals(model3.cars)
par(mfrow = c(1, 2))
acf(res.model3, lag.max = 36)
pacf(res.model3, lag.max = 36)

model4.cars <- arima(log.cars.ts, order = c(2, 1, 1), seasonal = list(order = c(1, 1, 2), period = 12), method = "ML")
coeftest(model4.cars)

res.model4 <- residuals(model4.cars)
par(mfrow = c(1, 2))
acf(res.model4, lag.max = 36)
pacf(res.model4, lag.max = 36)

model5.cars <- arima(log.cars.ts, order = c(2, 1, 2), seasonal = list(order = c(1, 1, 2), period = 12), method = "ML")
coeftest(model5.cars)

res.model5 <- residuals(model5.cars)
par(mfrow = c(1, 2))
acf(res.model5, lag.max = 36)
pacf(res.model5, lag.max = 36)

# Define scoring function for AIC and BIC
sort.score <- function(x, score = c("bic", "aic")) {
  if (score == "aic") {
    x[with(x, order(AIC)), ]
  } else if (score == "bic") {
    x[with(x, order(BIC)), ]
  } else {
    warning('score = "x" only accepts valid arguments ("aic", "bic")')
  }
}

# Calculate AIC and BIC for models
sc.AIC <- AIC(model1.cars, model2.cars, model3.cars, model4.cars, model5.cars)
sc.BIC <- AIC(model1.cars, model2.cars, model3.cars, model4.cars, model5.cars, k = log(length(cars.ts)))

# Sort models by AIC and BIC
sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")

# Function to analyze residuals
residual.analysis <- function(model, std = TRUE) {
  library(TSA)
  library(FitAR)
  if (std) {
    res.model <- rstandard(model)
  } else {
    res.model <- residuals(model)
  }
  par(mfrow = c(3, 2))
  plot(res.model, type = 'o', ylab = 'Standardized residuals', main = "Time series plot of standardized residuals")
  abline(h = 0)
  hist(res.model, main = "Histogram of standardized residuals")
  qqnorm(res.model, main = "QQ plot of standardized residuals")
  qqline(res.model, col = 2)
  acf(res.model, main = "ACF of standardized residuals")
  print(shapiro.test(res.model))
  LBQPlot(res.model, lag.max = length(model$residuals) - 1, StartLag = 1, k = 0, SquaredQ = FALSE)
}

# Analyze residuals models
residual.analysis(model = model2.cars)


# Forecast using SARIMA model
sarima.for(log.cars.ts, 120, 0, 1, 3, 1, 1, 2, 12)
