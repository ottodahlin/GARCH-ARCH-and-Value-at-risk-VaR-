#############################################################
## GARCH/ARCH. Value-at-Risk (VaR) - Otto Dahlin.
#############################################################

library(tidyverse)
library(FE)
library(psych)
library(xtable)
library(stargazer)
library(tseries)
library(gridExtra)
library(PerformanceAnalytics)
library(devtools)
library(FinTS)
library(AnalyzeTS)
library(lmtest)
library(tseries)
library(urca)
library(dynlm)
library(sandwich)
library(forecast)


##############################################################################
# daily dow jones data in log returns:
DJ <- DJ_d$r_Dow_Jones
DJ <- daily.dj
DJ

ts.plot(DJ, main="Daily Log returns of Daily Dow Jones index", ylab="Log Returns")
ts.plot(DJ^2, main="DAily Squared Log Returns of DAily Dow Jones", ylab="Squared Log Returns")


# TS object:
ts.DJ <- ts(DJ, start=c(1915, 1), end=c(1990, 2), frequency=365) #daily
ts.DJ.sq <- ts(DJ^2, start=c(1986, 1), end=c(1997, 31), frequency=365) # Squared log returns daily

par(mfrow=c(2,1))
plot(ts.DJ, main="Daily Log Returns DJ", ylab=" Log Returns")
plot(ts.DJ.sq , main="Daily Squared Log Returns DJ", ylab=" Squared Log Returns")

DJ.desc = describe(DJ)
DJ.sq = describe(DJ.sq)
DJ.desc
DJ.sq

DJ.desc.both = round(as.data.frame(rbind(DJ.desc = DJ.desc[2:12], DJ.sq = DJ.sq[2:12])), 6)
xtable(DJ.desc.both) # both DJ and DJ squared.


# QQ-plots:
par(mfrow=c(1,2))
qqnorm(DJ, main="Normal QQ-plot: Log Returns Daily DJ"); qqline(DJ.d, col="red")
# as extra below.... just to see squared..
qqnorm(DJ.sq, main="Normal QQ-plot: Squared Log returns Daily DJ"); qqline(DJ.sq, col="purple")

chart.Histogram(DJ, xlim= c(-0.08,0.08), main="Log Returns Daily DJ", 
                methods=c('add.density', 'add.normal'),
                colorset= c('blue', 'green', 'red'))
# normal curve in Red.
# density curve in green.


# Normality test:
jarque.bera.test(DJ) # reject H0 of normality
jarque.bera.test(DJ^2) # reject h0 of normality

# ACF:
acf(DJ)
acf(DJ^2) # severe autocorrelation in squared returns.This is a sign that we should test for ARCH-effects 

# Testing for ARCH-effects with Ljung-Box test here first. ARCH-LM test below:
Box.test(DJ^2,  type="Ljung-Box") 
# reject H0 of "no autocorrelations".
# Accept H0 of that there IS autoccorelation in squared log returns at 5% sign level.
Box.test(DJ, type="Ljung-Box")

# Testing for ARCH-effects in ARCH-LM test using package "FinTS"


ArchTest(DJ, lags = 5)
ArchTest(DJ^2, lags = 5)
# Reject H0 - ARCH-effects could be found in squared log returns.
# There is presence of ARCH-effects in data at 5% sign level


# BOX JENKIN-Approach: 1) Identification 2)Estimation of Modelspecification 3) Diagnostics
# Evaluate Stationarity
# Selection of the differencing level (d)
# Selection of AR level (p)
# Selection of MA level (q)

# Null hypothesis: Data is non-stationary
adf.test(DJ) # stationary
adf.test(DJ^2) # stationary

par(mfrow=c(2,1))
acf(DJ, main="ACF -  Log Returns Daily Dow Jones") # Q=3
pacf(DJ, main="PACF - Log Returns Daily Dow Jones" ) # P=2, 2 clearly significant where as maybe another 2 (maybe  P=3/4?) trivial

acf(DJ^2)
pacf(DJ^2)

# MODEL ESTIMATION:
auto.arima(DJ, trace=TRUE, ic="bic") # BEst Model: ARMA(2,3) suggests same model as mechanically estimated above.

# BEST MODEL: ARMA(2,3)
model1 <- arima(DJ, order=c(2,0,3))

# we do no differencing. the log returns series is stationary. confirmed as well as confirmed graphically.

#DIAGNOSTICS CHECKING of residuals of chosen model1 specification that is ARMA(2,3) for daily DJ.
model1$residuals %>% ggtsdisplay(plot.type = 'hist', lag.max=25)


################################ Diagnostics checking in Residuals ################################

# DIAGNOSTICS:

# Testing autocorrelation in residuals of model1 (PQ= 2,3). LB-test:
Box.test(resid(model1),  type="Ljung-Box") 
# fail to reject H0 of "No autocorrelation". Accept H0 so we have no autocorrelation in residuals
# in model specification ARMA(2,3)


jarque.bera.test(DJ) # reject H0 of normality

# Squared residuals from estimated ARMA(2,3)
model1.resid <-model1$residuals
acf(model1.resid^2, main="ACF: Estimated residuals of model ARMA(2,3) in Squared Log Returns DJ")
# As we can see, the ACF of squared residuals show many severe autocorrelation.
# There are indeed ARCH effects and that we should model volatility.


# fitting a GARCH to obtain the parameters and fitting parameters.

fit1 <- garchFit(formula= ~ arma(2,3)+ garch(1,1), data=DJ, cond.dist=c("norm"), trace=FALSE)
summary(fit1)

fit2 <- garchFit(formula= ~ arma(2,3)+ garch(1,2), data=DJ,cond.dist=c("norm"), trace=FALSE)
summary(fit2)

fit3 <- garchFit(formula= ~ arma(2,3)+ garch(2,1), data=DJ,cond.dist=c("norm"), trace=FALSE)
summary(fit3)

fit4 <- garchFit(formula= ~ arma(2,3)+ garch(2,2), data=DJ,cond.dist=c("norm"), trace=FALSE)
summary(fit4)


acf(fit1@residuals)


stand.res.fit1 <- fit1@residuals/fit1@sigma.t # standardized residuals ARMA(2,3)-GARCH(1,1)
stand.res.fit2 <- fit2@residuals/fit2@sigma.t # standardized residuals ARMA(2,3)-GARCH(1,2)
stand.res.fit3 <- fit3@residuals/fit3@sigma.t # standardized residuals ARMA(2,3)-GARCH(2,1)
stand.res.fit4 <- fit4@residuals/fit4@sigma.t # standardized residuals ARMA(2,3)-GARCH(2,2)


# Residuals of GARCH(1,1) - "fit1"
acf(stand.res.fit1)
acf(stand.res.fit1^2)

# GARCH (1,2)
acf(stand.res.fit2)
acf(stand.res.fit1^2)

# Garch(2,1)
acf(stand.res.fit3)
acf(stand.res.fit1^2)

# GARCH(2,2)
acf(stand.res.fit4)
acf(stand.res.fit1^4)

# Discrete white noise in all above.

# @fitted - fitted values
# @residuals - residuals
# @h.t - estimated variance
# @sigma.t - estimated standard deviation.

# standardized residuals - residuals divided by their standard deviation should be a white noise.
# Also their squares should be white noise

summary(fit1) # standardized residuals tests: Normality Test, Ljung-Box for standardized residuals and their squares.

# According to IC -BIC, ARMA(2,3)-GARCH(1,2) is the "best" specificaiton i.e. ("fit2") model from above.
# least generated value chosen
##########################################################################################################

# a) plotting original returns, the standardized residuals and conditional standard deviation
par(mfrow=c(3,1))
plot(ts.DJ, type="l", main="Daily DJ Log Returns - Original series", ylab="Log Returns") # original
stand.res.fit2 <- fit2@residuals/fit2@sigma.t #(ARMA 2,3) and GARCH (1,1) - standardizing with conditional standard deviation (sigma.t)
plot(stand.res.fit2, type="l", main="Plot- Standardized Residuals from ARMA(2,3)-GARCH(1,2)",
     ylab="Standardized Residuals") # standardized residuals where "sigma.t" is the conditional standard deviation.
plot(fit2@sigma.t, type="l", main="Conditional standard deviation estimated by ARMA(2,3)-GARCH(1,2)", 
     ylab="Conditional Standard Deviation")
# conditional standard deviation for model "fit2" i.e. ARMA(2,3) and GARCH(1,1)

# Standardized and UNstandardized residuals:
unstand.res.fit2 <- fit2@residuals # UNstandardized residuals

par(mfrow=c(1,2))
plot(stand.res.fit2, type="l", main="Standardized Residuals")
plot(unstand.res.fit2, type="l", main="Unstandardized Residuals")


#Property of standardized residuals is that they should be approximately normally distributed.
# To assess this compare the histogram of the standardized residuals with the standard normal distribution.

hist(stand.res.fit2, xlim=c(-7,7))

chart.Histogram(stand.res.fit2, xlim= c(-5,5), main="Distribution Standardized Residuals", 
                methods=c('add.density', 'add.normal'),
                colorset= c('blue', 'green'))

# discrete white noise below:
acf(stand.res.fit2)
acf(stand.res.fit2^2)
##################################################################################################

# evaluating autocorrelation functions of the standardized residuals and squared standardized residuals

# Assuming we for now use model "fit1" i.e. ARMA(2,3)-GARCH(1,2):
plot(fit2)
10 # acf of standardized residuals for model "fit2"
11 # acf of squared standardzied residuals for model "fit2"
13 # qq-plot of standardized residuals Normal Distribution
# trivial change....

#############################################################################

# Computing portmanteau statistics based on standardized residuals
stand.res.fit2 # standardized residuals from model "fit2" - ARMA(2,3)-GARCH(1,2)

Box.test(stand.res.fit2 , type="Ljung-Box", lag=6)

# lags have been altered....

###########################################################################

# d) Testing for (G)ARCH effects of higher orders by applying ARCH-LM test
# we use the "ArchTest" ARCH-LM test from "FinTS" package in R.

# ARCH-LM a portmanteau test with "q" lags tests whether
# there are ARCH effects at lags from 1 up to "q". It tests the joint
# significance of coefficents, alpha1, alpha2....alphaq in the equation.
# same results when testing for higher orders of ARCH LM test.

ArchTest(DJ, lags = 1)
ArchTest(DJ^2, lags = 5) # have tested for higher orders.. same results.
# we reject H0 despite chosen lag length.
#############################################################################

##################################################################################

# 
stand.res.fit2
qqnorm(stand.res.fit2,
       ylab="standardized residuals",
       main="QQ plot Standardized Residuals Model Fit1:  ARMA(2,3)-GARCH(1,1)")
qqline(stand.res.fit2, col="red")

jarque.bera.test(stand.res.fit2)
#####################################################################################
# T-DISTRIBUTION:

# 4) Re-estiamting model with T-distribution for all GARCH specifications. Using "std" as cond.dist
# form "fGarch" package. Estimating same models with T-distribution

library(fGarch)
fit1.T <- garchFit(formula= ~ arma(2,3)+ garch(1,1), data=DJ, cond.dist=c("std"), trace=FALSE)
summary(fit1.T)

fit2.T <- garchFit(formula= ~ arma(2,3)+ garch(1,2), data=DJ,cond.dist=c("std"), trace=FALSE)
summary(fit2.T)

fit3.T <- garchFit(formula= ~ arma(2,3)+ garch(2,1), data=DJ,cond.dist=c("std"), trace=FALSE)
summary(fit3.T)

fit4.T <- garchFit(formula= ~ arma(2,3)+ garch(2,2), data=DJ,cond.dist=c("std"), trace=FALSE)
summary(fit4.T)

# Judged from BIC then ARMA(2,3)-GARCH(1,2) is the most appropriate once again.

# STANDARDIZED RESIDUALS:
stand.res.fit1.T <- fit1.T@residuals/fit1.T@sigma.t # standardized residuals ARMA(2,3)-GARCH(1,1)
stand.res.fit2.T <- fit2.T@residuals/fit2.T@sigma.t # standardized residuals ARMA(2,3)-GARCH(1,2) # most preferred.
stand.res.fit3.T <- fit3.T@residuals/fit3.T@sigma.t # standardized residuals ARMA(2,3)-GARCH(2,1)
stand.res.fit4.T <- fit4.T@residuals/fit4.T@sigma.t # standardized residuals ARMA(2,3)-GARCH(2,2)

# Residuals of GARCH(1,2)
acf(stand.res.fit1.T)
acf(stand.res.fit1.T^2)

# GARCH (1,2)
acf(stand.res.fit2.T)
acf(stand.res.fit1.T^2)

# Garch(2,1)
acf(stand.res.fit3.T)
acf(stand.res.fit1.T^2)

# GARCH(2,2)
acf(stand.res.fit4.T)
acf(stand.res.fit1.T^4)

# acf and pacf of squared standardized residuals and standardized residuals of preferred model spec.
plot(fit2.T)
par(mfrow=c(1,2), main="hello")
10
11
13

# Standardized residuals from T-distribution: QQ-plot
plot(fit2.T) # choose nr 13 = "QQ-plot of standardized residuals"
13 # very big difference in tails between normal and t-distribution

# normal distribution QQ-plot
plot(fit2)
13

# Yes T-distribution is severly more appropriate.


###########################################################################################################

# Redoing above analysis for monthly DJ returns (Higher aggregation levels) - Converting to Monthly data

DJ <- daily.dj

# TS object: Converting from daily to Monthly series...
ts.DJ.m <- ts(DJ, start=c(1915, 1), end=c(1990, 2), frequency=12) #monthly
ts.plot(ts.DJ.m) # daily
class(ts.DJ.m)

DJ.m.df <- as.data.frame(ts.DJ.m)
DJ.m.df
DJ.m <- as.numeric(DJ.m.df$x)
ts.plot(DJ.m, main="Aggregated - Log Returns Monthly DJ series (from Daily To Monthly)", ylab="Log returns") # monthly Daily Dow jones aggregated from daily to monthly data.

# Normality test:
jarque.bera.test(DJ.m) # reject H0 of normality
jarque.bera.test(DJ.m^2) # reject h0 of normality

# ACF:
acf(DJ.m)
acf(DJ.m^2) # severe autocorrelation in squared returns.This is a sign that we should test for ARCH-effects 

# Testing for ARCH-effects with Ljung-Box (Portmanteau statistic) test here first. ARCH-LM test below:
Box.test(DJ.m^2,  type="Ljung-Box")  # squared monthly returns.
# reject H0 of "no autocorrelations".
# Accept H0 of that there IS autoccorelation in squared log returns at 5% sign level.

Box.test(DJ.m, type="Ljung-Box") # no autocorrelation in just log monthly DJ returns

# Testing for arch effects in Monthly DJ series
ArchTest(DJ.m, lags = 5)
ArchTest(DJ.m^2, lags = 1)
# Reject H0 - ARCH-effects could be found in squared log returns of monthly DJ returns.
# There is presence of ARCH-effects in data at 5% sign level in both just monthly DJ and squared monthly DJ returns
# despite lag length selection.


# BOX JENKIN-Approach: 1) Identification 2)Estimation of Modelspecification 3) Diagnostics
# Evaluate Stationarity
# Selection of the differencing level (d)
# Selection of AR level (p)
# Selection of MA level (q)

# Null hypothesis: Data is non-stationary
adf.test(DJ.m) # stationary
adf.test(DJ.m^2) # stationary

par(mfrow=c(2,1))
Acf(DJ.m, main="ACF -  Log Returns Monthly Dow Jones") 
pacf(DJ.m, main="PACF - Log Returns Monthly Dow Jones" ) 

acf(DJ.m^2)
pacf(DJ.m^2)

# MODEL ESTIMATION: disagreee with auto.arima function.... proposes other model.
auto.arima(DJ.m, trace=TRUE, ic="bic") 

# BEST MODEL: ARMA(3,2)
model1.m <- arima(DJ.m, order=c(4,0,1))

# we do no differencing. the log returns series is stationary. confirmed as well as confirmed graphically.

################################ Diagnostics checking in Residuals ################################

# DIAGNOSTICS:

# Testing autocorrelation in residuals of model1 . LB-test:
Box.test(resid(model1.m),  type="Ljung-Box") 
# fail to reject H0 of "No autocorrelation". Accept H0 so we have NO autocorrelation in residuals
# in model specification ARMA(4,1)

par(mfrow=c(2,2))
plot(residuals(model1.m), main ="Residualsplot | ARMA (4,1")
hist(residuals(model1.m),xlim= c(-0.10,0.10),  main ="Residual Distribution ARMA (4,1)", col = "gray", xlab = "Residual")
acf(residuals(model1.m), main = "ACF residuals | ARMA( 4,1)")
pacf(residuals(model1.m), main = "PACF residuals | ARMA (4,1)")

# acf and pacf indicates of relatively discrete white noise definitely.

jarque.bera.test(DJ.m) # reject H0 of normality

# Squared residuals from estimated ARMA(2,3)
model1.m.resid <-model1$residuals
acf(model1.m.resid^2, main="ACF: Estimated residuals of model ARMA(4,1) in Squared Log Returns  Monthly DJ")
# As we can see, the ACF of squared residuals show many severe autocorrelation.
# There are indeed ARCH effects and that we should fit a volatility model.


# fitting a GARCH to obtain the parameters and fitting parameters.

fit1.m <- garchFit(formula= ~ arma(4,1)+ garch(1,1), data=DJ.m, cond.dist=c("norm"), trace=FALSE)
summary(fit1.m) # BIC: -5.956915  # BEST PROPOSED MODEL

fit2.m <- garchFit(formula= ~ arma(4,1)+ garch(1,2), data=DJ.m,cond.dist=c("norm"), trace=FALSE)
summary(fit2.m) # BIC: -5.948244

fit3.m <- garchFit(formula= ~ arma(4,1)+ garch(2,1), data=DJ.m,cond.dist=c("norm"), trace=FALSE)
summary(fit3.m) # BIC: -5.942985

fit4.m <- garchFit(formula= ~ arma(4,1)+ garch(2,2), data=DJ.m,cond.dist=c("norm"), trace=FALSE)
summary(fit4.m) # -5.940700 

acf(fit1.m@residuals) # white noise..

stand.res.fit1.m <- fit1.m@residuals/fit1.m@sigma.t # standardized residuals ARMA(2,3)-GARCH(1,1)
stand.res.fit2.m <- fit2.m@residuals/fit2.m@sigma.t # standardized residuals ARMA(2,3)-GARCH(1,2)
stand.res.fit3.m <- fit3.m@residuals/fit3.m@sigma.t # standardized residuals ARMA(2,3)-GARCH(2,1)
stand.res.fit4.m <- fit4.m@residuals/fit4.m@sigma.t # standardized residuals ARMA(2,3)-GARCH(2,2)


# Residuals of GARCH(1,1) - "fit1 Montly DJ series"
acf(stand.res.fit1.m)
acf(stand.res.fit1.m^2)

# GARCH (1,2)
acf(stand.res.fit2.m)
acf(stand.res.fit2.m^2)

# Garch(2,1)
acf(stand.res.fit3.m)
acf(stand.res.fit3.m^2)

# GARCH(2,2)
acf(stand.res.fit4.m)
acf(stand.res.fit4.m^2)

# Discrete white noise in all above.

# @fitted - fitted values
# @residuals - residuals
# @h.t - estimated variance
# @sigma.t - estimated standard deviation.

# standardized residuals - residuals divided by their standard deviation should be a white noise.
# Also their squares should be white noise

summary(fit1.m) # standardized residuals tests: Normality Test, Ljung-Box for standardized residuals and their squares.

# According to IC -BIC, ARMA(4,1)-GARCH(1,1) is the "best" specificaiton 
# least generated value chosen
##########################################################################################################

# 2) ARMA(4,1)-GARCH(1,1)

# a) plotting original returns, the standardized residuals and conditional standard deviation
par(mfrow=c(3,1))
plot(ts.DJ.m, type="l", main="Daily DJ Log Returns - Original series", ylab="Log Returns") # original
stand.res.fit1.m <- fit1.m@residuals/fit1.m@sigma.t #(ARMA 4,1) and GARCH (1,1) - standardizing with conditional standard deviation (sigma.t)
plot(stand.res.fit1.m, type="l", main="Plot- Standardized Residuals from ARMA(4,1)-GARCH(1,1)",
     ylab="Standardized Residuals") # standardized residuals where "sigma.t" is the conditional standard deviation.
plot(fit1.m@sigma.t, type="l", main="Conditional standard deviation estimated by ARMA(4,1)-GARCH(1,1)", 
     ylab="Conditional Standard Deviation")


# Standardized and UNstandardized residuals:
unstand.res.fit1.m <- fit1.m@residuals # UNstandardized residuals
plot(unstand.res.fit1.m , type="l", main="unstandardized residuals")


par(mfrow=c(1,2))
plot(stand.res.fit1.m, type="l", main="Standardized Residuals")
plot(unstand.res.fit1.m, type="l", main="Unstandardized Residuals")


#Property of standardized residuals is that they should be approximately normally distributed.
# To assess this compare the histogram of the standardized residuals with the standard normal distribution.

hist(stand.res.fit1.m, xlim=c(-7,7), main="Histogram Standardized REsiduals")

chart.Histogram(stand.res.fit1.m, xlim= c(-5,5), main="Distribution Standardized Residuals", 
                methods=c('add.density', 'add.normal'),
                colorset= c('blue', 'green'))

# discrete white noise below:
acf(stand.res.fit1.m, main="Standardized Residuals DJ monthly series")
acf(stand.res.fit1.m^2, main="Squared Standardized Residuals DJ monthly Series")
##################################################################################################

# b) evaluating autocorrelation functions of the standardized residuals and squared standardized residuals

# Assuming we for now use model "fit1" i.e. ARMA(4,1)-GARCH(1,1):
plot(fit1.m)
10 # acf of standardized residuals for model "fit2"
11 # acf of squared standardzied residuals for model "fit2"
13 # qq-plot of standardized residuals Normal Distribution


#############################################################################

# c) Computing portmanteau statistics based on standardized residuals
stand.res.fit1.m # standardized residuals from model from ARMA(4,1)-GARCH(1,1)

# tried different lag lengths...
Box.test(stand.res.fit1.m , type="Ljung-Box", lag=3)

###########################################################################

# d) Testing for (G)ARCH effects of higher orders by applying ARCH-LM test
# we use the "ArchTest" ARCH-LM test from "FinTS" package in R.

# ARCH-LM a portmanteau test with "q" lags tests whether
# there are ARCH effects at lags from 1 up to "q". It tests the joint
# significance of coefficents, alpha1, alpha2....alphaq in the equation.
# same results when testing for higher orders of ARCH LM test.


#############################################################################

# e)  using information criteria - see summary()

# For information criterion model judgements see PDF /LATEX report.


# T-DISTRIBUTION: ARMA(4,1)-GARCH(1,1)

# 4) Re-estiamting model with T-distribution for all GARCH specifications. Using "std" as cond.dist
# form "fGarch" package. Estimating same models with T-distribution

fit1.T.m <- garchFit(formula= ~ arma(4,1)+ garch(1,1), data=DJ, cond.dist=c("std"), trace=FALSE)
summary(fit1.T.m) # BIC: 

fit2.T.m <- garchFit(formula= ~ arma(4,1)+ garch(1,2), data=DJ,cond.dist=c("std"), trace=FALSE)
summary(fit2.T.m) # BIC: 

fit3.T.m <- garchFit(formula= ~ arma(4,1)+ garch(2,1), data=DJ,cond.dist=c("std"), trace=FALSE)
summary(fit3.T.m) # BIC: 

fit4.T.m <- garchFit(formula= ~ arma(4,1)+ garch(2,2), data=DJ,cond.dist=c("std"), trace=FALSE)
summary(fit4.T.m) # BIC: 

# Judged from BIC then ARMA(4,1)-GARCH(1,2) is the most appropriate once again.

# STANDARDIZED RESIDUALS:
stand.res.fit1.T.m <- fit1.T.m@residuals/fit1.T.m@sigma.t # standardized residuals ARMA(4,1)-GARCH(1,1)
stand.res.fit2.T.m <- fit2.T.m@residuals/fit2.T.m@sigma.t # standardized residuals ARMA(4,1)-GARCH(1,2) # most preferred.
stand.res.fit3.T.m <- fit3.T.m@residuals/fit3.T.m@sigma.t # standardized residuals ARMA(4,1)-GARCH(2,1)
stand.res.fit4.T.m <- fit4.T.m@residuals/fit4.T.m@sigma.t # standardized residuals ARMA(4,1)-GARCH(2,2)

# Residuals of GARCH(1,2)
acf(stand.res.fit1.T.m)
acf(stand.res.fit1.T.m^2)

# GARCH (1,2)
acf(stand.res.fit2.T.m)
acf(stand.res.fit2.T.m^2)

# Garch(2,1)
acf(stand.res.fit3.T.m)
acf(stand.res.fit3.T.m^2)

# GARCH(2,2)
acf(stand.res.fit4.T.m)
acf(stand.res.fit4.T.m^4)

# acf and pacf of squared standardized residuals and standardized residuals of preferred model spec.
plot(fit2.T.m)
par(mfrow=c(1,2))
10
11
13

# Standardized residuals from T-distribution: QQ-plot
plot(fit2.T.m) # choose nr 13 = "QQ-plot of standardized residuals"
13 # very big difference in tails between normal and t-distribution

# normal distribution QQ-plot
plot(fit2.m)
13

# Yes T-distribution is severly more appropriate. Big difference in tails...
# Normal dist has more extreme outliers...

#################################################################################################
#################################################################################################
#################################################################################################



# PART B) Value-at-Risk


## daily ##
daily <- DJ_d$r_Dow_Jones
daily <- DJ
DJ
hist(daily)

mu <- mean(daily)  
sigma <- sqrt(var(daily))


horizon <- c(1,5,10) # forecast horizons: 1,5,10 days.
##################################################################

# 1) - Unconditional moments of a normal distribution

## Normal Daily, 0.01
x <- qnorm(0.01)
x
# h=1
VaR_0.01_h1 <- exp(qnorm(0.01) * sigma + mu)-1
VaR_0.01_h1

# h=5
VaR_0.01_h5 <- exp(qnorm(0.01) * sqrt((sigma^2*5))+ 5*mu) -1
VaR_0.01_h5

#h=10
VaR_0.01_h10 <- exp(qnorm(0.01) * sqrt((sigma^2*10))+ 10*mu) -1
VaR_0.01_h10

Normal_0.01 <- rbind(VaR_0.01_h1, VaR_0.01_h5, VaR_0.01_h10 )

#################################################################

# Normal Daily,  0.05, h=1,5,10 days forward
z <- qnorm(0.05)
z
# h=1
VaR_0.05_h1 <- exp(qnorm(0.05) * sigma + mu)-1
VaR_0.05_h1

# h=5
VaR_0.05_h5 <- exp(qnorm(0.05) * sqrt((sigma^2*5))+ 5*mu) -1
VaR_0.05_h5

# h=10
VaR_0.05_h10 <- exp(qnorm(0.05) * sqrt((sigma^2*10))+ 10*mu) -1
VaR_0.05_h10

Normal_0.05 <- rbind(VaR_0.05_h1, VaR_0.05_h5,VaR_0.05_h10 )

####################################################################

# Normal Daily,  0.005, h=1,5,10 days forward
p <- qnorm(0.005)
p
# h=1
VaR_0.005_h1 <- exp(qnorm(0.005) * sigma + mu)-1
VaR_0.005_h1

# h=5
VaR_0.005_h5 <- exp(qnorm(0.005) * sqrt((sigma^2*5))+ 5*mu) -1
VaR_0.005_h5

# h=10
VaR_0.005_h10 <- exp(qnorm(0.005) * sqrt((sigma^2*10))+ 10*mu) -1
VaR_0.005_h10

Normal_0.005 <- rbind(VaR_0.005_h1, VaR_0.005_h5,  VaR_0.005_h10 )
####################################################################
# Normal Daily,  0.1 with  h=1,5,10 days forward
q <- qnorm(0.1)
q
# h=1
VaR_0.1_h1 <- exp(qnorm(0.1) * sigma + mu)-1
VaR_0.1_h1

# h=5
VaR_0.1_h5 <- exp(qnorm(0.1) * sqrt((sigma^2*5))+ 5*mu) -1
VaR_0.1_h5

# h=10
VaR_0.1_h10 <- exp(qnorm(0.1) * sqrt((sigma^2*10))+ 10*mu) -1
VaR_0.1_h10

Normal_0.1 <- rbind(VaR_0.1_h1, VaR_0.1_h5,  VaR_0.1_h10 )
####################################################################

# h=1
VaR_0.001_h1 <- exp(qnorm(0.001) * sigma + mu)-1
VaR_0.001_h1

# h=5
VaR_0.001_h5 <- exp(qnorm(0.001) * sqrt((sigma^2*5))+ 5*mu) -1
VaR_0.001_h5

# h=10
VaR_0.001_h10 <- exp(qnorm(0.001) * sqrt((sigma^2*10))+ 10*mu) -1
VaR_0.001_h10

Normal_0.001 <- rbind(VaR_0.001_h1, VaR_0.001_h5,  VaR_0.001_h10 )

#####################################################################

# Weekly DJ
weekly <- DJ_w$r_close
weekly
hist(weekly)
mu.w <- mean(weekly)
sigma.w <- sqrt(var(weekly))

# Normal Weekly, 0.01 with h=1,5,10 days

VaR_0.01_h1.weekly <- exp(qnorm(0.01) * sigma.w + mu.w)-1
VaR_0.01_h1.weekly

# h=5
VaR_0.01_h5.weekly <- exp(qnorm(0.01) * sqrt((sigma.w^2*5))+ 5*mu.w) -1
VaR_0.01_h5.weekly

#h=10
VaR_0.01_h10.weekly <- exp(qnorm(0.01) * sqrt((sigma.w^2*10))+ 10*mu.w) -1
VaR_0.01_h10.weekly

Normal_0.01.weekly <- rbind(VaR_0.01_h1.weekly, VaR_0.01_h5.weekly, VaR_0.01_h10.weekly)

#######################################################################################

# h=1, 0.05 - normal weekly
VaR_0.05_h1.weekly <- exp(qnorm(0.05) * sigma.w + mu.w)-1
VaR_0.05_h1.weekly

# h=5
VaR_0.05_h5.weekly <- exp(qnorm(0.05) * sqrt((sigma.w^2*5))+ 5*mu.w) -1
VaR_0.05_h5.weekly

# h=10
VaR_0.05_h10.weekly <- exp(qnorm(0.05) * sqrt((sigma.w^2*10))+ 10*mu.w) -1
VaR_0.05_h10.weekly

Normal_0.05.weekly <- rbind(VaR_0.05_h1.weekly, VaR_0.05_h5.weekly,VaR_0.05_h10.weekly )
##############################################################################################


# Normal weekly,  0.005, h=1,5,10 days forward

# h=1
VaR_0.005_h1.weekly <- exp(qnorm(0.005) * sigma.w + mu.w)-1
VaR_0.005_h1.weekly

# h=5
VaR_0.005_h5.weekly <- exp(qnorm(0.005) * sqrt((sigma.w^2*5))+ 5*mu.w) -1
VaR_0.005_h5.weekly

# h=10
VaR_0.005_h10.weekly <- exp(qnorm(0.005) * sqrt((sigma.w^2*10))+ 10*mu.w) -1
VaR_0.005_h10.weekly

Normal_0.005.weekly <- rbind(VaR_0.005_h1.weekly, VaR_0.005_h5.weekly,  VaR_0.005_h10.weekly )

####################################################################
# Normal Weekly,  0.1 with  h=1,5,10 days forward

# h=1
VaR_0.1_h1.weekly <- exp(qnorm(0.1) * sigma.w + mu.w)-1
VaR_0.1_h1.weekly

# h=5
VaR_0.1_h5.weekly <- exp(qnorm(0.1) * sqrt((sigma.w^2*5))+ 5*mu.w) -1
VaR_0.1_h5.weekly

# h=10
VaR_0.1_h10.weekly <- exp(qnorm(0.1) * sqrt((sigma.w^2*10))+ 10*mu.w) -1
VaR_0.1_h10.weekly

Normal_0.1.weekly <- rbind(VaR_0.1_h1.weekly, VaR_0.1_h5.weekly,  VaR_0.1_h10.weekly )
####################################################################
# Normal Weekly,  0.1 with  h=1,5,10 days forward

# h=1
VaR_0.001_h1.weekly <- exp(qnorm(0.001) * sigma.w + mu.w)-1
VaR_0.001_h1.weekly

# h=5
VaR_0.001_h5.weekly <- exp(qnorm(0.001) * sqrt((sigma.w^2*5))+ 5*mu.w) -1
VaR_0.001_h5.weekly

# h=10
VaR_0.001_h10.weekly <- exp(qnorm(0.001) * sqrt((sigma.w^2*10))+ 10*mu.w) -1
VaR_0.001_h10.weekly

Normal_0.001.weekly <- rbind(VaR_0.001_h1.weekly, VaR_0.001_h5.weekly,  VaR_0.001_h10.weekly )


####################################################################
####################################################################

# 2) Unconditional moments of a T-distribution, Daily DJ

# N-1 degrees of freedom, DAily data
class(daily)
df.daily <- data.frame(daily)
daily.T <- df.daily[-18839,] # 18838 observations that is n-1 degrees of freedom
length(daily.T)

# obs, n=18838
mu.T <- mean(daily.T)
sigma.T <- sqrt(var(daily.T))


df=18838 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.01_h1.T <- exp(qt(0.01, df) * sigma.T + mu.T)-1
VaR_0.01_h1.T

?qt
# h=5
VaR_0.01_h5.T <- exp(qt(0.01, df)* sqrt((sigma.T^2*5)) + 5*mu.T) - 1
VaR_0.01_h5.T

# h=10
VaR_0.01_h10.T <- exp(qt(0.01, df)* sqrt((sigma.T^2*10)) + 10*mu.T) - 1
VaR_0.01_h10.T

Tdist_0.01.daily <- rbind(VaR_0.01_h1.T, VaR_0.01_h5.T,VaR_0.01_h10.T )
############################################################################

# 0.1
df=18838 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.1_h1.T <- exp(qt(0.1, df) * sigma.T + mu.T)-1
VaR_0.1_h1.T

# h=5
VaR_0.1_h5.T <- exp(qt(0.1, df)* sqrt((sigma.T^2*5)) + 5*mu.T) - 1
VaR_0.1_h5.T

# h=10
VaR_0.1_h10.T <- exp(qt(0.1, df)* sqrt((sigma.T^2*10)) + 10*mu.T) - 1
VaR_0.1_h10.T

Tdist_0.1.daily <- rbind(VaR_0.1_h1.T, VaR_0.1_h5.T, VaR_0.1_h10.T  )
############################################################################

# 0.05
df=18838 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.05_h1.T <- exp(qt(0.05, df) * sigma.T + mu.T)-1
VaR_0.05_h1.T

# h=5
VaR_0.05_h5.T <- exp(qt(0.05, df)* sqrt((sigma.T^2*5)) + 5*mu.T) - 1
VaR_0.05_h5.T

# h=10
VaR_0.05_h10.T <- exp(qt(0.05, df)* sqrt((sigma.T^2*10)) + 10*mu.T) - 1
VaR_0.05_h10.T

Tdist_0.05.daily <- rbind(VaR_0.05_h1.T, VaR_0.05_h5.T, VaR_0.05_h10.T )
############################################################################

# 0.005
df=18838 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.005_h1.T <- exp(qt(0.005, df) * sigma.T + mu.T)-1
VaR_0.005_h1.T

# h=5
VaR_0.005_h5.T <- exp(qt(0.005, df)* sqrt((sigma.T^2*5)) + 5*mu.T) - 1
VaR_0.005_h5.T

# h=10
VaR_0.005_h10.T <- exp(qt(0.005, df)* sqrt((sigma.T^2*10)) + 10*mu.T) - 1
VaR_0.005_h10.T

Tdist_0.005.daily <- rbind(VaR_0.005_h1.T, VaR_0.005_h5.T, VaR_0.005_h10.T)

############################################################################

# 0.001
df=18838 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.001_h1.T <- exp(qt(0.001, df) * sigma.T + mu.T)-1
VaR_0.001_h1.T

# h=5
VaR_0.001_h5.T <- exp(qt(0.001, df)* sqrt((sigma.T^2*5)) + 5*mu.T) - 1
VaR_0.001_h5.T

# h=10
VaR_0.001_h10.T <- exp(qt(0.001, df)* sqrt((sigma.T^2*10)) + 10*mu.T) - 1
VaR_0.001_h10.T

Tdist_0.001.daily <- rbind(VaR_0.001_h1.T, VaR_0.001_h5.T, VaR_0.001_h10.T)


############################################################################


# T distribution Weekly DJ
# N-1 degrees of freedom, Weekly Data
str(weekly) # 4686 obs.

df.weekly <- data.frame(weekly)
str(df.weekly)
weekly.T <- df.weekly[-4686,] # 4685 observations that is n-1 degrees of freedom
str(weekly.T)

# obs, n=18838
mu.T.weekly <- mean(weekly.T)
sigma.T.weekly <- sqrt(var(weekly.T))

# 0.1
df=4685 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.1_h1.T.weekly <- exp(qt(0.1, df) * sigma.T.weekly + mu.T.weekly)-1
VaR_0.1_h1.T.weekly

# h=5
VaR_0.1_h5.T.weekly <- exp(qt(0.1, df)* sqrt((sigma.T.weekly^2*5)) + 5*mu.T.weekly) - 1
VaR_0.1_h5.T.weekly

# h=10
VaR_0.1_h10.T.weekly <- exp(qt(0.1, df)* sqrt((sigma.T.weekly^2*10)) + 10*mu.T.weekly) - 1
VaR_0.1_h10.T.weekly

Tdist_0.1.weekly <- rbind(VaR_0.1_h1.T.weekly, VaR_0.1_h5.T.weekly,VaR_0.1_h10.T.weekly )
############################################################################


# 0.05
df=4685 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.05_h1.T.weekly <- exp(qt(0.05, df) * sigma.T.weekly + mu.T.weekly)-1
VaR_0.05_h1.T.weekly

# h=5
VaR_0.05_h5.T.weekly <- exp(qt(0.05, df)* sqrt((sigma.T.weekly^2*5)) + 5*mu.T.weekly) - 1
VaR_0.05_h5.T.weekly

# h=10
VaR_0.05_h10.T.weekly <- exp(qt(0.05, df)* sqrt((sigma.T.weekly^2*10)) + 10*mu.T.weekly) - 1
VaR_0.05_h10.T.weekly

Tdist_0.05.weekly <- rbind(VaR_0.05_h1.T.weekly,VaR_0.05_h5.T.weekly, VaR_0.05_h10.T.weekly )
############################################################################

# 0.01
df=4685 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.01_h1.T.weekly <- exp(qt(0.01, df) * sigma.T.weekly + mu.T.weekly)-1
VaR_0.01_h1.T.weekly

# h=5
VaR_0.01_h5.T.weekly <- exp(qt(0.01, df)* sqrt((sigma.T.weekly^2*5)) + 5*mu.T.weekly) - 1
VaR_0.01_h5.T.weekly

# h=10
VaR_0.01_h10.T.weekly <- exp(qt(0.01, df)* sqrt((sigma.T.weekly^2*10)) + 10*mu.T.weekly) - 1
VaR_0.01_h10.T.weekly

Tdist_0.01.weekly <- rbind(VaR_0.01_h1.T.weekly, VaR_0.01_h5.T.weekly,VaR_0.01_h10.T.weekly)
#######################################################################################

# 0.005
df=4685 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.005_h1.T.weekly <- exp(qt(0.005, df) * sigma.T.weekly + mu.T.weekly)-1
VaR_0.005_h1.T.weekly

# h=5
VaR_0.005_h5.T.weekly <- exp(qt(0.005, df)* sqrt((sigma.T.weekly^2*5)) + 5*mu.T.weekly) - 1
VaR_0.005_h5.T.weekly

# h=10
VaR_0.005_h10.T.weekly <- exp(qt(0.005, df)* sqrt((sigma.T.weekly^2*10)) + 10*mu.T.weekly) - 1
VaR_0.005_h10.T.weekly

Tdist_0.005.weekly <- rbind(VaR_0.005_h1.T.weekly,VaR_0.005_h5.T.weekly, VaR_0.005_h10.T.weekly )

#######################################################################################

# 0.001
df=4685 # n-1 degrees of freedom
length(daily.T)
# h=1
VaR_0.001_h1.T.weekly <- exp(qt(0.001, df) * sigma.T.weekly + mu.T.weekly)-1
VaR_0.001_h1.T.weekly

# h=5
VaR_0.001_h5.T.weekly <- exp(qt(0.001, df)* sqrt((sigma.T.weekly^2*5)) + 5*mu.T.weekly) - 1
VaR_0.001_h5.T.weekly

# h=10
VaR_0.001_h10.T.weekly <- exp(qt(0.001, df)* sqrt((sigma.T.weekly^2*10)) + 10*mu.T.weekly) - 1
VaR_0.001_h10.T.weekly

Tdist_0.001.weekly <- rbind(VaR_0.001_h1.T.weekly, VaR_0.001_h1.T.weekly, VaR_0.001_h10.T.weekly )


####################################################################################################
####################################################################################################

# 3) A gaussian GARCH(1,1) with constant conditional mean function
library(fGarch)

mean1 <- garchFit(formula= ~ garch(1,1), data=daily, cond.dist=c("norm"), trace=FALSE)
summary(mean1)
mean1

mean1
?predict
mu.garch11.norm <- predict(mean1, 1)[1]
mu.garch11.norm

sigma.garch11.norm <- predict(mean1, 1)[3] 
sigma.garch11.norm 

## Normal Daily, 0.1
# h=1
VaR_0.1_h1.garch11 <- exp(qnorm(0.1)* sigma.garch11.norm  + mu.garch11.norm) - 1
VaR_0.1_h1.garch11

#h=5
VaR_0.1_h5.garch11 <- exp(qnorm(0.1)* sqrt((sigma.garch11.norm^2*5)) + 5*mu.garch11.norm) - 1
VaR_0.1_h5.garch11

#h=10
VaR_0.1_h10.garch11 <- exp(qnorm(0.1)* sqrt((sigma.garch11.norm^2*10)) + 10*mu.garch11.norm) - 1
VaR_0.1_h10.garch11

normal_0.1.garch11 <- rbind(VaR_0.1_h1.garch11, VaR_0.1_h5.garch11,VaR_0.1_h10.garch11  )
######################################################################################################


## Normal Daily, 0.05
# h=1
VaR_0.05_h1.garch11 <- exp(qnorm(0.05)* sigma.garch11.norm  + mu.garch11.norm) - 1
VaR_0.05_h1.garch11

#h=5
VaR_0.05_h5.garch11 <- exp(qnorm(0.05)* sqrt((sigma.garch11.norm^2*5)) + 5*mu.garch11.norm) - 1
VaR_0.05_h5.garch11

#h=10
VaR_0.05_h10.garch11 <- exp(qnorm(0.05)* sqrt((sigma.garch11.norm^2*10)) + 10*mu.garch11.norm) - 1
VaR_0.05_h10.garch11

normal_0.05.garch11 <- rbind(VaR_0.05_h1.garch11,  VaR_0.05_h5.garch11,VaR_0.05_h10.garch11 )
######################################################################################################

## Normal Daily, 0.01
# h=1
VaR_0.01_h1.garch11 <- exp(qnorm(0.01)* sigma.garch11.norm  + mu.garch11.norm) - 1
VaR_0.01_h1.garch11

#h=5
VaR_0.01_h5.garch11 <- exp(qnorm(0.01)* sqrt((sigma.garch11.norm^2*5)) + 5*mu.garch11.norm) - 1
VaR_0.01_h5.garch11

#h=10
VaR_0.01_h10.garch11 <- exp(qnorm(0.01)* sqrt((sigma.garch11.norm^2*10)) + 10*mu.garch11.norm) - 1
VaR_0.01_h10.garch11

normal_0.01.garch11 <- rbind(VaR_0.01_h1.garch11,  VaR_0.01_h5.garch11, VaR_0.01_h10.garch11)
######################################################################################################

## Normal Daily, 0.005
# h=1
VaR_0.005_h1.garch11 <- exp(qnorm(0.005)* sigma.garch11.norm  + mu.garch11.norm) - 1
VaR_0.005_h1.garch11

#h=5
VaR_0.005_h5.garch11 <- exp(qnorm(0.005)* sqrt((sigma.garch11.norm^2*5)) + 5*mu.garch11.norm) - 1
VaR_0.005_h5.garch11

#h=10
VaR_0.005_h10.garch11 <- exp(qnorm(0.05)* sqrt((sigma.garch11.norm^2*10)) + 10*mu.garch11.norm) - 1
VaR_0.05_h10.garch11

normal_0.005.garch11 <- rbind(VaR_0.005_h1.garch11, VaR_0.005_h5.garch11,VaR_0.05_h10.garch11  )


######################################################################################################

## Normal Daily, 0.001
# h=1
VaR_0.001_h1.garch11 <- exp(qnorm(0.001)* sigma.garch11.norm  + mu.garch11.norm) - 1
VaR_0.001_h1.garch11

#h=5
VaR_0.001_h5.garch11 <- exp(qnorm(0.001)* sqrt((sigma.garch11.norm^2*5)) + 5*mu.garch11.norm) - 1
VaR_0.001_h5.garch11

#h=10
VaR_0.001_h10.garch11 <- exp(qnorm(0.001)* sqrt((sigma.garch11.norm^2*10)) + 10*mu.garch11.norm) - 1
VaR_0.001_h10.garch11

normal_0.001.garch11 <- rbind(VaR_0.001_h1.garch11, VaR_0.001_h5.garch11,VaR_0.001_h10.garch11  )

######################################################################################################


###### SAME AS ABOVE THOUGH FOR WEEKLY DATA: NORMAL, Weekly data ############################

mean1.w <- garchFit(formula= ~ garch(1,1), data=weekly, cond.dist=c("norm"), trace=FALSE)
summary(mean1.w)
mean1.w

mu.garch11.norm.w <- predict(mean1.w, 1)[1]
mu.garch11.norm.w

sigma.garch11.norm.w <- predict(mean1.w, 1)[3] 
sigma.garch11.norm.w


## Normal weekly, 0.1
# h=1
VaR_0.1_h1.garch11.w <- exp(qnorm(0.1)* sigma.garch11.norm.w  + mu.garch11.norm.w) - 1
VaR_0.1_h1.garch11.w

#h=5
VaR_0.1_h5.garch11.w <- exp(qnorm(0.1)* sqrt((sigma.garch11.norm.w^2*5)) + 5*mu.garch11.norm.w) - 1
VaR_0.1_h5.garch11.w

#h=10
VaR_0.1_h10.garch11.w <- exp(qnorm(0.1)* sqrt((sigma.garch11.norm.w^2*10)) + 10*mu.garch11.norm.w) - 1
VaR_0.1_h10.garch11.w

normal_0.1.garch11.w <- rbind(VaR_0.1_h1.garch11.w, VaR_0.1_h5.garch11.w,VaR_0.1_h10.garch11.w  )
######################################################################################################


## Normal weekly, 0.05
# h=1
VaR_0.05_h1.garch11.w <- exp(qnorm(0.05)* sigma.garch11.norm.w  + mu.garch11.norm.w) - 1
VaR_0.05_h1.garch11.w

#h=5
VaR_0.05_h5.garch11.w <- exp(qnorm(0.05)* sqrt((sigma.garch11.norm.w^2*5)) + 5*mu.garch11.norm.w) - 1
VaR_0.05_h5.garch11.w

#h=10
VaR_0.05_h10.garch11.w <- exp(qnorm(0.05)* sqrt((sigma.garch11.norm.w^2*10)) + 10*mu.garch11.norm.w) - 1
VaR_0.05_h10.garch11.w

normal_0.05.garch11.w <- rbind(VaR_0.05_h1.garch11.w,  VaR_0.05_h5.garch11.w, VaR_0.05_h10.garch11.w )
######################################################################################################

## Normal weekly, 0.01
# h=1
VaR_0.01_h1.garch11.w <- exp(qnorm(0.01)* sigma.garch11.norm.w  + mu.garch11.norm.w) - 1
VaR_0.01_h1.garch11.w

#h=5
VaR_0.01_h5.garch11.w <- exp(qnorm(0.01)* sqrt((sigma.garch11.norm.w^2*5)) + 5*mu.garch11.norm.w) - 1
VaR_0.01_h5.garch11.w

#h=10
VaR_0.01_h10.garch11.w <- exp(qnorm(0.01)* sqrt((sigma.garch11.norm.w^2*10)) + 10*mu.garch11.norm.w) - 1
VaR_0.01_h10.garch11.w

normal_0.01.garch11.w <- rbind(VaR_0.01_h1.garch11.w,  VaR_0.01_h5.garch11.w, VaR_0.01_h10.garch11.w)
######################################################################################################

## Normal weekly, 0.005
# h=1
VaR_0.005_h1.garch11.w <- exp(qnorm(0.005)* sigma.garch11.norm.w  + mu.garch11.norm.w) - 1
VaR_0.005_h1.garch11.w

#h=5
VaR_0.005_h5.garch11.w <- exp(qnorm(0.005)* sqrt((sigma.garch11.norm.w^2*5)) + 5*mu.garch11.norm.w) - 1
VaR_0.005_h5.garch11.w

#h=10
VaR_0.005_h10.garch11.w <- exp(qnorm(0.05)* sqrt((sigma.garch11.norm.w^2*10)) + 10*mu.garch11.norm.w) - 1
VaR_0.05_h10.garch11.w

normal_0.005.garch11.w <- rbind(VaR_0.005_h1.garch11.w, VaR_0.005_h5.garch11.w, VaR_0.05_h10.garch11.w)

######################################################################################################

## normal weekly, 0.001
# h=1
VaR_0.001_h1.garch11.w <- exp(qnorm(0.001)* sigma.garch11.norm.w  + mu.garch11.norm.w) - 1
VaR_0.001_h1.garch11.w

#h=5
VaR_0.001_h5.garch11.w <- exp(qnorm(0.001)* sqrt((sigma.garch11.norm.w^2*5)) + 5*mu.garch11.norm.w) - 1
VaR_0.001_h5.garch11.w

#h=10
VaR_0.001_h10.garch11.w <- exp(qnorm(0.001)* sqrt((sigma.garch11.norm.w^2*10)) + 10*mu.garch11.norm.w) - 1
VaR_0.001_h10.garch11.w

normal_0.001.garch11.w <- rbind(VaR_0.001_h1.garch11.w, VaR_0.001_h5.garch11.w, VaR_0.001_h10.garch11.w)


######################################################################################################
######################################################################################################

# 4) A GARCH (1,1) model based on T-distribution with constant conditional mean function:

# GARCH T-dist Daily
mean1.T <- garchFit(formula = ~garch(1,1), data=daily, cond.dist="std", trace = F)
mean1.T
summary(mean1.T)

mu.garch11.T <- predict(mean1.T, 1)[1]
mu.garch11.T

sigma.garch11.T <- predict(mean1.T, 1)[3] 
sigma.garch11.T

str(daily.T) # 18838 obs n-1 df
str(weekly.T) # 4685 obs n-1 df


# DAILY DJ T-dist:

df=18838 # n-1 degrees of freedom
# h=1
VaR_0.1_h1.T.garch11 <- exp(qt(0.1, df) * sigma.garch11.T + mu.garch11.T)-1
VaR_0.1_h1.T.garch11

# h=5
VaR_0.1_h5.T.garch11 <- exp(qt(0.1, df)* sqrt((sigma.garch11.T^2*5)) + 5*mu.garch11.T) - 1
VaR_0.1_h5.T.garch11

# h=10
VaR_0.1_h10.T.garch11 <- exp(qt(0.1, df)* sqrt((sigma.garch11.T^2*10)) + 10*mu.garch11.T) - 1
VaR_0.1_h10.T.garch11

Tdist_0.1.daily.garch11 <- rbind(VaR_0.1_h1.T.garch11, VaR_0.1_h5.T.garch11 , VaR_0.1_h10.T.garch11)
####################################################################################################

df=18838 # n-1 degrees of freedom
# h=1
VaR_0.05_h1.T.garch11 <- exp(qt(0.05, df) * sigma.garch11.T + mu.garch11.T)-1
VaR_0.05_h1.T.garch11

# h=5
VaR_0.05_h5.T.garch11 <- exp(qt(0.05, df)* sqrt((sigma.garch11.T^2*5)) + 5*mu.garch11.T) - 1
VaR_0.05_h5.T.garch11

# h=10
VaR_0.05_h10.T.garch11 <- exp(qt(0.1, df)* sqrt((sigma.garch11.T^2*10)) + 10*mu.garch11.T) - 1
VaR_0.05_h10.T.garch11

Tdist_0.05.daily.garch11 <- rbind(VaR_0.05_h1.T.garch11, VaR_0.05_h5.T.garch11,  VaR_0.05_h10.T.garch11)
####################################################################################################


df=18838 # n-1 degrees of freedom
# h=1
VaR_0.01_h1.T.garch11 <- exp(qt(0.01, df) * sigma.garch11.T + mu.garch11.T)-1
VaR_0.01_h1.T.garch11

# h=5
VaR_0.01_h5.T.garch11 <- exp(qt(0.01, df)* sqrt((sigma.garch11.T^2*5)) + 5*mu.garch11.T) - 1
VaR_0.01_h5.T.garch11

# h=10
VaR_0.01_h10.T.garch11 <- exp(qt(0.01, df)* sqrt((sigma.garch11.T^2*10)) + 10*mu.garch11.T) - 1
VaR_0.01_h10.T.garch11

Tdist_0.01.daily.garch11 <- rbind(VaR_0.01_h1.T.garch11, VaR_0.01_h5.T.garch11,  VaR_0.01_h10.T.garch11)
####################################################################################################

df=18838 # n-1 degrees of freedom
# h=1
VaR_0.005_h1.T.garch11 <- exp(qt(0.005, df) * sigma.garch11.T + mu.garch11.T)-1
VaR_0.005_h1.T.garch11

# h=5
VaR_0.005_h5.T.garch11 <- exp(qt(0.005, df)* sqrt((sigma.garch11.T^2*5)) + 5*mu.garch11.T) - 1
VaR_0.005_h5.T.garch11

# h=10
VaR_0.005_h10.T.garch11 <- exp(qt(0.005, df)* sqrt((sigma.garch11.T^2*10)) + 10*mu.garch11.T) - 1
VaR_0.005_h10.T.garch11

Tdist_0.005.daily.garch11 <- rbind(VaR_0.005_h1.T.garch11, VaR_0.005_h5.T.garch11,  VaR_0.005_h10.T.garch11)

####################################################################################################

df=18838 # n-1 degrees of freedom
# h=1
VaR_0.001_h1.T.garch11 <- exp(qt(0.001, df) * sigma.garch11.T + mu.garch11.T)-1
VaR_0.001_h1.T.garch11

# h=5
VaR_0.001_h5.T.garch11 <- exp(qt(0.001, df)* sqrt((sigma.garch11.T^2*5)) + 5*mu.garch11.T) - 1
VaR_0.001_h5.T.garch11

# h=10
VaR_0.001_h10.T.garch11 <- exp(qt(0.001, df)* sqrt((sigma.garch11.T^2*10)) + 10*mu.garch11.T) - 1
VaR_0.001_h10.T.garch11

Tdist_0.001.daily.garch11 <- rbind(VaR_0.001_h1.T.garch11, VaR_0.001_h5.T.garch11,  VaR_0.001_h10.T.garch11)

####################################################################################################
####################################################################################################

# Same as above however with weekly DJ data.


# GARCH T-dist Daily
mean1.T.w <- garchFit(formula = ~garch(1,1), data=weekly.T, cond.dist="std", trace = F)
mean1.T.w
summary(mean1.T.w)

mu.garch11.T.w <- predict(mean1.T.w, 1)[1]
mu.garch11.T.w

sigma.garch11.T.w <- predict(mean1.T.w, 1)[3] 
sigma.garch11.T.w

# WEEKLY DJ T-dist:

str(weekly.T)
df=4685 # n-1 degrees of freedom
# h=1
VaR_0.1_h1.T.garch11.w <- exp(qt(0.1, df) * sigma.garch11.T.w + mu.garch11.T.w)-1
VaR_0.1_h1.T.garch11.w

# h=5
VaR_0.1_h5.T.garch11.w <- exp(qt(0.1, df)* sqrt((sigma.garch11.T.w^2*5)) + 5*mu.garch11.T.w) - 1
VaR_0.1_h5.T.garch11.w

# h=10
VaR_0.1_h10.T.garch11.w <- exp(qt(0.1, df)* sqrt((sigma.garch11.T.w^2*10)) + 10*mu.garch11.T.w) - 1
VaR_0.1_h10.T.garch11.w

Tdist_0.1.daily.garch11.w <- rbind(VaR_0.1_h1.T.garch11.w, VaR_0.1_h5.T.garch11.w,VaR_0.1_h10.T.garch11.w)
####################################################################################################

df=4685# n-1 degrees of freedom
# h=1
VaR_0.05_h1.T.garch11.w <- exp(qt(0.05, df) * sigma.garch11.T.w + mu.garch11.T.w)-1
VaR_0.05_h1.T.garch11.w

# h=5
VaR_0.05_h5.T.garch11.w <- exp(qt(0.05, df)* sqrt((sigma.garch11.T.w^2*5)) + 5*mu.garch11.T.w) - 1
VaR_0.05_h5.T.garch11.w

# h=10
VaR_0.05_h10.T.garch11.w <- exp(qt(0.1, df)* sqrt((sigma.garch11.T.w^2*10)) + 10*mu.garch11.T.w) - 1
VaR_0.05_h10.T.garch11.w

Tdist_0.05.daily.garch11.w <- rbind(VaR_0.05_h1.T.garch11.w, VaR_0.05_h5.T.garch11.w,  VaR_0.05_h10.T.garch11.w)
####################################################################################################


df=4685 # n-1 degrees of freedom
# h=1
VaR_0.01_h1.T.garch11.w <- exp(qt(0.01, df) * sigma.garch11.T.w + mu.garch11.T.w)-1
VaR_0.01_h1.T.garch11.w

# h=5
VaR_0.01_h5.T.garch11.w <- exp(qt(0.01, df)* sqrt((sigma.garch11.T.w^2*5)) + 5*mu.garch11.T.w) - 1
VaR_0.01_h5.T.garch11.w

# h=10
VaR_0.01_h10.T.garch11.w <- exp(qt(0.01, df)* sqrt((sigma.garch11.T.w^2*10)) + 10*mu.garch11.T.w) - 1
VaR_0.01_h10.T.garch11.w

Tdist_0.01.daily.garch11.w <- rbind(VaR_0.01_h1.T.garch11.w, VaR_0.01_h5.T.garch11.w,  VaR_0.01_h10.T.garch11.w)
####################################################################################################

df=4685 # n-1 degrees of freedom
# h=1
VaR_0.005_h1.T.garch11.w <- exp(qt(0.005, df) * sigma.garch11.T.w + mu.garch11.T.w)-1
VaR_0.005_h1.T.garch11.w

# h=5
VaR_0.005_h5.T.garch11.w <- exp(qt(0.005, df)* sqrt((sigma.garch11.T.w^2*5)) + 5*mu.garch11.T.w) - 1
VaR_0.005_h5.T.garch11.w

# h=10
VaR_0.005_h10.T.garch11.w <- exp(qt(0.005, df)* sqrt((sigma.garch11.T.w^2*10)) + 10*mu.garch11.T.w) - 1
VaR_0.005_h10.T.garch11.w

Tdist_0.005.daily.garch11.w <- rbind(VaR_0.005_h1.T.garch11.w, VaR_0.005_h5.T.garch11.w,  VaR_0.005_h10.T.garch11.w)

####################################################################################################

df=4685 # n-1 degrees of freedom
# h=1
VaR_0.001_h1.T.garch11.w <- exp(qt(0.001, df) * sigma.garch11.T.w + mu.garch11.T.w)-1
VaR_0.001_h1.T.garch11.w

# h=5
VaR_0.001_h5.T.garch11.w <- exp(qt(0.001, df)* sqrt((sigma.garch11.T.w^2*5)) + 5*mu.garch11.T.w) - 1
VaR_0.001_h5.T.garch11.w

# h=10
VaR_0.001_h10.T.garch11.w <- exp(qt(0.001, df)* sqrt((sigma.garch11.T.w^2*10)) + 10*mu.garch11.T.w) - 1
VaR_0.001_h10.T.garch11.w

Tdist_0.001.daily.garch11.w <- rbind(VaR_0.001_h1.T.garch11.w, VaR_0.001_h5.T.garch11.w,  VaR_0.001_h10.T.garch11.w)

####################################################################################################
####################################################################################################


# PLOTS: 


## Normal Daily DJ Data:
Normal_0.01 <- rbind(VaR_0.01_h1, VaR_0.01_h5, VaR_0.01_h10 )
Normal_0.05 <- rbind(VaR_0.05_h1, VaR_0.05_h5,VaR_0.05_h10 )
Normal_0.005 <- rbind(VaR_0.005_h1, VaR_0.005_h5,  VaR_0.005_h10 )
Normal_0.1 <- rbind(VaR_0.1_h1, VaR_0.1_h5,  VaR_0.1_h10 )
Normal_0.001 <- rbind(VaR_0.001_h1, VaR_0.001_h5,  VaR_0.001_h10 )


# plotting Daily DJ normal dist.
plot(horizon, Normal_0.01, type="o", col="Green", pch="o", lty=1, ylim=c(-0.15,0), xlab="Horizon", ylab="Loss in Percentage (%)", main = "Daily DJ Returns | Unconditional Normal Distribution" )

points(horizon, Normal_0.05, col="Purple", pch="*")
lines(horizon, Normal_0.05, col="Purple",lty=2)

points(horizon, Normal_0.005, col="Red",pch="+")
lines(horizon, Normal_0.005, col="Red", lty=3)

points(horizon, Normal_0.1 , col="Blue", pch="*")
lines(horizon, Normal_0.1 , col="Blue",lty=2)

points(horizon, Normal_0.001, col="black", pch="o")
lines(horizon, Normal_0.001, col="black",lty=2)

legend("bottomleft",legend=c("0.01", "0.05","0.005", "0.1", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1)

#########################################################################################################################################


# Normal Weekly DJ data:
# Unconditional weekly DJ normal
Normal_0.01.weekly <- rbind(VaR_0.01_h1.weekly, VaR_0.01_h5.weekly, VaR_0.01_h10.weekly)
Normal_0.05.weekly <- rbind(VaR_0.05_h1.weekly, VaR_0.05_h5.weekly,VaR_0.05_h10.weekly )
Normal_0.005.weekly <- rbind(VaR_0.005_h1.weekly, VaR_0.005_h5.weekly,  VaR_0.005_h10.weekly )
Normal_0.1.weekly <- rbind(VaR_0.1_h1.weekly, VaR_0.1_h5.weekly,  VaR_0.1_h10.weekly )
Normal_0.001.weekly <- rbind(VaR_0.001_h1.weekly, VaR_0.001_h5.weekly,  VaR_0.001_h10.weekly )

plot(horizon, Normal_0.01.weekly, type="o", col="Green", pch="o", lty=1, ylim=c(-0.25,0), xlab="Horizon", ylab="Loss in Percentage (%)", main = "Weekly DJ Returns | Unconditional Normal Distribution" )

points(horizon, Normal_0.05.weekly, col="Purple", pch="*")
lines(horizon, Normal_0.05.weekly, col="Purple",lty=2)

points(horizon, Normal_0.005.weekly, col="Red",pch="+")
lines(horizon, Normal_0.005.weekly, col="Red", lty=3)

points(horizon, Normal_0.1.weekly , col="Blue", pch="*")
lines(horizon, Normal_0.1.weekly , col="Blue",lty=2)

points(horizon, Normal_0.001.weekly, col="black", pch="o")
lines(horizon, Normal_0.001.weekly, col="black",lty=2)

legend("bottomleft",legend=c("0.1", "0.05","0.005", "0.1", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1)
#########################################################################################################################################


# Unconditional Daily DJ T-distribution
Tdist_0.01.daily <- rbind(VaR_0.01_h1.T, VaR_0.01_h5.T,VaR_0.01_h10.T )
Tdist_0.1.daily <- rbind(VaR_0.1_h1.T, VaR_0.1_h5.T, VaR_0.1_h10.T  )
Tdist_0.05.daily <- rbind(VaR_0.05_h1.T, VaR_0.05_h5.T, VaR_0.05_h10.T )
Tdist_0.005.daily <- rbind(VaR_0.005_h1.T, VaR_0.005_h5.T, VaR_0.005_h10.T)
Tdist_0.001.daily <- rbind(VaR_0.001_h1.T, VaR_0.001_h5.T, VaR_0.001_h10)


plot(horizon, Tdist_0.01.daily, type="o", col="Green", pch="o", lty=1, ylim=c(-0.15,0), xlab="Horizon", ylab="Loss in Percentage (%)",
     main = "Daily DJ Returns | Unconditional T-Distribution" )

points(horizon, Tdist_0.1.daily, col="Purple", pch="*")
lines(horizon, Tdist_0.1.daily, col="Purple",lty=2)

points(horizon, Tdist_0.05.daily, col="Red",pch="+")
lines(horizon, Tdist_0.05.daily, col="Red", lty=3)

points(horizon, Tdist_0.005.daily , col="Blue", pch="*")
lines(horizon, Tdist_0.005.daily , col="Blue",lty=2)

points(horizon, Tdist_0.001.daily, col="black", pch="o")
lines(horizon, Tdist_0.001.daily, col="black",lty=2)

legend("bottomleft",legend=c("0.01", "0.1","0.05", "0.005", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1)

#########################################################################################################################################


# Unconditional Weekly DJ T-distribution
Tdist_0.1.weekly <- rbind(VaR_0.1_h1.T.weekly, VaR_0.1_h5.T.weekly,VaR_0.1_h10.T.weekly)
Tdist_0.05.weekly <- rbind(VaR_0.05_h1.T.weekly,VaR_0.05_h5.T.weekly, VaR_0.05_h10.T.weekly)
Tdist_0.01.weekly <- rbind(VaR_0.01_h1.T.weekly, VaR_0.01_h5.T.weekly,VaR_0.01_h10.T.weekly)
Tdist_0.005.weekly <- rbind(VaR_0.005_h1.T.weekly,VaR_0.005_h5.T.weekly, VaR_0.005_h10.T.weekly)
Tdist_0.001.weekly <- rbind(VaR_0.001_h1.T.weekly, VaR_0.001_h1.T.weekly, VaR_0.001_h10.T.weekly )

plot(horizon, Tdist_0.1.weekly, type="o", col="Green", pch="o", lty=1, ylim=c(-0.30,0), xlab="Horizon", ylab="Loss in Percentage (%)",
     main = "Weekly DJ Returns | Unconditional T-Distribution" )

points(horizon, Tdist_0.05.weekly, col="Purple", pch="*")
lines(horizon, Tdist_0.05.weekly, col="Purple",lty=2)

points(horizon, Tdist_0.01.weekly, col="Red",pch="+")
lines(horizon, Tdist_0.01.weekly, col="Red", lty=3)

points(horizon, Tdist_0.005.weekly , col="Blue", pch="*")
lines(horizon, Tdist_0.005.weekly , col="Blue",lty=2)

points(horizon, Tdist_0.001.weekly, col="black", pch="o")
lines(horizon, Tdist_0.001.weekly, col="black",lty=2)

legend("bottomleft",legend=c("0.1", "0.05","0.01", "0.005", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1)                           



####################################################################################################################################
# Gaussian GARCH (1,1) with costant conditional mean function: daily DJ
normal_0.1.garch11 <- rbind(VaR_0.1_h1.garch11, VaR_0.1_h5.garch11,VaR_0.1_h10.garch11)
normal_0.05.garch11 <- rbind(VaR_0.05_h1.garch11,  VaR_0.05_h5.garch11,VaR_0.05_h10.garch11)
normal_0.01.garch11 <- rbind(VaR_0.01_h1.garch11,  VaR_0.01_h5.garch11, VaR_0.01_h10.garch11)
normal_0.005.garch11 <- rbind(VaR_0.005_h1.garch11, VaR_0.005_h5.garch11,VaR_0.005_h10.garch11)
normal_0.001.garch11 <- rbind(VaR_0.001_h1.garch11, VaR_0.001_h5.garch11,VaR_0.001_h10.garch11)

normal_0.1.garch11$standardDeviation
normal_0.05.garch11$standardDeviation
normal_0.01.garch11$standardDeviation
normal_0.005.garch11$standardDeviation
normal_0.001.garch11$standardDeviation

plot(horizon, normal_0.1.garch11$standardDeviation, type="o", col="Green", pch="o", lty=1, ylim=c(-0.17,0), xlab="Horizon", ylab="Loss in Percentage (%)",
     main = "Daily DJ Returns | Gaussian GARCH(1,1) " )

points(horizon, normal_0.05.garch11$standardDeviation, col="Purple", pch="*")
lines(horizon, normal_0.05.garch11$standardDeviation, col="Purple",lty=2)

points(horizon, normal_0.01.garch11$standardDeviation, col="Red",pch="+")
lines(horizon, normal_0.01.garch11$standardDeviation, col="Red", lty=3)

points(horizon, normal_0.005.garch11$standardDeviation , col="Blue", pch="*")
lines(horizon, normal_0.005.garch11$standardDeviation , col="Blue",lty=2)

points(horizon, normal_0.001.garch11$standardDeviation, col="black", pch="o")
lines(horizon, normal_0.001.garch11$standardDeviation, col="black",lty=2)

legend("bottomleft",legend=c("0.1", "0.05","0.01", "0.005", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1)


#############################################################################################
# weekly daily DJ data - gaussian garch (1,1)

normal_0.1.garch11.w <- rbind(VaR_0.1_h1.garch11.w, VaR_0.1_h5.garch11.w, VaR_0.1_h10.garch11.w)
normal_0.05.garch11.w <- rbind(VaR_0.05_h1.garch11.w,  VaR_0.05_h5.garch11.w, VaR_0.05_h10.garch11.w)
normal_0.01.garch11.w <- rbind(VaR_0.01_h1.garch11.w,  VaR_0.01_h5.garch11.w, VaR_0.01_h10.garch11.w)
normal_0.005.garch11.w <- rbind(VaR_0.005_h1.garch11.w, VaR_0.005_h5.garch11.w, VaR_0.005_h10.garch11.w)
normal_0.001.garch11.w <- rbind(VaR_0.001_h1.garch11.w, VaR_0.001_h5.garch11.w, VaR_0.001_h10.garch11.w)

plot(horizon, normal_0.1.garch11.w$standardDeviation, type="o", col="Green", pch="o", lty=1, ylim=c(-0.25,0), xlab="Horizon", ylab="Loss in Percentage (%)",
     main = "Weekly DJ Returns | Gaussian GARCH(1,1) " )

points(horizon, normal_0.05.garch11.w$standardDeviation, col="Purple", pch="*")
lines(horizon, normal_0.05.garch11.w$standardDeviation, col="Purple",lty=2)

points(horizon, normal_0.01.garch11.w$standardDeviation, col="Red",pch="+")
lines(horizon, normal_0.01.garch11.w$standardDeviation, col="Red", lty=3)

points(horizon, normal_0.005.garch11.w$standardDeviation , col="Blue", pch="*")
lines(horizon, normal_0.005.garch11.w$standardDeviation , col="Blue",lty=2)

points(horizon, normal_0.001.garch11.w$standardDeviation, col="black", pch="o")
lines(horizon, normal_0.001.garch11.w$standardDeviation, col="black",lty=2)

legend("bottomleft",legend=c("0.1", "0.05","0.01", "0.005", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1, cex=1)


#############################################################################################
#############################################################################################

# GARCH(1,1) based on T-distribution
Tdist_0.1.daily.garch11 <- rbind(VaR_0.1_h1.T.garch11, VaR_0.1_h5.T.garch11 , VaR_0.1_h10.T.garch11)
Tdist_0.05.daily.garch11 <- rbind(VaR_0.05_h1.T.garch11, VaR_0.05_h5.T.garch11,  VaR_0.05_h10.T.garch11)
Tdist_0.01.daily.garch11 <- rbind(VaR_0.01_h1.T.garch11, VaR_0.01_h5.T.garch11,  VaR_0.01_h10.T.garch11)
Tdist_0.005.daily.garch11 <- rbind(VaR_0.005_h1.T.garch11, VaR_0.005_h5.T.garch11,  VaR_0.005_h10.T.garch11)
Tdist_0.001.daily.garch11 <- rbind(VaR_0.001_h1.T.garch11, VaR_0.001_h5.T.garch11,  VaR_0.001_h10.T.garch11)


plot(horizon, Tdist_0.1.daily.garch11$standardDeviation, type="o", col="Green", pch="o", lty=1, ylim=c(-0.15,0), xlab="Horizon", ylab="Loss in Percentage (%)",
     main = "Daily DJ Returns | GARCH(1,1) T-Distribution " )

points(horizon, Tdist_0.05.daily.garch11$standardDeviation, col="Purple", pch="*")
lines(horizon, Tdist_0.05.daily.garch11$standardDeviation, col="Purple",lty=2)

points(horizon, Tdist_0.01.daily.garch11$standardDeviation, col="Red",pch="+")
lines(horizon, Tdist_0.01.daily.garch11$standardDeviation, col="Red", lty=3)

points(horizon, Tdist_0.005.daily.garch11$standardDeviation , col="Blue", pch="*")
lines(horizon, Tdist_0.005.daily.garch11$standardDeviation , col="Blue",lty=2)

points(horizon, Tdist_0.001.daily.garch11$standardDeviation, col="black", pch="o")
lines(horizon, Tdist_0.001.daily.garch11$standardDeviation, col="black",lty=2)

legend("bottomleft",legend=c("0.1", "0.05","0.01", "0.005", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1, cex=1)


#################

# GARCH(1,1) T-distribution:

Tdist_0.1.daily.garch11.w <- rbind(VaR_0.1_h1.T.garch11.w, VaR_0.1_h5.T.garch11.w,VaR_0.1_h10.T.garch11.w)
Tdist_0.05.daily.garch11.w <- rbind(VaR_0.05_h1.T.garch11.w, VaR_0.05_h5.T.garch11.w,  VaR_0.05_h10.T.garch11.w)
Tdist_0.01.daily.garch11.w <- rbind(VaR_0.01_h1.T.garch11.w, VaR_0.01_h5.T.garch11.w,  VaR_0.01_h10.T.garch11.w)
Tdist_0.005.daily.garch11.w <- rbind(VaR_0.005_h1.T.garch11.w, VaR_0.005_h5.T.garch11.w,  VaR_0.005_h10.T.garch11.w)
Tdist_0.001.daily.garch11.w <- rbind(VaR_0.001_h1.T.garch11.w, VaR_0.001_h5.T.garch11.w,  VaR_0.001_h10.T.garch11.w)


plot(horizon, Tdist_0.1.daily.garch11.w$standardDeviation, type="o", col="Green", pch="o", lty=1, ylim=c(-0.25,0), xlab="Horizon", ylab="Loss in Percentage (%)",
     main = "Weekly DJ Returns | GARCH(1,1) T-Distribution " )

points(horizon, Tdist_0.05.daily.garch11.w$standardDeviation, col="Purple", pch="*")
lines(horizon, Tdist_0.05.daily.garch11.w$standardDeviation, col="Purple",lty=2)

points(horizon, Tdist_0.01.daily.garch11.w$standardDeviation, col="Red",pch="+")
lines(horizon, Tdist_0.01.daily.garch11.w$standardDeviation, col="Red", lty=3)

points(horizon, Tdist_0.005.daily.garch11.w$standardDeviation , col="Blue", pch="*")
lines(horizon, Tdist_0.005.daily.garch11.w$standardDeviation , col="Blue",lty=2)

points(horizon, Tdist_0.001.daily.garch11.w$standardDeviation, col="black", pch="o")
lines(horizon, Tdist_0.001.daily.garch11.w$standardDeviation, col="black",lty=2)

legend("bottomleft",legend=c("0.1", "0.05","0.01", "0.005", "0.001"), col=c("Green","Purple","Red", "Blue", "Black"),
       pch=c("o","","+", "", "o"),lty=c(1,2,3), ncol=1, cex=1)


######################################################################################################
######################################################################################################
# END
######################################################################################################
######################################################################################################



