# Interest Rate Risk
# 1

library(quantmod)
library(lubridate)
library(PerformanceAnalytics)
n <- c("SPY","IEF","SHY")
p_list1 <- lapply(n, function(s) get(getSymbols(s, from = "2002-07-01")) )
p_list2 <- lapply(p_list1, function(z) z[,6])
r_list <- lapply(p_list2, function(x) na.omit(log(x/lag(x)) ) )
# yearly return
ry_list <- lapply(r_list, function(v) apply.yearly(v, sum))
ry_list

# 2
# monthly return
rm_list <- lapply(r_list, function(v) apply.monthly(v, sum))
r <- na.omit(Reduce(merge, rm_list) )

# ETF adjusted performance 
charts.PerformanceSummary(r, main= "ETF Risk-Return")

m <- apply(r, 2, mean)*12
# SPY.Adjusted IEF.Adjusted SHY.Adjusted 
# 0.08568132   0.04864408   0.02010095 

s <- apply(r, 2 , sd)*sqrt(12)
# SPY.Adjusted IEF.Adjusted SHY.Adjusted 
# 0.13967881   0.06347702   0.01316692 
plot(m-s)

# Pearson correlation coefficients
cor(r)
# SPY.Adjusted IEF.Adjusted SHY.Adjusted
# SPY.Adjusted    1.0000000   -0.2976884   -0.3184548
# IEF.Adjusted   -0.2976884    1.0000000    0.7592255
# SHY.Adjusted   -0.3184548    0.7592255    1.0000000

# 3
fed <- get(getSymbols("FEDFUNDS", src= "FRED", from= "2000-01-01") )
head(fed)

fed$delta <- fed$FEDFUNDS - lag(fed$FEDFUNDS)

a <- date(r[1:10,1])
ceiling_date(a, "m")
fix <- ceiling_date(date(r), "m")

data1 <- data.frame(Date= fix, r)
data2 <- data.frame(Date= date(fed), fed/100)
data3 <- merge(data1, data2)
head(data3)
length(data3)

y1 <- cumsum(data3$SPY.Adjusted)
y2 <- cumsum(data3$IEF.Adjusted)
y3 <- cumsum(data3$SHY.Adjusted)
y4 <- data3$FEDFUNDS
x <- data3$Date
plot(y1~x, type="l", col= 1 , lwd = 2)
lines(y2~x, col= 2, lwd = 2)
lines(y3~x, col= 3, lwd = 2)
lines(y4~x, col= 4, lwd = 2)
par(new= TRUE)
plot(y4~x, type= "l", col = 5, axes= FALSE, labels=FALSE)
plot(y4~x, type= "l", col = 5, axis(4), ylab = "FED Rates") 

lm(SPY.Adjusted~delta, data= data3)
# beta=  6.292874 
lm(IEF.Adjusted~delta, data= data3)
# beta = -1.528488  
lm(SHY.Adjusted~delta, data= data3)
# beta = -0.983706 


# Managing Linear Risk
# 1.
library(lubridate)
b <- ("JPY=X")
s <- get(getSymbols(b, from = "2019-01-01", to = "2019-10-31"))
s <- s$`JPY=X.Adjusted`

plot(s, main = "JPY/USD")

#JPY/USD
s <- 1/s
plot(s, main = "USD/JPY")

# current date
t0 <- date("2019-10-31") 
# deliver of amount date
t1 <- date("2020-05-01")
# payment date
t2 <- date("2020-06-16") 

#spot rate
s0 <- last(s)

r <- na.omit(log(s/lag(s)) )
plot(r)
m <- mean(r)
# [1] 4.219405e-05
sd1 <- sd(r)
# [1] 0.00390532

sigma <- sqrt(252)*sd1
# [1] 0.06199503

theta <- m*252 + (sigma^2)/2
# [1] 0.01255459

# Simulation
s0 <- as.numeric(s0)
# [1] 0.009192106

drift <- (theta - (sigma^2)/2)
volat <- sigma

days <- as.numeric(t1 - t0)
trade_days <- round(days*252/365, 0)
# [1] 126
tau <- days/365
# [1] 0.5013699

N <- 1000
B_tau <- rnorm(N,0,sqrt(tau))
s1 <- s0*exp( drift*tau + volat*B_tau )
mean(s1)
# [1] 0.009255515

s0*exp(theta*tau)
# [1] 0.009250148

# VaR for the Unhedged
# 2 
q <- 125*(10^6)
# portfolio value
pl <- q*(s1-s0)
hist(pl)

mean(pl)
# [1] 7926.193

q*s0*(exp(theta*tau)-1)
# [1] 7255.267

var_pl <- mean(pl)- quantile (pl,0.01)
var_pl/1000
# 107.3918

# Unitary Hedge
# 3
rdiff <- theta
tau20 <- as.numeric(t2-t0)/365
tau21 <- as.numeric(t2-t1)/365
# f0 is fixed
f0 <- s0*exp(rdiff*tau20)
f1 <- s1*exp(rdiff*tau21)
pl <- q*(s1-s0) + q*(f0-f1)
hist(pl)

cor(f1-f0,s1-s0)
# [1] 1

mean(pl)
# [1] 7254.204
var_pl <- mean(pl)- quantile (pl,0.01)
var_pl/1000
# 0.1952814 

# future contract expiring in September 2020
tau21 <- as.numeric(t2-date("2020-09-01"))/365
b1 <- s1-f1
b2 <- s0-f0
mean(b2>b1)
# [1] 0

# Hedging using ETFs
# 4
# a
p <- get(getSymbols("IEF") )
p <- p$IEF.Adjusted
head(p)
head(s)
p1 <- na.omit(merge(s,p) )
head(p1)
s <- s0
f <- 112.3420
q <- 125* (10^6)
qf <- 1
head(p1)
r1 <- na.omit(log( p1/lag(p1)) )
cor(r1)
#                   JPY.X.Adjusted IEF.Adjusted
# JPY.X.Adjusted      1.0000000    0.1270229
# IEF.Adjusted        0.1270229    1.0000000

u <- lm(r1$JPY.X.Adjusted~ r1$IEF.Adjusted)
# beta: 1.489e-01

summary(u)
# Multiple R-squared:  0.01613,	Adjusted R-squared:  0.01138 

N <- -(1.489e-01)*(s0*q/f0*qf)
# [1] -18466470


