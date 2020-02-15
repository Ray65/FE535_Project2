library(quantmod)
library(lubridate)
library(plotrix)
library(gsubfn)
library(PerformanceAnalytics)
library(data.table)

# Problem 1
# 1, 2
n <- c("AMZN","SPY")
# getting the first day equity data
dataeq <- list()
for(s in n){
  cat(s,"\n")
  s1 <- lapply(n, function(s) get(getSymbols(s, from= as.Date("2020-02-10"), to= as.Date("2020-02-12"))) )
  dataeq1 <- lapply(s1, function(z) z[,6])
}
dataeq1

# adding VIX
# library(qrmdata)
# data(VIX)
# get(getSymbols("^VIX", src="yahoo", from= as.Date("2020-02-10"), to= as.Date("2020-02-12")))
# dataeq1 <- c(dataeq1, VIX[1,6])
# dataeq1

# getting the first day option data
# option data in feb 10 for two matrutrities of 17th of April and 19th of July
#data1op1amzn <- getOptionChain("AMZN", src="yahoo", Exp = "2020-04-17")
#data1op2amzn <- getOptionChain("AMZN", src="yahoo", Exp = "2020-06-19")
#data1opamzn <- as.data.table(getOptionChain("AMZN", src="yahoo", NULL))
#write.csv(data1op1amzn$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1amzncall.csv")
#write.csv(data1op1amzn$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1amznput.csv")
#write.csv(data1op2amzn$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2amzncall.csv")
#write.csv(data1op2amzn$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2amznput.csv")
#write.csv(data1opamzn, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1opamzn.csv")

#data1op1spy <- getOptionChain("SPY", src="yahoo", Exp = "2020-04-17")
#data1op2spy <- getOptionChain("SPY", src="yahoo", Exp = "2020-06-19")
#data1opspy <- as.data.table(getOptionChain("SPY", src="yahoo", NULL))
#write.csv(data1op1spy$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1spycall.csv")
#write.csv(data1op1spy$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1spyput.csv")
#write.csv(data1op2spy$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2spycall.csv")
#write.csv(data1op2spy$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2spyput.csv")
#write.csv(data1opspy, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1opspy.csv")

# getting the data for 15th of April and 17th of July
#data1op1vix <- getOptionChain("^VIX", src="yahoo", Exp = "2020-04-15")
#data1op2vix <- getOptionChain("^VIX", src="yahoo", Exp = "2020-06-17")
#data1opvix <- as.data.table(getOptionChain("^VIX", src="yahoo", NULL))
#write.csv(data1op1vix$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1vixcall.csv")
#write.csv(data1op1vix$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1vixput.csv")
#write.csv(data1op2vix$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2vixcall.csv")
#write.csv(data1op2vix$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2vixput.csv")

# indicating type of the option for first day
# option_chain1 <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1amzncall.csv")
# option_chain1$type <- "call"
# write.csv(option_chain1, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optionchain1amzn1call.csv")
# option_chain <-read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1amznput.csv")
#option_chain$type <- "put"
# option_chain2 <-read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2amzncall.csv")
# option_chain2$type <- "call"
# write.csv(option_chain2, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optionchain2amzn1call.csv")
# option_chain <-read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2amznput.csv")
#option_chain$type <- "put"

option_chain1 <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1spycall.csv")
option_chain1$type <- "call"
write.csv(option_chain1, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optionchain1spy1call.csv")
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1spyput.csv")
# option_chain$type <- "put"
option_chain2 <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2spycall.csv")
option_chain2$type <- "call"
write.csv(option_chain2, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optionchain2spy1call.csv")
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2spyput.csv")
# option_chain$type <- "put"

# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1vixcall.csv")
# option_chain$type <- "call"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op1vixput.csv")
# option_chain$type <- "put"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2vixcall.csv")
# option_chain$type <- "call"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data1op2vixput.csv")
# option_chain$type <- "put"

# getting the second day option data
# option data in feb 10 for two matrutrities of 17th of April and 19th of July
# data2op1amzn <- getOptionChain("AMZN", src="yahoo", Exp = "2020-04-17")
# data2op2amzn <- getOptionChain("AMZN", src="yahoo", Exp = "2020-06-19")
# data2opamzn <- as.data.table(getOptionChain("AMZN", src="yahoo", NULL))
# write.csv(data2op1amzn$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1amzncall.csv")
# write.csv(data2op1amzn$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1amznput.csv")
# write.csv(data2op2amzn$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2amzncall.csv")
# write.csv(data2op2amzn$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2amznput.csv")

# data2op1spy <- getOptionChain("SPY", src="yahoo", Exp = "2020-04-17")
# data2op2spy <- getOptionChain("SPY", src="yahoo", Exp = "2020-06-19")
# data2opspy <- as.data.table(getOptionChain("SPY", src="yahoo", NULL))
# write.csv(data2op1spy$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1spycall.csv")
# write.csv(data2op1spy$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1spyput.csv")
# write.csv(data2op2spy$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2spycall.csv")
# write.csv(data2op2spy$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2spyput.csv")

# getting the data for 15th of April and 17th of July
# data2op1vix <- getOptionChain("^VIX", src="yahoo", Exp = "2020-04-15")
# data2op2vix <- getOptionChain("^VIX", src="yahoo", Exp = "2020-06-17")
# data2opvix <- as.data.table(getOptionChain("^VIX", src="yahoo", NULL))
# write.csv(data2op1vix$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1vixcall.csv")
# write.csv(data2op1vix$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1vixput.csv")
# write.csv(data2op2vix$calls, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2vixcall.csv")
# write.csv(data2op2vix$puts, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2vixput.csv")

# indicating type of the option for second day
# option_chain1 <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1amzncall.csv")
# option_chain1$type <- "call"
# write.csv(option_chain1, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optionchain1amzn1call.csv")
# option_chain <-read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1amznput.csv")
# option_chain$type <- "put"
# option_chain2 <-read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2amzncall.csv")
# option_chain2$type <- "call"
# write.csv(option_chain2, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optionchain2amzn1call.csv")
# option_chain <-read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2amznput.csv")
#option_chain$type <- "put"

# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1spycall.csv")
# option_chain$type <- "call"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1spyput.csv")
# option_chain$type <- "put"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2spycall.csv")
# option_chain$type <- "call"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2spyput.csv")
# option_chain$type <- "put"

# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1vixcall.csv")
# option_chain$type <- "call"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op1vixput.csv")
# option_chain$type <- "put"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2vixcall.csv")
# option_chain$type <- "call"
# option_chain <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/data2op2vixput.csv")
# option_chain$type <- "put"

# Part 2
# 2.5
# BS Call
s <- 2079.28
k <- 720
sig <- 0.8572
r <- 0.0159
tau <- 0.333
bsc <- function(s,k,sig,r,tau){
  d2 <- (log(s/k)+(r-(sig^2)/2)*tau)/(sqrt(tau)*sig)
  d1 <- d2+sig*sqrt(tau)
  c_price <- s*pnorm(d1)- k*exp(-r*tau)*pnorm(d2)
  return(c_price)
}
bsc(s,k,sig,r,tau)
t_seq <- seq(0, tau , length=10000)
p_time <- sapply(t_seq, function(x) bsc(s,k,sig,r,x))
plot(p_time~t_seq, type='l')

# BS Put
k <- 25
sig <- 0.5
bsp <- function(s,k,sig,r,tau){
  d2 <- (log(s/k)+(r-(sig^2)/2)*tau)/(sqrt(tau)*sig)
  d1 <- d2+sig*sqrt(tau)
  p_price <- k*exp(-r*tau)*pnorm(-d2)-s*pnorm(-d1)
  return(p_price)
}
bsp(s,k,sig,r,tau)
t_seq <- seq(0, tau , length=10000)
p_time <- sapply(t_seq, function(x) bsp(s,k,sig,r,x))
plot(p_time~t_seq, type='l', xlab= "time", )

#BS with if
bs <- function(s,k,r,sig,tau,type){
  d2 <- (log(s/k)+(r-(sig^2)/2)*tau)/(sqrt(tau)*sig)
  d1 <- d2+sig*sqrt(tau)
  if(type=="call"){
    c_price <- s*pnorm(d1)- k*exp(-r*tau)*pnorm(d2)
    return(c_price)
  }
  if(type== "put"){
    p_price <- k*exp(-r*tau)*pnorm(-d2)-s*pnorm(-d1)
    return(p_price)
  }
}

bs(2079.28,720,0.0159,0.8572,0.333,"call")
# [1] 1366.354

# 2.6 and 2.8
# merging for amazon with different expiries in day 1 and reading it
option_chainamzn1call <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/day2 12.csv",
                                 header= TRUE, sep= ',')
# calculating days to maturity
option_chainamzn1call$daystoexpiry <- as.Date(option_chainamzn1call$Expiry,
                                              "%m/%d/%y")- Sys.Date()
# getting the average of bid and ask
option_chainamzn1call$premium <- (option_chainamzn1call$Bid+option_chainamzn1call$Ask)/2
head(option_chainamzn1call)
dim(option_chainamzn1call)
write.csv(option_chainamzn1call,"/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/day2 1212.csv")
# bisection method
# function written with help from Saeed Rahman, Feb,19, 2017
bisection <-function(S, K, t, r, type, option_price, max.iter=100000 ,tolerance=.0001){
  sigma.upper <- 2
  sigma.lower <- 0.001
  sigma.mid <- .5
  count <- 0
  fun.mid <- bs(s=S,k=K,r=r,sig=sigma.mid,tau=t,type=type)- option_price
  start.time <- Sys.time()
  while(abs(fun.mid) > tolerance && count<max.iter){
    fun.upper <- bs(s=S,k=K,tau=t,r=r,sig=sigma.upper,type=type)-option_price
    fun.lower <- bs(s=S,k=K,t=t,r=r,sig=sigma.lower,type=type)-option_price
    fun.mid <- bs(s=S,k=K,t=t,r=r,sig=sigma.mid,type=type)-option_price
    if(fun.mid*fun.lower < 0){
      sigma.upper <-sigma.mid
      sigma.mid <- (sigma.upper + sigma.lower)/2
    }else{
      sigma.lower<- sigma.mid
      sigma.mid <- (sigma.lower + sigma.upper)/2
    }
    count <- count + 1
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if(count>=max.iter){
    return(list(NA,time.taken,count))
  }else{
    return(list(sigma.mid,time.taken,count))
  }
}

# bisection model with dummy variables to calculate implied volatility
iv <- bisection(S=100,K=100,t=30/252,r=.05,type='call',option_price = 3.051184)
print(paste("Implied Volatility:",iv[1],"Time taken for calculations:",iv[2], "seconds",
            "Number of iterations:",iv[3]))
# [1] "Implied Volatility= 0.199972595214844 Time taken for 
# calculations= 0.000211000442504883 seconds Number of iterations= 14"

# running bisection root model with sample data
bisectionsam <- function(x){
  #stock_df <- as.data.frame(getSymbols(symbol, src="yahoo" ,from = as.Date("2020-02-10"), env = NULL))
  option_chain <- option_chainamzn1call
  #head(option_chain$daystoexpiry)
  #daystoexpiry
  libor <- 0.0159
  iv <- {}
  original_iv <- {}
  optionname <- {}
  strike <- {}
  daystoexpiry <- {}
  time_taken <- 0
  iterations <- 0
  for(i in 1:nrow(option_chain)){
      #print(i)
      bisection1 <- bisection(
        S <- x,
        K <- as.numeric(option_chain[i,"Strike"]),
        t <- as.numeric(option_chain[i, "daystoexpiry"])/252,
        r <- libor,
        type <- ifelse((option_chain[i,"type"]=="call"),"call","put"),
        option_price <- as.numeric(option_chain[i,"premium"]))
      
      iv <- append(iv, as.numeric(bisection1[1]))
      if(!is.na(bisection1[1])){
        time_taken <- as.numeric(bisection1[2])+time_taken
        iterations <- as.numeric(bisection1[3])+iterations
      }
      strike <- append(strike, as.numeric(option_chain[i,"Strike"]))
      optionname <- append(optionname, paste(option_chain[i,"Strike"], "-", option_chain[i,"type"],"Expiring on:",
                                             option_chain[i, "Expiry"]))

  }
  #length(option_chain$daystoexpiry)
  option_chain_df <- data.frame(option_chain$daystoexpiry,optionname,iv,strike)
  #option_chain_df <- data.frame(optionname,iv,strike)
  names(option_chain_df) <- c("daystoexpiry", "variable", "ImpliedVolatility", "Strike")
  time_taken <- time_taken/as.numeric(colSums(!is.na(option_chain_df))[3])
  iterations <- iterations/as.numeric(colSums(!is.na(option_chain_df))[3])
  return(list(option_chain_df, time_taken, iterations))
}

options.data <- {}
# AMZN price
iv.bisection <- bisectionsam(2133.91)
options.data <- iv.bisection[1]
options.data <- data.frame(options.data)
options.data <- options.data[complete.cases(options.data$ImpliedVolatility),]
head(options.data)
tail(options.data)

#options.data <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/plotspy2.csv")

 #plot(options.data$Strike, options.data$ImpliedVolatility, main='Implied Volatility for 2 months SPY call option',
     #xlab='Strike', ylab = 'Implied Volatility', col = "purple", type = 'p', lwd = 2)
#lines(options.data$Strike, options.data$ImpliedVolatility, col = "green", type = 'p')
#legend("topright", c("April Maturity", "June Maturity"), col = c( "purple", "green"),
       #pch = c(1.5,1.5,1.5), ncol = 1, cex = 0.6)

# write.csv(options.data, "/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/optiondataspy.csv")

# 3D plot
library(knitr)
library(rgl)
knit_hooks$set(webgl = hook_webgl)

plot3d(x=options.data$daystoexpiry,
       y=options.data$ImpliedVolatility,
       z=options.data$Strike,
       col = rainbow(1000))

# 2.7
# Newton Method
bsvega <- function(S,r,sig,K,t,type){
  d1 <- (log(S/K)+(r+0.5*sig^2)*t)/(sig*sqrt(t))
  return(S*sqrt(t)*dnorm(d1))
}
findvol <- function(option_price,S,r,K,t,type){
  maxIterations <- 100
  precision <- 1.0e-10
  # initial value
  sig <- 0.2
  for(i in 1:maxIterations){
    value <- bs(s=S,k=K,r=r,sig=sig,tau=t,type=type)
    vega <- bsvega(S,r,sig,K,t,type)
    if(vega==0){
      sig <- 1e-300
      return(sig)
    }
    diff <- option_price-value
    ifelse(abs(diff)<precision, return (sig),
           sig <- sig+diff/vega)
  }
  return(sig)
}
findvol(253.400, 2133.91, 0.0159, 1905, 0.5190, "call")
# [1] 0.1311012


# newtonsam <- function(x){
#   option_chain <- option_chainamzn1call
#   #head(option_chain$daystoexpiry)
#   #daystoexpiry
#   libor <- 0.0159
#   iv <- {}
#   original_iv <- {}
#   optionname <- {}
#   strike <- {}
#   daystoexpiry <- {}
#   for(i in 1){
#     print(i)
#     newton1 <- findvol(
#       S <- x,
#       K <- as.numeric(option_chain[i,"Strike"]),
#       t <- as.numeric(option_chain[i, "daystoexpiry"])/252,
#       r <- libor,
#       type <- ifelse((option_chain[i,"type"]=="call"),"call","put"),
#       option_price <- as.numeric(option_chain[i,"premium"]))
#   
#     iv <- append(iv, as.numeric(newton1))
#     strike <- append(strike, as.numeric(option_chain[i,"Strike"]))
#     optionname <- append(optionname, paste(option_chain[i,"Strike"], "-", option_chain[i,"type"],"Expiring on:",
#                                            option_chain[i, "Expiry"]))
#     
#   }
#   #length(option_chain$daystoexpiry)
#   option_chain_df <- data.frame(option_chain$daystoexpiry,optionname,iv,strike)
#   #option_chain_df <- data.frame(optionname,iv,strike)
#   names(option_chain_df) <- c("daystoexpiry", "variable", "ImpliedVolatility", "Strike")
#   return(option_chain_df)
# }
# 
# options.data <- {}
# iv.newton <- newtonsam(334.68)
# options.data <- iv.newton
# options.data <- data.frame(options.data)
# options.data <- options.data[complete.cases(options.data$ImpliedVolatility),]
# head(options.data)
# tail(options.data)

# 2.9
# put-call parity
# finding the put price
pcp <- function(c, s, k, r, t){
  p <- c-s+k*exp(-r*t)
  return(p)
}
pcp(13.725, 2133.91, 2700, 0.0159, 0.5190)
# [1] 557.626

# 2.11
# Greeks
# Analytic Delta
delta.an <- function(s,k,r,v,t,type){
  d1 <- (log(s/k)+(r+((v^2)/2)*t))/(v*sqrt(t))
  if(type == "call"){
    delta <- pnorm(d1)
  }else{
    delta <- pnorm(-d1)
  }
  return(delta)
}
delta.an(2133.91, 2700,0.0159,0.2307,0.5190,"call")
# [1] 0.1080491

# Numerical Delta
delta.num <- function(s,k,r,v,t,type){
  p1 <- bs(s,k,r,v,t,type)
  A <- s*0.1
  s1 <- s+A
  p2 <- bs(s1,k,r,v,t,type)
  delta <- 1/A*(p2-p1)
  return(delta)
}
delta.num(2133.91,2700,0.0159,0.2307,0.5190,"call")
# [1] 0.1639589

# Analytic Gamma
gamma.an <- function(s,k,r,v,t,type){
  d1 <- (log(s/k)+(r+((v^2/2)*t))/(v*sqrt(t)))
  gamma <- dnorm(d1)/(v*s*sqrt(t))
  return(gamma)
}
gamma.an(2133.91,2700,0.0159,0.2307,0.5190,"call")
# [1] 0.001123075

# Numerical Gamma
gamma.num <- function(s,k,r,v,t,type){
  p1 <- bs(s,k,r,v,t,type)
  A <- s*0.1
  s1 <- s+A
  s2 <- s-A
  p2 <- bs(s1,k,r,v,t,type)
  p3 <- bs(s2,k,r,v,t,type)
  gamma <- 1/(A^2)*(p2-(2*p1)+p3)
  return(gamma)
}
gamma.num(2133.91,2700,0.0159,0.2307,0.5190,"call")
# [1] 0.000494686

# Analytic Vega
vega.an <- function(s,k,r,v,t,type){
  d1 <- (log(s/k)+(r+((v^2/2)*t))/(v*sqrt(t)))
  vega <- s*dnorm(d1)*sqrt(t)
  return(vega)
}
vega.an(2133.91,2700,0.0159,0.2307,0.5190,"call")
# [1] 0.000494686

# Numerical Vega
vega.num <- function(s,k,r,v,t,type){
  p1 <- bs(s,k,r,v,t,type)
  A <- v*0.01
  v1 <- v+A
  p2 <- bs(s,k,r,v1,t,type)
  vega <- 1/A*(p2-p1)
  return(vega)
}
vega.num(2133.91,2700,0.0159,0.2307,0.5190,"call")
# [1] 271.7922

# 2.12
# bssam <- function(x){
#   option_chain1 <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/amzn12.csv")
#   option_chain2 <- read.csv("/Users/Apple/Desktop/MS FE/Semester 2/621- Computational Finance/Assignment 1/day2 1212.csv")
#   libor <- 0.0159
#   iv <- {}
#   original_iv <- {}
#   optionname <- {}
#   strike <- {}
#   daystoexpiry <- {}
#   #head(option_chain$daystoexpiry)
#   #daystoexpiry
#   
#   for(i in 1:nrow(option_chain2)){
#     #print(i)
#     bs1 <- bs(
#       s <- x,
#       k <- as.numeric(option_chain2[i,"Strike"]),
#       tau <- as.numeric(option_chain2[i, "daystoexpiry"])/252,
#       r <- libor,
#       type <- ifelse((option_chain2[i,"type"]=="call"),"call","put"),
#       sig <- as.numeric(option_chain1[i, "ImpliedVolatility"]))
#     
#     iv <- append(iv, as.numeric(bs1))
#     
#     strike <- append(strike, as.numeric(option_chain1[i,"Strike"]))
#     optionname <- append(optionname, paste(option_chain1[i,"variable"]))
#     
#   }
#   #length(option_chain$daystoexpiry)
#   option_chain_df <- data.frame(option_chain2$daystoexpiry,optionname,iv,strike)
#   #option_chain_df <- data.frame(optionname,iv,strike)
#   names(option_chain_df) <- c("daystoexpiry", "variable", "optionprice", "Strike")
#   
#   return(option_chain_df)
# }
# 
# options.data <- {}
# # AMZN price
# iv.bs <- bssam(2150.80)
# options.data <- iv.bs
# options.data <- data.frame(options.data)
# options.data <- options.data[complete.cases(options.data$optionprice),]
# head(options.data)
# tail(options.data)


# Part 3
# 3.1
# define the f function
f <- function(x){
  if(x==0){
    return(1)
  }else{
    y <- sin(x)
    return(y/x)
  }
}

# Trapezoidal approximation
trpz <- function(u,d,n){
  sumud <- sum(f(u),f(d))
  t <- (d-u)/n
  for(m in 1:n-1){
    c <- u+(m*t)
    sumud <- sum(sumud, 2*f(c))
  }
  g <- sumud*t/2
  return(g)
}
trpz(-1000000, 1000000, 4000000)
# [1] 3.141591

# Simpson approximation
smps <- function(u,d,n){
  sumud1 <- 0
  t <- (d-u)/n
  for(m in 1:n-1){
    c1 <- u+ (m*t)
    c2 <- u+ ((m+1)*t)
    c3 <- (c1+c2)/2
    h <- sum(sumud1,f(c1),(4*f(c3)+f(c2)))
    sumud1 <- sum(h)
  }
  k <- sumud1*(t/6)
  return(k)
}
smps(-1000000, 1000000, 4000000)
# [1] 3.141591

# 3.2
# increasing a
ferror1 <- function(n){
  errortrpz <- c()
  errorsmps <- c()
  for (i in c(30,300,3000,30000)){
    errortrpz <- c(errortrpz, (trpz(-i, i, n)-pi))
    errorsmps <- c(errorsmps, (smps(-i, i, n)-pi))
  }
  errordf <- data.frame(errortrpz, errorsmps)
  return(errordf)
}
ferror1(5000)
#     errortrpz     errorsmps
# 1 -0.0084746364 -0.0080795735
# 2 -0.0002305296  0.0001695228
# 3  0.0006580811  0.0006508952
# 4  6.2820452473 -2.0952367896

# increasing N
ferror2 <- function(i){
  errortrpz <- c()
  errorsmps <- c()
  for (n in c(10,100,1000,10000)){
    errortrpz <- c(errortrpz, (trpz(-i, i, n)-pi))
    errorsmps <- c(errorsmps, (smps(-i, i, n)-pi))
  }
  errordf <- data.frame(errortrpz, errorsmps)
  return(errordf)
}
ferror2(5000)
#     errortrpz     errorsmps
# 1 998.851394109  3.288246e+02
# 2  94.242211676  2.724154e+01
# 3   6.281303029 -2.094148e+00
# 4  -0.000254131 -6.181055e-05

# 3.3
# Trapezoidal steps
i <- 10000
for (m in 1:10000){
  diff <- abs(trpz(-i, i, m+1)-trpz(-i, i, m))
  if(diff<0.00001){
    print(m)
    break
  }
}
# 3432

# Simpson's steps
i <- 10000
for (m in 1:10000){
  diff <- abs(smps(-i, i, m+1)-smps(-i, i, m))
  if(diff<0.00001){
    print(m)
    break
  }
}
# 3324

# Problem 4
# b
# defining f1
f1 <- function(x,y){
  return(x*y)
}

# defining f2
f2 <- function(x,y){
  return(exp(x+y))
}

# Trepzoidal approximation
# function 1
trpz1 <- function(a,b,c,d,n,m){
  sum1 <- 0
  tx <- (b-a)/n
  ty <- (d-c)/m
  for(k in 0:n-1){
    cx <- sum(a,(k*tx))
    cxx <- sum(cx,tx)
    for(i in 0:m-1){
      cy <- sum(c,(i*ty))
      cyy <- sum(cy,ty)
      term1 <- sum(f1(cx,cy),f1(cx,cyy),f1(cxx,cy),f1(cxx,cyy))
      term2 <- 2*sum(f1((cx+cxx)/2,cy),f1((cx+cxx)/2,cyy),
                     f1(cx,(cy+cyy)/2),f1(cxx,(cy+cyy)/2))
      term3 <- 4*f1((cx+cxx)/2,(cy+cyy)/2)
      sum1 <- sum(sum1,term1,term2,term3)
    }
  }
  return(tx*ty*sum1/16)
}
trpz1(0,1,0,3,4000,6000)
# [1] 2.25

# function 2
trpz2 <- function(a,b,c,d,n,m){
  sum2 <- 0
  tx <- (b-a)/n
  ty <- (d-c)/m
  for(k in 0:n-1){
    cx <- sum(a,(k*tx))
    cxx <- sum(cx,tx)
    for(i in 0:m-1){
      cy <- sum(c,(i*ty))
      cyy <- sum(cy,ty)
      term1 <- sum(f2(cx,cy),f2(cx,cyy),f2(cxx,cy),f2(cxx,cyy))
      term2 <- 2*sum(f2((cx+cxx)/2,cy),f2((cx+cxx)/2,cyy),
                     f2(cx,(cy+cyy)/2),f2(cxx,(cy+cyy)/2))
      term3 <- 4*f2((cx+cxx)/2,(cy+cyy)/2)
      sum2 <- sum(sum2,term1,term2,term3)
    }
  }
  return(tx*ty*sum2/16)
}
trpz2(0,1,0,3,4000,6000)
# [1] 32.79996

# error
# increasing n and m
ferror1 <- function(i,j){
  errortrpz1 <- c()
  for (n in c(10,100,1000)){
    for (m in c(20,200,2000)){
      errortrpz1 <- c(errortrpz1, (trpz1(-i,i,-j,j,n,m)-2.25))
    }
  }
  return(errortrpz1)
}

# increasing i and j
ferror1 <- function(n,m){
  errortrpz1 <- c()
  for (i in c(10,100,1000)){
    for (j in c(20,200,2000)){
      errortrpz1 <- c(errortrpz1, (trpz1(-i,i,-j,j,n,m)-2.25))
    }
  }
  return(errortrpz1)
}


