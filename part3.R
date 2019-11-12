#value initialization
r <- 0.0175
sig <- 0.2
N <- 10^(4.5)
S <- 100
k_seq <- 1:5
Z <- rnorm(N)

#part 1 - True calculation
F_true <- S*exp(r*(k_seq))

#part 2 - MC simulation
f_sim <- function(k) {
  R_N <- (r-(sig^2)/2)*k + sig*sqrt(k)*Z
  S_k <- S*exp(R_N)
  return(mean(S_k))
}
F_sim <- sapply(k_seq, f_sim)

#boxplot
boxplot(F_sim)
summary(boxplot(F_sim))


plot(F_true, type = "o", col = "blue", lty = 1, ylab = "Prices", xlab = "Years")
points(F_sim, col = "red", pch = "*")
lines(F_sim, col = "red", lty = 2)
legend(1, 109, legend=c("True Value", "MC Simulated Value"),
       col=c("blue", "red"), pch = c("o","*"), lty=1:2, cex=0.8)

#part 3
#Asset price
spot <- 100
r <- 0.0175
k <- 1
futureVal <- spot*exp(r*k)
futureVal_disc <- spot*exp(-(r*k))
#Forward price
N <- 10^(8)
Z <- rnorm(N)
R_N <- (r-(sig^2)/2)*k + sig*sqrt(k)*Z
S_k <- spot*exp(R_N)
S_k_mean <- mean(S*exp(R_N))
#Future price of asset 1 year from now = 101.7654
#According to simulation, the price of the forward contract can go up or down (i.e greater than 101.7654 or less) - use strategy accordingly
#two sim values of forward: 101.6492, 101.9684

#VaR calculation????
var5 <- quantile(S_k, 0.05)
plot(density(S_k))
points(futureVal, col = "red", pch = "o")
var <- ecdf(S_k)(futureVal)

s0 <- 100
r <- 0.0175
k <- 1
sig <- 0.2

N <- 10^(4)
Z <- rnorm(N)
R_N <- (r-(sig^2)/2)*k + sig*sqrt(k)*Z
S_k <- s0*exp(R_N)
S_k_mean <- mean(s0*exp(R_N))

f0 <- s0 * exp(r*k) #correct
s1 <- s0 * exp(r*k) #correct

f1 <- S_k_mean

pl <- (s1-s0) + (f0-f1)
mean(pl)
var_pl <- mean(pl) - quantile(pl, 0.1)

#part 4
#Calculating forward price minus $0.25
k <- 1
N <- 10^(4)
Z <- rnorm(N)
R_N <- (r-(sig^2)/2)*k + sig*sqrt(k)*Z

S_k <- spot*exp(R_N)
S_k_mean <- mean(spot*exp(R_N))
S_k_meanReduc <- S_k_mean - 0.25

#meanReduc <- c(S_k_meanReduc2, S_k_meanReduc3, S_k_meanReduc)

#Asset price
futureVal <- spot*exp(r*k)  #stays same
futureVal_disc <- spot*exp(-(r*k))  #decreases

#Plotting
plot(S_k_meanReduc, col = "black", lty = 1, pch = "o", ylim = c(95, 105), ylab = "Price", xlab = "Years")
points(futureVal, col = "blue", pch = "o")
lines(futureVal, col = "blue", lty = 2)
points(futureVal_disc, col = "red", pch = "o")
lines(futureVal_disc, col = "red", lty = 2)

legend(0.6, 105, legend=c("Price of stock increases with time", "Forward contract market price reduced by $0.25", "Price of stock decreases with time"),
       col=c("blue", "black", "red"), pch = "o", cex=0.8)

#################xxxxxxxxxxxxxxxxx#######################
#Creating price path
GMB2 <- function(n) {
  S <- 100
  dt <- 1/252
  mu <- 0.0175
  sig <- 0.2
  for(i in 1:252) {
    dR <- rnorm(1,
                dt*(mu - 0.5*sig^2),
                sig*sqrt(dt) )
    S_dt <- S[i]*exp(dR)
    S <- c(S,S_dt)
  }
  return(S)
}
r <- 0.0175
sig <- 0.2
N <- 10^(2)
S <- 100
k_seq <- 1:252
Z <- rnorm(N)




f_sim2 <- function(k) {
  R_N <- (r-(sig^2)/2)*(k/252) + sig*sqrt(k/252)*Z
  S_k <- S*exp(R_N)
  return(S_k)
}
F_sim2 <- sapply(k_seq, f_sim2)
plot(F_sim2)
#################xxxxxxxxxxxxxxxxx#######################



