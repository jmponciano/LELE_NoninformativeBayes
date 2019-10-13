# R Code for Figure 1 and Jeffrey's prior

library(boot)
B = 100000
p1 = runif(B)
psi1 = rnorm(B,0,10)
# Effect of transformation on the prior distribution
par(mfrow=c(2,2))
plot(density(p1), main="Uniform prior on probability scale",xlim=c(0,1),xlab="")
plot(density(log(p1/(1-p1))), main="Implied prior on Logit scale",xlab="")
plot(density(psi1), main="Noninformative prior on Logit scale",,xlab="")
plot(density(inv.logit(psi1)), main="Implied prior on probability scale",xlim=c(0,1),xlab="")

par(mfrow=c(1,1))
p2 =rbeta(100000,0.5,0.5)
plot(density(p2,from=0,to=1), main="Jeffrey's prior on probability of success",xlim=c(0,1),xlab="")
