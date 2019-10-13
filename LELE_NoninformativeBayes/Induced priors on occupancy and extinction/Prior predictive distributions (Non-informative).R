# Induced priors on parameters of interest (functions of parameters): Prior predictive distributions
# This program will produce figures 5 and 6. 

# Occupancy probability: Priors on the natural parameters affect the prior predictive distribution.

# Uniform priors on the original scale. This is also affected by the number of visits.

k = 2
B = 10000
p1 = runif(B,0,1)
psi1 = runif(B,0,1)
parms.mat = cbind(p1,psi1)
occupancy.pred = (((1-p1)^k)*psi1)/(((1-p1)^k)*psi1 + (1-psi1))

plot(density(occupancy.pred,from=0,to=1),ylim=c(0,4),main="",xlab="")
# hist(occupancy.pred)
summary(occupancy.pred)

# Non-informative prior on the Logit scale
k = 2
B = 10000
library(boot)
Y1 = rnorm(B,0,10)
Y2 = rnorm(B,0,10)
p1 = inv.logit(Y1)
psi1 = inv.logit(Y2)
parms.mat = cbind(p1,psi1)
occupancy.pred = (((1-p1)^k)*psi1)/(((1-p1)^k)*psi1 + (1-psi1))
par(new=T)
plot(density(occupancy.pred,from=0,to=1),ylim=c(0,4),col="red",main="Prior predictive occupancy",xlab="Probability")
# hist(occupancy.pred)
summary(occupancy.pred)






# PVA and induced priors (Dennis et al 1991, equation 15)

B = 10000
mu1 = rlnorm(B,0,10)
sigma1 = rlnorm(B,0,10)
xd = 10
PseudoExtinct = exp(-2*mu1*xd/(sigma1^2))
plot(density(PseudoExtinct,from=0,to=1),ylim=c(0,10),col="blue",main="Prior Extinction probability",xlab="",ylab="")
# hist(PseudoExtinct)
summary(PseudoExtinct)

B = 10000
mu1 = runif(B,0,10)
sigma1 = runif(B,0,10)
xd = 10
PseudoExtinct = exp(-2*mu1*xd/(sigma1^2))
par(new=T)
plot(density(PseudoExtinct,from=0,to=1),ylim=c(0,10),col="red",main="",xlab="Probability of quasi-extinction",ylab="Density")

summary(PseudoExtinct)


