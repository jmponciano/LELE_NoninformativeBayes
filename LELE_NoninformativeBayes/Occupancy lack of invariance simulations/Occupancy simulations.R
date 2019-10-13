#  Simulations for showing non-invariance of the occupancy model inference to different non-informative priors

source("Source functions for occupancy model non-variance.R")

psi = 0.3   # Occupancy probability
p =0.3      # detection probability

N = 30   # number of sites
K = 2      # number of visit
S = 100
SimOcc.out = SimOcc.fn(N=N,K=K,p=p,psi=psi,S=S)

apply(SimOcc.out,2,summary)
summary(SimOcc.out[,3]/SimOcc.out[,2])
summary(SimOcc.out[,3]/SimOcc.out[,1])
summary(SimOcc.out[,2]/SimOcc.out[,1])


# Plot for Jeffrey's prior on 'p'

x = seq(0,1,0.001)
hx = dbeta(x,0.5,0.5)
plot(x,hx,type="l",main="Jeffreys prior for p",xlab="p",ylab="Density")

# Transformation of priors

lp = rnorm(1000000,0,0.01)
plot(density(lp))
p.tran = exp(lp)/(1+exp(lp))
plot(density(p.tran))

p1 = runif(1000000,0,1)
lp.tran = log(p1/(1-p1))
plot(density(lp.tran))
prec.lp = var(lp.tran)

