 # This program is for the analysis of Kit fox data to conduct the PPI under Bayesian and Frequentist methods. This program will give the results used in Table 1 and will draw Figure 2 comparing PPI under likelihood inference and non-informative Bayesian inference. You can choose your own non-informative priors in the non-informative Bayesian inference to see the effect of the lack of invariance.
 
 # May need to change the library to access the source program.
 
 source("Source PPI Bayes Fisher.R")
 
 
# This is the Kit Fox data along with SEs.
kitfox = c(117,220,161,164,135,166,131,117,46,88,190,363,133)
ses = c(11.5,16.6,15.9,14.3,13.4,14.7,13.9,13.0,9.2,9.2,15.9,23.5,11)



# Plot of SEs vs abundance to see if Poisson error is okay or not.

plot(kitfox,ses^2,ylab="Sample variance",xlab="Abundance")
title("Poisson sampling error")
lines(abline(0,1))

# **********************************************************************

# Bayesian PPI computation

# **********************************************************************

dat1 <- list(N=kitfox,y1=log(kitfox[1]),T=length(kitfox), Q.future = 10)

PoissonAB.BayesPPI <- jags.fit(dat1, c("a","K","sigma.e","x[2:23]"),PoissonAB.BayesPPImodel, n.update=5000, n.chains=5,n.iter = 10000)

# These are simple diagnostic checks for the MCMC convergence etc. and also for looking at the estimates. 

BayesABPPI.summary= summary(PoissonAB.BayesPPI)
plot(PoissonAB.BayesPPI)
gelman.diag(PoissonAB.BayesPPI,multivariate=F)
BayesABPPI.summary

PoissonRK.BayesPPI <- jags.fit(dat1, c("a","K","sigma.e","x[2:23]"),PoissonRK.BayesPPImodel, n.update=5000, n.chains=5,n.iter = 10000)
BayesRKPPI.summary= summary(PoissonRK.BayesPPI)
plot(PoissonRK.BayesPPI)
gelman.diag(PoissonRK.BayesPPI,multivariate=F)
BayesRKPPI.summary

StatesAB.est = exp(BayesABPPI.summary$quantiles[-(1:3),])
StatesRK.est = exp(BayesRKPPI.summary$quantiles[-(1:3),])


#   *********************************************************************************************

#                    Frequentist approach using data cloning

#  **********************************************************************************************


N = kitfox
y = log(kitfox)

# It is important to use dcdim and data.matrix commands to make sure cloning is done appropriately.
# We will first do simple estimation to make sure the estimates are similar and so are the precision matrices. 


dat = list(Nkk=dcdim(data.matrix(N)),ykk=y[1],T=length(N), kk=1)

PoissonAB.DC <- dc.fit(dat, c("la","lK","lsigma.e"),PoissonAB.DCmodel, n.clones = c(1,20), multiply = "kk", unchanged = "T", n.update=5000, n.chains=5,n.iter = 10000)

PoissonRK.DC <- dc.fit(dat, c("la","lK","lsigma.e"),PoissonRK.DCmodel, n.clones = c(1,20),multiply = "kk", unchanged = "T", n.update=5000, n.chains=5,n.iter = 10000)

# Do the diagnostic checks

plot(PoissonAB.DC)
plot(PoissonRK.DC)

gelman.diag(PoissonAB.DC,multivariate=F)
gelman.diag(PoissonRK.DC,multivariate=F)

mu.parmsAB = summary(PoissonAB.DC)[[1]][,1]
prec.parmsAB = solve(vcov(PoissonAB.DC))

mu.parmsRK = summary(PoissonRK.DC)[[1]][,1]
prec.parmsRK = solve(vcov(PoissonRK.DC))

mu.parmsAB-mu.parmsRK 
vcov(PoissonAB.DC)/vcov(PoissonRK.DC)


# MLE Parameter estimates under the two parameterizations

exp(mu.parmsAB)   # ML estimates under AB model
exp(mu.parmsRK)   # ML estimates under RK model

# Now we will compute the PPI under the two parameterizations


dat1 <- list(N=kitfox,y1=log(kitfox[1]),T=length(kitfox), Q.future = 10, mu.parms=mu.parmsAB, prec.parms = prec.parmsAB)
PoissonAB.DCPPI <- jags.fit(dat1, c("x[2:23]"),PoissonAB.DCPPImodel, n.update=5000, n.chains=5,n.iter = 10000)

dat1 <- list(N=kitfox,y1=log(kitfox[1]),T=length(kitfox), Q.future = 10, mu.parms=mu.parmsRK, prec.parms = prec.parmsRK)
PoissonRK.DCPPI <- jags.fit(dat1, c("x[2:23]"),PoissonRK.DCPPImodel, n.update=5000, n.chains=5,n.iter = 10000)




DC_AB_PPI.summary= summary(PoissonAB.DCPPI)
plot(PoissonAB.DCPPI)
gelman.diag(PoissonAB.DCPPI)
DCStatesAB.est = exp(DC_AB_PPI.summary$quantiles[,])

DC_RK_PPI.summary= summary(PoissonRK.DCPPI)
plot(PoissonRK.DCPPI)
gelman.diag(PoissonRK.DCPPI)
DCStatesRK.est = exp(DC_RK_PPI.summary$quantiles[,])

cbind(DCStatesAB.est[-(1:12),1],DCStatesRK.est[-(1:12),1])

# These are essentially identical to each other; hence parameterization invariant. 

cbind(DCStatesAB.est[-(1:12),1],DCStatesRK.est[-(1:12),1],StatesAB.est[-(1:12),1],StatesRK.est[-(1:12),1])


# Here is the final plot showing invariance of the likelihood and non-invariance of the Bayesian approach. We need to add legend to this.

plot(DCStatesRK.est[-(1:12),1],cex=0.1,type="l",ylim=c(0,100), lty=1,xlab="Time",ylab="PPI")
par(new=T)
plot(DCStatesAB.est[-(1:12),1],cex=0.1,type="l",ylim=c(0,100),lty=2,xlab="Time",ylab="PPI")
par(new=T)
plot(StatesAB.est[-(1:12),1],cex=0.1,type="l",ylim=c(0,100),lty=3,,xlab="Time",ylab="PPI")
par(new=T)
plot(StatesRK.est[-(1:12),1],cex=0.1,type="l",ylim=c(0,100),lty=4,xlab="Time",ylab="PPI")
title("Population Prediction Intervals under two parameterizations")
legend(7,95,c("BayesAB","BayesAK","LikeAB","LikeAK"), lty = c(3,4,2,1))


