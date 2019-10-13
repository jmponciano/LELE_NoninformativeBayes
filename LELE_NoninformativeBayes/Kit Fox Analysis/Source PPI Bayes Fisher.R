#  This is a source program for PPI Bayes Frequentist program. 

library(dclone)


# ***********************************************************

   #     BAYES PPI COMPUTATIONS
   
# ***********************************************************


# We do the Bayesian PPI computations. These are quite easy by modifying the 
# basic model to predict the future states. 

PoissonRK.BayesPPImodel<-function(){
	
                     x[1] <- y1
                     for (j in 2:T) {
                         x[j] ~ dnorm(mu[j], prcx)
			             mu[j] <- a*(1- (min(exp(x[j-1]),10000)/K)) +  x[j-1]
                         N[j] ~ dpois(min(exp(x[j]),10000))
}

for ( j1 in (T+1):(T+Q.future)){
	
	x[j1] ~ dnorm(mu[j1], prcx)
    mu[j1] <- a*(1- (min(exp(x[j1-1]),10000)/K)) +  x[j1-1]
}

# Prior specification

prcx <- 1/(sigma.e^2)
a ~ dnorm(0,0.1)
K ~ dgamma(0.01,0.01)
sigma.e ~ dlnorm(0,0.1)

}

PoissonAB.BayesPPImodel<-function(){
	
                     x[1] <- y1
                     for (j in 2:T) {
                         x[j] ~ dnorm(mu[j], prcx)
                         mu[j] <- a - b*min(exp(x[j-1]),10000) +  x[j-1]
			            N[j] ~ dpois(min(exp(x[j]),10000))
}

for ( j1 in (T+1):(T+Q.future) ){
	
	x[j1] ~ dnorm(mu[j1], prcx)
    mu[j1] <-  a - b*min(exp(x[j1-1]),10000) +  x[j1-1]
}

# Prior specification

K <- a/b
prcx <- 1/(sigma.e^2)
a ~ dnorm(0,0.1)
b ~ dunif(0.001,1)
sigma.e ~ dlnorm(0,0.1)

}



# ****************************************************************************

#         DATA CLONING PPI COMPUTATION

# ****************************************************************************

# Now write the PPI program for DC analysis. We have to substitute the priors by the 
# asymptotic covariance matrix. Otherwise the program is same as the Bayesian one. We need to parameterize the model so that parameters are on unrestricted range and multivariate normality assumption is likely to hold. 

PoissonAB.DCmodel<-function(){
	for (i in 1:kk){
                     x[1,i] <- ykk[i]
                     for (j in 2:T) {
                         x[j,i] ~ dnorm(mu[j,i], prcx)
                         mu[j,i] <- a - b*min(exp(x[j-1,i]),10000) +  x[j-1,i]
			            
                         Nkk[j,i] ~ dpois(min(exp(x[j,i]),10000))
}
}

# Prior specification (identical to the Bayes PPI model)

K <- a/b
prcx <- 1/(sigma.e^2)
a ~ dlnorm(0,0.1)
b ~ dunif(0.001,1)
sigma.e ~ dlnorm(0,0.1)

# We use the transformation so that asymptotic Normal distribution is appropriate. We only monitor these transformed parameters for convergence. 

lsigma.e <- log(sigma.e)
la <- log(a)
# lb <- log(b)
lK <- log(K)
}


PoissonAB.DCPPImodel<-function(){
	
                     x[1] <- y1
                     for (j in 2:T) {
                         x[j] ~ dnorm(mu[j], prcx)
                         mu[j] <- a - b*min(exp(x[j-1]),10000) +  x[j-1]
			            N[j] ~ dpois(min(exp(x[j]),10000))
}

for ( j1 in (T+1):(T+Q.future) ){
	
	x[j1] ~ dnorm(mu[j1], prcx)
    mu[j1] <- a - b*min(exp(x[j1-1]),10000) +  x[j1-1]
}

a <- exp(parms[2])
b <- exp(parms[2])/exp(parms[1])
sigma.e <- exp(parms[3])

prcx <- 1/(sigma.e*sigma.e)
parms ~ dmnorm(mu.parms,prec.parms)

}


# Now we will do the same analysis using the (r,K) parameterization. 

PoissonRK.DCmodel<-function(){
	
	for (i in 1:kk){
		
                     x[1,i] <- ykk[i]
                     for (j in 2:T) {
                         x[j,i] ~ dnorm(mu[j,i], prcx)
			             mu[j,i] <- a*(1- (min(exp(x[j-1,i]),10000)/K)) +  x[j-1,i]
                         Nkk[j,i] ~ dpois(min(exp(x[j,i]),10000))
}
}

# Prior specification (identical to the BayesRK PPI model)

prcx <- 1/(sigma.e^2)
a ~ dlnorm(0,0.1)
K ~ dgamma(0.01,0.01)
sigma.e ~ dlnorm(0,0.1)

# We use the transformation so that asymptotic Normal distribution is appropriate. We only monitor these transformed parameters for convergence. We do not need to transform 'a' because in this model parameterization, it is not constrained to be positive. 

lsigma.e <- log(sigma.e)
lK <- log(K)
la <- log(a)
}


PoissonRK.DCPPImodel<-function(){
	
                     x[1] <- y1
                     for (j in 2:T) {
                         x[j] ~ dnorm(mu[j], prcx)
			             mu[j] <- a*(1- (min(exp(x[j-1]),10000)/K)) +  x[j-1]
                         N[j] ~ dpois(min(exp(x[j]),10000))
}

for ( j1 in (T+1):(T+Q.future)){
	
	x[j1] ~ dnorm(mu[j1], prcx)
    mu[j1] <- a*(1- (min(exp(x[j1-1]),10000)/K)) +  x[j1-1]
}

a <- exp(parms[2])
K <- exp(parms[1])
sigma.e <- exp(parms[3])
prcx <- 1/(sigma.e*sigma.e)
parms ~ dmnorm(mu.parms,prec.parms)
}


# *********************************************************************************************
