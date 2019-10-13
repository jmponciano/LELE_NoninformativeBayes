# Source file for studying the non-invariance of the estimators and predictors for occupancy studies. 

# Data generation function


OccData.fn = function(N,K,p,psi){
	
	Y = matrix(0,N,K)
	O = Y
	for (i in 1:N){
		Y[i] = rbinom(1,1,psi)
		for (j in 1:K){
			O[i,j] = rbinom(1,1,p)*Y[i]
		}
	}
	return(O)
}


# Model files for different parameterization

library(dclone)

Occ_p.model = function(){
	
	for (i in 1:N){
		Y[i]~dbin(psi,1)
		for (j in 1:K){
			O[i,j]~dbin(p,Y[i])
		}
	}

# Occupancy probability
	
num1 <- pow((1-p),K)*psi
den1 <- num1 + (1-psi)
OccProb <- num1/den1

# This part is only for the American Toad analysis. We compute the occupancy rate. 

OccRate <- (10 + (17*OccProb))/27

	# Non-informative priors
	p ~ dunif(0,1)
	psi ~ dunif(0,1)

}


Occ_Lodds.model = function(){
	
	for (i in 1:N){
		Y[i]~dbin(psi,1)
		for (j in 1:K){
			O[i,j]~dbin(p,Y[i])
		}
	}

# Occupancy probability
	
num1 <- pow((1-p),K)*psi
den1 <- num1 + (1-psi)
OccProb <- num1/den1

# This part is only for the American Toad analysis. We compute the occupancy rate. 

OccRate <- (10 + (17*OccProb))/27


	p <- exp(lp)/(1+exp(lp))
	psi <- exp(lpsi)/(1+exp(lpsi))
	
	lp ~ dnorm(0,0.01)
	lpsi ~ dnorm(0,0.01)
}


# Simulation function


SimOcc.fn = function(N,K,p,psi,S){
	
out = matrix(0,S,7)
	
for (s in 1:S){
	
OccData = OccData.fn(N,K,p,psi)
Occ.data = list(O=dcdim(data.matrix(OccData)),N=N,K=K)

Occ.fit = jags.fit(Occ.data,c("p","psi"),Occ_p.model,inits = list(Y=as.vector(apply(OccData,1,max))))
Occ.Loddsfit = jags.fit(Occ.data,c("p","psi"),Occ_Lodds.model,inits = list(Y=as.vector(apply(OccData,1,max))))

tmp1 = summary(Occ.fit)
tmp2 = summary(Occ.Loddsfit)

# Probability of presence when observed to be unoccupied

numT = ((1-p)^K)*psi
denT = (1-psi) + numT
OccProbT = numT/denT

p1 = tmp1[[1]][1,1]
psi1 = tmp1[[1]][2,1]

p2 = tmp2[[1]][1,1]
psi2 = tmp2[[1]][2,1]

num1 = ((1-p1)^K)*psi1
den1 = num1 + (1-psi1)

num2 = ((1-p2)^K)*psi2
den2 = num2 + (1-psi2)

OccProb1 = num1/den1
OccProb2 = num2/den2

out[s,] = c(OccProbT,OccProb1,OccProb2,p1,p2,psi1,psi2)
}

return(out)

}


