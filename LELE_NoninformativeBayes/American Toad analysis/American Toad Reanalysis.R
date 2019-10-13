# Reanalysis of American Toad data from MacKenzie et al. (2002, Ecology). We take a subset with 3 visits and fit only the constant probability model. This code will lead to results in Table 3 and Figure 3.



# Read the data file

AmToad = read.csv("AmericanToad.csv")

AmToad = AmToad[-27,]

tmp = apply(AmToad,1,max)

sum(tmp)/length(tmp)   # Naive estimate of occupancy


# Now let us fit the occupancy model using two different parameterizations and see what happens.

source("Source functions for occupancy model non-variance.R")

OccData = AmToad
N = nrow(OccData)
K = ncol(OccData)
Occ.data = list(O=dcdim(data.matrix(OccData)),N=N,K=K)

Occ.fit = jags.fit(Occ.data,c("p","psi","OccProb","OccRate"),Occ_p.model,inits = list(Y=as.vector(apply(OccData,1,max))), n.adapt=10000, n.update = 30000, n.iter=30000)

Occ.Loddsfit = jags.fit(Occ.data,c("p","psi","OccProb","OccRate"),Occ_Lodds.model,inits = list(Y=as.vector(apply(OccData,1,max))),n.adapt=10000, n.update = 30000, n.iter=30000)

# Tables are based on the following results (Due to Monte Carlo variation the numbers will be slightly different each time the program is run)

tmp1 = summary(Occ.fit)
tmp2 = summary(Occ.Loddsfit)

# Values reported in the table are obtained from the following output. The actual numbers for each run of the program will be slightly different because of the Monte Carlo variation. 

tmp1
tmp2

# Plot the posterior distributions for the two occupancy rates 

try1 = c(Occ.fit[[1]][,2],Occ.fit[[2]][,2],Occ.fit[[3]][,2])
try2 = c(Occ.Loddsfit[[1]][,2],Occ.Loddsfit[[2]][,2],Occ.Loddsfit[[3]][,2])

plot(density(try2),ylim = c(0,8),xlim=c(0,1), xlab="Occupancy rate",ylab="Belief",main="",lty=2)
par(new=TRUE)
plot(density(try1),ylim = c(0,8),xlim=c(0,1), xlab="Occupancy rate",ylab="Belief",main="",lty=1)
title("Posterior distributions for occupancy rate")
legend(0.01,7,c("Probability","Logit"),lty=c(1,2))

# Let us do the likelihood based analysis (estimation and prediction) to see if this is invariant to the choice of the parameterization.

Occ.data = list(O=data.matrix(OccData),N=N,K=K)
inits = list(Y=as.vector(apply(OccData,1,max)))

ifun = function(model,n.clones){
	list(Y=dclone(as.vector(apply(OccData,1,max)),n.clones))
}


Pfit.DC = dc.fit(Occ.data,c("p","psi"),Occ_p.model,inits = inits, n.adapt=10000, n.update = 30000, n.iter=30000,n.clones = c(1,20),unchanged="K",multiply="N",initsfun=ifun)
summary(Pfit.DC)

Loddsfit.DC = dc.fit(Occ.data,c("p","psi"),Occ_Lodds.model,inits = inits, n.adapt=10000, n.update = 30000, n.iter=30000,n.clones = c(1,5),unchanged="K",multiply="N",initsfun=ifun)
summary(Loddsfit.DC)

