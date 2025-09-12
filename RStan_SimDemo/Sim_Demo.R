##############################################################################
###--- R code for simulating and fitting log-RT data with hierarchical &---###
###--- 					non-hierarchical Models							---###
##############################################################################
### Code: Udo Boehm
### Email: u.bohm@rug.nl

### Sample means and sds of log-transformed RTs of simulated pps from distributions fitted to real lexical decision data
### Simulate pps, sample RTs, fit hier. model to the data ###

library(rstan)
#load the estimated group-level distribution parameters estimated from the Wagenmakers et al., 2008 data
load("/data/p264576/LexDec/LexDatParams8.RData")


### Simulation parameters
params <- list()
params$mu1.m <- mu.m
params$delta <- c(0)
params$sigma.m <- sigma.m
params$sigma.s <- sigma.s
params$mu2.m <- params$mu1.m+params$delta*params$sigma.m
params$k <- 5 #trials per pp
params$N <- 30 #pps per group
params$NSim <- 20 #number of repetitions of the simulation to run - running more than 10 will take a while

### SIMULATION ###

### Simulate data, fit hierarchical model to the data, store summary stats of the posterior samples
for(iSim in 1:params$NSim){
	### Generate similuated data
	
	# Sample pp means and sds
	params$mu2.m <- params$mu1.m + idelta*params$sigma.m
	params$M1 <- rnorm(iN, params$mu1.m, params$sigma.m)
	params$M2 <- rnorm(iN, params$mu2.m, params$sigma.m)

	params$SIGMA1 <- abs(rnorm(iN, 0, params$sigma.s))
	params$SIGMA2 <- abs(rnorm(iN, 0, params$sigma.s))
				
	# Sample log-transformed RT data
	simdat <- NULL
				
	for(isubj in 1:iN){
		x <- rnorm(ik, params$M1[isubj], params$SIGMA1[isubj])
		y <- rnorm(ik, params$M2[isubj], params$SIGMA2[isubj])
		
		simdat <- rbind(simdat, data.frame(subj = isubj, x = x, y = y, m1=params$M1[isubj], m2=params$M2[isubj], s1=params$SIGMA1[isubj], s2=params$SIGMA2[isubj]))
	}
	
	# Get the data in the right format
	rt <- array(NA, c(2,ik,iN))
	rt[1,,] <- as.matrix(unstack(simdat, x~subj))
	rt[2,,] <- as.matrix(unstack(simdat, y~subj))
	k <- ik
	n <- iN
	dat <- list(rt=rt, k=k, n=n)#"gr", "subj",
	
	### Fit the hierarchical model to the simulated data
	
	# Parameters for which we want posterior samples
	mod.params <- c("m","s","sigma_m","delta","mu_m","sigma_s")

	# Initialise parameters
	inits <- function(){
		list("mu_m1" = runif(1,5,7), "delta" = rnorm(1,0,1), "sigma_m" = runif(1,.001,4.999), "sigma_s" = runif(1,.001,9.999), "m"=cbind(runif(n,5.5,6.5),runif(n,5.5,6.5)),"s"=matrix(runif(2*n,.001,.999),ncol=2,byrow=F))
	}

	# Select hierarchical Bayesian model
	model.file <- "full model.stan"
	
	# Obtain posterior samples
	ft <- stan(file = model.file, data = dat, pars=mod.params, chains=3, warmup=2000, thin=4, iter=2e4, init=inits, cores=3)

	# Select summary statistic of the posterior samples
	hier.out <- list()
	hier.out$mean <- list(lapply(extract(ft, pars=c("sigma_m", "delta","sigma_s")), mean), "mu_m"=apply(extract(ft, "mu_m")$mu_m,2,mean), "m"=apply(extract(ft, pars=c("m"))$m, c(2,3), mean), "s"=apply(extract(ft, pars=c("s"))$s, c(2,3), mean))
	hier.out$median <- list(lapply(extract(ft, pars=c("sigma_m", "delta","sigma_s")), median), "mu_m"=apply(extract(ft, "mu_m")$mu_m,2,median), "m"=apply(extract(ft, pars=c("m"))$m, c(2,3), median), "s"=apply(extract(ft, pars=c("s"))$s, c(2,3), median))
	# Store posterior samples for delta
	delta.hier.post <- extract(ft, pars="delta")[[1]]
	
	# Store the simulated data from the generating model
	od <- cbind(tapply(simdat$m1, simdat$subj, unique), tapply(simdat$m2, simdat$subj, unique))

	### Fit the non-hierarchical model to the simulated data
	
	# Get the data in the right format
	m.rt <- array(NA, c(iN,2))
	m.rt[,1] <- tapply(simdat$x, simdat$subj, mean)
	m.rt[,2] <- tapply(simdat$y, simdat$subj, mean)
	n <- iN
	dat <- list(m_rt=m.rt, n=n)

	mod.params <- c("mu_m","delta","sigma_m")

	# Initialise parameters
	inits <- function(){
		list("mu_m1" = runif(1,5,7), "delta" = rnorm(1,0,1), "sigma_m" = runif(1,.001,4.999))
	}

	model.file <- "non-hier model.stan"

	ft <- stan(file=model.file, data=dat, pars=mod.params, chains=3, warmup=500, thin=1, iter=5000, init=inits)

	delta.parthier.post <- extract(ft, pars="delta")$delta
	save("params","simdat","od","hier.out","delta.hier.post","delta.parthier.post",file="/out/S5-demo.RData")
}
