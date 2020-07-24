########################################
#Gopher Tortoise - Fort Stewart 	       #
#Spatial mark-recapture analysis       #
#E. Hunter, July 2018        	       #
########################################

rm(list=ls())

ELIZABETHOFFICE = F
ELIZABETH = T

###################################
##########  SET WORKING DIRECTORY


if(ELIZABETHOFFICE) 
	rootDir <- "C:\\Users\\elizabethhunter\\Dropbox\\GSU_Research\\GopherTortoise\\Data\\ReceivedData\\FtStewart\\"
if(ELIZABETH) rootDir <- "~/Dropbox/GSU_Research/GopherTortoise/Data/ReceivedData/FtStewart/"

ScriptDir <- paste(rootDir,"Rcode",sep="")
DataDir <- paste(rootDir,"Data",sep="")
FiguresDir <- paste(rootDir,"Figures",sep="")
BUGSDir <- paste(rootDir,"BUGS",sep="")

setwd(DataDir)

getwd()

##################################
##########  LOAD PACKAGES

library(runjags)
library(MuMIn)
library(MASS)

###################################
##########   READ IN DATA

setwd(DataDir)

#Data: capture histories, covariates, etc.
source("GT_data_20200708.txt")


###############
# BUGS MODEL
###############

##START THE FILE FOR JAGS

BUGSfilename <- "GT_BUGS_v17.txt"
setwd(BUGSDir)

cat("
    
    model{
    
#############
# PRIORS
#############
      
      ###########
      # PROBABILITY OF CAPTURE IN SURVEYED AREA
      ###########
      
      p.male.eff ~ dunif(-3,3)    # effect of sex on mean capture probability
	  p.scl.dist.eff ~ dunif(-3,3)  #interactive effect of size and distance
	  p.effort.eff ~ dunif(-5,5)   #effect of effort
	  p.tavg.eff ~ dunif(-3,3)   #effect of temperature
	  p.bucket.eff ~ dunif(-5,5)  #effect of bucket sampling
	  p.dist.eff ~ dunif(-3,3)   #effect of distance
	  sigma ~ dunif(0, 4)  #effective strip half-width
	  sigma2 <- sigma*sigma

      
      
      #### mean capture probability
    for(s in 1:nsites){

	  p0[s] ~ dunif(0,1)                 # mean/intercept detection prob within the surveyed area
      logit.p0[s] <- log(p0[s]/(1-p0[s]))   # convert to logit
	  p.burn.eff[s] ~ dunif(-3,3)  #site-specific effect of years-since-burn on capture (interacts with distance)

	for(ind in 1:ninds[s]){  # loop through (data augmented) individuals
        for(period in first.period[s]:last.period[s]){ 
          for(session in 1:nsessions[period]){
		    dist[s,ind,period,session] <- exp(-1*pow(rd.dist[s,ind,period,session],2) / (2*sigma2)) 
		    logit(thisp0[s,ind,period,session]) <- logit.p0[s] + p.male.eff*is.male[s,ind] 
									+ p.scl.dist.eff*scl[s,ind]*dist[s,ind,period,session] 
									+ p.effort.eff*effort[s,period,session]
									+ p.bucket.eff*bucket[period]
									+ p.tavg.eff*tavg[period,session]
									+ p.burn.eff[s]*burn.loc[s,ind,period,session]*dist[s,ind,period,session] 
									+ p.dist.eff*dist[s,ind,period,session]
									
          }
        }
      }
	}
      

      ###########
      # PROBABILITY OF RECRUITMENT
      ###########
      
      # Dirichlet prior for entrance probabilities (following Royle's S-A formulation)
     

	for(s in 1:nsites){
       init.entranceprob[s] ~ dgamma(1,1)          #first entrance probability is fundamentally different from the others   
					# probability of entering the population at time 1 is estimated separately 
     	 beta[s,first.period[s]] <- init.entranceprob[s]
     	 gamma[s,first.period[s]] <- beta[s,first.period[s]]/sum(beta[s,first.period[s]:last.period[s]]) 
     	 for(period in (first.period[s]+1):last.period[s]){
      	  beta[s,period] ~ dgamma(interval[(period-1)],1)     #close to equal probability of entrance each occasion
     		  gamma[s,period] <- beta[s,period]/sum(beta[s,first.period[s]:last.period[s]])    # probability of entrance
      }

      # convert to conditional entrance probs
      
      cprob[s,first.period[s]] <- gamma[s,first.period[s]]
      for(period in (first.period[s]+1):last.period[s]){
        cprob[s,period]<-gamma[s,period]/(1-sum(gamma[s,first.period[s]:(period-1)]))
      }
	}#end site loop


	###########
      # PROBABILITY OF SURVIVING/AGING 
      ###########

      phi.male.eff ~ dunif(-3,3)    # effect of sex 
	  phi.scl.eff ~ dunif(-3,3)	# effect of size


	for(s in 1:nsites){
      	phi0[s] ~ dunif(0.1,1)                
      	phi0.logit[s] <- log(phi0[s]/(1-phi0[s])) 

		phi.burn.eff[s] ~ dunif(-3,3)	# site-specific effect of years-since-burn

		for(period in (first.period[s]+1):last.period[s]){
		for(ind in 1:ninds[s]){ 
			logit(phi[s,ind,period]) <- phi0.logit[s] + phi.male.eff*is.male[s,ind] 
								+ phi.scl.eff*scl[s,ind] + phi.burn.eff[s]*ysb[s,period]
	}
	}
	}
      
      
      ###########
      # PROBABILITY OF BEING REAL
      ###########
      
	for(s in 1:nsites){
      	psi[s] ~ dunif(0,1)   #data augmentation parameter
	}   
      
      
      #############
      # DEAL WITH MISSING DATA
      #############
      
	for(s in 1:nsites){
     		prob.male[s] ~ dunif(0,1)    #probability of an individual being male
      
      	for(ind in 1:ninds[s]){
      	 is.male[s,ind] ~ dbern(prob.male[s])
      	}

		scl.prec[s] ~ dgamma(0.1,0.1)
		scl.sd[s] <- pow((1/scl.prec[s]),0.5)
		for(ind in 1:ninds[s]){
			scl[s,ind] ~ dnorm(0, scl.prec[s])
		}

		#Distance from road (and location-specific burning data - missing data for 2019)
		rd.shape[s] ~ dunif(0.2, 3)
		rd.rate[s] ~ dunif(1, 10)
		mu.burn.loc[s] ~ dunif(0,5)
		burn.loc.prec[s] ~ dgamma(0.1, 0.1)
		for(ind in 1:ninds[s]){
			for(period in first.period[s]:last.period[s]){
			for(session in 1:nsessions[period]){
				rd.dist[s,ind,period,session] ~ dgamma(rd.shape[s], rd.rate[s])
				burn.loc[s,ind,period,session] ~ dnorm(mu.burn.loc[s], burn.loc.prec[s])
			}
			}
		}

		#Years-since-burn (missing data for 2019)
		mu.ysb[s] ~ dunif(0,5)
		ysb.prec[s] ~ dgamma(0.1, 0.1)
		for(period in first.period[s]:last.period[s]){
			ysb[s,period] ~ dnorm(mu.ysb[s], ysb.prec[s])
		}

	}


      
#############
# LIKELIHOOD
#############
      for(s in 1:nsites){
      for(ind in 1:ninds[s]){  # loop through augmented individuals
        thissex[s,ind] <- 2-is.male[s,ind]   # male is 1, female is 2
        z[s,ind] ~ dbern(psi[s])    # is this individual real or fake? (data augmentation)
        in.pop.now[s,ind,first.period[s]] ~ dbern(gamma[s,first.period[s]])   # initial entrance probability...  
        
        for(session in 1:nsessions[first.period[s]]){
		#p=probability of catching an individual at a given point given it is in the population (alive), in the study area, and real (first PP)
          p[s,ind,first.period[s],session] <- thisp0[s,ind,first.period[s],session] * z[s,ind] * in.pop.now[s,ind,first.period[s]]  
          y[s,ind,first.period[s],session] ~ dbern(p[s,ind,first.period[s],session])
        }
        
        recruitable[s,ind,first.period[s]] <- 1-in.pop.now[s,ind,first.period[s]]     # 1 if the indiv is not yet in the population 
        recruited[s,ind,first.period[s]] <- in.pop.now[s,ind,first.period[s]]
        
        for(period in (first.period[s]+1):last.period[s]){
          survival.part[s,ind,period] <- pow(phi[s,ind,period],interval[(period-1)]) * in.pop.now[s,ind,(period-1)]
          recruitment.part[s,ind,period] <- cprob[s,period]*recruitable[s,ind,(period-1)]
          expected.inpop[s,ind,period] <- survival.part[s,ind,period] + recruitment.part[s,ind,period]  # either it is still in the pop or it just entered
          in.pop.now[s,ind,period] ~ dbern(expected.inpop[s,ind,period])
          recruitable[s,ind,period] <- recruitable[s,ind,(period-1)] * (1-in.pop.now[s,ind,period])       # is it still not yet in the study population?
          recruited[s,ind,period] <- (1-in.pop.now[s,ind,(period-1)]) * in.pop.now[s,ind,period]   # was it recruited this year?
          
          for(session in 1:nsessions[period]){
		#probability of catching an individual at a given point given it is in the population, in the study area, and real (after first PP)
            p[s,ind,period,session] <- thisp0[s,ind,period,session] * z[s,ind] * in.pop.now[s,ind,period] 
            y[s,ind,period,session] ~ dbern(p[s,ind,period,session])                           
          }
          captured_this_period[s,ind,period] <- step(sum(y[s,ind,period,1:nsessions[period]])-1)  # was it captured at least once?
        }    
      }
      }
      
      
      #############
      # DERIVED TERMS
      #############
      
	#N and per capita recruitment
     for(s in 1:nsites){
      for(period in first.period[s]:last.period[s]){
        N[s,period] <- inprod(in.pop.now[s,1:ninds[s],period],z[s,1:ninds[s]])
        N.recruited[s,period] <- inprod(recruited[s,1:ninds[s],period],z[s,1:ninds[s]]) / N[s,period]
      }
	}

    
    }   ## end BUGS model
    
    ",file=BUGSfilename)






######################
# PREPARE DATA FOR BUGS
######################


data.for.bugs <- list(
	y = y,
	effort = effort,
	bucket = bucket,
	ninds = ninds,
	is.male = is.male,
	scl = scl,
	tavg = tavg,
	nsessions = nsessions,
	first.period = first.period,
	last.period = last.period,
	rd.dist = rd.dist,
	ysb = ysb,
	burn.loc = burn.loc,
	nsites = nsites,
	interval = interval
)


#Initial alive (1 if real, 0 if augmented)
#Initial in.pop.now (1 for all individuals)
max.inds <- max(ninds)
nprimaryG <- length(nsessions)
init.z <- array(NA, dim=c(nsites, max.inds))
init.in.pop.now <- array(NA, dim=c(nsites, max.inds, nprimaryG))
for(s in 1:nsites){
	init.z[s,] <- c(rep(1, times=ninds[s]/3), rep(0, times=ninds[s]-(ninds[s]/3)), rep(NA, times=max.inds-ninds[s]))
	for(i in 1:ninds[s]){
	for(p in first.period[s]:last.period[s]){
		init.in.pop.now[s,i,p] <- 1
}}}

initz.bugs<-function(){
  list(
    z=init.z,
    p0=runif(nsitesG,0.1,0.15),
    p.male.eff=runif(1,-0.1,0.1),
    p.effort.eff=runif(1,0.01,0.1),
    #p.scl.eff=runif(1,-0.1,0.1),
    p.scl.dist.eff=runif(1,-0.1,0.1),
    p.tavg.eff=runif(1,-0.1,0.1),
    p.bucket.eff=runif(1,-0.1,0.1),
    p.burn.eff=runif(nsitesG, -0.1, 0.1),
    p.dist.eff=runif(1, -0.1, 0.1),
    
    phi0=runif(nsitesG,0.9,0.99),
    #phi.year.prec=runif(nsitesG,0.01, 0.05),
    phi.male.eff=runif(1,-0.1,0.1),
    phi.scl.eff=runif(1,-0.1,0.1),
    phi.burn.eff=runif(nsitesG,-0.1,0.1),
    psi = runif(nsitesG,0.5,0.9),
    init.entranceprob = runif(nsitesG,3,4),

    sigma = runif(1, 1, 2),
    rd.shape = runif(nsitesG, 0.5, 2),
    rd.rate = runif(nsitesG, 3, 7),
    
    ##  initialize every indiv as being in the population?  
    in.pop.now = init.in.pop.now

    prob.male = c(0.5, 0.5, 0.5),  

    mu.ysb = c(1.62, 1.68, 1.69),
    ysb.prec = c(1.59, 1.18, 0.67),

    mu.burn.loc = c(1.62, 1.68, 1.69),
    burn.loc.prec = c(1.59, 1.18, 0.67)
  )
}


#Running JAGS
setwd(BUGSDir)
system.time(
  mod<-run.jags(
    model=BUGSfilename,
    monitor=c("z", "p0", "p.male.eff", "p.effort.eff", "p.bucket.eff", "p.burn.eff", "p.scl.dist.eff", "p.tavg.eff", "p.dist.eff","phi0", "phi.male.eff", "phi.scl.eff", "phi.burn.eff", "psi", "init.entranceprob", "sigma", "rd.shape", "rd.rate", "in.pop.now", "prob.male", "mu.ysb", "ysb.prec", "mu.burn.loc", "burn.loc.prec"),
    data=data.for.bugs,
    n.chains = 3,
    inits=initz.bugs,
    burnin = 20000,  
    sample=5000,  
    thin=5,
    method="parallel"
    #clearWD=FALSE
  )
)

#Examine posteriors
Mod.mcmc <- mod$mcmc
Mod.mcmc.list <- mcmc.list(Mod.mcmc)
heidel.diag(Mod.mcmc.list)
gelman.diag(Mod.mcmc.list)
plot(Mod.mcmc, ask=T)


#####################################
#ANALYZING JAGS MODEL OUTPUT (post hoc analyses)
DataDir <- paste(rootDir,"Data/JAGS_output_20200708",sep="")
DataDir2 <- paste(rootDir,"Data", sep="")

#Model was analyzed using a high powered computing cluster, so read in output files:
setwd(DataDir)
idx1 <- read.csv("output_20200708_chain1index.csv", header=FALSE)
chn1 <- read.csv("output_20200708_chain1chain1.csv", header=FALSE)
chn2 <- read.csv("output_20200708_chain2chain1.csv", header=FALSE)
chn3 <- read.csv("output_20200708_chain3chain1.csv", header=FALSE)

#For each row in the "idx", want to get the values from the 3 chains for the values between the two indices
out.list <- list()

for(i in 1:nrow(idx1)){
	var.id <- as.character(idx1[i,1])
	temp.idx <- as.numeric(idx1[i,2:3])
	temp.chn1 <- chn1[temp.idx[1]:temp.idx[2],2]
	temp.chn2 <- chn2[temp.idx[1]:temp.idx[2],2]
	temp.chn3 <- chn3[temp.idx[1]:temp.idx[2],2]
	var.out <- as.data.frame(cbind(temp.chn1, temp.chn2, temp.chn3))
	names(var.out) <- c(paste(var.id, 1, sep="-"), paste(var.id, 2, sep="-"), paste(var.id, 3, sep="-"))
	out.list[[i]] <- var.out
}



#Pull out p effects and get histograms
p.0 <- out.list[[1]]
p0.1 <- c(p.0[,1], p.0[,2], p.0[,3])
hist(p0.1)
quantile(p0.1, probs=c(0.025, 0.5, 0.975))

p.0 <- out.list[[2]]
p0.2 <- c(p.0[,1], p.0[,2], p.0[,3])
hist(p0.2)
quantile(p0.2, probs=c(0.025, 0.5, 0.975))

p.0 <- out.list[[3]]
p0.3 <- c(p.0[,1], p.0[,2], p.0[,3])
hist(p0.3)
quantile(p0.3, probs=c(0.025, 0.5, 0.975))

p.male <- out.list[[4]]
p.male <- c(p.male[,1], p.male[,2], p.male[,3])
hist(p.male)
quantile(p.male, probs=c(0.025, 0.5, 0.975))

p.effort <- out.list[[5]]
p.effort <- c(p.effort[,1], p.effort[,2], p.effort[,3])
hist(p.effort)
quantile(p.effort, probs=c(0.025, 0.5, 0.975))

p.bucket <- out.list[[6]]
p.bucket <- c(p.bucket[,1], p.bucket[,2], p.bucket[,3])
hist(p.bucket)
quantile(p.bucket, probs=c(0.025, 0.5, 0.975))

p.scl.dist <- out.list[[10]]
p.scl.dist <- c(p.scl.dist[,1], p.scl.dist[,2], p.scl.dist[,3])
hist(p.scl.dist)
quantile(p.scl.dist, probs=c(0.025, 0.5, 0.975))

p.tavg <- out.list[[14]]
p.tavg <- c(p.tavg[,1], p.tavg[,2], p.tavg[,3])
hist(p.tavg)
quantile(p.tavg, probs=c(0.025, 0.5, 0.975))

p.burn <- out.list[[15]]
p.burn1 <- c(p.burn[,1], p.burn[,2], p.burn[,3])
hist(p.burn1)
quantile(p.burn1, probs=c(0.025, 0.5, 0.975))

p.burn <- out.list[[16]]
p.burn2 <- c(p.burn[,1], p.burn[,2], p.burn[,3])
hist(p.burn2)
quantile(p.burn2, probs=c(0.025, 0.5, 0.975))

p.burn <- out.list[[17]]
p.burn3 <- c(p.burn[,1], p.burn[,2], p.burn[,3])
hist(p.burn3)
quantile(p.burn3, probs=c(0.025, 0.5, 0.975))

p.dist <- out.list[[18]]  
p.dist <- c(p.dist[,1], p.dist[,2], p.dist[,3])
hist(p.dist)
quantile(p.dist, probs=c(0.025, 0.5, 0.975))

sigma <- out.list[[128]]
sigma <- c(sigma[,1], sigma[,2], sigma[,3])
hist(sigma)
quantile(sigma, probs=c(0.025, 0.5, 0.975))

prob.male1 <- out.list[[7]]
prob.male1 <- c(prob.male1[,1], prob.male1[,2], prob.male1[,3])
hist(prob.male1)
quantile(prob.male1, probs=c(0.025, 0.5, 0.975))

prob.male2 <- out.list[[8]]
prob.male2 <- c(prob.male2[,1], prob.male2[,2], prob.male2[,3])
hist(prob.male2)
quantile(prob.male2, probs=c(0.025, 0.5, 0.975))

prob.male3 <- out.list[[9]]
prob.male3 <- c(prob.male3[,1], prob.male3[,2], prob.male3[,3])
hist(prob.male3)
quantile(prob.male3, probs=c(0.025, 0.5, 0.975))

#Phi effects
phi0.1 <- out.list[[19]]
phi0.1 <- c(phi0.1[,1], phi0.1[,2], phi0.1[,3])
hist(phi0.1)
quantile(phi0.1, probs=c(0.025, 0.5, 0.975))

phi0.2 <- out.list[[20]]
phi0.2 <- c(phi0.2[,1], phi0.2[,2], phi0.2[,3])
hist(phi0.2)
quantile(phi0.2, probs=c(0.025, 0.5, 0.975))

phi0.3 <- out.list[[21]]
phi0.3 <- c(phi0.3[,1], phi0.3[,2], phi0.3[,3])
hist(phi0.3)
quantile(phi0.3, probs=c(0.025, 0.5, 0.975))

phi.male <- out.list[[22]]
phi.male <- c(phi.male[,1], phi.male[,2], phi.male[,3])
hist(phi.male)
quantile(phi.male, probs=c(0.025, 0.5, 0.975))

phi.scl <- out.list[[23]]
phi.scl <- c(phi.scl[,1], phi.scl[,2], phi.scl[,3])
hist(phi.scl)
quantile(phi.scl, probs=c(0.025, 0.5, 0.975))

phi.burn1 <- out.list[[24]]
phi.burn1 <- c(phi.burn1[,1], phi.burn1[,2], phi.burn1[,3])
hist(phi.burn1)
quantile(phi.burn1, probs=c(0.025, 0.5, 0.975))

phi.burn2 <- out.list[[25]]
phi.burn2 <- c(phi.burn2[,1], phi.burn2[,2], phi.burn2[,3])
hist(phi.burn2)
quantile(phi.burn2, probs=c(0.025, 0.5, 0.975))

phi.burn3 <- out.list[[26]]
phi.burn3 <- c(phi.burn3[,1], phi.burn3[,2], phi.burn3[,3])
hist(phi.burn3)
quantile(phi.burn3, probs=c(0.025, 0.5, 0.975))


#Check psi's (to make sure augmenting enough)
psi1 <- out.list[[73]]
psi1 <- c(psi1[,1], psi1[,2], psi1[,3])
hist(psi1)


#Pull out N's and get CIs
N.idx <- idx1[27:72,]
N.out <- list()
for(i in 1:46){
	N.out[[i]] <- out.list[[i+26]]
}

Ns <- NULL
for(i in 1:46){
	temp <- N.out[[i]]
	temp2 <- c(temp[,1], temp[,2], temp[,3])
	Ns <- cbind(Ns, temp2)
}
Ns <- as.data.frame(Ns)
names(Ns) <- idx1$V1[27:72]
#Credible intervals
N.ci.up <- apply(Ns, 2, function(x) quantile(x, probs=0.975))
N.ci.md <- apply(Ns, 2, function(x) quantile(x, probs=0.5))
N.ci.lo <- apply(Ns, 2, function(x) quantile(x, probs=0.025))
N.ci <- rbind(N.ci.up, N.ci.md, N.ci.lo)



#Pull out Nrec's and get CIs (immigration and recruitment)
Nrec.idx <- idx1[76:121,]
Nrec.out <- list()
for(i in 1:46){
	Nrec.out[[i]] <- out.list[[i+75]]
}

Nrecs <- NULL
for(i in 1:46){
	temp <- Nrec.out[[i]]
	temp2 <- c(temp[,1], temp[,2], temp[,3])
	Nrecs <- cbind(Nrecs, temp2)
}
Nrecs <- as.data.frame(Nrecs)
names(Nrecs) <- idx1$V1[76:121]
#Credible intervals
Nrec.ci.up <- apply(Nrecs, 2, function(x) quantile(x, probs=0.975))
Nrec.ci.md <- apply(Nrecs, 2, function(x) quantile(x, probs=0.5))
Nrec.ci.lo <- apply(Nrecs, 2, function(x) quantile(x, probs=0.025))
Nrec.ci <- rbind(Nrec.ci.up, Nrec.ci.md, Nrec.ci.lo)


#Get population size in format for analysis (columns: site, year, median abundance)
out <- data.frame(N.ci.md, N.ci.up, N.ci.lo, Nrec.ci.md, Nrec.ci.up, Nrec.ci.lo)
out$site <- substr(names(Ns), 3, 3)
out$period <- ifelse(nchar(names(Ns))==6, substr(names(Ns), 5, 5), substr(names(Ns), 5, 6))
yrs <- c(1994:2008, 2010:2011, 2014:2016, 2018:2019)
yrs.periods <- as.data.frame(cbind(yrs, 1:22))
names(yrs.periods) <- c("yrs", "prds")
out$yr <- yrs.periods$yrs[match(out$period, yrs.periods$prds)]
out$site.nm <- ifelse(out$site==1, "Fs", ifelse(out$site==2, "Es", "E12"))
out$zone.yr <- paste(out$site.nm, out$yr, sep="-")



###########
#Read in burn data
setwd(DataDir2)
burn <- read.csv("yrs_burn_pls1.csv", header=TRUE)  #Not weighted by point density
burn$zone.yr <- paste(burn$Zone, burn$Year, sep="-")
#Years 2010, 2014, and 2018 come after a gap in tort data collection: ysb should be averaged w/ previous year (2014, previous 2 years)
burn$ysb.tort <- NULL
for(i in 1:nrow(burn)){
	burn$ysb.tort[i] <- ifelse(burn$Year[i]==2010 | burn$Year[i]==2018, mean(c(burn$Yrs_since_mn[i], burn$Yrs_since_mn[i-1])), ifelse(burn$Year[i]==2014, mean(c(burn$Yrs_since_mn[i], burn$Yrs_since_mn[i-1], burn$Yrs_since_mn[i-2])), burn$Yrs_since_mn[i]))
}



#Weighted linear regression (weighted by CIs) accounting for serial autocorrelation
out.burn <- merge(out, burn, by="zone.yr")
#Add in the 95% CI for N (weights are 1/size of CI)
out.burn$N.weight <- 1/(out.burn$N.ci.up - out.burn$N.ci.lo)

w <- out.burn$N.weight[out.burn$Zone=="Fs"]
y <- out.burn$N.ci.md[out.burn$Zone=="Fs"]
x <- out.burn$ysb.tort[out.burn$Zone=="Fs"]
t <- out.burn$yr[out.burn$Zone=="Fs"]
n.fs.ar.w <- gls(y ~ x, weights=varFunc(~w), correlation=corCAR1(form = ~1 | t))
summary(n.fs.ar.w)

w <- out.burn$N.weight[out.burn$Zone=="Es"]
y <- out.burn$N.ci.md[out.burn$Zone=="Es"]
x <- out.burn$ysb.tort[out.burn$Zone=="Es"]
t <- out.burn$yr[out.burn$Zone=="Es"]
n.es.ar.w <- gls(y ~ x, weights=varFunc(~w), correlation=corCAR1(form = ~1 | t))
summary(n.es.ar.w)


#Immigration/recruitment (no serial autocorrelation for immigration/recruitment)
#Remove first period for each zone (equivalent to N)
out.burn.rec <- out.burn[out.burn$zone.yr != "Es-1994",]
out.burn.rec <- out.burn.rec[out.burn.rec$zone.yr != "Fs-1996",]
out.burn.rec <- out.burn.rec[out.burn.rec$zone.yr != "E12-2001",]
out.burn.rec$Nrec.weight <- 1/(out.burn.rec$Nrec.ci.up - out.burn.rec$Nrec.ci.lo)

nrec.fs.w <- rlm(out.burn.rec$Nrec.ci.md[out.burn.rec$Zone=="Fs"] ~ out.burn.rec$Yrs_since_mn[out.burn.rec$Zone=="Fs"], weights=out.burn.rec$Nrec.weight[out.burn.rec$Zone=="Fs"], wt.method="inv.var")
summary(nrec.fs.w)

nrec.es.w <- rlm(out.burn.rec$Nrec.ci.md[out.burn.rec$Zone=="Es"] ~ out.burn.rec$Yrs_since_mn[out.burn.rec$Zone=="Es"], weights=out.burn.rec$Nrec.weight[out.burn.rec$Zone=="Es"], wt.method="inv.var")
summary(nrec.es.w)





#
