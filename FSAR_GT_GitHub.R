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
source("GT_data_20201127.txt")



###############
# BUGS MODEL
###############

##START THE FILE FOR JAGS

BUGSfilename <- "GT_BUGS_v23.txt"
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
	p.scl.eff ~ dunif(-3,3)
	p.scl.dist.eff ~ dunif(-3,3)
	p.effort.eff ~ dunif(-5,5)
	p.tavg.eff ~ dunif(-3,3)
	p.bucket.eff ~ dunif(-5,5)
	p.burn.dist.eff ~ dunif(-3,3)
	p.burn.eff ~ dunif(-3,3)
	p.dist.eff ~ dunif(-5,5)
	sigma ~ dunif(0, 4)  #drop off in detection with distance
	sigma2 <- sigma*sigma

      
      
      #### mean capture probability
    for(s in 1:nsites){

	p0[s] ~ dunif(0,1)                 # mean/intercept detection prob within the surveyed area
      logit.p0[s] <- log(p0[s]/(1-p0[s]))   # convert to logit

	for(ind in 1:ninds[s]){  # loop through (data augmented) individuals
        for(period in first.period[s]:last.period[s]){ 
          for(session in 1:nsessions[period]){
		dist[s,ind,period,session] <- exp(-1*pow(rd.dist[s,ind,period,session],2) / (2*sigma2)) 
		logit(thisp0[s,ind,period,session]) <- logit.p0[s] + p.male.eff*is.male[s,ind] 
									+ p.scl.dist.eff*scl[s,ind]*dist[s,ind,period,session] #Size interacts w/ distance
									+ p.scl.eff*scl[s,ind]  #Size also has a main effect (larger tortoises still easier to see at 0 distance)
									+ p.effort.eff*effort[s,period,session]
									+ p.bucket.eff*bucket[period]
									+ p.tavg.eff*tavg[period,session]
									+ p.burn.dist.eff*burn.loc[s,ind,period,session]*dist[s,ind,period,session] #Burn interacts w/ distance
									+ p.burn.eff*burn.loc[s,ind,period,session]
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
     		  gamma[s,period] <- beta[s,period]/sum(beta[s,first.period[s]:last.period[s]])    # probability of entrance...
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

		phi.burn.site.eff[s] ~ dunif(-3,3)	# site-specific effect of years-since-burn

		for(period in (first.period[s]+1):last.period[s]){
		for(ind in 1:ninds[s]){ 
			logit(phi[s,ind,period]) <- phi0.logit[s] + phi.male.eff*is.male[s,ind] 
								+ phi.scl.eff*scl[s,ind] + phi.burn.site.eff[s]*ysb[s,period]
	}
	}
	}
      
      
      ###########
      # PROBABILITY OF BEING REAL
      ###########
      
	for(s in 1:nsites){
      	psi[s] ~ dunif(0,1)   
	}   
      
      
      #############
      # DEAL WITH MISSING DATA
      #############
      
	for(s in 1:nsites){
     		prob.male[s] ~ dunif(0,1)    
      
      	for(ind in 1:ninds[s]){
      	 is.male[s,ind] ~ dbern(prob.male[s])
      	}

		scl.prec[s] ~ dgamma(0.1,0.1)
		scl.sd[s] <- pow((1/scl.prec[s]),0.5)
		for(ind in 1:ninds[s]){
			scl[s,ind] ~ dnorm(0, scl.prec[s])
		}

		#Distance from road (and location specific burning data)
		rd.shape[s] ~ dunif(0.2, 3)
		rd.rate[s] ~ dunif(1, 10)
		for(period in first.period[s]:last.period[s]){
			mu.burn.loc[s,period] ~ dunif(0,5)  #Changed in v20 to be indexed to period as well as site
			burn.loc.prec[s,period] ~ dgamma(0.1, 0.1)
			for(ind in 1:ninds[s]){
				for(session in 1:nsessions[period]){
					rd.dist[s,ind,period,session] ~ dgamma(rd.shape[s], rd.rate[s])
					burn.loc[s,ind,period,session] ~ dnorm(mu.burn.loc[s,period], burn.loc.prec[s,period])
				}
			}
		}

		#Years-since-burn (missing data for 2019, 2020)
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

	#Slope and intercept of ysb effect on N and N.rec (not 2019 or 2020, so do E12 separately)
	for(s in 1:2){
		temp.N[s,1:2] <- inverse(t(X[first.period[s]:(last.period[s]-2),,s]) %*% X[first.period[s]:(last.period[s]-2),,s]) 
			%*% t(X[first.period[s]:(last.period[s]-2),,s]) %*% N[s,first.period[s]:(last.period[s]-2)]
		ysb.N.int[s] <- temp.N[s,1]
		ysb.N.eff[s] <- temp.N[s,2]
		
		#Not first year for Nrec (all individuals are recruited in first year)
		temp.Nrec[s,1:2] <- inverse(t(X[(first.period[s]+1):(last.period[s]-2),,s]) %*% X[(first.period[s]+1):(last.period[s]-2),,s]) 
			%*% t(X[(first.period[s]+1):(last.period[s]-2),,s]) %*% N.recruited[s,(first.period[s]+1):(last.period[s]-2)]
		ysb.Nrec.int[s] <- temp.Nrec[s,1]
		ysb.Nrec.eff[s] <- temp.Nrec[s,2]
	}
	#Now s=3 (E12)
		temp.N[3,1:2] <- inverse(t(X[first.period[3]:last.period[3],,3]) %*% X[first.period[3]:last.period[3],,3]) 
			%*% t(X[first.period[3]:last.period[3],,3]) %*% N[3,first.period[3]:last.period[3]]
		ysb.N.int[3] <- temp.N[3,1]
		ysb.N.eff[3] <- temp.N[3,2]
		
		#Not first year for Nrec (all individuals are recruited in first year)
		temp.Nrec[3,1:2] <- inverse(t(X[(first.period[3]+1):last.period[3],,3]) %*% X[(first.period[3]+1):last.period[3],,3]) 
			%*% t(X[(first.period[3]+1):last.period[3],,3]) %*% N.recruited[3,(first.period[3]+1):last.period[3]]
		ysb.Nrec.int[3] <- temp.Nrec[3,1]
		ysb.Nrec.eff[3] <- temp.Nrec[3,2]


    
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
	X = X,
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
    p.scl.eff=runif(1,-0.1,0.1),
    p.scl.dist.eff=runif(1,-0.1,0.1),
    p.tavg.eff=runif(1,-0.1,0.1),
    p.bucket.eff=runif(1,-0.1,0.1),
    p.burn.dist.eff=runif(1, -0.1, 0.1),
    p.burn.eff=runif(1, -0.1, 0.1),
    p.dist.eff=runif(1, -0.1, 0.1),
    
    phi0=runif(nsitesG,0.9,0.99),
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
    monitor=c("z", "p0", "p.male.eff", "p.effort.eff", "p.bucket.eff", "p.burn.eff", "p.scl.eff", "p.scl.dist.eff", "p.tavg.eff", "p.dist.eff", "p.burn.dist.eff","phi0", "phi.male.eff", "phi.scl.eff", "phi.burn.eff", "psi", "init.entranceprob", "sigma", "rd.shape", "rd.rate", "in.pop.now", "prob.male", "mu.ysb", "ysb.prec", "mu.burn.loc", "burn.loc.prec"),
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



#
