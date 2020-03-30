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
source("GT_data_20191127.txt")


###############
# BUGS MODEL
###############

##START THE FILE FOR JAGS

BUGSfilename <- "GT_BUGS_v8.txt"
setwd(BUGSDir)

cat("
    
    model{
    
#############
# PRIORS
#############
      
      ###########
      # PROBABILITY OF CAPTURE IN SURVEYED AREA
      ###########
      
      p0 ~ dunif(0,1)                 # mean/intercept detection prob within the surveyed area
      logit.p0 <- log(p0/(1-p0))   # convert to logit
      p.male.eff ~ dunif(-3,3)    # effect of sex on mean capture probability
	p.scl.eff ~ dunif(-3,3)
	p.effort.eff ~ dunif(-5,5)
	p.tavg.eff ~ dunif(-3,3)
	p.bucket.eff ~ dunif(-5,5)
	sigma ~ dunif(0, 4)  #drop off in detection with distance
	sigma2 <- sigma*sigma

      
      #### add random effect for period to soak up any additional variance in capture probability b/w sessions...?
      
      #### mean capture probability
    for(s in 1:nsites){
	for(ind in 1:ninds[s]){  # loop through (data augmented) individuals
        for(period in first.period[s]:last.period[s]){
		dist[s,ind,period] <- exp(-1*pow(rd.dist[s,ind,period],2) / (2*sigma2))  
          for(session in 1:nsessions[period]){
		mu.g[s,ind,period,session] <- 1/(1+exp(-1*logit.p0 + p.male.eff*is.male[s,ind] 
									+ p.scl.eff*scl[s,ind] 
									+ p.effort.eff*effort[s,period,session]
									+ p.bucket.eff*bucket[period]
									+ p.tavg.eff*tavg[period,session]
									))
            thisp0[s,ind,period,session] <- mu.g[s,ind,period,session] * dist[s,ind,period]
									
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

      phi.male.eff ~ dunif(-3,3)    # effect of sex on mean capture probability
	phi.scl.eff ~ dunif(-3,3)

	for(s in 1:nsites){
      	phi0[s] ~ dunif(0.1,1)                
      	phi0.logit[s] <- log(phi0[s]/(1-phi0[s])) 

		phi.year.prec[s] ~ dgamma(0.01,0.01)
		phi.year.sd[s] <- pow((1/phi.year.prec[s]),0.5)
		#phi[s,,first.period[s]] <- 1  #Version 14: changed from phi[s,1] <- 1	
		for(period in (first.period[s]+1):last.period[s]){
			phi.year.eff[s,period] ~ dnorm(0,phi.year.prec[s])
		for(ind in 1:ninds[s]){ #Version 14: added in individual level effects on survival (sex, size)
			logit(phi[s,ind,period]) <- phi0.logit[s] + phi.year.eff[s,period] + phi.male.eff*is.male[s,ind] + phi.scl.eff*scl[s,ind] 
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

		#Distance from road
		rd.shape[s] ~ dunif(0.2, 3)
		rd.rate[s] ~ dunif(1, 10)
		for(ind in 1:ninds[s]){
			for(period in first.period[s]:last.period[s]){
				rd.dist[s,ind,period] ~ dgamma(rd.shape[s], rd.rate[s])
			}
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
      
     for(s in 1:nsites){
      for(period in first.period[s]:last.period[s]){
        N[s,period] <- inprod(in.pop.now[s,1:ninds[s],period],z[s,1:ninds[s]])
        N.recruited[s,period] <- inprod(recruited[s,1:ninds[s],period],z[s,1:ninds[s]])
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
    p0=runif(1,0.1,0.15),
    p.male.eff=runif(1,-0.1,0.1),
    p.effort.eff=runif(1,0.01,0.1),
    p.scl.eff=runif(1,-0.1,0.1),
    p.tavg.eff=runif(1,-0.1,0.1),
    p.bucket.eff=runif(1,-0.1,0.1),
    
    phi0=runif(nsites,0.9,0.99),  
    phi.year.prec=runif(nsites,0.01, 0.05),
    phi.male.eff=runif(1,-0.1,0.1),
    phi.scl.eff=runif(1,-0.1,0.1),
    psi = runif(nsites,0.5,0.9),
    init.entranceprob = runif(nsites,3,4),

    sigma = runif(1, 1, 2),
    rd.shape = runif(nsites, 0.5, 2),
    rd.rate = runif(nsites, 3, 7),
    
    ##  initialize every indiv as being in the population
    in.pop.now = init.in.pop.now,

    prob.male = c(0.5, 0.5, 0.5)   

  )
}


#Running JAGS
setwd(BUGSDir)
system.time(
  mod<-run.jags(
    model=BUGSfilename,
    monitor=c("p0","p.male.eff","prob.male", "p.effort.eff", "p.scl.eff", "scl.sd",	"p.tavg.eff",
              "phi0", "phi", "phi.year.sd",
              "psi","sigma", "rd.shape", "rd.rate", "N", "N.recruited",
              "gamma","init.entranceprob"),
    data=data.for.bugs,
    n.chains = 3,
    inits=initz.bugs,
    burnin = 10000,  
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
#ANALYZING JAGS MODEL OUTPUT
DataDir <- paste(rootDir,"Data/JAGS_output_20191127",sep="")
DataDir2 <- paste(rootDir,"Data", sep="")

#Model was analyzed using a high powered computing cluster, so read in output files:
setwd(DataDir)
idx1 <- read.csv("output_20191127_chain1index.csv", header=FALSE)
chn1 <- read.csv("output_20191127_chain1chain1.csv", header=FALSE)
chn2 <- read.csv("output_20191127_chain2chain1.csv", header=FALSE)
chn3 <- read.csv("output_20191127_chain3chain1.csv", header=FALSE)

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



#Pull out phi year effects and get CIs
phi.idx <- idx1[16:58,]
phi.out <- list()
for(i in 1:43){
	phi.out[[i]] <- out.list[[i+15]]
}

phis <- NULL
for(i in 1:43){
	temp <- phi.out[[i]]
	temp2 <- c(temp[,1], temp[,2], temp[,3])
	phis <- cbind(phis, temp2)
}
phis <- as.data.frame(phis)
names(phis) <- idx1$V1[16:58]

#Back calculate to 0-1 prob range, add in phi0 (don't need to add in scl.eff because centered on 0, and just report female, not male phis)
phi0.idx <- idx1[13:15,]
phi0.out <- list()
for(i in 1:3){
	phi0.out[[i]] <- out.list[[i+12]]
}

phi0s <- NULL
for(i in 1:3){
	temp <- phi0.out[[i]]
	temp2 <- c(temp[,1], temp[,2], temp[,3])
	phi0s <- cbind(phi0s, temp2)
}
phi0s <- as.data.frame(phi0s)
names(phi0s) <- idx1$V1[13:15]
temp <- apply(phi0s, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
1/(1+exp(-1*temp))

phis.c <- NULL
for(i in 1:43){
	temp <- phis[,i]
	temp.site <- as.numeric(substr(names(phis)[i], 14, 14))
	temp.phi <- 1/(1 + exp(-1*(temp + phi0s[temp.site])))
	temp.phi <- unlist(temp.phi)
	phis.c <- cbind(phis.c, temp.phi)
}
phis.c <- as.data.frame(phis.c)
names(phis.c) <- names(phis)

#Credible intervals
phi.ci.up <- apply(phis.c, 2, function(x) quantile(x, probs=0.975))
phi.ci.md <- apply(phis.c, 2, function(x) quantile(x, probs=0.5))
phi.ci.lo <- apply(phis.c, 2, function(x) quantile(x, probs=0.025))
phi.ci <- rbind(phi.ci.up, phi.ci.md, phi.ci.lo)



#Pull out N's and get CIs
N.idx <- idx1[61:106,]
N.out <- list()
for(i in 1:46){
	N.out[[i]] <- out.list[[i+60]]
}

Ns <- NULL
for(i in 1:46){
	temp <- N.out[[i]]
	temp2 <- c(temp[,1], temp[,2], temp[,3])
	Ns <- cbind(Ns, temp2)
}
Ns <- as.data.frame(Ns)
names(Ns) <- idx1$V1[61:106]
#Credible intervals
N.ci.up <- apply(Ns, 2, function(x) quantile(x, probs=0.975))
N.ci.md <- apply(Ns, 2, function(x) quantile(x, probs=0.5))
N.ci.lo <- apply(Ns, 2, function(x) quantile(x, probs=0.025))
N.ci <- rbind(N.ci.up, N.ci.md, N.ci.lo)

x <- c(1:46)
plot(N.ci.md ~ x, ylim=c(0,250))
segments(x0=x, x1=x, y0=N.ci.lo, y1=N.ci.up)



#Pull out N recruited and get CIs
Nrec.idx <- idx1[110:155,]
Nrec.out <- list()
for(i in 1:46){
	Nrec.out[[i]] <- out.list[[i+109]]
}

Nrecs <- NULL
for(i in 1:46){
	temp <- Nrec.out[[i]]
	temp2 <- c(temp[,1], temp[,2], temp[,3])
	Nrecs <- cbind(Nrecs, temp2)
}
Nrecs <- as.data.frame(Nrecs)
names(Nrecs) <- idx1$V1[110:155]
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


#Read in prescribed burn data
setwd(DataDir2)
burn <- read.csv("yrs_burn_pls1.csv", header=TRUE)  #Not weighted by point density
burn$zone.yr <- paste(burn$Zone, burn$Year, sep="-")

#Dataset for N and N.rec
out.burn <- merge(out, burn, by="zone.yr")
#Remove first period for each zone
out.burn <- out.burn[out.burn$zone.yr != "Es-1994",]
out.burn <- out.burn[out.burn$zone.yr != "Fs-1996",]
out.burn <- out.burn[out.burn$zone.yr != "E12-2001",]
#Add in the 95% CI for N (weights are 1/size of CI)
out.burn$N.weight <- 1/(out.burn$N.ci.up - out.burn$N.ci.lo)
out.burn$Nrec.weight <- 1/(out.burn$Nrec.ci.up - out.burn$Nrec.ci.lo)

#Different dataset for phi because of fewer estimated years for phi
#Remove E12 because only 2 years
out.phi <- data.frame(phi.ci.md, phi.ci.up, phi.ci.lo)
out.phi$site <- substr(names(phis), 14, 14)
out.phi$period <- ifelse(nchar(names(phis))==17, substr(names(phis), 16, 16), substr(names(phis), 16, 17))
yrs <- c(1994:2008, 2010:2011, 2014:2016, 2018:2019)
yrs.periods <- as.data.frame(cbind(yrs, 1:22))
names(yrs.periods) <- c("yrs", "prds")
out.phi$yr <- yrs.periods$yrs[match(out.phi$period, yrs.periods$prds)]
out.phi$site.nm <- ifelse(out.phi$site==1, "Fs", ifelse(out.phi$site==2, "Es", "E12"))
out.phi$zone.yr <- paste(out.phi$site.nm, out.phi$yr, sep="-")
out.phi <- merge(out.phi, burn, by="zone.yr")
out.t <- out.phi[out.phi$Zone!="E12",]
out.t <- out.t[out.t$zone.yr!="Es-1994" & out.t$zone.yr!="Fs-1994" & out.t$zone.yr!="Fs-1997",]
#Add in the 95% CI for phi (weights are 1/size of CI)
out.t$phi.weight <- 1/(out.t$phi.ci.up - out.t$phi.ci.lo)

#######
#Population variables vs. years-since-burn
#Weighted linear regressions

#N - population size
n.fs.w <- rlm(out.burn$N.ci.md[out.burn$Zone=="Fs"] ~ out.burn$Yrs_since_mn[out.burn$Zone=="Fs"], weights=out.burn$N.weight[out.burn$Zone=="Fs"], wt.method="inv.var")
n.es.w <- rlm(out.burn$N.ci.md[out.burn$Zone=="Es"] ~ out.burn$Yrs_since_mn[out.burn$Zone=="Es"], weights=out.burn$N.weight[out.burn$Zone=="Es"], wt.method="inv.var")

#Phi - mortality/emigration (change to "probability of leaving")
#Remove years 2010 and 2014 because they are years following a gap in surveying (phi is averaged across 2 years but burning is not)
phi.fs.w <- rlm(1-out.t$phi.ci.md[out.t$Zone=="Fs" & out.t$yr!=2014 & out.t$yr!=2010] ~ out.t$Yrs_since_mn[out.t$Zone=="Fs" & out.t$yr!=2014 & out.t$yr!=2010], weights=out.t$phi.weight[out.t$Zone=="Fs" & out.t$yr!=2014 & out.t$yr!=2010], wt.method="inv.var")
phi.es.w <- rlm(1-out.t$phi.ci.md[out.t$Zone=="Es" & out.t$yr!=2014 & out.t$yr!=2010] ~ out.t$Yrs_since_mn[out.t$Zone=="Es" & out.t$yr!=2014 & out.t$yr!=2010], weights=out.t$phi.weight[out.t$Zone=="Es" & out.t$yr!=2014 & out.t$yr!=2010], wt.method="inv.var")

#Nrec - recruitment/immigration
#Remove years 2010 and 2014 because they are years following a gap in surveying (Nrec is averaged across 2 years but burning is not)
nrec.fs.w <- rlm(out.burn$Nrec.ci.md[out.burn$Zone=="Fs" & out.burn$yr!=2014 & out.burn$yr!=2010]/out.burn$N.ci.md[out.burn$Zone=="Fs" & out.burn$yr!=2014 & out.burn$yr!=2010] ~ out.burn$Yrs_since_mn[out.burn$Zone=="Fs" & out.burn$yr!=2014 & out.burn$yr!=2010], weights=out.burn$Nrec.weight[out.burn$Zone=="Fs" & out.burn$yr!=2014 & out.burn$yr!=2010], wt.method="inv.var")
nrec.es.w <- rlm(out.burn$Nrec.ci.md[out.burn$Zone=="Es" & out.burn$yr!=2014 & out.burn$yr!=2010]/out.burn$N.ci.md[out.burn$Zone=="Es" & out.burn$yr!=2014 & out.burn$yr!=2010] ~ out.burn$Yrs_since_mn[out.burn$Zone=="Es" & out.burn$yr!=2014 & out.burn$yr!=2010], weights=out.burn$Nrec.weight[out.burn$Zone=="Es" & out.burn$yr!=2014 & out.burn$yr!=2010], wt.method="inv.var")






#