## Simulating case studies for MSE
## Started August 2015
## adapted by Jessica Walsh from Working Group simulations (by Iago Mosquiera, Coilin Minto, Chato Osio)

library(dplyr)
library(plyr)
library(FLBRP)		
library(FLAssess)	
library(FLash)

source('analysis/functions.R')

## VARS
#need to run this for each harvest dynamics - OW and ED03 manually

set.seed(1234)
nyears <- 60 # Max. number of years of simulations
nyear_proj <- 20 # number of years projected forward in management
iters <- 600 # number of replicates bump up to 1000 if needed, but 600 ensured convergence
ed = "OW" #one-way trip = OW. Set to "ED03" for bio-economic coupled effort dynamics 

#+++++++++++++++++++++++++++++++++++++
#Life history parameters 
# details about these parameters and references are in excel sheet data-raw/sp_params.xlsx
# sce: Scenarios list {{{
sce <- list(
	LH=list(
		# BC Bocaccio rockfish
		BC =list(
			par=FLPar(linf=75.9, sl=1.5, sr=2000, s=0.662, a50 = 5.7, v=47268, k=0.184, a=7.36E-06, 
								b= 3.114, bg= 3.114),
			range=c(min=1, max=48, minfbar=2, maxfbar=48),
			ID=0.7, rho = 0.426, margSD = 0.778, sigmaC = 0.2, UR = 0.2, 
			resilience = "low", species_cat = "Miscellaneous demersal fishes"),
		# SA Pacific sardine 
		SA =list(
			par=FLPar(linf=23.46, sl=2, sr=2000, s=0.8, a50 = 1.5, v=1008430,  k = 0.386, a=7.524E-06, 
								b= 3.233, bg= 3.233),
			range=c(min=1, max=15, minfbar=1, maxfbar=15),
			ID=0.1, rho = 0.435, margSD = 0.766, sigmaC = 0.2, UR = 0.2,
			resilience = "medium", species_cat = "Herrings, sardines, anchovies"),
		# PS Petrale sole 
		PS =list(
			par=FLPar(linf=54.31, sl=3.25, sr=2000, s=0.86, a50 = 10, v=32426, k = 0.13, a=2.083E-06, 
								b= 3.474, bg= 3.474),
			range=c(min=1, max=20, minfbar=3, maxfbar=20),
			ID=0.8, rho = 0.437, margSD = 0.636, sigmaC = 0.2, UR = 0.2,
			resilience = "medium", species_cat = "Flounders, halibuts, soles"),
		# SJ Skipjack 
		#note: k growth estimate is reduced by 80% to make mortality more realistic
		SJ =list(
			par=FLPar(linf=81, sl=0.4, sr=5000, s=0.8, a50 = 2, v=69931650, k = 0.44, a=5.529E-06, 
								b= 3.336, bg= 3.336),
			range=c(min=1, max=12, minfbar=2, maxfbar=12),
			ID=0.1, rho = 0.466, margSD = 0.777, sigmaC = 0.5, UR = 0.5,
			resilience = "medium", species_cat = "Tunas, bonitos, billfishes"),
		# BF Bullet and Frigate tuna
		#note: k growth estimate is reduced by 80% to make mortality more realistic
		BF =list(
			par=FLPar(linf=48.2, sl = 0.5, sr=2000, s=0.9, a50 = 1.5, v=4657840, a50=1.5, k=0.5, a=6E-6, 
								b = 3.194, bg = 3.194),
			range=c(min=1, max=5, minfbar=1, maxfbar=5),
			ID=0.1, rho = 0.466, margSD = 0.777, sigmaC = 0.5, UR = 0.5,
			resilience = "medium", species_cat = "Tunas, bonitos, billfishes"),
		# CR Corvina reina
		CR =list(
			par=FLPar(linf=122, sl=1.5, sr=2000, s=0.6, v=1165620, a50=6, k=0.17, a=2.4E-5, 
								b = 2.824, bg = 2.824),
			range=c(min=1, max=18, minfbar=3, maxfbar=18),
			ID=0.1, rho = 0.466, margSD = 0.777, sigmaC = 0.5, UR = 0.5,
			resilience = "low", species_cat = "Miscellaneous coastal fishes")
	),
	# Effort/F dynamics, x value: RC, ED0, ED0.3, OW
	ED=list(ED03=0.3, OW=0.8)
	#ED=list(ED03=0.3, OW=0.80, RC=0.80)
	)

saveRDS(sce, file="data-raw/sim_life_history_settings.rds")

# RUN for sims
# LH
#lh = "BC"
for(lh in names(sce$LH)) {
	sims <- list()
	par <- gislasim(sce$LH[[lh]]$par)
	
	#set mortality constant at k*1.5 (Jensen, from review paper)
	brp <- lh(par, range=sce$LH[[lh]]$range, fnM= function(par, len) 1.5 * par["k"] * len/len)

	#plot(1:12,m(brp), main = "SJ")
	#plot(fbar(brp),ypr(brp),main="SJ")
	
	# Autocorrelation of recruitment residuals
	srres=ar1rlnorm(rho=sce$LH[[lh]]$rho, iters=iters, years=1:nyears, margSD=sce$LH[[lh]]$margSD)
	
	# Initial depletion
	#stk <- setupStock(brp, iniBiomass=sce$LH[[lh]]$par["v"] * (1-sce$LH[[lh]]$ID), nyears) 
	stk <- setupStock(brp, iniBiomass=sce$LH[[lh]]$par["v"], nyears, iters = iters) # no initial depletion
	
	# ED
	#for(ed in names(sce$ED)) {
		stock <- switch(ed, 
										# one way trip
										"OW" =oneWayTrip(stk, fmax=refpts(brp)['crash', 'harvest']*sce$ED[[ed]],
																		sr=list(model='bevholt', params=params(brp)), years=2:nyears, nyear_proj = nyear_proj,
																		start_yr_f = 1, srres=srres),
										##Note - for BF tuna, we made fcrash = 3 because FLbrp wasn't calculating this value in the ref.pts. 
										#this simulation BF_OW was done manually, after other simulations were constructed. 
										# "OW" =oneWayTrip(stk, fmax=3*sce$ED[[ed]], 
										# 								 sr=list(model='bevholt', params=params(brp)), years=2:nyears, nyear_proj = nyear_proj, 
										# 								 start_yr_f = 1, srres=srres),
										
										## roller coaster
										# "RC1"=rollerCoaster(stk, fmax=refpts(brp)['crash', 'harvest']*sce$ED[[ed]], 
										# 									 fmsy=refpts(brp)['msy', 'harvest'], years=2:nyears, up=0.1, down=0.05,
										# 									 sr=list(model='bevholt', params=params(brp)), srres=srres),
										
										# # RC2
										# "RC"=rollerCoaster2(stk, fmax=refpts(brp)['crash', 'harvest']*sce$ED[[ed]], 
										# 										fmsy=refpts(brp)['msy', 'harvest'], years=2:nyears, upy=25, top=5, downy=30,
										# 										sr=list(model='bevholt', params=params(brp)), srres=srres),
										# effort dynamics
										# "ED0"=effortDynamics(stk, bmsy=c(refpts(brp)['msy', 'biomass']), iters=iters,
										# 										 sr=list(model='bevholt', params=params(brp)), years=2:nyears, xp=sce$ED[[ed]], 
										# 										 srres=srres),
										 "ED03"=effortDynamics(stk, bmsy=c(refpts(brp)['msy', 'biomass']), iters=iters, rsd = 0,
																					sr=list(model='bevholt', params=params(brp)), years=2:nyears, xp=sce$ED[[ed]], a = 0.5,
																					srres=srres)
		)
		
		# NAME
		name <- paste(lh)
		name(stock) <- name
		desc(stock) <- paste(name)
		
		# SIM
		# adding in underreporting if required (currently not using)
		sims[[name]] <- list(lh=par, code=name, stock=stock,
												 brp=brp, 
												 catch= catch(stock),
												 #catchUR=catch(stock)*(1-sce$LH[[lh]]$UR), 
												 effdynamic= ed,
												 resilience=sce$LH[[lh]]$resilience,
												 species_cat=sce$LH[[lh]]$species_cat,
												 margSD=sce$LH[[lh]]$margSD,
												 rho=sce$LH[[lh]]$rho)
		
		print(paste(name, ed))
		
		# Error in C: 20 or 50% CV {{{
		sims[[lh]]$catchE <- FLQuant(aperm(apply(sims[[lh]]$catch, 1:5, function(x) rlnorm(iters, log(x) -sce$LH[[lh]]$sigmaC^2/2, sce$LH[[lh]]$sigmaC)),
																			 c(2,3,4,5,6,1)), dimnames=dimnames(sims[[lh]]$catch))
		
		save(sims, file=paste("data-raw/sims-", lh,"-",ed,".RData", sep=""))
#	}
	# save SRRES file
	save(srres, file=paste("data-raw/srres/srres-",lh,".RData", sep=""))
}

