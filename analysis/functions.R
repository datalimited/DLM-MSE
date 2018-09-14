# functions.R

# Copyright 2003-2012 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC (modified for MSE project by Jessica Walsh)

# setupStock {{{
setupStock <- function(brp, iniBiomass, nyears, iters) {
	
	# find corresponding F
	#idx <- max(which(stock(brp) <= iniBiomass)[1], 2)
	#idx <- max(which(ssb(brp) <= iniBiomass)[1], 2) #original
	idx <- which(stock(brp) <= as.numeric(iniBiomass))[1]
	
	# create from FLBRP
	stk <- as(brp, 'FLStock')[, idx]

	# change 'year' dimnames
	dimnames(stk) <- list(year=1)

	# expand
	stk <- stf(stk, nyears-1, 1)

	# get total biomass recomputed
	stock(stk) <- computeStock(stk)

	# add iters
	stk <- propagate(stk, iters)

	return(stk)
} # }}}

# effortDynamics {{{
effortDynamics <- function(stk, bmsy, sr, years=2:dims(stk)$maxyear, xp, a= 1, iters=1, rsd = 0,
	srres=rnorm(iters, FLQuant(0, dimnames=list(year=years)), 0)) {
for (year in years) {
		har <- fbar(stk)[,year-1]
		#bio <- sum(stock.n(stk)[,year-1] * stock.wt(stk)[,year-1] * catch.sel(brp)) # exploited biomass, but need to keep it consitent with bmsy
		bio <- stock(stk)[,year-1]
	
		# Target F dynamics
		# F_{y} = F_{y-1} * (B_{y-1} / (BMSY * a)) ^ x
		eff <- har * (bio / (bmsy * a)) ^ xp
		future_bau_iters <- c(eff) * rlnorm(1 * iters, meanlog = 0 - (rsd^2/2), sdlog=rsd) 
		#note this implementation error only works for 1 iteration. 
		#if need more iterations you need to add a target array
		
		# fwdControl
		fctl <- fwdControl(
			data.frame(year=(FLCore::as.data.frame(stk@catch[,year]) %>% dplyr::select(year))[1,], 
								 quantity='f', 
								 val=c(future_bau_iters[1])))

		# fwd
		stk <- fwd(stk, ctrl = fctl, sr=sr, 
							 sr.residuals=srres[,year], 
							 sr.residuals.mult=TRUE)
	}
	#cat("\n")
	return(stk)
} # }}}

# oneWayTrip {{{
oneWayTrip <- function(stk, sr, fmax=refpts(brp)['crash', 'harvest']*0.80,
	years=2:dims(stk)$maxyear, srres=rnorm(iters, FLQuant(0, dimnames=list(year=years)), 0),
	nyear_proj, start_yr_f, rsd = 0) {
	# rsd only works if there is 1 iteration
	# check if it still works with simulation 01-sims.R. 
	
	# limits
	#start_yr_f is the index of the year that f0 starts from. e.g. in a sim from year 1-10, start_yr_f = 1. 
	# If extending harvest dynamic forward from an existing time series of 20 yrs, start_yr_f could equal 20
	f0 <- c(fbar(stk)[1,start_yr_f,,,,1])
	#f0 <- c(fbar(stk)[,1])
	fmax <- c(fmax)
	rate <- exp((log(fmax) - log(f0)) / (length(years)+nyear_proj))

	# linear trend
	f <-rate^(0:length(years)) *f0
	future_bau_iters <- c(f) * rlnorm(length(f), meanlog = 0 - (rsd^2/2), sdlog=rsd) 

	# fwdControl
	fctl <- fwdControl(data.frame(year=years, 
																quantity='f', 
																val=future_bau_iters[-1]))

	# fwd
	stk <- fwd(stk, fctl, sr=sr, sr.residuals=srres, sr.residuals.mult=TRUE)
	
	return(stk)
} # }}}


# oneWayQuickTrip {{{
oneWayQuickTrip <- function(stk, sr, fmax=refpts(brp)['crash', 'harvest']*0.80,
														rate=1.20, years=2:dims(stk)$maxyear,
														srres=rnorm(iters, FLQuant(0, dimnames=list(year=years)), 0)) {
	
	# limits
	f0 <- c(fbar(stk)[,1])
	fmax <- c(fmax)
	
	# linear trend
	f <- rate^(0:length(years))*f0
	f[f > fmax] <- fmax
	
	
	# fwdControl
	fctl <- fwdControl(data.frame(year=years, quantity='f', val=f[-1]))
	
	# fwd
	stk <- fwd(stk, fctl, sr=sr, sr.residuals=srres, sr.residuals.mult=TRUE)
	
	return(stk)
} # }}}


# rollerCoaster {{{
rollerCoaster <- function(stk, sr, fmax=refpts(brp)['crash', 'harvest']*0.80,
	fmsy=refpts(brp)['msy', 'harvest'], years=2:dims(stk)$maxyear, up=0.05, down=0.04,
	srres=rnorm(iters, FLQuant(0, dimnames=list(year=years)), 0)) {

	# F
	f <- rep(NA, length(years))

	# limits
	f0 <- c(fbar(stk)[,1])
	fmax <- c(fmax)

	# linear trend
	rateup <- log(fmax/f0) / log(1 + up)
	fup <- f0 * ((1 + up) ^ (0:ceiling(rateup)))
	lfup <- length(fup)
	f[1:lfup] <- fup

	# at the top
	f[lfup:(lfup+5)] <- fup[lfup]

	# coming down!
	ratedo <- log(fmsy/f[length(fup)+5]) / log(1 + down)
	lfdo <- length(f) - (lfup +6) + 1
	fdo <- f[lfup+5] * ((1 + down) ^ seq(0, ceiling(ratedo), length=lfdo))
	f[(lfup+6):length(f)] <- fdo[1:lfdo]
	
	# fwdControl
	fctl <- fwdControl(data.frame(year=years, quantity='f', val=f))

	# fwd
	stk <- fwd(stk, fctl, sr=sr, sr.residuals=srres, sr.residuals.mult=TRUE)
	
	return(stk)
} # }}}



# rollerCoaster2 {{{
rollerCoaster2 <- function(stk, sr, fmax=refpts(brp)['crash', 'harvest']*0.80,
													 fmsy=refpts(brp)['msy', 'harvest'], years=2:dims(stk)$maxyear, upy=25, top=5,
													 downy=(dims(stk)$maxyear)-(upy+top),
													 srres=rnorm(iters, FLQuant(0, dimnames=list(year=years)), 0)) {
	
	# F
	f <- rep(NA, length(years))
	
	# limits
	f0 <- c(fbar(stk)[,1])
	fmax <- c(fmax)
	
	# linear trend up: 1:upy
	fup <- seq(f0, fmax, length=upy)
	f[1:upy] <- fup
	
	# at the top: upy+1:upy+6
	f[(upy+1):(upy+top-1)] <- fmax
	
	# coming down!
	fdo <- seq(fmax, fmsy, length=downy+1)[-1]
	f[(upy+top):(upy+top+downy-1)] <- fdo
	
	# fwdControl
	fctl <- fwdControl(data.frame(year=years, quantity='f', val=f))
	
	# fwd
	stk <- fwd(stk, fctl, sr=sr, sr.residuals=srres, sr.residuals.mult=TRUE)
	
	return(stk)
} # }}}


# ar1rlnorm {{{
#From Thorson et al. 2014 paper about residiuals and recruitment.
#But errors in formulas 16&17. This is the correct formula. 
#haven't discarded negative values. should I, as described in Thorsons paper? 
ar1rlnorm <- function(rho, years, iters=1, margSD=0.6) {
 
	n_t <- length(years)
	rhosq <- rho ^ 2
	
	res <- matrix(rnorm(n_t*iters, mean=0, sd=margSD), nrow=n_t, ncol=iters)
	res <- apply(res, 2, function(x) {
		for(t in 2:n_t)
		x[t] = rho*x[t-1] + sqrt(1-rhosq)*x[t]
		return(exp(x - margSD^2/2))
				})
	
		return(FLQuant(array(res, dim=c(1,n_t,1,1,1,iters)),
									 dimnames=list(year=years, iter=seq(1, iters))))
	
}

# }}}

# HCR 4010 Rule {{{
HCR_4010 <- function(yaxis_msy, HCR_lim, HCR_trig, HCR_curr){
	##y axis is either Hmsy for HR HCR or Msy for Catch HCR
	## each argument should be either 'True' or 'estimated' 
	## HCR_lim, HCR_trig & HCR_curr should all be in the same units, e.g. Biomass or BBmsy
	intercept <- - yaxis_msy * HCR_lim / (HCR_trig - HCR_lim)
	slope <- yaxis_msy / (HCR_trig - HCR_lim)
	
	hcr <- ifelse(HCR_curr <= HCR_lim, 0,
								ifelse(HCR_curr > HCR_lim & HCR_curr < HCR_trig,
											 intercept + slope * HCR_curr, yaxis_msy))
	return(hcr)
}

#}}}


# Stepwise HCR {{{
HCR_stepwise <- function(yaxis_msy, HCR_lim_sw, HCR_trig_sw, HCR_curr){
	##y axis is either Hmsy for HR HCR or Msy for Catch HCR
	## each argument should be either 'True' or 'estimated' 
	## HCR_lim_sw, HCR_trig_sw & HCR_curr should all be in the same units, e.g. Biomass or BBmsy
	hcr <-	ifelse (HCR_curr >= HCR_trig_sw, yaxis_msy,
								 ifelse (HCR_curr < HCR_trig_sw & HCR_curr >= HCR_lim_sw, 
								 				yaxis_msy *0.5, 0 ))
	return(hcr)
}

#}}}