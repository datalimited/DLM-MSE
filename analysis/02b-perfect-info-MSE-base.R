#---------------------------------------------
##4. Conduct MSE using perfect knowledge ####
#---------------------------------------------
#Projects forward using TAC HCR, HR HCR or no management with true values - not using data-limited
#code written by Jessica Walsh, Coilin Minto, Ernesto Jardim etc.
#combined code from MSE tutorial on FLR course.

#note - effort dynamics of 'true' business as usual is HD0.3 or OW0.8

library(FLash)
library(FLBRP)
library(FLa4a)
library(FLAssess)
library(ggplotFL)
library(dplyr)

source('analysis/functions.R')

#### Setup for testing ####
# ts = 20
# sp = "SA"
# j = 2
# niters = 1
# nyear_proj = 20
# rsd = 0.1
# ed = "OW"
# buffer <- 0.75 #or buffer of 80% Hmsy and 80% msy based on buffer used for 
# filename <- paste0(sp,"_TS",ts,"_",ed,"_rep",j)
# load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
# srres <- srres_all[[sp]]
# sr =list(model='bevholt', params=params(sims[[sp]]$brp))
# stock_full <- FLCore::iter(sims[[1]][["stock"]],j)


#### Function ####
HCR_projection_true <- function(stock_full, stock_preassess, sp, ed, ts, rep, nyear_proj, buffer, niters, rsd, srres, sr) {
	
	#load estimated b.bmsy data from DL preassessment
	bbmsy_pre_all <- stock_preassess$bbmsy_pre_COM

	#set stock to work with
	#stock_full <- FLCore::iter(sims[[1]][["stock"]],rep)
	final_yr_sims <- FLCore::dims(stock_full)$maxyear
	stock_pre <- FLCore::trim(stock_full, year=((length(stock_full@catch)-(ts-1)):length(stock_full@catch)))

	true_pre_df <- filter(bbmsy_pre_all, method == "TRUE")
	
	#set.seed(1234)
	
	#------------------
	## Calculated true parameters ####
	#------------------
	#HCRs for Simulated - TRUE data
	#sr <- list(model='bevholt', params = FLCore::params(sims[[1]]$brp)) 
	#true_refpts <- brp(FLBRP(stock_pre, sr = sr))@refpts 
	true_refpts <- FLBRP::brp(FLBRP::FLBRP(stock_full, sr = sr))@refpts
	
	## CALCULATE THE HARVEST RATIO using effort as proxy for 'f'
	Hmsy_true <- an(true_refpts["msy", "harvest"]) 
	Ymsy_true <- an(true_refpts["msy", "yield"])
	K_true <- an(true_refpts["virgin","biomass"]) 
	Bmsy_true <- an(true_refpts["msy", "biomass"])
	BBmsy_curr_true <- an(median(stock(stock_pre)[,c((ts-4):ts),,,,]/ Bmsy_true))
	
	#SSBmsy_true <- an(true_refpts["msy", "ssb"]) 
	#SSBK_true <- an(true_refpts["virgin","ssb"]) 
	
	# Calculate current harvest rate - divide all catch.sims by biomass/year  
	# catch observation error not included here, dependng on DL assessment in 03-stock-assess-dl
	harvest_rate_true <- an(catch(stock_pre))/an(stock(stock_pre)) #instantaneous harvest rate per year 
	harvest_rate_true_curr <- median(harvest_rate_true[(length(harvest_rate_true)-4):length(harvest_rate_true)])
	catch_true_curr <- median(an(catch(stock_pre))[(ts-4):ts])
	

	#------------------
	## 40-10 HCR set up ####
	#------------------
	# HCR based on harvest rate (effort) or based on catch (output - TAC) 
	# Buffer set above 
	
	##Harvest control rule function
	## B 40 (Btrigger) and B10 (Blimit)
	BBmsy_trig <- 0.8 #40% K or 80% BMSY
	BBmsy_lim <- 0.2 	#10% K or 20% BMSY
	
	#Calculating 4010 HCRs - TRUE
	hr_4010_true <- HCR_4010(yaxis_msy = Hmsy_true*buffer, HCR_trig = BBmsy_trig, HCR_lim = BBmsy_lim,
													 HCR_curr = BBmsy_curr_true)
	
	catch_4010_true <- HCR_4010(yaxis_msy = Ymsy_true*buffer, HCR_trig = BBmsy_trig, HCR_lim = BBmsy_lim,
															HCR_curr = BBmsy_curr_true)
	
	#------------------
	## Stepwise HCR set up ####
	#----- -------------
	#Stepwise rule using TAC (Fish at msy until B/Bmsy = 1. When 0.5 < B/Bmsy < 1, fish at 0.5MSY. 
	#When B/Bmsy < 0.5, fishing stops).
	BBmsy_lim_sw <- 0.6
	BBmsy_trig_sw <- 1.1
	
	#Calculating stepwise HCRs - TRUE
	hr_sw_true <- HCR_stepwise(yaxis_msy = Hmsy_true*buffer, HCR_trig_sw = BBmsy_trig_sw, HCR_lim_sw = BBmsy_lim_sw,
														 HCR_curr = BBmsy_curr_true)
	
	catch_sw_true <- HCR_stepwise(yaxis_msy = Ymsy_true*buffer, HCR_trig_sw = BBmsy_trig_sw, HCR_lim_sw = BBmsy_lim_sw,
																HCR_curr = BBmsy_curr_true)	 
	
	param_comp <- data.frame(Model = c("estimated", "true"),
													 # B_curr = c(B_curr_est, B_curr_true), 
													 Bmsy = c(NA, Bmsy_true), 
													 BBmsy_current = c(NA,BBmsy_curr_true),
													 K = c(NA, K_true), 
													 Fmsy= c(NA, Hmsy_true),
													 Ymsy=c(NA,Ymsy_true),
													 HR_curr = c(NA,harvest_rate_true_curr),
													 TAC_4010 = c(NA,catch_4010_true),
													 HR_4010 = c(NA,hr_4010_true),
													 TAC_sw = c(NA,catch_sw_true),
													 HR_sw = c(NA, hr_sw_true))
	param_comp
	#----------------------------
	#Start MSE Projections ####
	#----------------------------
	
	#-------------------------------
	# Settings ####
	#-------------------------------
	
	# new_trgtArray <- array(NA, dim=c(nyear_proj,3,niters), 
	#                        dimnames = list((final_yr_sims+1):(final_yr_sims+nyear_proj), c("min","val","max"),
	#                                        iter=1:niters))
	
	#standardise SR.residual error across simulations 
	#srres=ar1rlnorm(rho=sims[[sp]]$rho, iters=niters, years=(final_yr_sims+1):(nyear_proj+final_yr_sims), margSD=sims[[sp]]$margSD)
	#srres <- FLQuant(rlnorm(niters*nyear_proj, 0 - rsd^2/2, sd=rsd), dimnames = list(year=(final_yr_sims+1):(nyear_proj+final_yr_sims), iter=1:niters))
	
	#-------------------------------
	#FOR 40-10 Catch HCR (TAC) ####
	#-------------------------------
	
	stock_proj <- stock_pre
	val <- catch_4010_true
	
	last.yr <- final_yr_sims
	message(paste(sp, "TAC 4010", rep))
	
	## Extend stock an additional year
	stock_proj <- stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	#fix implementation error - if HCR is 0, then how to add implementation error
	future_catch_iters <- val * rlnorm(nyear_proj * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd)   
	#new_trgtArray[,"val",] <- future_catch_iters
	
	# using TAC - output HCR
	fwd.ctrl <- fwdControl(
		data.frame(year = (last.yr+1):(last.yr + nyear_proj),
							 val = future_catch_iters,
							 #val = val, 
							 quantity = "catch") #,
		#trgtArray = new_trgtArray
	)
	
	stock_proj <- fwd(stock_proj,
										ctrl = fwd.ctrl,
										sr = sr,
										#sr.residuals = srres[1,ac(ts+(n)),1,1,1,1:niters], 
										sr.residuals = srres[1,,1,1,1,rep],
										sr.residuals.mult = TRUE)
	
	stock_proj_TAC4010 <- stock_proj
	
	#--------------------------------------------
	#FOR 40-10 Effort - Harvest Ratio input HCR ####
	#--------------------------------------------
	stock_proj <- stock_pre
	val <- hr_4010_true/harvest_rate_true_curr
	
	# new_trgtArray <- array(NA, dim=c(nyear_proj,3,niters), 
	# 											 dimnames = list((final_yr_sims+1):(final_yr_sims+nyear_proj), c("min","val","max"),
	# 											 iter=1:niters))
	
	
	last.yr <- final_yr_sims
	message(paste(sp, "HR 4010",rep))
	
	## Extend stock for management period 
	stock_proj <- stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	future_harvest_iters <- val * rlnorm(nyear_proj * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd)
	# new_trgtArray[,"val",] <- future_harvest_iters
	
	## set the forward projection using Effort - Input HCR
	fwd.ctrl <- fwdControl(
		data.frame(year = (last.yr+1):(last.yr + nyear_proj),
							 #val = val,
							 val = future_harvest_iters,
							 quantity="f",
							 rel.year=last.yr) #,
		#trgtArray=new_trgtArray
	)
	
	stock_proj <- fwd(stock_proj,
										ctrl = fwd.ctrl,
										sr = sr,
										#sr.residuals = srres[1,ac(ts+(n)),1,1,1,1:niters],
										sr.residuals = srres[1,,1,1,1,rep],
										sr.residuals.mult = TRUE)
	
	stock_proj_HR4010 <- stock_proj
	
	#--------------------------------------------
	#Stepwise HCR - Catch TAC
	#--------------------------------------------
	
	stock_proj <- stock_pre
	val <- catch_sw_true
	
	last.yr <- final_yr_sims
	message(paste(sp, "TAC sw", rep))
	
	## Extend stock an additional year
	stock_proj <-stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	#fix implementation error - if HCR is 0, then how to add implementation error
	future_catch_iters <- val * rlnorm(nyear_proj * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd)   
	#new_trgtArray[,"val",] <- future_catch_iters
	
	# using TAC - output HCR
	fwd.ctrl <- fwdControl(
		data.frame(year = (last.yr+1):(last.yr + nyear_proj),
							 val = future_catch_iters,
							 #val = val, 
							 quantity = "catch") #,
		#trgtArray = new_trgtArray
	)
	
	stock_proj <- fwd(stock_proj,
										ctrl = fwd.ctrl,
										sr = sr,
										#sr.residuals = srres[1,ac(ts+(n)),1,1,1,1:niters], 
										sr.residuals = srres[1,,1,1,1,rep],
										sr.residuals.mult = TRUE)
	
	stock_proj_TACsw <- stock_proj
	
	#--------------------------------------------
	#Stepwise HCR - Harvest rate 
	#--------------------------------------------
	
	stock_proj <- stock_pre
	val <- hr_sw_true/harvest_rate_true_curr
	
	# new_trgtArray <- array(NA, dim=c(nyear_proj,3,niters), 
	# 											 dimnames = list((final_yr_sims+1):(final_yr_sims+nyear_proj), c("min","val","max"),
	# 											 iter=1:niters))
	
	last.yr <- final_yr_sims
	message(paste(sp, "HR sw",rep))
	
	## Extend stock for management period 
	stock_proj <- stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	future_harvest_iters <- val * rlnorm(nyear_proj * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd)
	# new_trgtArray[,"val",] <- future_harvest_iters
	
	## set the forward projection using Effort - Input HCR
	fwd.ctrl <- fwdControl(
		data.frame(year = (last.yr+1):(last.yr + nyear_proj),
							 #val = val,
							 val = future_harvest_iters,
							 quantity="f",
							 rel.year=last.yr) #,
		#trgtArray=new_trgtArray
	)
	
	stock_proj <- fwd(stock_proj,
										ctrl = fwd.ctrl,
										sr = sr,
										#sr.residuals = srres[1,ac(ts+(n)),1,1,1,1:niters],
										sr.residuals = srres[1,,1,1,1,rep],
										sr.residuals.mult = TRUE)
	
	stock_proj_HRsw <- stock_proj
	
	#--------------------------------------------
	#Business as usual scenario ####
	#--------------------------------------------
	##taken from effort dynamics function in "analysis/functions.R"
	message(paste(sp, "BAU", ed, rep))
	
	stock_proj <- stock_pre
	stock_proj <- stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	years_proj =	(final_yr_sims+1):(final_yr_sims+nyear_proj)
	years_proj_idx = (1+ts):(nyear_proj + ts)
	
	srres_bau = FLCore::expand(srres,year = (final_yr_sims-ts+1):(final_yr_sims+nyear_proj))
	
	# ED
	stock_proj <- switch(ed,
											 # effort dynamics
											 "ED03"=effortDynamics(stk = stock_proj, bmsy=Bmsy_true, iters=niters, years= years_proj_idx,
											 											#sr=list(model='bevholt', params=params(sims[[sp]]$brp)),
											 											sr=sr,
											 											xp=0.3, a = 0.5, srres=srres_bau[1,,1,1,1,rep], rsd = rsd),
											 # one way trip
											 "OW" =oneWayTrip(stk = stock_proj, fmax=true_refpts['crash', 'harvest']*0.8,
											 								 #sr=list(model='bevholt', params=params(sims[[sp]]$brp)), 
											 								 sr = sr,
											 								 years=years_proj, srres=srres_bau[1,,1,1,1,rep], nyear_proj, start_yr_f = ts, rsd = rsd)#,
											 #note: start_f_yrs = 20 = ts. 
											 # # roller coaster
											 # "RC1"=rollerCoaster(stk, fmax=true_refpts['crash', 'harvest']*0.8,
											 # 										fmsy=refpts(brp)['msy', 'harvest'], years=years_proj, up=0.1, down=0.05,
											 # 										sr=list(model='bevholt', params=params(sims[[sp]]$brp)), srres=srres_bau[1,,1,1,1,rep])
	)
	
	
	stock_proj_BAU <- stock_proj
	
	#--------------------------------------------
	#Fishing at Fmsy for yield comparison ####
	#--------------------------------------------
	##Rather than comparing performance of yield objective relative to MSY,
	# we project stock forward fishing at Fmsy (true optimal yield)
	# this is because most stocks are already below MSY before management starts
	# with no implementation error 
	
	stock_proj <- stock_pre
	val <- Hmsy_true
	
	last.yr <- final_yr_sims
	message(paste(sp, "Fmsy true",rep))
	
	## Extend stock for management period 
	stock_proj <- stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	#future_harvest_iters <- val * rlnorm(nyear_proj * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd)
	
	## set the forward projection using Effort - Input HCR
	fwd.ctrl <- fwdControl(
		data.frame(year = (last.yr+1):(last.yr + nyear_proj),
							 val = val,
							 quantity="f")
	)
	
	stock_proj <- fwd(stock_proj,
										ctrl = fwd.ctrl,
										sr = sr,
										sr.residuals = srres[1,,1,1,1,rep],
										sr.residuals.mult = TRUE)
	
	stock_proj_Fmsy <- stock_proj
	
	#compile results
	list(stock_proj_TAC4010 = stock_proj_TAC4010, stock_proj_HR4010 = stock_proj_HR4010,
			 stock_proj_TACsw =stock_proj_TACsw, stock_proj_HRsw =stock_proj_HRsw,
			 stock_proj_BAU=stock_proj_BAU, stock_proj_Fmsy = stock_proj_Fmsy, 
			 param_comp = param_comp)
}


# #--------------------------------------------
# #Plotting results ####
# #--------------------------------------------
#results.vec <- FLStocks('Effort 4010 HCR'=stock_proj_HR4010,'TAC 4010 HCR'=stock_proj_TAC4010, 'Effort sw HCR'=stock_proj_HRsw, 'TAC sw HCR'=stock_proj_TACsw, 'BAU' = stock_proj_BAU)
#results.vec <- FLStocks('Effort HCR'=mse_runs$stock_proj_HR,'TAC HCR'=mse_runs$stock_proj_TAC, "BAU"=mse_runs$stock_proj_BAU)
# 
# # 
# # pdf(paste0("plots/projected_stock","-",sp,"-TS",ts,"-rep",rep,".pdf"),height=4,width=6)
# plot(results.vec) + geom_vline(xintercept = final_yr_sims, linetype = "dashed")
# # dev.off()
# # 
