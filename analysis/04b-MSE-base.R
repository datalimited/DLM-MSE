#---------------------------------------------
##4. Conduct MSE ####
#---------------------------------------------
#Projects forward using TAC HCR, HR HCR or no management
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
# rep = j
# niters = 1
# nyear_proj = 20
# rsd = 0.1
# ed = "OW"
# ur = 0
# buffer <- 0.75 #or buffer of 80% Hmsy and 80% msy based on buffer used for
# filename <- paste0(sp,"_TS",ts,"_",ed,"_rep",j)
# load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
# sr =list(model='bevholt', params=params(sims[[sp]]$brp))
# stock_full <- FLCore::iter(sims[[1]][["stock"]],j)
# catchE_est <- FLCore::iter(sims[[1]]$catchE,j) * (1-ur)
# stock_preassess <- readRDS(file=paste0("data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
# srres <- srres_all[[sp]]

# ##Estimated ratios of TBmsy/TB0 from Thorson et al. 2012 Spawning biomass reference points paper CJFAS - estimates for each Order
# Bmsy_B0_ratio_proxy = data.frame(rep(c('BC','SA','PS','SJ','BF','CR')),
# 											 c(0.572, 0.413, 0.410, 0.311, 0.311, 0.311))
# 
# colnames(Bmsy_B0_ratio_proxy)<-c('Species', 'Bmsy/B0_ratio')


#### Function ####
HCR_projection <- function(stock_full, stock_preassess, catchE_est, sp, ed, ts, rep, buffer, nyear_proj, niters, rsd, srres, sr) {
	
	#load estimated b.bmsy data from DL preassessment
	bbmsy_pre_all <- stock_preassess$bbmsy_pre_COM
	#cmsy_pre_fit <- stock_preassess$cmsy_pre_fit
	#comsir_pre_fit <- stock_preassess$comsir_pre_fit
	
	#set stock to work with
	#stock_full <- FLCore::iter(sims[[1]][["stock"]],rep)
	final_yr_sims <- FLCore::dims(stock_full)$maxyear
	stock_pre <- FLCore::trim(stock_full, year=((length(stock_full@catch)-(ts-1)):length(stock_full@catch)))

	#cmsy_pre_df <- filter(bbmsy_pre_all, method == "CMSY")
	ens_pre_df <- filter(bbmsy_pre_all, method == "Ensemble")
	#true_pre_df <- filter(bbmsy_pre_all, method == "TRUE") # this is if we rerun 3-dlm-assessment.R code again. there was an error in how we calculated Bmsy. if we rerun this code, the error is now fixed, and we can use this true_pre_df again. 
	# but for now, use this code below (lines 61-74)
	
	#sr <- list(model='bevholt', params = FLCore::params(sims[[1]]$brp)) 
	true_refpts <- FLBRP::brp(FLBRP::FLBRP(stock_full, sr = sr))@refpts
	
	bmsy_true <- true_refpts["msy", "biomass"]
	
	## true b/bmsy 
	#use bbmsy_pre_all dataframe instead of stock_pre

	true_df_biomass <- FLCore::as.data.frame(stock(stock_pre)) %>%
		dplyr::select(year, data) %>%
		dplyr::rename(biomass = data)
	
	true_pre_df <- true_df_biomass %>%
		mutate("BBmsy" = biomass / bmsy_true ) %>%
		dplyr::select(Year = year, BBmsy)
	true_pre_df$method <- "TRUE"
	
	# ggplot(bbmsy_pre_all, aes(Year, BBmsy, colour = method)) +
	# 	geom_line()
	
	#set.seed(1234)
	
	#------------------
	## Calculated true parameters ####
	#------------------
	#HCRs for Simulated - TRUE data

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
	## Calculated estimated parameters ####
	#------------------
	#calculate current harvest rate and biomass, K, and harvest rate and biomass at MSY
	# using average mean of cmsy and comsir estimates. 
	
	Hmsy_est <- mean(c(stock_preassess$com_pre_variables$cmsy_Hmsy_est,stock_preassess$com_pre_variables$comsir_Hmsy_est), na.rm = TRUE)
	Ymsy_est <- mean(c(stock_preassess$com_pre_variables$cmsy_Ymsy_est,stock_preassess$com_pre_variables$comsir_Ymsy_est), na.rm = TRUE)
	K_est <-  mean(c(stock_preassess$com_pre_variables$cmsy_K_est,stock_preassess$com_pre_variables$comsir_K_est), na.rm= TRUE)
	##add a stop or move on function here if Hmsy or Ymsy are NAs
	
	#BBmsy_curr_est - estimated as the last ensemble estimate
	# the last value is actually the mean over the past 5 years, (hence the 2 NAs)
	# as the ensemble uses the mean of the last 5 years for each model to calculate its bbmsy estimate
	BBmsy_curr_est <- ens_pre_df$BBmsy[(length(ens_pre_df$BBmsy)-2)]
	
	# calculate current harvest rate estimated - divide all catch.sims by biomass/year
	# catch observation error included here, dependng on DL assessment in 03-stock-assess-dl
	#catchE_est <- FLCore::iter(sims[[1]]$catchE,rep)
	catchE_est <- FLCore::trim(catchE_est,year = (length(catchE_est)-(ts-1)):length(catchE_est))
	
	#calculate current hr estimated - mean from cmsy and comsir
	#cmsy current harvest rate
	cmsy_hr_est <- stock_preassess$com_pre_variables$cmsy_hr_est
	cmsy_hr_est_curr <- if(!is.na(cmsy_hr_est[1])) median(cmsy_hr_est[(length(cmsy_hr_est)-4):length(cmsy_hr_est)],na.rm = TRUE) else NA
	
	# #comsir current harvest rate calculation
	comsir_hr_est <- stock_preassess$com_pre_variables$comsir_hr_est
	comsir_hr_est_curr <- if(!is.na(comsir_hr_est[1])) median(comsir_hr_est[(length(comsir_hr_est)-4):length(comsir_hr_est)],na.rm = TRUE) else NA
	
	harvest_rate_est_curr <- mean(c(cmsy_hr_est_curr, comsir_hr_est_curr), na.rm = TRUE)
	
	#------------------
	## 40-10 HCR set up ####
	#------------------
	# HCR based on harvest rate (effort) or based on catch (output - TAC) 
	#buffer <- 0.75 #buffer set above: 75% or 80% Hmsy and msy based on buffer used for 
	# setting ABC acceptable biological catch from OFL overfishing limit.
	
	##Harvest control rule function
	## B 40 (Btrigger) and B10 (Blimit)
	BBmsy_trig <- 0.8 #40% K or 80% BMSY
	BBmsy_lim <- 0.2 	#10% K or 20% BMSY
	
	#Calculating 4010 HCRs - TRUE
	hr_4010_true <- HCR_4010(yaxis_msy = Hmsy_true*buffer, HCR_trig = BBmsy_trig, HCR_lim = BBmsy_lim,
													 HCR_curr = BBmsy_curr_true)
	
	catch_4010_true <- HCR_4010(yaxis_msy = Ymsy_true*buffer, HCR_trig = BBmsy_trig, HCR_lim = BBmsy_lim,
															HCR_curr = BBmsy_curr_true)
	
	
	#Calculating 4010HCRs for ensemble COM estimated data
	
	hr_4010_est <- HCR_4010(yaxis_msy = Hmsy_est*buffer, HCR_trig = BBmsy_trig, HCR_lim = BBmsy_lim,
													HCR_curr = BBmsy_curr_est)
	
	catch_4010_est <- HCR_4010(yaxis_msy = Ymsy_est*buffer, HCR_trig = BBmsy_trig, HCR_lim = BBmsy_lim,
														 HCR_curr = BBmsy_curr_est)
	
	
	##Using U (catch/biomass) and Baranov equation to solve for F
	# solveBarF <-  function (F, M, U )
	# {
	# 	Y <- - F - M - log ( 1 - (F+M)*U/F )
	# 	Y
	# }
	# 
	# Fseq <- seq ( from = 0.2, to = 10, length = 10 )
	# Y <- solveBarF ( F = Fseq, M = 0.5, U = 0.2 )
	# 
	# uniroot ( f = solveBarF, interval = c(0.2,0.73),
	# 					M = 0.5, U = 0.2 )
	
	#or we should just convert U into F. 
	# U = 1-exp(-F)or 
	# F = -ln (1-U)
	
	
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
	
	
	#Calculating stepwise HCRs for ensemble COM estimated data
	hr_sw_est <- HCR_stepwise(yaxis_msy = Hmsy_est*buffer, HCR_trig_sw = BBmsy_trig_sw, HCR_lim_sw = BBmsy_lim_sw,
														HCR_curr = BBmsy_curr_est)
	
	catch_sw_est <- HCR_stepwise(yaxis_msy = Ymsy_est*buffer, HCR_trig_sw = BBmsy_trig_sw, HCR_lim_sw = BBmsy_lim_sw,
															 HCR_curr = BBmsy_curr_est)
	
	
	#------------------
	#check estimates ####
	#----- -------------
	
	param_comp <- data.frame(Model = c("estimated", "true"),
													 # B_curr = c(B_curr_est, B_curr_true), 
													 Bmsy = c(NA, Bmsy_true), 
													 BBmsy_current = c(BBmsy_curr_est,BBmsy_curr_true),
													 K = c(K_est, K_true), 
													 Fmsy= c(Hmsy_est, Hmsy_true),
													 Ymsy=c(Ymsy_est,Ymsy_true),
													 HR_curr = c(harvest_rate_est_curr,harvest_rate_true_curr),
													 TAC_4010 = c(catch_4010_est,catch_4010_true),
													 HR_4010 = c(hr_4010_est,hr_4010_true),
													 TAC_sw = c(catch_sw_est,catch_sw_true),
													 HR_sw = c(hr_sw_est, hr_sw_true))
	param_comp
	
	HCR_df <- data.frame(Model = c("estimated","estimated","true", "true","estimated","estimated","true", "true"), 
											 HCR = c("HR_4010","TAC_4010","HR_4010","TAC_4010","HR_sw", "TAC_sw","HR_sw","TAC_sw"),
											 val = c(hr_4010_est, catch_4010_est, hr_4010_true, catch_4010_true,
											 				hr_sw_est, catch_sw_est, hr_sw_true, catch_sw_true))
	HCR_df
	
	dl_model_ests <- data.frame(Metric = c("Hmsy", "Ymsy","HR_curr"), 
															CMSY = c(stock_preassess$com_pre_variables$cmsy_Hmsy_est, stock_preassess$com_pre_variables$cmsy_Ymsy_est,cmsy_hr_est_curr),
															COMSIR = c(stock_preassess$com_pre_variables$comsir_Hmsy_est, stock_preassess$com_pre_variables$comsir_Ymsy_est,comsir_hr_est_curr), 
															mean = c(Hmsy_est, Ymsy_est,harvest_rate_est_curr))
	dl_model_ests
	
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
	# srres=ar1rlnorm(rho=sims[[sp]]$rho, iters=niters, years=(final_yr_sims+1):(nyear_proj+final_yr_sims), margSD=sims[[sp]]$margSD)
	#srres now is created before looping to make sure all iterations are different (see 2b-mse-run.R code)
	
	#-------------------------------
	#FOR 40-10 Catch HCR (TAC) ####
	#-------------------------------
	
	stock_proj <- stock_pre
	val <- catch_4010_est
	
	last.yr <- final_yr_sims
	message(paste(sp, "TAC 4010", rep))
	
	## Extend stock an additional year
	stock_proj <- stf(stock_proj, nyears = nyear_proj)
	stock_proj <- propagate(stock_proj, niters)
	
	#fix implementation error - if HCR is 0, then how can we add implementation error?
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
	val <- hr_4010_est/harvest_rate_est_curr
	
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
	val <- catch_sw_est
	
	last.yr <- final_yr_sims
	message(paste(sp, "TAC sw", rep))
	
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
	
	stock_proj_TACsw <- stock_proj
	
	#--------------------------------------------
	#Stepwise HCR - Harvest rate 
	#--------------------------------------------
	
	stock_proj <- stock_pre
	val <- hr_sw_est/harvest_rate_est_curr
	
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
											 											sr = sr,
											 											xp=0.3, a = 0.5, srres=srres_bau[1,,1,1,1,rep], rsd = rsd),
											 # one way trip
											 "OW" =oneWayTrip(stk = stock_proj, fmax=true_refpts['crash', 'harvest']*0.8,
											 								 #sr=list(model='bevholt', params=params(sims[[sp]]$brp)), 
											 								 sr = sr,
											 								 years=years_proj, srres=srres_bau[1,,1,1,1,rep], nyear_proj, start_yr_f = ts, rsd=rsd)#,
											 #note: start_f_yrs = 20 = ts. 
											 # # roller coaster
											 # "RC1"=rollerCoaster(stk, fmax=true_refpts['crash', 'harvest']*0.8,
											 # 										fmsy=refpts(brp)['msy', 'harvest'], years=years_proj, up=0.1, down=0.05,
											 # 										sr=list(model='bevholt', params=params(sims[[sp]]$brp)), srres=srres_bau[1,,1,1,1,rep])
	)
	
	message(paste(sp, "BAU", ed, rep))
	
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
#mse_runs <-
list(stock_proj_TAC4010 = stock_proj_TAC4010, stock_proj_HR4010 = stock_proj_HR4010,
			 stock_proj_TACsw =stock_proj_TACsw, stock_proj_HRsw =stock_proj_HRsw, stock_proj_Fmsy = stock_proj_Fmsy,
			 stock_proj_BAU=stock_proj_BAU, HCR_compar = HCR_df, param_comp = param_comp, 
		   dl_model_ests = dl_model_ests, true_refpts = true_refpts)
}


# #--------------------------------------------
####Spare code
# #--------------------------------------------

##Catch HCR
# When looping projections through each year - not needed
# for(n in 1:nyear_proj){
# 	last.yr <- dims(stock_proj)$maxyear
# 	message(paste(sp, "TAC Working on projected year:", n))
# 	
# 	## Extend stock an additional year
# 	stock_proj <-stf(stock_proj, nyears = 1)
# 	stock_proj <- propagate(stock_proj, niters)
# 	
# 	#fix implementation error - if HCR is 0, then how to add implementation error
# 	future_catch_iters <- val * rlnorm(1 * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd)   
# 	new_trgtArray[,"val",] <- future_catch_iters
# 	
# 	# using TAC - output HCR
# 	fwd.ctrl <- fwdControl(
# 		data.frame(year = last.yr + 1,
# 							 val = val, 
# 							 quantity = "catch"),
# 		trgtArray = new_trgtArray
# 	)
# 	
# 	stock_proj <- fwd(stock_proj,
# 										ctrl = fwd.ctrl,
# 										sr = sr,
# 										sr.residuals = srres[1,ac(ts+(n)),1,1,1,1:niters], 
# 										sr.residuals.mult = TRUE)
# }


##### For Harvest Rate HCR
# for(n in 1:nyear_proj){
#   last.yr <- dims(stock_proj)$maxyear
#   message(paste(sp,"HR Working on projected year:", n))
#   
#   ## Extend stock an additional year
#   stock_proj <- stf(stock_proj, nyears = 1)
#   stock_proj <- propagate(stock_proj, niters)
#   
#   future_harvest_iters <- val * rlnorm(1 * niters, meanlog = 0 - (rsd^2/2), sdlog=rsd) 
#   new_trgtArray[,"val",] <- future_harvest_iters
#   
#   ## set the forward projection using Effort - Input HCR
#   fwd.ctrl <- fwdControl(
#         data.frame(year = last.yr + 1,
#                   val = val,  
#                   quantity="f"),
#                   rel.year=dims(stock_proj)$maxyear,     #check if this is right data (origianlly sim.stock)
#         					rel.year=last.yr,    #or should this relative year remain as TS (i.e. 60)      
#         					trgtArray=new_trgtArray)
#         
#   
#   stock_proj <- fwd(stock_proj,
#                     ctrl = fwd.ctrl,
#                     sr = sr,
#                     sr.residuals = srres[1,ac(ts+(n)),1,1,1,1:niters],
#   									sr.residuals.mult = TRUE)
# }

# #--------------------------------------------
# #Plotting results ####
# #--------------------------------------------
#results.vec <- FLStocks('Effort 4010 HCR'=stock_proj_HR4010,'TAC 4010 HCR'=stock_proj_TAC4010, 'Effort sw HCR'=stock_proj_HRsw, 'TAC sw HCR'=stock_proj_TACsw, 'BAU' = stock_proj_BAU)
#results.vec <- FLStocks('Effort HCR'=mse_runs$stock_proj_HR,'TAC HCR'=mse_runs$stock_proj_TAC, "BAU"=mse_runs$stock_proj_BAU)
# 
# # 
# # pdf(paste0("plots/projected_stock","-",sp,"-TS",ts,"-rep",rep,".pdf"),height=4,width=6)
##plot(results.vec) + geom_vline(xintercept = final_yr_sims, linetype = "dashed")
# # dev.off()
# # 
