#---
#title: "06-performance-metrics-UR"
#author: "Jessica Walsh"
#date: "4 July 2017"
#---

#setwd("~/MSE_project")

library(ggplot2)
library(FLash)
library(FLa4a)
library(FLAssess)
library(dplyr)
library(gridExtra)


sp_lookup = data.frame(c('BC','SA','PS','SJ','CR'),
											 c('Rockfish','Sardine','Sole',
											 	'Tuna', 'Corvina'))
colnames(sp_lookup)<-c('Symbol', 'Name')

species <- factor(levels = c("BC","PS","SJ","SA","CR"))

sce_lookup = data.frame(c("ED03","ED03","OW","OW"),
												c("5", "20", "5","20"),
												c("ED03_05yr", "ED03_20yr","OW_05yr", "OW_20yr"))
colnames(sce_lookup) <- c("ED", "REA_YR", "scenario")

scenario <- factor(c("ED03_05yr", "ED03_20yr","OW_05yr", "OW_20yr"))

##parameters
#iters = 1
it = 1

ts = 20					#TS	#length of catch time series before first assessment
#proj_yr = 20 		# num years to project forward
length_sim = 60 # num years in simulations

SP <- c("BC","PS","CR","SA","SJ")
ED <- c("ED03","OW")
REA_YR = c("5", "20")		# time after first assessment until second assessment
reps <- 1:600
#TS <- c(20, 60)
HCR <- c("HR4010","TAC4010","HRsw","TACsw","BAU")
#REA = c("base", "ens_plus")
BUF <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
UR = "yes"

sce <- readRDS(file="data-raw/sim_life_history_settings.rds")
UR_all = list()
for(lh in names(sce$LH)) {
	UR_temp = sce$LH[[lh]]$UR
	name <- paste(lh)
	UR_all[[name]] <- UR_temp
}
rm(UR_temp)

## make_metric_dataframe
# ##Create datafile that with all performance metrics across scenarios 
#last created 28/1/2018
#note this takes an hour or 2 to run

## need to redo the data compiling if MSE has been run again

require(doParallel)
registerDoParallel(cores = 3)
getDoParWorkers()

plyr::l_ply(SP, .paropts = list(.export=c("ED", "REA_YR", "UR", "reps", "HCR", "BUF",
																					"it", "ts", "length_sim", "UR_all"),
																.packages = c("FLCore", "dplyr", "gridExtra")),
						.fun = function(sp) {
							
							metrics_UR <- data.frame(matrix(vector(), 0, 19, dimnames=list(c(), c("sp", "ts", "ed", "ur", "j","rea_yr", "buffer", "ctrl","BBmsy_final", "prop_25bmsy", "prop_10bmsy","change_bbmsy", "FFmsy_final", "prop_abovefmsy","change_ffmsy", "catch_mean", "catch_relFmsy_mean", "rel_catch_tot", "catch_sd"))), stringsAsFactors=F)
							
							HCR_target_vals_UR <- data.frame(matrix(vector(),0, 17,dimnames = list(c(),c("filename_buf", "sp", "ts", "ed", "ur", "j","buffer", "Model","Bmsy","BBmsy_current","Fmsy","Ymsy","HR_curr","TAC_4010","HR_4010","TAC_sw", "HR_sw"))), stringsAsFactors=F)
							
							ur <- switch(UR,
													 "yes" = UR_all[[sp]],
													 "no" = 0)
							
							ref_pts <- list()
							for (ed in ED) {
								#for (sp in SP) {
								load(paste0("data-raw/sims-", sp,"-",ed,".RData"))
								
								for (buffer in BUF) {
									buff_name <- buffer * 100
									
									for (j in reps) {
										filename_buf <- paste0(sp,"_TS",ts,"_",ed, "_UR", ur*100, "_rep",j,"_buf",buff_name)
										mse_runs <- readRDS(file=paste0("data-generated/mse_results/UR/buffer",buff_name,"/",ed,"/mse_projection_",filename_buf,".rds"))
										
										##collect param_comps for each filename
										params <- data.frame(filename_buf, sp, ts, ed, ur, j, buffer, mse_runs$param_comp)
										
										Bmsy_sim <- dplyr::filter(params, Model == "true") %>% dplyr::select(Bmsy) %>% an()
										Fmsy_sim <- dplyr::filter(params, Model == "true") %>% dplyr::select(Fmsy) %>% an()
										Msy_sim <- dplyr::filter(params, Model == "true") %>% dplyr::select(Ymsy) %>% an()
										
										HCR_target_vals_UR <- bind_rows(HCR_target_vals_UR,params)
										
										
										## Collect each performance metric for all scenarios
										for (rea_yr in REA_YR) {
											
											mgmt_yrs = (length_sim+1):(length_sim+an(rea_yr))
											years = (length_sim-(ts-1)):(length_sim+an(rea_yr))
											
											for (ctrl in HCR) {
												
												stock_proj <- eval(parse(text = paste0("mse_runs$stock_proj_", ctrl)))
												stock_proj <- FLCore::iter(stock_proj,it) %>%
													FLCore::trim(year=years)
												
												
												## B/Bmsy status at end of management period (mean of last 3 years)
												bbmsy <- FLQuant(stock(stock_proj)/Bmsy_sim) %>%
													FLCore::as.data.frame() %>%
													dplyr::select(year, data)
												
												BBmsy_final <- mean(bbmsy$data[(length(bbmsy$data)-2):length(bbmsy$data)],na.rm=TRUE)
												
												## Prop of years that have a projected B/Bmsy below 25%Bmsy
												prop_25bmsy <- sum(bbmsy$data[bbmsy$year==mgmt_yrs] < 0.25) / an(rea_yr) * 100
												
												## Prop of years that have a projected B/Bmsy below 10%Bmsy
												prop_10bmsy <- sum(bbmsy$data[bbmsy$year==mgmt_yrs] < 0.10) / an(rea_yr) * 100
												
												##Change in B/Bmsy since start of management
												change_bbmsy <- bbmsy$data[bbmsy$year == max(years)] - bbmsy$data[bbmsy$year == (min(mgmt_yrs)-1)]  
												
												## F/Fmsy status at end of management period (mean of last 3 years)
												ffmsy <- FLQuant(fbar(stock_proj)/Fmsy_sim) %>%
													FLCore::as.data.frame() %>%
													dplyr::select(year, data)
												
												FFmsy_final <- mean(ffmsy$data[(length(ffmsy$data)-2):length(ffmsy$data)],na.rm=TRUE)
												
												## F/Fmsy status (proportion of years where F/Fmsy is above Fmsy)
												prop_abovefmsy <- sum(ffmsy$data[ffmsy$year==mgmt_yrs] > 1)/ an(rea_yr) * 100
												
												##Change in F/Fmsy since start of management 
												# (or should this be mean last 3 years - mean last 3 years before mgmt?)
												change_ffmsy <- ffmsy$data[ffmsy$year == max(years)]- ffmsy$data[ffmsy$year == (min(mgmt_yrs)-1)] 
												
												## Mean catch (across multiple iterations) per year in management period, relative to MSY
												catch_msy <- FLQuant(catch(stock_proj)/Msy_sim) %>%
													FLCore::as.data.frame() %>%
													dplyr::select(year, data)
												catch_mean <- mean(catch_msy$data[catch_msy$year==mgmt_yrs], na.rm=TRUE)
												
												## Mean catch per year in management period, relative to fishing at Fmsy
												#finding yield across management period if project forward with Fmsy
												catch_Fmsy <- catch(mse_runs$stock_proj_Fmsy) %>% 
													FLCore::trim(year=years) 
												
												catch_relFmsy <- FLQuant(catch(stock_proj)/catch_Fmsy) %>%
													FLCore::as.data.frame() %>%
													dplyr::select(year, data)
												catch_relFmsy_mean <- mean(catch_relFmsy$data[catch_relFmsy$year==mgmt_yrs], na.rm=TRUE)
												
												## Total catch across all years in management period (per repitition), relative to virgin biomass
													rel_catch_vB <- FLQuant(catch(stock_proj)/as.numeric(sims[[sp]]$lh["v"])) %>%
													FLCore::as.data.frame() %>%
												 dplyr::select(year, data)
												 rel_catch_tot <- sum(rel_catch_vB$data[rel_catch_vB$year == mgmt_yrs], na.rm=TRUE)
												
												## Annual variation in catches (standard deviation for each repition) over management period relative to msy
												catch_sd <- sd(catch_msy$data[catch_msy$year==mgmt_yrs],	na.rm =TRUE)
												
												out <- data.frame(sp, ts, ed, ur, j, rea_yr, buffer, ctrl, BBmsy_final, prop_25bmsy, prop_10bmsy, change_bbmsy, FFmsy_final, prop_abovefmsy, change_ffmsy, catch_mean, catch_relFmsy_mean, rel_catch_tot, catch_sd)
												
												metrics_UR <- bind_rows(metrics_UR,out)
												rm(out)
												
											}
										}
										message(filename_buf)
									}
								}
							}
							saveRDS(metrics_UR, file=paste0("data-generated/performance_metrics_dat_UR_", sp,".rds"))
							saveRDS(HCR_target_vals_UR, file=paste0("data-generated/HCR_target_vals_dat_UR_", sp, ".rds"))
						}, .parallel = TRUE)

metrics_BC <- readRDS(paste0("data-generated/performance_metrics_dat_UR_BC.rds"))
metrics_CR <- readRDS(paste0("data-generated/performance_metrics_dat_UR_CR.rds"))
metrics_SA <- readRDS(paste0("data-generated/performance_metrics_dat_UR_SA.rds"))
metrics_SJ <- readRDS(paste0("data-generated/performance_metrics_dat_UR_SJ.rds"))
metrics_PS <- readRDS(paste0("data-generated/performance_metrics_dat_UR_PS.rds"))

metrics_full <- bind_rows(metrics_BC, metrics_CR) %>%
	bind_rows(metrics_SA) %>% bind_rows(metrics_SJ) %>%
	bind_rows(metrics_PS)

saveRDS(metrics_full, paste0("data-generated/performance_metrics_dat_UR.rds"))
metrics <- metrics_full

HCR_target_BC <- readRDS(paste0("data-generated/HCR_target_vals_dat_UR_BC.rds"))
HCR_target_CR <- readRDS(paste0("data-generated/HCR_target_vals_dat_UR_CR.rds"))
HCR_target_SA <- readRDS(paste0("data-generated/HCR_target_vals_dat_UR_SA.rds"))
HCR_target_SJ <- readRDS(paste0("data-generated/HCR_target_vals_dat_UR_SJ.rds"))
HCR_target_PS <- readRDS(paste0("data-generated/HCR_target_vals_dat_UR_PS.rds"))

HCR_target_full <- bind_rows(HCR_target_BC, HCR_target_CR) %>%
	bind_rows(HCR_target_SA) %>% bind_rows(HCR_target_SJ) %>%
	bind_rows(HCR_target_PS)

saveRDS(HCR_target_full, paste0("data-generated/HCR_target_vals_dat_UR.rds"))
HCR_target_vals <- HCR_target_full

