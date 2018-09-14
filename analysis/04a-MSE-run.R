# 4 Select Stock to conduct MSE

# Authors Jessica Walsh, Sean Anderson and Coilin Minto, Ernesto Jardim
# began 17 Aug 2015

## set scenario 
#SP <- "CR"			
#ED <- "OW"
reps <- 1:600
BUF = c(0.3,0.4)
UR <- "no"
	
SP <- c("BC","SA","PS","SJ","CR")
#TS <- c(20, 60)
ED <- c("ED03", "OW")
#BUF <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0)
#UR <- c("yes","no")

rsd_impl = 0.1 
it = 1
iters = 1

ts = 20					#TS	#length of catch time series before first assessment
proj_yr = 20 		# num years to project forward
rea_year = 20		# time after first assessment until second assessment
length_sim = 60 # num years in simulations

set.seed(1234)
source("analysis/04b-MSE-base.R")

#standardise SR.residual error across simulations 
#same across HCRs and OW/ED03, different across iterations and species
sce <- readRDS(file="data-raw/sim_life_history_settings.rds")
srres_all = list()
for(lh in names(sce$LH)) {
	srres_temp = ar1rlnorm(rho=sce$LH[[lh]]$rho, iters=1000, years=(length_sim+1):(proj_yr+length_sim), margSD=sce$LH[[lh]]$margSD)
	name <- paste(lh)
	srres_all[[name]] <- srres_temp
}
rm(srres_temp)

#add underreporting if UR = yes, at different levels for different species
UR_all = list()
for(lh in names(sce$LH)) {
	UR_temp = sce$LH[[lh]]$UR
	name <- paste(lh)
	UR_all[[name]] <- UR_temp
}
rm(UR_temp)

#### Stock projections ####
require(doParallel)
registerDoParallel(cores = 3)
getDoParWorkers()

# looping by Sp
# seems to be more efficient than looping by iteration
start.time <- Sys.time()

plyr::l_ply(SP, .paropts = list(.export=c("reps", "ED", "ts", "UR", "BUF", "proj_yr","iters", "rsd_impl",
																						"HCR_projection","HCR_4010","HCR_stepwise",
																						"effortDynamics","oneWayTrip","ar1rlnorm", "srres_all", "UR_all"),
																	.packages = c("FLCore","FLash","FLBRP","FLa4a",
																								"FLAssess", "dplyr")),
						.fun = function(sp) {
							
							#### Projecting forward ####
							srres <- srres_all[[sp]]
							
							for (buffer in BUF) {
							buff_name <- buffer * 100
							
							for (ed in ED) {
								#for (sp in SP) {
								for (u in UR) {
									if(u == "no") {
										load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
										sr =list(model='bevholt', params=params(sims[[sp]]$brp))
										
										for (j in reps) {
											filename <- paste0(sp,"_TS",ts,"_",ed,"_rep",j)
											stock_full <- FLCore::iter(sims[[1]][["stock"]],j)
											catchE_est <- FLCore::iter(sims[[1]]$catchE,j) 
											
											preassess_stock <- readRDS(file=paste0("data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
											
											mse_runs <- HCR_projection(stock_full = stock_full, catchE_est = catchE_est, stock_preassess = preassess_stock, 
																								 sp = sp, ed = ed, ts = ts, rep = j, nyear_proj = proj_yr, buffer = buffer,
																								 niters = iters, rsd = rsd_impl, srres = srres, sr = sr)
											
											filename_buf <- paste0(sp,"_TS",ts,"_",ed,"_rep",j,"_buf",buff_name)
											saveRDS(mse_runs, file=paste0("data-generated/mse_results/buffer",buff_name,"/", ed, "/mse_projection_",filename_buf,".rds"))
										}
									}
								#   if (u == "yes") {
								# 		ur <- switch(u,
								# 								 "yes" = UR_all[[sp]],
								# 								 "no" = 0)
								# 		
								# 			load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
								# 			sr =list(model='bevholt', params=params(sims[[sp]]$brp))
								# 			
								# 			for (j in reps) {
								# 			filename <- paste0(sp,"_TS",ts,"_",ed,"_UR", ur*100, "_rep",j)
								# 				stock_full <- FLCore::iter(sims[[1]][["stock"]],j)
								# 				catchE_est <- FLCore::iter(sims[[1]]$catchE,j) * (1-ur)
								# 				
								# 				preassess_stock <- readRDS(file=paste0("data-generated/preassess_results/UR/",ed,"/stock_preassess_",filename,".rds"))
								# 				
								# 				mse_runs <- HCR_projection(stock_full = stock_full, catchE_est = catchE_est, stock_preassess = preassess_stock, 
								# 																	 sp = sp, ed = ed, ts = ts, rep = j, nyear_proj = proj_yr, buffer = buffer,
								# 																	 niters = iters, rsd = rsd_impl, srres = srres, sr = sr)
								# 				
								# 				filename_buf <- paste0(sp,"_TS",ts,"_",ed, "_UR", ur*100, "_rep",j,"_buf",buff_name)
								# 				saveRDS(mse_runs, file=paste0("data-generated/mse_results/UR/buffer",buff_name,"/", ed, "/mse_projection_",filename_buf,".rds"))
								# 				}
								#   }   
								}     #UR
							# }			#sp
								}    #ed
							}     #buffer
						}, .parallel = TRUE)

beepr::beep()
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


# # #looping through iterations, ED and SP
# #slower
# plyr::l_ply(reps, .paropts = list(.export=c("SP", "ED", "ts","buffer", "proj_yr","iters", "rsd_impl",
# 																						"HCR_projection","HCR_4010","HCR_stepwise",
# 																						"effortDynamics","oneWayTrip","ar1rlnorm", "srres_all"),
# 																	.packages = c("FLCore","FLash","FLBRP","FLa4a",
# 																								"FLAssess", "dplyr")),
# 						.fun = function(j) {
# 							
# 							#### Projecting forward ####
# 							#for (buffer in BUF) {
# 								buff_name <- buffer * 100
# 								
# 								for (ed in ED) {
# 									for (sp in SP) {
# 										#for (j in reps) {
# 										load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
# 										filename <- paste0(sp,"_TS",ts,"_",ed,"_rep",j)
# 										srres <- srres_all[[sp]]
# 										sr =list(model='bevholt', params=params(sims[[sp]]$brp))
# 										stock_full <- FLCore::iter(sims[[1]][["stock"]],j)
# 										
# 										preassess_stock <- readRDS(file=paste0("data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
# 										
# 										srres <- srres_all[[sp]]
# 										
# 										mse_runs <- HCR_projection(stock_full = stock_full, stock_preassess = preassess_stock, sp = sp,
# 																							 ed = ed, ts = ts, rep = j, nyear_proj = proj_yr, buffer = buffer,
# 																							 niters = iters, rsd = rsd_impl, srres = srres, sr = sr)
# 										
# 										filename_buf <- paste0(sp,"_TS",ts,"_",ed,"_rep",j,"_buf",buff_name)
# 										saveRDS(mse_runs, file=paste0("data-generated/mse_results/buffer",buff_name,"/", ed, "/mse_projection_",filename_buf,".rds"))
# 									}
# 									#}
# 								}
# 							#}
# 						}, .parallel = TRUE)
# beepr::beep()
# 
