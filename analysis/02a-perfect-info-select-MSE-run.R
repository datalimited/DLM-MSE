# 2 Select Stock to conducts MSE for perfect information - fine tuning HCR buffers
# Authors Jessica Walsh, Sean Anderson and Coilin Minto
# began 03 Sept 2016
# major edits 08 Jan 2018


#ED = "OW"
reps <-1:600
buffer = 0.5 #or test buffer with loop in BUF used for setting ABC acceptable biological catch from OFL overfishing limit.

SP <- c("BC","SA","PS","SJ","CR")
#TS <- c(20, 60)
ED = c("ED03", "OW")
#BUF = c(0.3, 0.4, 0.5, 0.6, 0.70, 0.80, 0.90, 1.0)

rsd_impl = 0.1 
it = 1
iters = 1

ts = 20					#TS	#length of catch time series before first assessment
proj_yr = 20 		# num years to project forward
rea_year = 20		# time after first assessment until second assessment
length_sim = 60 # num years in simulations

set.seed(1234)
source("analysis/02b-perfect-info-MSE-base.R")

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

#### Stock preassessments ####
require(doParallel)
registerDoParallel(cores = 3)
getDoParWorkers()

start.time <- Sys.time()

#looping by sp
plyr::l_ply(SP, .paropts = list(.export=c("reps", "ED", "ts","buffer","proj_yr","iters", "rsd_impl",
																					"HCR_projection_true","HCR_4010","HCR_stepwise",
																					"effortDynamics","oneWayTrip","ar1rlnorm","srres_all"),
																.packages = c("FLCore","FLash","FLBRP","FLa4a",
																							"FLAssess", "dplyr")),
						.fun = function(sp) {
							
							#### Projecting forward ####
								srres <- srres_all[[sp]]
							
								# for (buffer in BUF) {
								buff_name <- buffer * 100
								
								for (ed in ED) {
									load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
									sr =list(model='bevholt', params=params(sims[[sp]]$brp))
								
									for (j in reps) {
										filename <- paste0(sp,"_TS",ts,"_",ed,"_rep",j)
										stock_full <- FLCore::iter(sims[[1]][["stock"]],j)
										preassess_stock <- readRDS(file=paste0("data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
										
										mse_runs <- HCR_projection_true(stock_full = stock_full, stock_preassess = preassess_stock, sp = sp,
																										ed = ed, ts = ts, rep = j, nyear_proj = proj_yr, buffer = buffer,
																										niters = iters, rsd = rsd_impl, srres=srres, sr = sr)
										
										filename_buf <- paste0(sp,"_TS",ts,"_",ed,"_rep",j,"_buf",buff_name)
										saveRDS(mse_runs, file=paste0("data-generated/mse_results/true/buffer",buff_name,"/", ed, "/mse_projection_",filename_buf,"_true.rds"))
										#saveRDS(mse_runs, file=paste0("~/scratch/data-generated/mse_results/buffer",buff_name,"/", ed, "/mse_projection_",filename_buf,"_true.rds")
									}
								}
						#	}
						}, .parallel = TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

beepr::beep()

