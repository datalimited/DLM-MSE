# 2 Select Stock to conduct Stock preassessment
# Authors Jessica Walsh, Sean Anderson and Coilin Minto
# began 17 Aug 2015

## set scenario

reps <- 1:600
SP <- c("BC","SA","PS","SJ","CR")
#TS <- c(20, 60)
ED = c("ED03","OW")
#UR = c("no","yes")
UR = "yes"

# rsd_impl = 0.1
# it = 1
# iters = 1

ts = 20					#TS	#length of catch time series before first assessment
proj_yr = 20 		# num years to project forward
#rea_year = 20		# time after first assessment until second assessment
length_sim = 60 # num years in simulations

source("analysis/03b-dl-preassess-base.R")

set.seed(1234)

sce <- readRDS(file="data-raw/sim_life_history_settings.rds")
UR_all = list()
for(lh in names(sce$LH)) {
	UR_temp = sce$LH[[lh]]$UR
	name <- paste(lh)
	UR_all[[name]] <- UR_temp
}
rm(UR_temp)

#### Stock preassessments using Catch-only models ####
require(doParallel)
registerDoParallel(cores=3)
getDoParWorkers()

start.time <- Sys.time()

plyr::l_ply(SP, .paropts = list(.export=c("reps", "ED", "UR", "ts","stock_assess", "UR_all"),
																.packages = c("plyr","dplyr","datalimited","fishensembles",
																							"FLCore", "FLBRP", "randomForest","zoo",
																							"reshape2","splitstackshape")),
						.fun = function(sp) {
							
							for (ed in ED) {
								for (j in reps) {
									# for (sp in SP) {
									for (u in UR) {
										ur <- switch(u,
											"yes" = UR_all[[sp]],
											"no" = 0)
										
										filename <- paste0(sp,"_TS",ts,"_",ed,"_UR", ur*100, "_rep",j)
										message(filename)
										load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
										
										stock_preassess <- stock_assess(sims=sims, ts=ts, sp= sp, rep=j, ur= ur,
																										resilience = sims[[1]]$resilience,
																										species_cat = sims[[1]]$species_cat)
										
										
										filename_saves <- switch(as.character(u), 
										"no" = paste0("data-generated/preassess_results/", ed, "/stock_preassess_",filename,".rds"),
										"yes" = paste0("data-generated/preassess_results/UR/", ed, "/stock_preassess_",filename,".rds")
																						 )
											
										saveRDS(stock_preassess, file= filename_saves)
										#saveRDS(stock_preassess, file=paste0("~/scratch/data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
									}
									# }
								}
							}
						}, .parallel = TRUE)

end.time <- Sys.time()
end.time - start.time

# ##reps looping
# plyr::l_ply(reps, .paropts = list(.export=c("SP", "ED", "ts","stock_assess"),
# 																	.packages = c("plyr","dplyr","datalimited","fishensembles",
# 																								"FLCore","FLBRP","randomForest","zoo",
# 																								"reshape2","splitstackshape")),
# 						.fun = function(j) {
# 							
# 							for (ed in ED) {
# 								#	for (j in reps) {
# 								for (sp in SP) {
#										for (ur in UR) {
# 									filename <- paste0(sp,"_TS",ts,"_",ed,"_rep",j)
# 									message(filename)
# 									load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
# 									
# 									stock_preassess <- stock_assess(sims=sims, ts=ts, ur = ur, sp= sp, rep=j,
# 																									resilience = sims[[1]]$resilience,
# 																									species_cat = sims[[1]]$species_cat)
# 									
# 									saveRDS(stock_preassess, file=paste0("data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
# 									#saveRDS(stock_preassess, file=paste0("~/scratch/data-generated/preassess_results/",ed,"/stock_preassess_",filename,".rds"))
# 								}
#									}
# 								#	}
# 							}
# 						}, .parallel = TRUE)
