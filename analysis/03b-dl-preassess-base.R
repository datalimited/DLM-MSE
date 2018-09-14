## PERFORM Data-limited ASSESSMENT
##code writen by Jessica Walsh, Coilin Minto, Sean Anderson

library(plyr)
library(dplyr)	
library(datalimited)
library(FLCore)
library(FLBRP)
library(ggplotFL)
library(fishensembles)
library(randomForest)
library(zoo)
library(splitstackshape)
library(reshape2)

# sp = "SJ"
# ts = 20
# rep = 2
# ed = "ED03"
# ur = 0
# load(paste0("data-raw/sims-", sp,"-",ed,".Rdata"))
# species_cat <- sims[[sp]]$species_cat
# resilience <- sims[[sp]]$resilience

stock_assess <- function(sims, sp, ts, rep, ur, resilience, species_cat) {
	
## has multiple iterations, extract the first one and for correct time series length
stock_full <- FLCore::iter(sims[[1]]$stock, rep)
stock <- FLCore::trim(stock_full, year=((length(stock_full@catch)-(ts-1)):length(stock_full@catch)))

## using true catch here 
#catch_sim <- catch(stock)

##using catch with observation error and bias
catchE <- FLCore::iter(sims[[sp]][["catchE"]], rep) * (1 -ur)
catch_sim <- trim(catchE, year=(length(stock_full@catch)-(ts-1)):length(stock_full@catch))

####Assess stock using 3 catch-only models####

##CMSY
cmsy_pre <- cmsy(yr = an(dimnames(catch_sim)$year), ct = an(catch_sim), reps = 3e5, start_r = resilience(resilience))
#saveRDS(cmsy_pre,file="../data-generated/pre-fits-cmsy.rds")

fake_data <- data.frame(bbmsy_q2.5 = NA, bbmsy_q25 = NA, bbmsy_q50 = NA, 
												bbmsy_q75 = NA, bbmsy_q97.5 = NA)

blank_data <- fake_data[rep(1:nrow(fake_data),times = ts),] %>%
	mutate(year= (length(stock_full@catch)-(ts-1)):length(stock_full@catch))

# blank_data <- splitstackshape:: expandRows(
# 	cbind(fake_data, count = ts), "count") %>%
# 	mutate(year= (length(stock_full@catch)-(ts-1)):length(stock_full@catch))

cmsy_bbmsy <- tryCatch({
 	bbmsy_out <- dplyr::select(cmsy_pre$bbmsy, -bbmsy_mean, -bbmsy_sd, -catch)
 	bbmsy_out}, error = function(e) blank_data)
cmsy_bbmsy$method <- "CMSY"

##COM-SIR
comsir_pre <- comsir(ct = an(catch_sim), yr =an(dimnames(catch_sim)$year) , nsim = 3e6,
										 start_r = resilience(resilience),
										 n_posterior = 5e3)

comsir_bbmsy <- tryCatch({
	if (comsir_pre$msd > 1) {
	bbmsy_out = blank_data
	bbmsy_out
} else {
	bbmsy_out <- dplyr::select(comsir_pre$bbmsy, -bbmsy_mean, -bbmsy_sd, -catch)
	bbmsy_out
}
	}, error = function(e) blank_data)
comsir_bbmsy$method <- "COMSIR"

#saveRDS(comsir_pre, file = "../data-generated/pre-fits-comsir.rds")

##mPRM 
#predicting bbmsy using RAM fits
ram <- dplyr::inner_join(datalimited::ram_ts, 
												 datalimited::spp_categories, by = "scientificname")
ram_prm_dat <- plyr::ddply(ram, "stockid", function(x) {
	format_prm(year = x$year, catch = x$catch, bbmsy = x$bbmsy_ram,
						 species_cat = x$spp_category[1L])
})
m <- fit_prm(ram_prm_dat)

#fitting mPRM to sims
mprm_sim_format <- format_prm(year = an(dimnames(catch_sim)$year), 
													catch = an(catch_sim), 
													bbmsy = rep(NA, length(an(catch_sim))),
					 species_cat = species_cat)

		mprm_pre <- predict_prm(mprm_sim_format, model = m, ci = TRUE)
		mprm_pre$year <- mprm_sim_format$year
		mprm_pre$method <- "mPRM"
#saveRDS(mprm_pre, file = "data-generated/bbmsy-pre-mprm.rds")

rm(mprm_sim_format)

####Fitting ensemble####

## make data frame of 3 catch only models responses
cmsy_pre_bbmsy <- dplyr:: select(cmsy_bbmsy, year, bbmsy_q50, method)  
comsir_pre_bbmsy <- dplyr::select(comsir_bbmsy, year, bbmsy_q50, method)
mprm_pre_bbmsy <- dplyr::select(mprm_pre, year, bbmsy_q50,method)

pre_COM_fit <- bind_rows(cmsy_pre_bbmsy, comsir_pre_bbmsy) %>%
	bind_rows(mprm_pre_bbmsy)

##Calculate Spectral frequencies
# currently, cutting stocks with fewer than 10 data points
cutoff_yr = 10
# using AR as smoother, empirical didn't seem to confer much more benefit
train_spec_mat <- function(x, freq_vec = 1/c(5, 20)) {
	if(length(x) >= cutoff_yr) {
		sp <- spec.ar(x/max(x), plot = FALSE)
		# approximate at fixed frequencies - necessary as series of different length
		approx(x = sp$freq, y = sp$spec, xout = freq_vec) %>% as.data.frame
	} else {
		data.frame(x=NA, y=NA, xout=NA)
	}
}

cdat <- FLCore::as.data.frame(catch_sim) %>%
	dplyr::select(year, data) %>%
	dplyr::rename(catch = data)

spec <- train_spec_mat(cdat$catch) %>%
	dplyr::rename(spec_freq = x, spec_dens = y) %>%
	as.data.frame()

nrow(spec)
spec <- na.omit(spec)
nrow(spec)

spec_wide <- spec %>%
	mutate(spec_freq = paste0("spec_freq_", spec_freq), stock = sp) %>%
	reshape2::dcast(stock ~ spec_freq, value.var = "spec_dens")

##Calculate mean bbmsy over 5 years based on DL estimates
# mean_bbmsy <- function(dat, years_window = 5L) {
# 	# chunk of data must have columns: b_bmsy_true, b_bmsy_est
# 	# message(paste(unique(dat$stockid), unique(dat$iter), sep = "-"))
# 		if (nrow(dat) > years_window) { 
# 			.n <- nrow(dat)
# 			i <- seq(.n-(years_window-1), .n)
# 			bbmsy_true_mean = mean(dat$b_bmsy_true[i])
# 			bbmsy_est_mean = mean(dat$b_bmsy_est[i])
# 			ytrue <- dat$b_bmsy_true[i]
# 			yest <- dat$b_bmsy_est[i]
# 			data.frame(bbmsy_true_mean, bbmsy_est_mean)
# 		}
# }

pre_COM_fits <- arrange(pre_COM_fit, method, year) %>%
	dplyr::mutate(b_bmsy_true =  as.character(NA)) %>%
	dplyr::rename(b_bmsy_est = bbmsy_q50)
pre_COM_fits$stock <- eval(sp) 
pre_COM_fits <-	arrange(pre_COM_fits, stock, method, year)

##for only end 5 years
#library("doParallel")
#registerDoParallel(cores = 4)

# for predicting bbmsy using ensemble for mean of past 5 years
# pre_COM_sum <- plyr::ddply(pre_COM_fits, c("method","stock"), .parallel = TRUE, .fun = mean_bbmsy)
# #why does it get warning? is this important? 
# pre_COM_sum$bbmsy_true_mean <- NULL
#  
# pre_COM_sum <- reshape2::dcast(pre_COM_sum, stock ~ method, value.var = "bbmsy_est_mean")
# pre_COM_sum <- inner_join(pre_COM_sum, spec_wide) %>%
# 	rename(Costello = mPRM)

##creating running means across 5 years
#all time series fit with ensemble.
pre_COM_yr <- reshape2::dcast(pre_COM_fits, stock+year ~method,value.var = "b_bmsy_est")
pre_COM_yr <- left_join(pre_COM_yr, spec_wide) %>%
		dplyr::rename(Costello = mPRM)

myrollmean <- function(x){
	c(NA, NA, rollapply(x, 5, mean, na.rm= TRUE), NA, NA)
}

pre_COM_sum_yr <- group_by(pre_COM_yr, stock) %>%
	arrange(year) %>%
	mutate(CMSY = myrollmean(CMSY),
				 COMSIR = myrollmean(COMSIR),
				 Costello = myrollmean(Costello)) %>%
	as.data.frame()

pre_COM_sum_yr$COMSIR[pre_COM_sum_yr$COMSIR == "NaN"] = NA  
pre_COM_sum_yr$CMSY[pre_COM_sum_yr$CMSY == "NaN"] = NA  

#### Training ensemble model to simulation data ####

#fitting pre-assessment ensemble model to case studies
if (!file.exists("data-generated/ensemble_training_model.rds")) {

ensemble_m <- fishensembles::make(ntree=1000, cores=2, formula = log(bbmsy_true_mean)~
                                    CMSY + COMSIR + Costello + spec_freq_0.05 + spec_freq_0.2)

ensemble_m_cmsy <- fishensembles::make(ntree=1000, cores=2, formula = log(bbmsy_true_mean)~
                                         COMSIR + Costello + spec_freq_0.05 + spec_freq_0.2)

ensemble_m_comsir <- fishensembles::make(ntree=1000, cores=2, formula = log(bbmsy_true_mean)~
                                           CMSY + Costello + spec_freq_0.05 + spec_freq_0.2)

ensemble_m_costello <- fishensembles::make(ntree=1000, cores=2, formula = log(bbmsy_true_mean)~
                                             CMSY + COMSIR +  spec_freq_0.05 + spec_freq_0.2)

ensemble_training <- list(ensemble_m = ensemble_m, ensemble_m_cmsy = ensemble_m_cmsy,
                            ensemble_m_comsir = ensemble_m_comsir,
                          ensemble_m_costello = ensemble_m_costello)
saveRDS(ensemble_training, file = "data-generated/ensemble_training_model.rds")
} else {
ensemble_training <- readRDS("data-generated/ensemble_training_model.rds")
  }

##Fit ensemble to our new data
# accounting for species where some models didn't converge

#Ensemble for all models
if(!is.na(pre_COM_sum_yr$COMSIR[3]) & !is.na(pre_COM_sum_yr$CMSY[3]) & !is.na(pre_COM_sum_yr$Costello[3]) ) {
	 pre_COM_sum_yr$ensemble <- exp(predict(ensemble_training$ensemble_m$model, newdata=pre_COM_sum_yr))
		pre_COM_sum_yr$ens_type = "all"
} else { 
	#ensemble for data missing CMSY
	if (is.na(pre_COM_sum_yr$CMSY[3]) & !is.na(pre_COM_sum_yr$COMSIR[3]) & !is.na(pre_COM_sum_yr$Costello[3])) { 
	pre_COM_sum_yr$ensemble <- exp(predict(ensemble_training$ensemble_m_cmsy$model, newdata=dplyr::select(pre_COM_sum_yr,-CMSY)))
	 pre_COM_sum_yr$ens_type = "nocmsy"
	} 
#ensemble for data missing COMSIR
	if (is.na(pre_COM_sum_yr$COMSIR[3]) & !is.na(pre_COM_sum_yr$CMSY[3]) & !is.na(pre_COM_sum_yr$Costello[3])) {
	pre_COM_sum_yr$ensemble <- exp(predict(ensemble_training$ensemble_m_comsir$model, newdata=dplyr::select(pre_COM_sum_yr,-COMSIR)))
	 pre_COM_sum_yr$ens_type = "nocomsir"
	} else {
		pre_COM_sum_yr$ensemble <- NA
		pre_COM_sum_yr$ens_type <- NA
	}
}
#saveRDS(pre_COM_sum_yr,file = "../data-generated/bbmsy_pre_fits_ens.rds")

#### -----------------------------------------------------
#### extract the true and estimated b/bmsy for plotting 
#### -----------------------------------------------------

## true bmsy
# refpts.sim <- sims[[1]]$brp@refpts #originally used this one, 
# but brp is the original FLBRP object before simulations, so not accurate for each specific iteration
sr <- list(model='bevholt', params = FLCore::params(sims[[1]]$brp)) 
refpts.sim <- FLBRP::brp(FLBRP::FLBRP(stock_full, sr = sr))@refpts

bmsy <- refpts.sim["msy", "biomass"]

## true b/bmsy 
true_df_biomass <- FLCore::as.data.frame(stock(stock)) %>%
	dplyr::select(year, data) %>%
	dplyr::rename(biomass = data)

true_df <- true_df_biomass %>%
	mutate("BBmsy" = biomass / bmsy ) %>%
	dplyr::select(Year = year, BBmsy)
true_df$method <- "TRUE"

#estimated bbmsy
est_ens <- data.frame(Year = pre_COM_sum_yr$year, 
											BBmsy = pre_COM_sum_yr$ensemble)
est_ens$method <- "Ensemble"

#joining estimated and true bbmsy to data-limited bbmsy estimates
pre_COM_fit <- dplyr::rename(pre_COM_fit, BBmsy = bbmsy_q50, Year = year) %>%
	bind_rows(true_df) %>%
	bind_rows(est_ens)

#saveRDS(pre_COM_fit, file="../data-generated/bbmsy-pre-COM.rds")

#### -----------------------------------------------------
#### extract parameters from  needed for next step - MSE
#### -----------------------------------------------------
#------------------
## Calculated estimated parameters ####
#------------------
#calculate current harvest rate and biomass, K, and harvest rate and biomass at MSY
# using average mean of cmsy and comsir estimates. 

if (!is.na(cmsy_bbmsy$bbmsy_q50[3])) {
	cmsy_Hmsy_est <- median(cmsy_pre$theta$r[cmsy_pre$theta$ell == 1]/2)
	cmsy_Ymsy_est <- median(cmsy_pre$msy[cmsy_pre$theta$ell == 1])
	cmsy_K_est <- median(cmsy_pre$theta$k[cmsy_pre$theta$ell == 1])
	cmsy_biomass <- apply(cmsy_pre$biomass, 2, median)
	cmsy_bmsy <- median(cmsy_pre$bmsy[cmsy_pre$theta$ell == 1])
	cmsy_hr_est <- cmsy_pre$bbmsy$catch / cmsy_biomass[-length(cmsy_biomass)]
	} else {
	cmsy_Hmsy_est <- NA
	cmsy_Ymsy_est <- NA
	cmsy_bmsy <- NA
	cmsy_hr_est <- NA
	cmsy_k_est <- NA
}

if(!is.na(comsir_bbmsy$bbmsy_q50[3])) {
	comsir_Hmsy_est <- median(comsir_pre$posterior$r/2)
	comsir_Ymsy_est <- median(comsir_pre$posterior$r*comsir_pre$posterior$k/4)
	comsir_K_est <- median(comsir_pre$posterior$k)
	comsir_output <- comsir_pre$quantities %>%
		arrange(yr) %>% dplyr::group_by(yr)
	comsir_summ <- dplyr::summarize(comsir_output, med_bio = median(predbio),med_ct = median(ct))
	comsir_hr_est <- comsir_summ$med_ct/comsir_summ$med_bio
	} else {
	comsir_Hmsy_est <- NA
	comsir_Ymsy_est <- NA
	comsir_K_est <- NA
	comsir_summ <- NA
	comsir_hr_est <- NA
}

com_estimated_values <- list(cmsy_Hmsy_est=cmsy_Hmsy_est, cmsy_Ymsy_est=cmsy_Ymsy_est, cmsy_K_est = cmsy_K_est,
														 cmsy_hr_est=cmsy_hr_est, cmsy_bmsy = cmsy_bmsy,
														 comsir_Hmsy_est=comsir_Hmsy_est, comsir_Ymsy_est=comsir_Ymsy_est, comsir_K_est = comsir_K_est,
														 comsir_summ=comsir_summ,comsir_hr_est=comsir_hr_est)


list(bbmsy_pre_COM = pre_COM_fit, cmsy_summary = cmsy_bbmsy, comsir_summary = comsir_bbmsy, 
		 mprm_summary = mprm_pre, com_pre_variables = com_estimated_values, ensemble = pre_COM_sum_yr$ens_type[[1]])

}


# ## quick plot 
# pdf(file=paste0("plots/preassess-vs-true-plot-",filename,".pdf", width = 4, height = 3))
# ggplot(pre_COM_fit, aes(x = Year, y = BBmsy)) +
#   geom_line(aes(colour = method)) +
#   geom_hline(yintercept = 1, linetype = "dashed")
# dev.off()

