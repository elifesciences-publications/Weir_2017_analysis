## Weir_fig_4D_4F_S4B_S4C_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

#### MITOCHONDRIAL DATA (FIGURE 4D AND S4A) ####
# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_2/fig_2C_S2_data.csv')

# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$bestradiol <- split_fnames$X4
input_data$timepoint <- as.numeric(substr(split_fnames$X6,2,3))
input_data$timepoint[input_data$timepoint != 0] <- input_data$timepoint[input_data$timepoint != 0] - 1
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log(input_data$yfp_mean)
# each timepoint represents 15 minutes
input_data$time_in_hr <- input_data$timepoint*0.25
# calculate means for YFP-15 and no YFP for bgrd subtraction
summ <- ddply(.data = input_data, .variables = 'substrate', summarize,
              mean_log_yfp = mean(log_yfp))
# generate background-subtracted YFP values by subtracting mean YFP density at 
# mitos in cells lacking YFP
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# summarize data for model fitting
plot_summ <- ddply(.data = input_data, .variables = c('substrate','timepoint','bestradiol'),
                   summarize,
                   bsub_yfp_mean = mean(bsub_yfp),
                   bsub_yfp_sd = sd(bsub_yfp),
                   log_bsub_yfp_mean = mean(log_bsub_yfp, na.rm = TRUE),
                   log_bsub_yfp_sd = sd(log_bsub_yfp, na.rm = TRUE),
                   nobs = length(bsub_yfp),
                   nobs_below_0 = length(log_bsub_yfp[is.na(log_bsub_yfp)]),
                   nobs_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd]),
                   frac_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd])/length(log_yfp))
plot_summ$log_bsub_yfp_se <- plot_summ$log_bsub_yfp_sd/sqrt(plot_summ$nobs-plot_summ$nobs_below_0)
plot_summ$bestradiol <- as.character(plot_summ$bestradiol)
plot_summ$bestradiol[plot_summ$bestradiol == 'minusbest'] <- '-MSP1'
plot_summ$bestradiol[plot_summ$bestradiol == 'plusbest'] <- '+MSP1'
plot_summ$bestradiol <- factor(plot_summ$bestradiol, levels = c('-MSP1','+MSP1'))
## preparing data for modeling ##
summ_for_mod <- subset(plot_summ, timepoint > 0 & substrate == 'deltac30')
# MSP1 expression took 30 minutes to occur after imaging began.
# We therefore to began analysis starting at the 2nd timepoint.
summ_for_mod$timepoint <- summ_for_mod$timepoint - 1
summ_for_mod$time_in_hr <- summ_for_mod$timepoint*0.25
for (s in unique(summ_for_mod$substrate)){
  for (b in unique(summ_for_mod$bestradiol)){
    summ_for_mod$norm_log_yfp[summ_for_mod$bestradiol == b] <- summ_for_mod$log_bsub_yfp_mean[summ_for_mod$bestradiol == b] - summ_for_mod$log_bsub_yfp_mean[summ_for_mod$bestradiol == b & summ_for_mod$timepoint == 0] 
  }
}
# enter the function to be optimized into R. See Methods for reference details.
min.2state.log.fun <- function(data, param){
  with(data, sum((log((param[1]-param[3])*exp(-(param[1]+param[2])*time_in_hr)+param[2]*exp(-param[3]*time_in_hr))-log(param[1]+param[2]-param[3])-norm_log_yfp)^2))
}
log_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())
# use optim to estimate starting values to fit 2-state model to data

require(nls2)
# fit two-state model using brute force
nls.2state.plus.log <- nls2(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                            data = subset(summ_for_mod, bestradiol == '+MSP1' & time_in_hr <= 1),
                            start = data.frame(a = c(-1,10), b = c(-1,10), d = c(-1,10)),
                            algorithm = 'brute-force', control = nls.control(maxiter = 100000))

nls.2state.minus.log <-nls2(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                            data = subset(summ_for_mod, bestradiol == '-MSP1' & time_in_hr <= 1),
                            start = data.frame(a = c(-1,10), b = c(-1,10), d = c(-1,10)),
                            algorithm = 'brute-force', control = nls.control(maxiter = 100000))
# fit linear models
lm.plus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(summ_for_mod, bestradiol == '+MSP1' & time_in_hr <= 1))
lm.minus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(summ_for_mod, bestradiol == '-MSP1' & time_in_hr <= 1))


# calculating RSSs for table
RSS.2state.plus <- sum(resid(nls.2state.plus.log)^2)
RSS.2state.minus <- sum(resid(nls.2state.minus.log)^2)
RSS.1state.plus <- sum(resid(lm.plus.log)^2)
RSS.1state.minus <- sum(resid(lm.minus.log)^2)

# calculating AICs as described in paper referenced
AIC.2state.plus <- 2*3 + 5*log(sum(resid(nls.2state.plus.log)^2)/5) + 2*3*4
AIC.2state.minus <- 2*3 + 5*log(sum(resid(nls.2state.minus.log)^2)/5) + 2*3*4
AIC.1state.plus <- 2 + 5*log(sum(resid(lm.plus.log)^2)/5) + 2*1*2/3
AIC.1state.minus <- 2 + 5*log(sum(resid(lm.minus.log)^2)/5) + 2*1*2/3

# p-value for 2-state being better (1-pval is probability shown in figure - standard "p this is not true")
pval.plus <- exp((AIC.1state.plus-AIC.2state.plus)/2)/(1+exp((AIC.1state.plus-AIC.2state.plus)/2))
pval.minus <- exp((AIC.1state.minus-AIC.2state.minus)/2)/(1+exp((AIC.1state.minus-AIC.2state.minus)/2))

# extract parameters from linear models for reporting in supplemental table
summary(lm.plus.log)
summary(lm.minus.log)

# function for plotting 2-state fits
C.2state.log.fun <- function(time_in_hr, params){
  log((params$a-params$d)*exp(-(params$a+params$b)*time_in_hr)+params$b*exp(-params$d*time_in_hr))-log(params$a+params$b-params$d)
}


ggplot(subset(summ_for_mod, time_in_hr <= 1),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  scale_color_manual(values = c('black','red')) +
  geom_abline(intercept = 0, slope = coef(lm.minus.log), color = 'black') +
  geom_abline(intercept = 0, slope = coef(lm.plus.log), color = 'red') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.minus.log))),
                color = 'black', linetype = 2) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.plus.log))),
                color = 'red', linetype = 2) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed mitochondrial YFP density\nnormalized to t = 0')
#ggsave('/path/to/fig/4D/plot.pdf', useDingbats = FALSE, device = cairo_pdf(width = 6, height = 5))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# plotting the fraction of mitochondria with undetectable YFP for S4A
ggplot(subset(summ_for_mod, time_in_hr <=2 & time_in_hr >= 0),
       aes(x = time_in_hr, y = frac_below_bg, color = bestradiol)) +
  geom_line() +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x = 'Time after MSP1 expression (hr)',
       y = 'Fraction of mitochondria without detectable YFP-Pex15')
#ggsave('/path/to/fig/S4A/plot.pdf', device = cairo_pdf(width = 6, height = 5))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

#### PEROXISOMAL DATA (FIGURE 4F AND S4B) ####

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_4/fig_4A_4F_S4B_data.csv')


# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$replicate <- split_fnames$X4
input_data$bestradiol <- split_fnames$X5
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
input_data$timepoint <- input_data$timepoint - 1
# calculate YFP density by dividing total YFP intensity by particle volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# log-transform because data is roughly log-normal
input_data$log_yfp <- log(input_data$yfp_mean)
# each timepoint is 15 minutes
input_data$time_in_hr <- input_data$timepoint*0.25
# calculate mean background in non-YFP-containing strain for bgrd subtraction
summ <- ddply(.data = input_data, .variables ='substrate', summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
input_data$bsub_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
# log-transform
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# summarize data for model fitting and plotting
# because object size change during the experiment artifactually complicated analyses,
# we only imaged peroxisomes that fit the same size as the distribution at t=0.
lg_plot_summ <- ddply(.data = subset(input_data, volume >150), .variables = c('substrate','time_in_hr','bestradiol','replicate'),
                      summarize,
                      log_bsub_yfp_mean = mean(log_bsub_yfp, na.rm = TRUE),
                      log_bsub_yfp_sd = sd(log_bsub_yfp, na.rm = TRUE),
                      nobs = length(bsub_yfp),
                      nobs_below_0 = length(log_bsub_yfp[is.na(log_bsub_yfp)]),
                      nobs_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd]),
                      frac_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd])/length(log_yfp))
lg_plot_summ$log_bsub_yfp_se <- lg_plot_summ$log_bsub_yfp_sd/sqrt(lg_plot_summ$nobs-lg_plot_summ$nobs_below_0)
lg_plot_summ$bestradiol <- as.character(lg_plot_summ$bestradiol)
lg_plot_summ$bestradiol[lg_plot_summ$bestradiol == 'minusbest'] <- '-MSP1'
lg_plot_summ$bestradiol[lg_plot_summ$bestradiol == 'plusbest'] <- '+MSP1'
lg_plot_summ$bestradiol <- factor(lg_plot_summ$bestradiol, levels = c('-MSP1','+MSP1'))
lg_plot_summ$log_postmean <- log(lg_plot_summ$bsub_yfp_mean)
# as with mito deltaC30 data, took 30 mins for MSP1 expression to begin after starting imaging.
lg_plot_summ$time_in_hr <- lg_plot_summ$time_in_hr - 0.5
# normalize to t=0 for the replicate
lg_plot_summ$norm_log_yfp <- NA
for (i in unique(lg_plot_summ$replicate)){
  for (b in unique(lg_plot_summ$bestradiol)){
    lg_plot_summ$norm_log_yfp[lg_plot_summ$replicate == i & lg_plot_summ$substrate == 'yfp15' & lg_plot_summ$bestradiol == b] <- lg_plot_summ$log_bsub_yfp_mean[lg_plot_summ$replicate == i & lg_plot_summ$substrate == 'yfp15' & lg_plot_summ$bestradiol == b] -  lg_plot_summ$log_bsub_yfp_mean[lg_plot_summ$replicate == i & lg_plot_summ$substrate == 'yfp15' & lg_plot_summ$time_in_hr == 0 & lg_plot_summ$bestradiol == b]
  }
}
ggplot(subset(lg_plot_summ, time_in_hr <= 2 & time_in_hr >= 0 & substrate == 'yfp15'),
       aes(x = time_in_hr, y = norm_log_yfp, color = replicate, linetype = bestradiol)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed peroxisomal YFP density\nnormalized to replicate minus b-estradiol t = 0')
#ggsave('/path/to/fig/S4B/plot.pdf', device = cairo_pdf(width = 7, height = 5), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

## model fitting. See mito section above and methods for details.##

min.2state.log.fun <- function(data, param){
  with(data, sum((log((param[1]-param[3])*exp(-(param[1]+param[2])*time_in_hr)+param[2]*exp(-param[3]*time_in_hr))-log(param[1]+param[2]-param[3])-norm_log_yfp)^2))
}
min.2state.lin.fun <- function(data, param){
  with(data, sum((log(((param[1]-param[3])*exp(-(param[1]+param[2])*time_in_hr)+param[2]*exp(-param[3]*time_in_hr))/(param[1]+param[2]-param[3]))-log(norm_lin_yfp))^2))
}

log_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())
lin_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())

for (b in unique(lg_plot_summ$bestradiol)){
  last_tpt = 1
  log_model_estimates = optim(par = c(2, 0.4, .6), min.2state.log.fun,
                              data = subset(lg_plot_summ, substrate == 'yfp15' & replicate == 'rep1' & bestradiol == b & time_in_hr <= last_tpt & time_in_hr >= 0))
  log_par_estimates <- rbind(log_par_estimates, data.frame(bestradiol = b, 
                                                           a = log_model_estimates$par[1], b = log_model_estimates$par[2], d = log_model_estimates$par[3], last_tpt = last_tpt))
  
}

nls.2state.plus.log <- nls(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                           data = subset(lg_plot_summ, bestradiol == '+MSP1' & time_in_hr <= last_tpt & time_in_hr >= 0 & replicate == 'rep1' & substrate == 'yfp15'),
                           start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == '+MSP1'],
                                        b = log_par_estimates$b[log_par_estimates$bestradiol == '+MSP1'],
                                        d = log_par_estimates$d[log_par_estimates$bestradiol == '+MSP1']))

nls.2state.minus.log <- nls(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                            data = subset(lg_plot_summ, bestradiol == '-MSP1' & time_in_hr <= last_tpt & time_in_hr >= 0 & replicate == 'rep1' & substrate == 'yfp15'),
                            start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == '-MSP1'],
                                         b = log_par_estimates$b[log_par_estimates$bestradiol == '-MSP1'],
                                         d = log_par_estimates$d[log_par_estimates$bestradiol == '-MSP1']))
# adding brute force to better fit the non-linear model for -MSP1
require(nls2)
nls.2state.minus.log <- nls2(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                             data = subset(lg_plot_summ, bestradiol == '-MSP1' & time_in_hr <= last_tpt & time_in_hr >= 0 & replicate == 'rep1' & substrate == 'yfp15'),
                             start = data.frame(a = c(-1,5), b = c(-1,5), d = c(-1,5)),
                             algorithm = 'brute-force', control = nls.control(maxiter = 100000))

lm.plus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(lg_plot_summ, bestradiol == '+MSP1' & time_in_hr <= last_tpt & time_in_hr >=0 & replicate == 'rep1' & substrate == 'yfp15'))
lm.minus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(lg_plot_summ, bestradiol == '-MSP1' & time_in_hr <= last_tpt & time_in_hr >= 0 & replicate == 'rep1' & substrate == 'yfp15'))

# calculating RSSs
RSS.2state.plus <- sum(resid(nls.2state.plus.log)^2)
RSS.2state.minus <- sum(resid(nls.2state.minus.log)^2)
RSS.1state.plus <- sum(resid(lm.plus.log)^2)
RSS.1state.minus <- sum(resid(lm.minus.log)^2)

# calculating AICs
n_tpts <- last_tpt*4 + 1
AIC.2state.plus <- 2*3 + n_tpts*log(sum(resid(nls.2state.plus.log)^2)/n_tpts) + 2*3*4/(n_tpts-4)
AIC.2state.minus <- 2*3 + n_tpts*log(sum(resid(nls.2state.minus.log)^2)/n_tpts) + 2*3*4/(n_tpts-4)
AIC.1state.plus <- 2 + n_tpts*log(sum(resid(lm.plus.log)^2)/n_tpts) + 2*1*2/(n_tpts-2)
AIC.1state.minus <- 2 + n_tpts*log(sum(resid(lm.minus.log)^2)/n_tpts) + 2*1*2/(n_tpts-2)

# p-value for 2-state being better
pval.plus <- 1/(1+exp((AIC.2state.plus-AIC.1state.plus)/2))
pval.minus <- exp((AIC.1state.minus-AIC.2state.minus)/2)/(1+exp((AIC.1state.minus-AIC.2state.minus)/2))

C.2state.log.fun <- function(time_in_hr, params){
  log((params$a-params$d)*exp(-(params$a+params$b)*time_in_hr)+params$b*exp(-params$d*time_in_hr))-log(params$a+params$b-params$d)
}

ggplot(subset(lg_plot_summ, substrate == 'yfp15' & replicate == 'rep1' & time_in_hr <= last_tpt & time_in_hr >=0),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  scale_color_manual(values = c('black','red')) +
  geom_abline(intercept = 0, slope = coef(lm.minus.log), color = 'black') +
  geom_abline(intercept = 0, slope = coef(lm.plus.log), color = 'red') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.minus.log))),
                color = 'black', linetype = 2) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.plus.log))),
                color = 'red', linetype = 2) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed peroxisomal YFP density\nnormalized to t = 0')
#ggsave('/path/to/fig/4F/plot.pdf', useDingbats = FALSE, device = cairo_pdf(width = 6, height = 5))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# generate a 2nd plot to use as an inset showing the diff with longer timepoints
nls.2state.plus.log.2hr <- nls(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                               data = subset(lg_plot_summ, bestradiol == '+MSP1' & time_in_hr <= 2 & time_in_hr >= 0 & replicate == 'rep1' & substrate == 'yfp15'),
                               start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == '+MSP1'],
                                            b = log_par_estimates$b[log_par_estimates$bestradiol == '+MSP1'],
                                            d = log_par_estimates$d[log_par_estimates$bestradiol == '+MSP1']))
nls.2state.minus.log.2hr <- nls2(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                                 data = subset(lg_plot_summ, bestradiol == '-MSP1' & time_in_hr <= 2 & time_in_hr >= 0 & replicate == 'rep1' & substrate == 'yfp15'),
                                 start = data.frame(a = c(-1,5), b = c(-1,5), d = c(-1,5)),
                                 algorithm = 'brute-force', control = nls.control(maxiter = 100000))

lm.plus.log.2hr <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(lg_plot_summ, bestradiol == '+MSP1' & time_in_hr <= 2 & time_in_hr >=0 & replicate == 'rep1' & substrate == 'yfp15'))
lm.minus.log.2hr <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(lg_plot_summ, bestradiol == '-MSP1' & time_in_hr <= 2 & time_in_hr >=0 & replicate == 'rep1' & substrate == 'yfp15'))


ggplot(subset(lg_plot_summ, substrate == 'yfp15' & replicate == 'rep1' & time_in_hr <= 2 & time_in_hr >=0),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  geom_abline(intercept = 0, slope = coef(lm.plus.log.2hr), color = 'red') +
  geom_abline(intercept = 0, slope = coef(lm.minus.log.2hr), color = 'black') +
  scale_x_continuous(breaks = c(0,1,2,3)) +
  scale_color_manual(values = c('black','red')) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.minus.log.2hr))),
                color = 'black', linetype = 2) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.plus.log.2hr))),
                color = 'red', linetype = 2) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed peroxisomal YFP density\nnormalized to t = 0')
#ggsave('/path/to/fig/4F/plot.pdf', useDingbats = FALSE, device = cairo_pdf(width = 6, height = 5))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
