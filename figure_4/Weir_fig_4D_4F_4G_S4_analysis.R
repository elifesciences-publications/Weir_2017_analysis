## Weir_fig_4D_4F_S4_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
require(minpack.lm)

#### MITOCHONDRIAL DATA (FIGURE 4D AND S4A) ####
# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_2/fig_2C_4D_data.csv')

# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))

input_data$substrate <- split_fnames$X3
input_data$bestradiol <- split_fnames$X4
input_data$timepoint <- as.numeric(substr(split_fnames$X6,2,3))
input_data$timepoint <- input_data$timepoint - 1
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
# calculate background mean and sd for tabulation
ctrl_mean <- mean(input_data$log_yfp[input_data$substrate == 'noyfp'])
ctrl_sd <- sd(input_data$log_yfp[input_data$substrate == 'noyfp'])
# summarize data for model fitting
summ_for_plot <- ddply(.data=input_data, .variables=c('substrate','bestradiol','time_in_hr'),
                       summarize,
                       mean_log_yfp=mean(log_bsub_yfp, na.rm=TRUE),
                       sd_log_yfp = sd(log_bsub_yfp, na.rm=TRUE), 
                       nobs = length(log_bsub_yfp),
                       nobs_below_0 = length(log_bsub_yfp[is.na(log_bsub_yfp)]),
                       nobs_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd]),
                       frac_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd])/length(log_yfp))
summ_for_plot$se_log_yfp <- summ_for_plot$sd_log_yfp/sqrt(summ_for_plot$nobs)

## preparing data for modeling ##
plot_summ <- subset(summ_for_plot, substrate == 'deltac30')
for (j in unique(plot_summ$bestradiol)){
  plot_summ$norm_log_yfp[plot_summ$bestradiol == j] <- plot_summ$mean_log_yfp[plot_summ$bestradiol == j] -  plot_summ$mean_log_yfp[plot_summ$time_in_hr == 0 & plot_summ$bestradiol == j]
}

# enter the function to be optimized into R. See Methods for reference details.
min.2state.log.fun <- function(data, param){
  with(data, sum((log((param[1]-param[3])*exp(-(param[1]+param[2])*time_in_hr)+param[2]*exp(-param[3]*time_in_hr))-log(param[1]+param[2]-param[3])-norm_log_yfp)^2))
}
# use optim to estimate starting values to fit 2-state model to data
log_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())
for (b in unique(plot_summ$bestradiol)){
  last_tpt = 0.75
  log_model_estimates = optim(par = c(2, 1, 1), min.2state.log.fun,
                              data = subset(plot_summ, bestradiol == b & time_in_hr <= last_tpt & time_in_hr >= 0))
  log_par_estimates <- rbind(log_par_estimates, data.frame(bestradiol = b, 
                                                           a = log_model_estimates$par[1], b = log_model_estimates$par[2], d = log_model_estimates$par[3], last_tpt = last_tpt))
  
}

require(nls2)
# fit two-state models
nls.2state.plus.log <- nlsLM(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                             data = subset(plot_summ, bestradiol == 'plusbest' & time_in_hr <= last_tpt & time_in_hr >= 0),
                             start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == 'plusbest'],
                                          b = log_par_estimates$b[log_par_estimates$bestradiol == 'plusbest'],
                                          d = log_par_estimates$d[log_par_estimates$bestradiol == 'plusbest']))

nls.2state.minus.log <- nlsLM(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                              data = subset(plot_summ, bestradiol == 'minusbest' & time_in_hr <= last_tpt & time_in_hr >= 0),
                              start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == 'minusbest'],
                                           b = log_par_estimates$b[log_par_estimates$bestradiol == 'minusbest'],
                                           d = log_par_estimates$d[log_par_estimates$bestradiol == 'minusbest']))

# fit linear models
lm.plus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(plot_summ, bestradiol == 'plusbest' & time_in_hr <= last_tpt & time_in_hr >=0))
lm.minus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(plot_summ, bestradiol == 'minusbest' & time_in_hr <= last_tpt & time_in_hr >= 0))


# calculating RSSs for table
RSS.2state.plus <- sum(resid(nls.2state.plus.log)^2)
RSS.2state.minus <- sum(resid(nls.2state.minus.log)^2)
RSS.1state.plus <- sum(resid(lm.plus.log)^2)
RSS.1state.minus <- sum(resid(lm.minus.log)^2)

# integrate area between the 1-state and 2-state models for each sample
abs.diff.fun <- function(x, params, time_in_hr){
  abs(time_in_hr*x - (log((params$a-params$d)*exp(-(params$a+params$b)*x)+params$b*exp(-params$d*x))-log(params$a+params$b-params$d)))
}

abs.diff.minus <- integrate(abs.diff.fun, 0,0.75, params=as.list(coef(nls.2state.minus.log)), time_in_hr=coef(lm.minus.log))$value/0.75
abs.diff.plus <- integrate(abs.diff.fun, 0, 0.75, params=as.list(coef(nls.2state.plus.log)), time_in_hr=coef(lm.plus.log))$value/0.75

# these values are used for figure 4G plotting later in this file.

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

# plot.
ggplot(subset(plot_summ, time_in_hr <= last_tpt & time_in_hr >=0),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = norm_log_yfp-se_log_yfp, ymax = norm_log_yfp+se_log_yfp), width = 0.1) +
  scale_color_manual(values = c('black','red')) +
  geom_abline(intercept = 0, slope = coef(lm.minus.log), color = 'black') +
  geom_abline(intercept = 0, slope = coef(lm.plus.log), color = 'red') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)) +
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
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_3/fig_3C_4F_data.csv')


# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$promoter <- split_fnames$X4
input_data$bestradiol <- split_fnames$X5
input_data$position <- as.character(split_fnames$X6)
input_data$position <- as.numeric(input_data$position)
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
input_data$timepoint <- input_data$timepoint - 1
input_data$yfp_mean <- input_data$yfp/input_data$volume
# I accidentially mislabeled the negative controls (positions 13-15)
# as pmsp1 +bestradiol. Rectifying this below.
input_data$substrate <- as.character(input_data$substrate)
input_data$substrate[input_data$position %in% c(13,14,15)] <- 'noyfp'
input_data$promoter[input_data$position %in% c(13,14,15)] <- 'pzev'
input_data$bestradiol[input_data$position %in% c(13,14,15)] <- 'minusbest'
input_data$substrate <- factor(input_data$substrate)
# calculate YFP density by dividing total YFP intensity by particle volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# log-transform because data is roughly log-normal
input_data$log_yfp <- log(input_data$yfp_mean)
# each timepoint is 15 minutes
input_data$time_in_hr <- input_data$timepoint*0.25
# calculate mean background in non-YFP-containing strain for bgrd subtraction
summ <- ddply(.data = input_data, .variables = c('substrate', 'time_in_hr'), summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# there's no variation in background over time, so I'll combine all of the no yfp ctrls.
summ <- ddply(.data = input_data, .variables = c('substrate'), summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# perform background subtraction
input_data$bsub_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# calculate control mean and standard deviation for tabulation later
ctrl_mean <- mean(input_data$log_yfp[input_data$substrate == 'noyfp'])
ctrl_sd <- sd(input_data$log_yfp[input_data$substrate == 'noyfp'])
# log-transform
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# summarize data for model fitting and plotting
plot_summ <- ddply(.data = input_data, .variables = c('substrate','time_in_hr','bestradiol','promoter'),
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
plot_summ$log_postmean <- log(plot_summ$bsub_yfp_mean)
# took 30 mins for MSP1 expression to begin after starting imaging.
plot_summ$time_in_hr <- plot_summ$time_in_hr - 0.5
plot_summ$norm_log_yfp <- NA
# normalize to t=0 for the replicate
for (i in unique(plot_summ$promoter)){
  for (j in unique(plot_summ$bestradiol)){
    plot_summ$norm_log_yfp[plot_summ$promoter == i & plot_summ$substrate == 'yfp15' & plot_summ$bestradiol == j] <- plot_summ$log_bsub_yfp_mean[plot_summ$promoter == i & plot_summ$substrate == 'yfp15' & plot_summ$bestradiol == j] -  plot_summ$log_bsub_yfp_mean[plot_summ$promoter == i & plot_summ$substrate == 'yfp15' & plot_summ$time_in_hr == 0 & plot_summ$bestradiol == j]
  }
}

## model fitting. See mito section above and methods for details.##

min.2state.log.fun <- function(data, param){
  with(data, sum((log((param[1]-param[3])*exp(-(param[1]+param[2])*time_in_hr)+param[2]*exp(-param[3]*time_in_hr))-log(param[1]+param[2]-param[3])-norm_log_yfp)^2))
}
min.2state.lin.fun <- function(data, param){
  with(data, sum((log(((param[1]-param[3])*exp(-(param[1]+param[2])*time_in_hr)+param[2]*exp(-param[3]*time_in_hr))/(param[1]+param[2]-param[3]))-log(norm_lin_yfp))^2))
}

log_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())

for (b in unique(plot_summ$bestradiol)){
  last_tpt = 2
  log_model_estimates = optim(par = c(1, 1, 0.5), min.2state.log.fun,
                              data = subset(plot_summ, substrate == 'yfp15' & promoter == 'pzev' & bestradiol == b & time_in_hr <= last_tpt & time_in_hr >= 0))
  log_par_estimates <- rbind(log_par_estimates, data.frame(bestradiol = b, 
                                                           a = log_model_estimates$par[1], b = log_model_estimates$par[2], d = log_model_estimates$par[3], last_tpt = last_tpt))
  
}

nls.2state.plus.log.lm <- nlsLM(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                                data = subset(plot_summ, bestradiol == 'plusbest' & time_in_hr <= last_tpt & time_in_hr >= 0 & promoter == 'pzev' & substrate == 'yfp15'),
                                start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == 'plusbest'],
                                             b = log_par_estimates$b[log_par_estimates$bestradiol == 'plusbest'],
                                             d = log_par_estimates$d[log_par_estimates$bestradiol == 'plusbest']))

nls.2state.minus.log.lm <- nlsLM(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                                 data = subset(plot_summ, bestradiol == 'minusbest' & time_in_hr <= last_tpt & time_in_hr >= 0 & promoter == 'pzev' & substrate == 'yfp15'),
                                 start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == 'minusbest'],
                                              b = log_par_estimates$b[log_par_estimates$bestradiol == 'minusbest'],
                                              d = log_par_estimates$d[log_par_estimates$bestradiol == 'minusbest']))

lm.plus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(plot_summ, bestradiol == 'plusbest' & time_in_hr <= last_tpt & time_in_hr >=0 & promoter == 'pzev' & substrate == 'yfp15'))
lm.minus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(plot_summ, bestradiol == 'minusbest' & time_in_hr <= last_tpt & time_in_hr >= 0 &promoter == 'pzev' & substrate == 'yfp15'))

# absolute value of difference between the two curves
abs.diff.fun <- function(x, params, time_in_hr){
  abs(time_in_hr*x - (log((params$a-params$d)*exp(-(params$a+params$b)*x)+params$b*exp(-params$d*x))-log(params$a+params$b-params$d)))
}

abs.diff.minus <- integrate(abs.diff.fun, 0, 2, params=as.list(coef(nls.2state.minus.log.lm)), time_in_hr=coef(lm.minus.log))$value/2
abs.diff.plus <- integrate(abs.diff.fun, 0, 2, params=as.list(coef(nls.2state.plus.log.lm)), time_in_hr=coef(lm.plus.log))$value/2

# plotting the four values generated from abs.diff integration, here and in the deltac30 experiment.
# values:
# FL minusbest: 0.0503 per hr
# FL plusbest: 0.293 per hr
# deltac30 minusbest: 0.0296 per hr
# deltac30 plusbest: 0.0388 per hr

for_abs_plt = data.frame(substrate = c('FL','FL','deltac30','deltac30'),
                         bestradiol = c('minusbest','plusbest','minusbest','plusbest'),
                         magnitude = c(0.0503, 0.293, 0.0296, 0.0388))

ggplot(for_abs_plt, aes(x=interaction(bestradiol,substrate), y=magnitude)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.background = element_rect(fill=NA, color = NA),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.line = element_line(color='black'))
#ggsave('/path/to/fig/4G/plot.pdf', device=cairo_pdf(width=3, height=5), useDingbats=FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# generate a 2nd plot to use as an inset showing the diff with longer timepoints

ggplot(subset(plot_summ, substrate == 'yfp15' & promoter == 'pzev' & time_in_hr <= last_tpt & time_in_hr >=0),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  scale_color_manual(values = c('black','red')) +
  geom_abline(intercept = 0, slope =coef(lm.minus.log), color = 'black') +
  geom_abline(intercept = 0, slope = coef(lm.plus.log), color = 'red') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)) +
  stat_function(fun = C.2state.log.fun,
                args=list(params = as.list(coef(nls.2state.minus.log.lm))),
                color = 'black', linetype = 2) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.plus.log.lm))),
                color = 'red', linetype = 2) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed mitochondrial YFP density\nnormalized to t = 0')
#ggsave('/path/to/4F/plot.pdf', useDingbats = FALSE, device = cairo_pdf(width = 6, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.


# repeating fitting with shorter timecourse
log_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())
lin_par_estimates <- data.frame(bestradiol = character(), a = numeric(),b = numeric(), d = numeric(), last_tpt = numeric())

for (b in unique(plot_summ$bestradiol)){
  last_tpt = 0.75
  log_model_estimates = optim(par = c(2.9, 0.5, .5), min.2state.log.fun,
                              data = subset(plot_summ, substrate == 'yfp15' & promoter == 'pzev' & bestradiol == b & time_in_hr <= last_tpt & time_in_hr >= 0))
  log_par_estimates <- rbind(log_par_estimates, data.frame(bestradiol = b, 
                                                           a = log_model_estimates$par[1], b = log_model_estimates$par[2], d = log_model_estimates$par[3], last_tpt = last_tpt))
  
}

nls.2state.45min.plus.log.3 <- nlsLM(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                                     data = subset(plot_summ, bestradiol == 'plusbest' & time_in_hr <= last_tpt & time_in_hr >= 0 & promoter == 'pzev' & substrate == 'yfp15'),
                                     start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == 'plusbest'],
                                                  b = log_par_estimates$b[log_par_estimates$bestradiol == 'plusbest'],
                                                  d = log_par_estimates$d[log_par_estimates$bestradiol == 'plusbest']))
nls.2state.45min.minus.log <- nlsLM(norm_log_yfp ~ log((a-d)*exp(-(a+b)*time_in_hr)+b*exp(-d*time_in_hr))-log(a+b-d),
                                    data = subset(plot_summ, bestradiol == 'minusbest' & time_in_hr <= last_tpt & time_in_hr >= 0 & promoter == 'pzev' & substrate == 'yfp15'),
                                    start = list(a = log_par_estimates$a[log_par_estimates$bestradiol == 'minusbest'],
                                                 b = log_par_estimates$b[log_par_estimates$bestradiol == 'minusbest'],
                                                 d = log_par_estimates$d[log_par_estimates$bestradiol == 'minusbest']))
# generating a 2nd plot to use as an inset showing the diff with shorter timepoints
lm.45min.plus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(plot_summ, bestradiol == 'plusbest' & time_in_hr <= 0.75 & time_in_hr >=0 & promoter == 'pzev' & substrate == 'yfp15'))
lm.45min.minus.log <- lm(norm_log_yfp ~ time_in_hr - 1, data = subset(plot_summ, bestradiol == 'minusbest' & time_in_hr <= 0.75 & time_in_hr >= 0 &promoter == 'pzev' & substrate == 'yfp15'))


ggplot(subset(plot_summ, substrate == 'yfp15' & promoter == 'pzev' & time_in_hr <= 0.75 & time_in_hr >=0),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  scale_color_manual(values = c('black','red')) +
  geom_abline(intercept = 0, slope = coef(lm.45min.minus.log), color = 'black') +
  geom_abline(intercept = 0, slope = coef(lm.45min.plus.log), color = 'red') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75)) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.45min.minus.log))),
                color = 'black', linetype = 2) +
  stat_function(fun = C.2state.log.fun, 
                args=list(params = as.list(coef(nls.2state.45min.plus.log.3))),
                color = 'red', linetype = 2) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed mitochondrial YFP density\nnormalized to t = 0')
#ggsave('/path/to/fig/S4B/plot.pdf', device=cairo_pdf(width=6, height=4), useDingbats=FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
