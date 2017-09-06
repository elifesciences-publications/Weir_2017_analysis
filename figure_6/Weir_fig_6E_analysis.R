## Weir_fig_6C_analysis.R ## 

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_6/Fig_6E_data.csv')
# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$pex3 <- split_fnames$X4
input_data$bestradiol <- split_fnames$X5
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
input_data$timepoint <- input_data$timepoint - 1
input_data <- subset(input_data, volume != 0)
input_data$time_in_hr<- input_data$timepoint*0.25
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log(input_data$yfp_mean)
# calculate means for bgrd subtraction
summ <- ddply(.data = input_data, .variables = c('substrate'), summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# generate background-subtracted YFP values by subtracting mean YFP density at 
# pexs in cells lacking YFP
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# summarize data for plotting

plot_summ <- ddply(.data = input_data, .variables = c('substrate','time_in_hr','bestradiol','pex3'),
                   summarize,
                   bsub_yfp_mean = mean(bsub_yfp),
                   bsub_yfp_sd = sd(bsub_yfp),
                   log_bsub_yfp_mean = mean(log_bsub_yfp, na.rm = TRUE),
                   log_bsub_yfp_sd = sd(log_bsub_yfp, na.rm = TRUE),
                   nobs = length(bsub_yfp),
                   nobs_below_0 = length(log_bsub_yfp[is.na(log_bsub_yfp)]))
plot_summ$time_in_hr <- plot_summ$time_in_hr - 0.25
# calculate standard error
plot_summ$log_bsub_yfp_se <- plot_summ$log_bsub_yfp_sd/sqrt(plot_summ$nobs-plot_summ$nobs_below_0)
# normalize to t=0
plot_summ$norm_log_yfp <- NA
for (c in unique(plot_summ$pex3)){
  for (s in unique(plot_summ$substrate)){
    plot_summ$norm_log_yfp[plot_summ$pex3 == c & plot_summ$substrate == s] <- plot_summ$log_bsub_yfp_mean[plot_summ$pex3 == c & plot_summ$substrate == s] -  plot_summ$log_bsub_yfp_mean[plot_summ$pex3 == c & plot_summ$substrate == s & plot_summ$time_in_hr == 0 & plot_summ$bestradiol == 'minusbest']
  }
}
# plot
ggplot(subset(plot_summ, substrate != 'noyfp' & time_in_hr >= 0 & time_in_hr <= 2),
       aes(x = time_in_hr, y = norm_log_yfp, color = bestradiol, linetype=pex3)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c('black','red')) +
  geom_errorbar(aes(ymin = norm_log_yfp-log_bsub_yfp_se, ymax = norm_log_yfp+log_bsub_yfp_se), width = 0.1) +
  theme(panel.background = element_rect(fill=NA, color = 'black'),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black')) +
  labs(x='Time after MSP1 expression (hr)', y = 'Log-transformed peroxisomal YFP density\nnormalized to minus b-estradiol t = 0')
#ggsave('/path/to/fig/6E/plot.pdf', device = cairo_pdf(width = 2, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.