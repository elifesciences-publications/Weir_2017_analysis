## Weir_fig_7E_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_7/Fig_7E_data.csv')

# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))

input_data$substrate <- split_fnames$X3
input_data$tmd <- split_fnames$X4
input_data$bestradiol <- split_fnames$X5
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
input_data$timepoint <- input_data$timepoint - 1
# each timepoint corresponds to 15 minutes
input_data$time_in_hr = 0.25*input_data$timepoint
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log(input_data$yfp_mean)
# calculate means for YFP-15 and no YFP for bgrd subtraction
summ <- ddply(.data = input_data, .variables ='substrate', summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# generate background-subtracted YFP values by subtracting mean YFP density at 
# mitos in cells lacking YFP
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# calculate background mean and sd for tabulation
ctrl_mean <- mean(input_data$log_yfp[input_data$substrate == 'noyfp'])
ctrl_sd <- sd(input_data$log_yfp[input_data$substrate == 'noyfp'])
# summarize data for model fitting
summ_for_plot <- ddply(.data=input_data, .variables=c('substrate','tmd','bestradiol','time_in_hr'),
                       summarize,
                       mean_log_yfp=mean(log_bsub_yfp, na.rm=TRUE),
                       sd_log_yfp = sd(log_bsub_yfp, na.rm=TRUE), 
                       nobs = length(log_bsub_yfp))
summ_for_plot$se_log_yfp <- summ_for_plot$sd_log_yfp/sqrt(summ_for_plot$nobs)
# normalize to sample t=0
for (t in unique(summ_for_plot$tmd)){
  summ_for_plot$norm_log_yfp[summ_for_plot$tmd == t] <- summ_for_plot$mean_log_yfp[summ_for_plot$tmd == t] - summ_for_plot$mean_log_yfp[summ_for_plot$tmd == t & summ_for_plot$bestradiol == 'minusbest' & summ_for_plot$time_in_hr == 0]
}
# MSP1 induction took 30 minutes after beginning treatment.
summ_for_plot$time_in_hr <- summ_for_plot$time_in_hr - 0.5
# plot
ggplot(subset(summ_for_plot, substrate == 'yfp15' & time_in_hr >= 0),
       aes(x = factor(time_in_hr), y = norm_log_yfp, color = tmd, linetype = bestradiol, group = interaction(tmd, bestradiol))) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = norm_log_yfp-se_log_yfp, ymax = norm_log_yfp+se_log_yfp)) +
  scale_color_manual(values = c('black', 'red')) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank())
#ggsave('/path/to/fig/7E/plot.pdf', device=cairo_pdf(width=6, height=4), useDingbats=FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
