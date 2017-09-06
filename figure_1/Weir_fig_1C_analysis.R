## Weir_fig_1C_analysis.R ##

# load dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_1/fig_1C_data.csv')
# generate YFP density by dividing total YFP intensity by volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate volume = 0 particles (artifact of image processing)
input_data <- subset(input_data, volume != 0)
# data is closer to log-normal than linear normal distribution, so log-transform
# for identifying means and sds of populations
input_data$log_yfp <- log(input_data$yfp_mean)
# calculate mean of background for background subtraction
summ <- ddply(.data = input_data, .variables = 'substrate', summarize,
              mean_log_yfp = mean(log_yfp))
# generate background subtracted values using no YFP control
input_data$bsub_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$bgnorm_yfp <- input_data$yfp_mean - min(input_data$yfp_mean)
input_data$bgnorm_yfp <- input_data$bgnorm_yfp/mean(input_data$bgnorm_yfp[input_data$substrate == 'noyfp'])

# plot for figure 1C
ggplot(subset(input_data, substrate == 'yfp15' & time_in_hr <= 0.5),
       aes(x = factor(time_in_hr), y = bgnorm_yfp, color=bestradiol)) +
  geom_violin(position = position_dodge(width = 0.75), width = 0.8) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_log10(breaks=c(1,2,3,4,5,6,7,8,9,
                         10,20,30,40,50,60,70,80,90,
                         100),
                labels=c(1,'','','','','','','','',
                         10,'','','','','','','','',
                         100),
                limits=c(1,100)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank())
#ggsave('~/path/to/Fig/1C/plot.pdf', device = cairo_pdf(width = 6, height = 5), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.