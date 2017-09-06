## Weir_fig_2C_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_2/Fig_2C_4D_data.csv')
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
summ <- ddply(.data = input_data, .variables ='substrate', summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# generate background-subtracted YFP values by subtracting mean YFP density at 
# mitos in cells lacking YFP
input_data$bsub_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# plot for Figure 2C
input_data$bgnorm_yfp <- input_data$yfp_mean - min(input_data$yfp_mean)
input_data$bgnorm_yfp <- input_data$bgnorm_yfp/mean(input_data$bgnorm_yfp[input_data$substrate == 'noyfp'])

ggplot(subset(input_data, substrate %in% c('noyfp','deltac30') & time_in_hr < 1),
       aes(x = factor(time_in_hr), y = bgnorm_yfp, color=bestradiol)) +
  geom_violin(position = position_dodge(width = 0.75), width = 0.8) + facet_grid(. ~ substrate) +
  scale_color_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  #coord_cartesian(ylim = c(1,500)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank())
#ggsave('/path/to/fig/2C/plot.pdf', device = cairo_pdf(width = 6, height = 5))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.