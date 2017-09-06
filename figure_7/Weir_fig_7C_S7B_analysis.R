## Weir_fig_7C_S7B_analysis.R ## 

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_7/Fig_7C_S7B_data.csv')
# extract relevant sample information from image filenames

split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))

input_data$substrate <- split_fnames$X1
input_data$MSP1_allele <- split_fnames$X3
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log(input_data$yfp_mean)
# calculate means for bgrd subtraction
summ <- ddply(.data = input_data, .variables = c('substrate', 'MSP1_allele'), summarize,
              mean_yfp = mean(yfp_mean))

# generate background-subtracted YFP values by subtracting mean YFP density at 
# pexs in cells lacking YFP
input_data$bsub_yfp <- input_data$yfp_mean - summ$mean_yfp[summ$substrate == 'mitocontrol']
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# summarize data for plotting by calculating fraction of activity in each mutant.
# see figure legend.
calc_resc <- ddply(.data = subset(input_data, substrate == 'deltac30'),
                   .variables = 'MSP1_allele', summarize,
                   mean_bsub = mean(bsub_yfp), bsub_se = sd(bsub_yfp)/sqrt(length(bsub_yfp)))
calc_resc$rescue <- NA
calc_resc$rescue_se <- NA
calc_resc$rescue <- (calc_resc$mean_bsub[calc_resc$MSP1_allele == 'deltamsp1'] - calc_resc$mean_bsub)/(calc_resc$mean_bsub[calc_resc$MSP1_allele == 'deltamsp1'] - calc_resc$mean_bsub[calc_resc$MSP1_allele == 'msp1wt'])
calc_resc$rescue_se <- (calc_resc$bsub_se[calc_resc$MSP1_allele == 'deltamsp1'] + calc_resc$bsub_se)/(calc_resc$mean_bsub[calc_resc$MSP1_allele == 'deltamsp1'] - calc_resc$mean_bsub[calc_resc$MSP1_allele == 'msp1wt'])

# plot data.
ggplot(subset(calc_resc, MSP1_allele %in% c('tom70n','pex22n')),
       aes(x = MSP1_allele, y = rescue, color = MSP1_allele)) +
  geom_bar(stat = 'identity', fill = 'white') +
  geom_errorbar(aes(ymin = rescue - rescue_se, ymax = rescue + rescue_se), width = 0.67) +
  coord_cartesian(ylim = c(-0.1,1)) +
  scale_color_manual(values = c('darkblue','darkgreen')) +
  geom_hline(yintercept = 0, size = 1) +
  scale_y_continuous() +
  theme(panel.background = element_rect(color = NA,fill = NA),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 270, vjust = 0.5),
        axis.line.y = element_line(color = 'black', linetype = 'solid', size = 1))
#ggsave('/path/to/fig/7C/plot.pdf', device = cairo_pdf(width = 2, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# plot violins for supplement
ggplot(subset(input_data, substrate == 'deltac30'), aes(x = MSP1_allele, y = bsub_yfp)) + geom_violin() +
  theme(panel.background = element_rect(color = 'black',fill = NA),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.position = 'none',
      axis.text = element_text(color = 'black'),
      axis.text.x = element_text(angle = 270, vjust = 0.5))
#ggsave('/path/to/fig/S7B/plot.pdf', device = cairo_pdf(width = 2, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
