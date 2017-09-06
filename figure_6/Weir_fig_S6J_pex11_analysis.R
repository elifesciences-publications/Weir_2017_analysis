## Weir_fig_S6J_pex11_analysis.R ## 

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_6/Fig_S6J_pex11_data.csv')
# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))

input_data$allele <- split_fnames$X4
input_data$tir <- split_fnames$X3
input_data$auxin <- split_fnames$X5
input_data$timepoint <- split_fnames$X1
input_data$timepoint <- as.character(input_data$timepoint)
input_data$timepoint[input_data$timepoint == '0hr'] <- 0
input_data$timepoint[input_data$timepoint == '0.5hr'] <- 1
input_data$timepoint[input_data$timepoint == '1hr'] <- 2
input_data$timepoint[input_data$timepoint == '1.5hr'] <- 3
input_data$timepoint[input_data$timepoint == '2hr'] <- 4
input_data$timepoint <- as.numeric(input_data$timepoint)
input_data$time_in_hr = 0.5*input_data$timepoint
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log(input_data$yfp_mean)
# calculate means for YFP-15 and no YFP for bgrd subtraction
summ <- ddply(.data = input_data, .variables ='tir', summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# generate background-subtracted YFP values by subtracting mean YFP density at 
# pexs in cells lacking YFP
input_data$bsub_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$tir == 'notft'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)

# summarize data for plotting
no_lo = subset(input_data, bsub_yfp > 10) # eliminate cells lacking fluorescence
summ_for_plot <- ddply(.data=no_lo, .variables=c('allele','auxin','tir','time_in_hr'),
                       summarize,
                       mean_log_yfp=mean(log_bsub_yfp, na.rm=TRUE),
                       sd_log_yfp = sd(log_bsub_yfp, na.rm=TRUE), 
                       nobs = length(log_bsub_yfp))
# calculate standard error
summ_for_plot$se_log_yfp = summ_for_plot$sd_log_yfp/sqrt(summ_for_plot$nobs)
# normalize to t=0 MSP1
summ_for_plot$norm_log_yfp = summ_for_plot$mean_log_yfp - summ_for_plot$mean_log_yfp[summ_for_plot$allele == 'msp1' & summ_for_plot$auxin == 'plusauxin' & summ_for_plot$tir == '2773' & summ_for_plot$time_in_hr == 0]
# plot.
ggplot(subset(summ_for_plot, tir != 'notft' & auxin == 'plusauxin' & time_in_hr %in% c(0,0.5,1,1.5,2)),
       aes(x = factor(time_in_hr), y = norm_log_yfp, color = allele, group = allele)) + 
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = norm_log_yfp-se_log_yfp, ymax = norm_log_yfp+se_log_yfp)) +
  scale_color_manual(values = c('black', 'red')) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank())
#ggsave('/path/to/fig/S6J_pex11/plot.pdf', device = cairo_pdf(width = 3, height = 3), useDingbats = F)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.