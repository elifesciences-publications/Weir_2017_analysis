# Weir_fig_1C_S1D_analysis.R

# load dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_1/fig_1C_S1D_data.csv')
# parse input image filenames to extract relevant sample properties
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$dox <- split_fnames$X4
input_data$bestradiol <- split_fnames$X5
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
input_data$timepoint[input_data$timepoint != 0] <- input_data$timepoint[input_data$timepoint != 0] - 1
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
input_data$norm_to_bg_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$norm_to_bg_yfp <- input_data$yfp_mean/exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
# plot for figure 1C
ggplot(subset(input_data, substrate == 'yfp15'),
       aes(x = factor(timepoint), color = bestradiol, y = bsub_yfp)) + 
  geom_violin(width = 0.8, position = position_dodge(width = 0.8)) + 
  scale_y_log10(breaks = c(10,20,30,40,50,60,70,80,90,
                           100,200,300,400,500,600,700,800,900,
                           1000,2000,3000,4000,5000,6000,7000,8000,9000,
                           10000,20000,30000),
                labels = c(10,'','','','','','','','',
                           100,'','','','','','','','',
                           1000,'','','','','','','','',
                           10000,'','')) +
  scale_color_manual(values = c('black','purple')) +
  theme(panel.background = element_rect(color = 'black',fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Time after MSP1 expression (hr)',
       y = 'YFP-Pex15 remaining at mitochondria\nnormalized to -MSP1 mean at time 0',
       title = 'Mitochondrial YFP-Pex15 turnover')
#ggsave('~/path/to/Fig/1C/plot.pdf', device = cairo_pdf(width = 6, height = 5), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# plot for figure 1-figure supplement 1D
input_data$minadded_yfp <- input_data$bsub_yfp - min(input_data$bsub_yfp) + 1
input_data$minadded_norm_yfp <- input_data$minadded_yfp/exp(mean(log(input_data$minadded_yfp[input_data$substrate == 'noyfp'])))
ggplot(subset(input_data, timepoint %in% c(0,1,2)),
       aes(x = factor(timepoint), color = bestradiol, linetype = substrate, y = minadded_norm_yfp)) + 
  geom_violin(width = 2, position = position_dodge(width = 1)) + 
  scale_y_log10(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                           1,2,3,4,5,6,7,8,9,
                           10,20,30,40,50,60,70,80,90,
                           100,200,300,400,500,600,700,800,900,1000),
                labels = c(0.1,'','','','','','','','',
                           1,'','','','','','','','',
                           10,'','','','','','','','',
                           100,'','','','','','','','',1000)) +
  scale_linetype_manual(values = c(3,1)) +
  scale_color_manual(values = c('black','purple')) +
  theme(panel.background = element_rect(color = 'black',fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = mean(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp']) + 3*sd(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp']), 
             color = 'red', size = 0.5, linetype = 2) +
  labs(x = 'Time after MSP1 expression (hr)',
       y = 'YFP-Pex15 remaining at mitochondria\nnormalized to autofluorescence',
       title = 'Mitochondrial YFP-Pex15 turnover')
#ggsave('/path/to/Fig/S1D/plot.pdf', useDingbats=FALSE, device = cairo_pdf(width = 7, height = 3))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
