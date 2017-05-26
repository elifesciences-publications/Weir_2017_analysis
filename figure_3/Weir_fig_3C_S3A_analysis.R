## Weir_fig_3C_S3A_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
# change the following path to point to your clone of the data.
input_data <- read_csv('~/Dropbox/2017_Weir_paper/Weir_2017_analysis/figure_3/fig_3C_S3A_data.csv')
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_3/fig_3C_S3A_data.csv')
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
# find mean intensity for background
summ <- ddply(.data = input_data, .variables = 'substrate', summarize,
              mean_log_yfp = mean(log_yfp))
# subtract mean of no-YFP sample to get background subtracted pixel intensities
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# each timepoint is 30 mins
input_data$time_in_hr <- input_data$timepoint*0.5
# plot for figure 3C
ggplot(subset(input_data, substrate == 'wt' & time_in_hr %in% c(0,1,2)), aes(x = factor(time_in_hr), y = bsub_yfp, color = bestradiol)) +
  geom_violin(position = position_dodge(width = 0.6), width = 0.7) +
  scale_color_manual(values = c('black','red')) +
  scale_y_log10(breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000),
                labels = c(1,'','','','','','','','',
                           10,'','','','','','','','',
                           100,'','','','','','','','',
                           1000,'','','','','','','','',
                           10000,'')) +
  theme(panel.background = element_rect(color = 'black', fill = NA),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black')) +
  labs(x = 'Time after MSP1 expression, hr',
       y = 'Peroxisomal YFP density, AU',
       title = 'Peroxisomal YFP-Pex15 turnover\nfollowing MSP1 expression')
#ggsave('/path/to/fig/3C/plot.pdf', useDingbats = FALSE, device = cairo_pdf(width = 6, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# plotting for supplemental figure 3A showing fold above background
input_data$minadded_yfp <- input_data$bsub_yfp - min(input_data$bsub_yfp) + 1
input_data$minadded_norm_yfp <- input_data$minadded_yfp/exp(mean(log(input_data$minadded_yfp[input_data$substrate == 'noyfp'])))
plot_b <- ggplot(subset(input_data, time_in_hr %in% c(0,1,2)),
                 aes(x = factor(time_in_hr), color = bestradiol,
                     y = minadded_norm_yfp, linetype = substrate)) + 
  geom_violin(position = position_dodge(width = 0.5), width = 1) + 
  scale_y_log10() +
  scale_linetype_manual(values = c(3,1)) +
  scale_color_manual(values = c('black','purple')) +
  theme(panel.background = element_rect(color = 'black',fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = mean(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp' & input_data$time_in_hr <=2]) + 3*sd(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp' & input_data$time_in_hr <=2]),
             color = 'red', linetype = 2, size = 0.5) +
  labs(x = 'Time after MSP1 expression (hr)',
       y = 'YFP-Pex15 remaining at peroxisomes\nnormalized to autofluorescence',
       title = 'Peroxisomal YFP-Pex15 turnover')
#ggsave('/path/to/fig/S3A/plot.pdf', device = cairo_pdf(width = 12, height = 3), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.