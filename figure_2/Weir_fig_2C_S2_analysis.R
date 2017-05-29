## Weir_fig_2C_S2_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
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
# plot for Figure 2C
ggplot(subset(input_data, substrate == 'deltac30' & time_in_hr %in% c(0,0.5,1)),
       aes(x = factor(time_in_hr), y = bsub_yfp, color = bestradiol)) +
  geom_violin(position = position_dodge(width = 0.75), width = 0.8) +
  scale_color_manual(values = c('black','purple')) +
  scale_y_log10(breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000),
                labels = c(1,'','','','','','','','',
                           10,'','','','','','','','',
                           100,'','','','','','','','',
                           1000,'','','','')) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'))
#ggsave('/path/to/fig/2C/plot.pdf', device = cairo_pdf(width = 6, height = 5))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# prepare data to plot for figure S2. 
# Because there were relatively few background values measured at each timepoint,
# they were combined to generate one background distribution.
input_data$minadded_yfp <- input_data$bsub_yfp - min(input_data$bsub_yfp) # for plotting background on log scale
# normalize to background mean
input_data$minadded_norm_yfp <- input_data$minadded_yfp/exp(mean(log(input_data$minadded_yfp[input_data$substrate == 'noyfp'])))
plot_b <- ggplot(subset(input_data, time_in_hr %in% c(0,0.5,1) & substrate == 'deltac30'),
                 aes(x = factor(time_in_hr), color = bestradiol, y = minadded_norm_yfp)) + 
  geom_violin() + 
  scale_y_log10(limits = c(0.5, 100), 
                breaks = c(0.5,0.6,0.7,0.8,0.9,1,
                           2,3,4,5,6,7,8,9,10,
                           20,30,40,50,60,70,80,90,100),
                labels = c('','','','','',1,
                           '','','','','','','','',10,
                           '','','','','','','','',100)) +
  scale_linetype_manual(values = c(1,3)) +
  scale_color_manual(values = c('black','purple')) +
  theme(panel.background = element_rect(color = 'black',fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = mean(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp' & input_data$time_in_hr <=2]) + 3*sd(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp' & input_data$time_in_hr <=1]),
             color = 'red', linetype = 2, size = 0.5) +
  labs(x = 'Time after MSP1 expression (hr)',
       y = 'YFP-Pex15 remaining at mitochondria\nnormalized to autofluorescence',
       title = 'Mitochondrial YFP-Pex15 turnover')     

plot_a <- ggplot(subset(input_data, time_in_hr <= 2 & substrate == 'noyfp'),
                 aes(color = bestradiol, x = substrate, y = minadded_norm_yfp)) + 
  geom_violin(color = 'black', linetype = 2) +
  scale_y_log10(limits = c(0.5, 100), 
                breaks = c(0.5,0.6,0.7,0.8,0.9,1,
                           2,3,4,5,6,7,8,9,10,
                           20,30,40,50,60,70,80,90,100),
                labels = c('','','','','',1,
                           '','','','','','','','',10,
                           '','','','','','','','',100)) +
  geom_hline(yintercept = mean(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp' & input_data$time_in_hr <=2]) + 3*sd(input_data$minadded_norm_yfp[input_data$substrate == 'noyfp' & input_data$time_in_hr <=1]),
             color = 'red', linetype = 2, size = 0.5) +
  theme(panel.background = element_rect(color = 'black',fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Time after MSP1 expression (hr)',
       y = 'YFP-Pex15 remaining at mitochondria\nnormalized to autofluorescence',
       title = 'Mitochondrial YFP-Pex15 turnover')  

grid.arrange(plot_a, plot_b, widths = c(1,3))
combined_plots <- arrangeGrob(plot_a, plot_b, widths = c(1,3))
#ggsave('/path/to/fig/S2/plot.pdf', plot = combined_plots, device = cairo_pdf(width = 6, height = 3), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

