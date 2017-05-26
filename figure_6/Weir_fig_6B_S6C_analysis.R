## Weir_fig_6B_S6C_analysis.R ## 

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_6/fig_6B_S6C_data.csv')
# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$msp1_allele <- split_fnames$X4
input_data$auxin <- split_fnames$X5
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
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
# pexs in cells lacking YFP
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
#fixing auxin labels
input_data$auxin <- as.character(input_data$auxin)
input_data$auxin[input_data$auxin == '.minusauxin'] <- 'minusauxin'
input_data$auxin[input_data$auxin == '.plusauxin'] <- 'plusauxin'
input_data$auxin <- factor(input_data$auxin, levels = c('minusauxin','plusauxin'))
# summarize data for plotting
plot_summ <- ddply(.data = input_data, .variables = c('substrate','timepoint','msp1_allele','auxin'),
                   summarize,
                   bsub_yfp_mean = mean(bsub_yfp),
                   bsub_yfp_sd = sd(bsub_yfp),
                   log_bsub_yfp_mean = mean(log_bsub_yfp, na.rm = TRUE),
                   log_bsub_yfp_sd = sd(log_bsub_yfp, na.rm = TRUE),
                   nobs = length(bsub_yfp),
                   nobs_below_0 = length(log_bsub_yfp[is.na(log_bsub_yfp)]),
                   nobs_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd]),
                   frac_below_bg = length(log_yfp[log_yfp<ctrl_mean+2*ctrl_sd])/length(log_yfp))
# calculate standard error
plot_summ$log_bsub_yfp_se <- plot_summ$log_bsub_yfp_sd/sqrt(plot_summ$nobs-plot_summ$nobs_below_0)
# add time in hour
plot_summ$time_in_hr <- plot_summ$timepoint * 0.25
plot_summ$norm_to_zero <- NA
# normalize to t = 0 for the same parent strain w/o auxin
# return to linear space for supplemental figure plotting
plot_summ$delogged_bsub_yfp_mean <- exp(plot_summ$log_bsub_yfp_mean)
plot_summ$delogged_ebar_min <- exp(plot_summ$log_bsub_yfp_mean - plot_summ$log_bsub_yfp_se)
plot_summ$delogged_ebar_max <- exp(plot_summ$log_bsub_yfp_mean + plot_summ$log_bsub_yfp_se)
for (m in plot_summ$msp1_allele){
  plot_summ$norm_to_zero[plot_summ$msp1_allele == m] <- plot_summ$log_bsub_yfp_mean[plot_summ$msp1_allele == m] - plot_summ$log_bsub_yfp_mean[plot_summ$msp1_allele == m & plot_summ$time_in_hr == 0 & plot_summ$auxin == 'minusauxin' & plot_summ$substrate == 'yfp15']
}

ggplot(subset(plot_summ, substrate == 'yfp15'),
       aes(x = time_in_hr, y = norm_to_zero, linetype = msp1_allele, color = auxin, group = interaction(auxin, msp1_allele))) +
  geom_point() + geom_line(size = 0.33) +
  geom_errorbar(aes(ymin = norm_to_zero - log_bsub_yfp_se,
                    ymax = norm_to_zero + log_bsub_yfp_se), width = 0.2, linetype = 1, size = 0.33) +
  scale_color_manual(values = c('black','red', 'green','blue')) +
  theme(panel.background = element_rect(color = 'black', fill = NA),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        legend.title = element_blank()) +
  labs(x = 'Time after addition of b-estradiol, hr',
       y = 'Log-transformed peroxisomal YFP density, AU',
       title = 'Peroxisomal YFP-Pex15 turnover\nfollowing MSP1 expression')
#ggsave('/path/to/fig/6B/plot.pdf', device = cairo_pdf(width = 6, height = 5), useDingbats = F)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# plot comparison of absolute starting levels for Figure 6-Figure Supplement 1C
ggplot(subset(plot_summ, time_in_hr == 0 & substrate == 'yfp15' & auxin == 'minusauxin'),
       aes(x = msp1_allele, y = delogged_bsub_yfp_mean)) +
  geom_bar(stat = 'identity', width = 0.75, color = 'black', fill = 'gray50') +
  geom_errorbar(aes(ymin = delogged_ebar_min, ymax = delogged_ebar_max), position = position_dodge(width = 0.75), width = 0.375) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5000)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = 'black'))
#ggsave('/path/to/fig/S6C/plot.pdf', device = cairo_pdf(width = 2, height = 4))