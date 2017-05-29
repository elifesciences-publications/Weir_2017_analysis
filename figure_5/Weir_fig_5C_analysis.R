## Weir_fig_5C_analysis.R ##

# change the following path to point to your clone of the data.
input_data <- read_csv('~/Dropbox/2017_Weir_paper/Weir_2017_analysis/figure_5/fig_5C_data.csv')

# generate fluorescence density by dividing intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
input_data$mch_mean <- input_data$mcherry/input_data$volume
# eliminate too-small and too-large particles
input_data <- subset(input_data, volume > 5 & volume < 3000)
# fluorescence density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log10(input_data$yfp_mean)
input_data$log_mch <- log10(input_data$mch_mean)
summ <- ddply(.data = input_data, .variables = 'substrate', summarize,
              mean_log_yfp = mean(log_yfp),
              sd_log_yfp = sd(log_yfp),
              mean_log_mch = mean(log_mch),
              sd_log_mch = sd(log_mch))
# eliminate particles with mean fluorescence less than 2 sd above the bgrd
input_data$yfp_mean[input_data$yfp_mean < 10^(summ$mean_log_yfp[summ$substrate == 'untagged'] + 2*summ$sd_log_yfp[summ$substrate == 'untagged'])] <- NA
input_data$mch_mean[input_data$mch_mean < 10^(summ$mean_log_mch[summ$substrate == 'untagged'] + 2*summ$sd_log_mch[summ$substrate == 'untagged'])] <- NA
# use only values with positive fluorescence 
# (highly mobile particles sometimes incorrectly measured)
input_data <- input_data[complete.cases(input_data), ]
# subtract background values
input_data$bsub_yfp <- input_data$yfp_mean - 10^(summ$mean_log_yfp[summ$substrate == 'untagged'])
input_data$log_bsub_yfp <- log10(input_data$bsub_yfp)
input_data$bsub_mch <- input_data$mch_mean - 10^(summ$mean_log_mch[summ$substrate == 'untagged'])
input_data$log_bsub_mch <- log10(input_data$bsub_mch)


#calculate tFT ratio by dividing bsub_mch by bsub_yfp
input_data$tFT_ratio <- input_data$bsub_mch/input_data$bsub_yfp
input_data$log_tFT_ratio <- log10(input_data$tFT_ratio)
# summarize data for plotting
f_summ <- ddply(input_data, .variables = c('substrate','timepoint','bestradiol'),
                summarize,
                log_yfp_mean = mean(log_bsub_yfp),
                log_yfp_sd = sd(log_bsub_yfp),
                log_mch_mean = mean(log_bsub_mch),
                log_mch_sd = sd(log_bsub_mch),
                log_tft_mean = mean(log_tFT_ratio),
                log_tft_sd = sd(log_tFT_ratio),
                nobs = length(log_tFT_ratio[!is.na(log_tFT_ratio)]))
f_summ <- subset(f_summ, substrate == 'tftpex15')
# generate standard error values
f_summ$log_yfp_se <- f_summ$log_yfp_sd/sqrt(f_summ$nobs)
f_summ$log_mch_se <- f_summ$log_mch_sd/sqrt(f_summ$nobs)
f_summ$log_tft_se <- f_summ$log_tft_sd/sqrt(f_summ$nobs)
# transform back to linear space
f_summ$lin_yfp_mean <- 10^f_summ$log_yfp_mean
f_summ$lin_mch_mean <- 10^f_summ$log_mch_mean
f_summ$lin_tft_mean <- 10^f_summ$log_tft_mean
# calculate values for error bar limits for plotting
f_summ$yfp_ebar_min <- 10^(f_summ$log_yfp_mean-f_summ$log_yfp_se)
f_summ$yfp_ebar_max <- 10^(f_summ$log_yfp_mean+f_summ$log_yfp_se)
f_summ$mch_ebar_min <- 10^(f_summ$log_mch_mean-f_summ$log_mch_se)
f_summ$mch_ebar_max <- 10^(f_summ$log_mch_mean+f_summ$log_mch_se)
f_summ$tft_ebar_min <- 10^(f_summ$log_tft_mean-f_summ$log_tft_se)
f_summ$tft_ebar_max <- 10^(f_summ$log_tft_mean+f_summ$log_tft_se)
f_summ$time_in_min <- f_summ$timepoint * 60
f_summ$norm_yfp <- NA
f_summ$norm_yfp_ebar_min <- NA
f_summ$norm_yfp_ebar_max <- NA
f_summ$norm_tft <- NA
f_summ$norm_tft_ebar_min <- NA
f_summ$norm_tft_ebar_max <- NA
# normalize all values to the 0 hour timepoint
for (t in unique(f_summ$bestradiol)){
  f_summ$norm_yfp[f_summ$bestradiol == t] <- f_summ$lin_yfp_mean[f_summ$bestradiol == t]/f_summ$lin_yfp_mean[f_summ$bestradiol == t & f_summ$timepoint == 0]
  f_summ$norm_yfp_ebar_min[f_summ$bestradiol == t] <- f_summ$yfp_ebar_min[f_summ$bestradiol == t]/f_summ$lin_yfp_mean[f_summ$bestradiol == t & f_summ$timepoint == 0]
  f_summ$norm_yfp_ebar_max[f_summ$bestradiol == t] <- f_summ$yfp_ebar_max[f_summ$bestradiol == t]/f_summ$lin_yfp_mean[f_summ$bestradiol == t & f_summ$timepoint == 0]
  f_summ$norm_tft[f_summ$bestradiol == t] <- f_summ$lin_tft_mean[f_summ$bestradiol == t]/f_summ$lin_tft_mean[f_summ$bestradiol == t & f_summ$timepoint == 0]
  f_summ$norm_tft_ebar_min[f_summ$bestradiol == t] <- f_summ$tft_ebar_min[f_summ$bestradiol == t]/f_summ$lin_tft_mean[f_summ$bestradiol == t & f_summ$timepoint == 0]
  f_summ$norm_tft_ebar_max[f_summ$bestradiol == t] <- f_summ$tft_ebar_max[f_summ$bestradiol == t]/f_summ$lin_tft_mean[f_summ$bestradiol == t & f_summ$timepoint == 0]
}

ggplot(f_summ, aes(x = time_in_min, y = norm_yfp, color = bestradiol, group = bestradiol)) +
  geom_point() + geom_line() +
  scale_x_continuous(limits = c(-25,265), breaks = c(0,60,120,180,240)) +
  scale_color_manual(values = c('black','purple')) +
  geom_errorbar(width = 50, aes(ymin = norm_yfp_ebar_min, ymax = norm_yfp_ebar_max)) +
  scale_y_log10(limits = c(0.1, 1.1), 
                breaks = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                labels = c('','','','','','','',0.1,'','','','','','','','',1)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', size = 0.25),
        legend.key = element_blank(),
        legend.title = element_blank())
#ggsave('/path/to/fig/5C/yfp_plot.pdf', device = cairo_pdf(width = 4, height = 3), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

ggplot(f_summ, aes(x = time_in_min, y = norm_tft, color = bestradiol, group = bestradiol)) +
  geom_point() + geom_line() +
  scale_x_continuous(limits = c(-25,265), breaks = c(0,60,120,180,240)) +
  scale_color_manual(values = c('black','purple')) +
  geom_errorbar(width = 50, aes(ymin = norm_tft_ebar_min, ymax = norm_tft_ebar_max)) +
  scale_y_log10(limits = c(0.4, 1.35), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                labels = c(0.1,'','','',0.5,'','','','',1)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', size = 0.25),
        legend.key = element_blank(),
        legend.title = element_blank())
#ggsave('/path/to/fig/5C/tft_plot.pdf', device = cairo_pdf(width = 4, height = 3), useDingbats = FALSE)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.