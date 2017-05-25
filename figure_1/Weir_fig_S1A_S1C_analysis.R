# Weir_fig_S1A_S1C_analysis.R
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following two paths to point to your clone of the data.
singletpt_input_data <- read_csv('/path/to/Weir_2017_analysis/figure_1/S1_singletpt_analysis_output.csv')
timelapse_input_data <- read_csv('/path/to/Weir_2017_analysis/figure_1/S1_timelapse_analysis_output.csv')
# extract relevant sample information from image filenames in singletpt_input_data
split_fnames <- strsplit(singletpt_input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
singletpt_input_data$substrate <- split_fnames$X3
# add additional treatment information to match timelapse
singletpt_input_data$bestradiol <- 'minusbest'
singletpt_input_data$timepoint <- 0
# extract relevant sample information from image filenames in timelapse_input_data
split_fnames <- strsplit(timelapse_input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
timelapse_input_data$substrate <- split_fnames$X4
timelapse_input_data$bestradiol <- split_fnames$X3
timelapse_input_data$timepoint <- as.numeric(substr(split_fnames$X6,2,3))
timelapse_input_data$timepoint <- timelapse_input_data$timepoint - 1 # first timepoint = 0 hour
# merge t=0 single timepoint data and timelapse data
input_data <- rbind(singletpt_input_data, timelapse_input_data)
# generate YFP density by dividing total YFP by object volumes
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# log-transform because fluorescence density is roughly log-normal
input_data$log_yfp <- log(input_data$yfp_mean)
# get mean and sd of data for plotting
summ <- ddply(.data = input_data, .variables = c('substrate','obj_channel'), summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
input_data$bsub_yfp <- NA
# subtract background values in each object channel (447 = mitos, 594 = pexs) using mean value 
# from strain with no MSP1YFP
for (c in unique(summ$obj_channel)){
  input_data$bsub_yfp[input_data$obj_channel == c] <- input_data$yfp_mean[input_data$obj_channel == c] - exp(summ$mean_log_yfp[summ$substrate == 'noyfp' & summ$obj_channel == c])
}
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
input_data$obj_channel <- as.character(input_data$obj_channel)
input_data$obj_channel[input_data$obj_channel == '447'] <- 'mitochondria'
input_data$obj_channel[input_data$obj_channel == '594'] <- 'peroxisomes'
# generate mean and sd for plotting
plot_summ <- ddply(input_data, .variables = c('substrate','timepoint','bestradiol','obj_channel'),
                   summarize, yfp_mean = mean(bsub_yfp),
                   yfp_sd = sd(bsub_yfp),
                   nobs = length(bsub_yfp))
# generate standard error
plot_summ$yfp_se <- plot_summ$yfp_sd/sqrt(plot_summ$nobs)
# normalize to endogenously expressed Msp1-YFP intensity levels
plot_summ$yfp_norm <- NA
plot_summ$yfp_norm_se <- NA
for (c in unique(plot_summ$obj_channel)){
  plot_summ$yfp_norm[plot_summ$obj_channel == c] <- plot_summ$yfp_mean[plot_summ$obj_channel == c]/plot_summ$yfp_mean[plot_summ$substrate == 'msp1yfp' & plot_summ$obj_channel == c]
  plot_summ$yfp_norm_se[plot_summ$obj_channel == c] <- plot_summ$yfp_se[plot_summ$obj_channel == c]/plot_summ$yfp_mean[plot_summ$substrate == 'msp1yfp' & plot_summ$obj_channel == c]
}
# plot data for figure S1C.
# because experiments never proceed beyond 4 hours for this paper, only plotted data
# for the first 4 hours of MSP1 expression.
ggplot(subset(plot_summ, substrate == 'zevmsp1yfp' & timepoint < 5),
       aes(x = timepoint, y = yfp_norm, color = bestradiol, group = bestradiol)) +
  facet_grid(obj_channel ~ .) + geom_line() + geom_point() +
  geom_hline(yintercept = 1, linetype = 2, color = 'gray50') +
  scale_color_manual(values = c('purple','black')) +
  geom_errorbar(aes(ymin = yfp_norm- yfp_norm_se, ymax = yfp_norm + yfp_norm_se), width = 0.5) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', size = 0.5),
        strip.background = element_blank(),
        legend.key = element_blank())
#ggsave('/path/to/fig/S1C/plot.pdf', device = cairo_pdf(width = 6, height = 3), useDingbats = F)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
plot_summ$bestradiol[plot_summ$bestradiol == 'minusbest'] <- 'nobest' # correcting inconsistent labeling
ggplot(subset(plot_summ, substrate %in% c('zevmsp1yfp','msp1yfp') & timepoint == 0 & bestradiol == 'nobest'),
       aes(x = obj_channel, y = yfp_norm, fill = substrate)) + geom_bar(color = 'black', stat = 'identity', position = position_dodge(width = 0.75), width = 0.75) +
  geom_errorbar(aes(ymin = yfp_norm-yfp_norm_se, ymax = yfp_norm + yfp_norm_se), position = position_dodge(width = 0.75), width = 0.25) +
  geom_hline(yintercept = 0, color = 'black', size = 0.5) +
  scale_fill_manual(values = c('white','gray50')) +
  theme(panel.background = element_rect(fill = NA, color = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', size = 0.5),
        strip.background = element_blank(),
        legend.key = element_blank(),
        axis.line.y = element_line(size = 0.5, color = 'black'))
#ggsave('/path/to/fig/S1A/plot.pdf', device = cairo_pdf(width = 6, height = 5), useDingbats = F)
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.