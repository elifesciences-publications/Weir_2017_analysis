## Weir_fig_S6A_analysis.R ## 

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_6/fig_S6A_data.csv')
# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$pex_allele <- split_fnames$X4
input_data$msp1_allele <- split_fnames$X5
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
input_data$log_yfp <- log(input_data$yfp_mean)
# calculate mean YFP and subtract background using non-fluorescent strain
summ <- ddply(input_data, .variables = c('substrate'),
              summarize, log_yfp_mean = mean(log_yfp),
              nobs = length(log_yfp))
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$log_yfp_mean[summ$substrate == 'noyft'])


# peroxisome size changes artifactually altered intensity due to convolution. 
# we normalized using a linear model.
for_plt <- subset(input_data, substrate == 'yft15')
wt_size_lm <- lm(log(bsub_yfp) ~ log(volume), data = subset(for_plt, substrate == 'yft15' & pex_allele == 'wt' & msp1_allele == 'msp1'))
for_plt$log_bsub_yfp <- log(for_plt$bsub_yfp)
for_plt <- subset(for_plt, !is.na(log_bsub_yfp))
for_plt$log_norm_yfp <- for_plt$log_bsub_yfp - log(for_plt$volume) * wt_size_lm$coefficients[2]
# summarize data for plotting
summ_for_plt <- ddply(for_plt, .variables = c('pex_allele','msp1_allele'),
                      summarize,
                      log_norm_yfp_mean = mean(log_norm_yfp),
                      log_norm_yfp_sd = sd(log_norm_yfp),
                      log_yfp_mean = mean(log_bsub_yfp),
                      log_yfp_sd = sd(log_bsub_yfp),
                      nobs = length(log_norm_yfp))
# calculate standard error values
summ_for_plt$log_norm_yfp_se <- summ_for_plt$log_norm_yfp_sd/sqrt(summ_for_plt$nobs)
summ_for_plt$log_yfp_se <- summ_for_plt$log_yfp_sd/sqrt(summ_for_plt$nobs)
# fix msp1 allele labels
summ_for_plt$msp1_allele <- factor(summ_for_plt$msp1_allele, levels = c('msp1','deltamsp1'))
# get error bar limits
summ_for_plt$log_ebar_min <- summ_for_plt$log_norm_yfp_mean-summ_for_plt$log_norm_yfp_se
summ_for_plt$log_ebar_max <- summ_for_plt$log_norm_yfp_mean+summ_for_plt$log_norm_yfp_se
# return to linear space
summ_for_plt$norm_yfp_mean <- exp(summ_for_plt$log_norm_yfp_mean)
summ_for_plt$yfp_mean <- exp(summ_for_plt$log_yfp_mean)
summ_for_plt$yfp_min <- exp(summ_for_plt$log_yfp_mean-summ_for_plt$log_yfp_se)
summ_for_plt$yfp_max <- exp(summ_for_plt$log_yfp_mean+summ_for_plt$log_yfp_se)
summ_for_plt$ebar_min <- exp(summ_for_plt$log_ebar_min)
summ_for_plt$ebar_max <- exp(summ_for_plt$log_ebar_max)
# fix pex allele labels
summ_for_plt$pex_allele <- factor(summ_for_plt$pex_allele, levels = c('wt','deltapex1','deltapex6'))
# normalize to wt pex for plotting
summ_for_plt$wtnorm_yfp <- summ_for_plt$yfp_mean/summ_for_plt$yfp_mean[summ_for_plt$pex_allele == 'wt' & summ_for_plt$msp1_allele == 'msp1']
summ_for_plt$wtnorm_yfp_min <- summ_for_plt$yfp_min/summ_for_plt$yfp_mean[summ_for_plt$pex_allele == 'wt' & summ_for_plt$msp1_allele == 'msp1']
summ_for_plt$wtnorm_yfp_max <- summ_for_plt$yfp_max/summ_for_plt$yfp_mean[summ_for_plt$pex_allele == 'wt' & summ_for_plt$msp1_allele == 'msp1']

# plot.
ggplot(subset(summ_for_plt, msp1_allele == 'msp1'), aes(x = pex_allele, y = wtnorm_yfp)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.8, color = 'black',fill = 'gray50') +
  geom_errorbar(aes(ymin = wtnorm_yfp_min,
                    ymax = wtnorm_yfp_max),
                position = position_dodge(width = 0.8), width = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,2), breaks = c(0,1,2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(size = 0.5, color = 'black'))
#ggsave('/path/to/fig/S6A/plot.pdf', device = cairo_pdf(width = 4, height = 3))

