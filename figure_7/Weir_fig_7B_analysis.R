## Weir_fig_7B_analysis.R ## 

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_7/Fig_7B_data.csv')
# extract relevant sample information from image filenames

split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))

input_data$MSP1_allele <- split_fnames$X3
# generate YFP density by dividing YFP intensity total by object volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# YFP density is roughly log-normal, so log-transform for mean calculation
# the mean background is 2170 with these scope settings.
input_data$bsub_yfp <- input_data$yfp_mean - 2170

# summarize data for plotting.
summ <- ddply(input_data, .variables = c('MSP1_allele', 'organelle'), summarize,
              yfp_mean = mean(bsub_yfp, na.rm = TRUE),
              yfp_se = sd(bsub_yfp, na.rm = TRUE)/sqrt(length(bsub_yfp)),
              yfp_sd = sd(bsub_yfp, na.rm = TRUE))
# normalize values to Msp1YFP wild type mean at each organelle.
summ$norm_yfp_mean <- NA
for (i in unique(summ$organelle)){
  summ$norm_yfp_mean[summ$organelle == i] <- summ$yfp_mean[summ$organelle == i]/summ$yfp_mean[summ$organelle == i & summ$MSP1_allele == 'msp1yfp']
  summ$norm_se[summ$organelle == i] <- summ$yfp_se[summ$organelle == i]/summ$yfp_mean[summ$organelle == i & summ$MSP1_allele == 'msp1yfp']
  summ$norm_sd[summ$organelle == i] <- summ$yfp_sd[summ$organelle == i]/summ$yfp_mean[summ$organelle == i & summ$MSP1_allele == 'msp1yfp']
}
summ$organelle <- factor(summ$organelle, levels = c('pex','mito'))

# plot data.
ggplot(subset(summ, MSP1_allele %in% c('pex22nmsp1cyfp','tom70nmsp1cyfp-2')), aes(x = MSP1_allele, y = norm_yfp_mean, color = MSP1_allele)) +
  geom_bar(stat = 'identity',width = 0.75, position = position_dodge(width = 0.85), fill = 'white') +
  facet_grid(organelle ~ .) +
  geom_errorbar(aes(ymin = norm_yfp_mean - norm_se, ymax = norm_yfp_mean + norm_se), position = position_dodge(width = 0.85), width = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.35)) +
  geom_hline(yintercept = 0, size = 0.5) +
  theme(panel.background = element_rect(color = NA,fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.line.y = element_line(color = 'black', size = 0.5, linetype = 'solid'),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  scale_color_manual(values = c('darkblue','darkgreen'))
#ggsave('/path/to/fig/7B/plot.pdf', device = cairo_pdf(width = 2, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.
