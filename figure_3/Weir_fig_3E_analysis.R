## Weir_fig_3E_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# because the intensities for mitochondrial deltaC30 and peroxisomal Wt Pex15
# vary dramatically, we used different laser power for these two experiments,
# which meant we needed to measure background separately. therefore, there are
# two different datasets for these two experiments with distinct background
# measurements.

# pex tFT-Pex15 Wt analysis and plotting
# change the next line to point at your clone of the data
input_data <- read_csv('path/to/Weir_2017_analysis/figure_3/fig_3E_pex_data.csv')

# eliminate too-large and too-small particles
input_data <- input_data[input_data$volume < 3000 & input_data$volume > 5, ]
# extract sample information from filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$MSP1_allele <- as.character(split_fnames[,4])
input_data$MSP1_allele[input_data$MSP1_allele == 'wt'] <- 'MSP1'
input_data$MSP1_allele <- factor(input_data$MSP1_allele, levels = c('MSP1','deltamsp1'))
input_data$PEX15_variant <- as.character(split_fnames[,3])
input_data$PEX15_variant[input_data$PEX15_variant == 'pex15wt'] <- 'Pex15 Wt'
input_data$PEX15_variant <- factor(input_data$PEX15_variant, levels = c('Pex15 Wt', 'ctrl'))
# get YFP density by dividing total YFP intensity by volume of particle
input_data$yfp_mean <- input_data$yfp/input_data$volume
summarized_df <- ddply(input_data, .variables = c('PEX15_variant'), summarize,
                       mean_yfp = mean(yfp_mean))
input_data$bsub_yfp <- input_data$yfp_mean - summarized_df$mean_yfp[summarized_df$PEX15_variant=='ctrl']
ggplot(subset(input_data, PEX15_variant %in% c('Pex15 Wt')), aes(x = PEX15_variant, y = bsub_yfp, color = MSP1_allele)) +
  geom_violin(position = position_dodge(width = 1)) +
  scale_color_manual('',values = c('black','red')) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 12),
        axis.ticks = element_line(color = 'black')) +
  labs(y = 'Peroxisomal sfYFP\nfluorescence density, AU')
#ggsave('/path/to/fig/3E/pex_plot.pdf', device = cairo_pdf(width = 3, height = 3))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.

# Plotting mitochondrial data
# change the next line to point to the directory where you cloned the data
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_3/fig_3E_mito_data.csv')
# extract sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$MSP1_allele <- split_fnames$X4
# get YFP density by dividing total YFP intensity by particle volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# subtract background using strain lacking tagged Pex15
summ <- ddply(.data = input_data, .variables = c('substrate', 'MSP1_allele'), summarize,
              mean_yfp = mean(yfp_mean))
input_data$bsub_yfp <- input_data$yfp_mean - summ$mean_yfp[summ$substrate=='no15']
# plot data
ggplot(subset(noctrls, bsub_yfp < 1000 & substrate == 'deltac30-4'), aes(x = substrate, y = norm_yfp, 
                                                                         color = MSP1_allele)) + 
  geom_violin(position = position_dodge(width = 0.7)) + 
  scale_color_manual(values = c('black','red')) + 
  theme(panel.background = element_rect(color = 'black',fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14)) +
  labs(y = 'Mitochondrial YFP density\nnormalized to msp1D mean')
#ggsave('/path/to/fig/3E/mito_plot.pdf', device = cairo_pdf(width = 3, height = 3))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.