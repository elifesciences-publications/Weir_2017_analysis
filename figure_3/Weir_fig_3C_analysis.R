## Weir_fig_3C_analysis.R ##

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
# change the following path to point to your clone of the data.
input_data <- read_csv('/path/to/Weir_2017_analysis/figure_3/Fig_3C_4F_data.csv')
# extract relevant sample information from image filenames
split_fnames <- strsplit(input_data$img, '_')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
input_data$substrate <- split_fnames$X3
input_data$promoter <- split_fnames$X4
input_data$bestradiol <- split_fnames$X5
input_data$position <- as.character(split_fnames$X6)
input_data$position <- as.numeric(input_data$position)
input_data$timepoint <- as.numeric(substr(split_fnames$X7,2,3))
input_data$timepoint <- input_data$timepoint - 1
input_data$yfp_mean <- input_data$yfp/input_data$volume
# I accidentially mislabeled the negative controls (positions 13-15)
# as pmsp1 +bestradiol. Rectifying this below.
input_data$substrate <- as.character(input_data$substrate)
input_data$substrate[input_data$position %in% c(13,14,15)] <- 'noyfp'
input_data$promoter[input_data$position %in% c(13,14,15)] <- 'pzev'
input_data$bestradiol[input_data$position %in% c(13,14,15)] <- 'minusbest'
input_data$substrate <- factor(input_data$substrate)
# calculate YFP density by dividing total YFP intensity by particle volume
input_data$yfp_mean <- input_data$yfp/input_data$volume
# eliminate artifactual objects with volume 0
input_data <- subset(input_data, volume != 0)
# log-transform because data is roughly log-normal
input_data$log_yfp <- log(input_data$yfp_mean)
# each timepoint is 15 minutes
input_data$time_in_hr <- input_data$timepoint*0.25
# calculate mean background in non-YFP-containing strain for bgrd subtraction
summ <- ddply(.data = input_data, .variables = c('substrate'), summarize,
              mean_log_yfp = mean(log_yfp),
              log_yfp_sd = sd(log_yfp))
# perform background subtraction
input_data$bsub_yfp <- NA
input_data$bsub_yfp <- input_data$yfp_mean - exp(summ$mean_log_yfp[summ$substrate == 'noyfp'])
input_data$log_bsub_yfp <- log(input_data$bsub_yfp)
# plot for figure 3C

ggplot(subset(input_data, substrate == 'yfp15' & time_in_hr %in% c(0, 1, 2)),
       aes(x = factor(time_in_hr), y = bsub_yfp, linetype=promoter, color = bestradiol)) +
  geom_violin(position = position_dodge(width = 0.75), width = 0.8) +
  scale_color_manual(values = c('black', 'red')) +
  scale_linetype_manual(values=c(2,1)) +
  scale_y_log10(breaks=c(1,2,3,4,5,6,7,8,9,
                         10,20,30,40,50,60,70,80,90,
                         100,200,300,400,500,600,700,800,900,
                         1000,2000,3000,4000,5000,6000,7000,8000,9000,
                         10000,20000),
                labels=c(1,'','','','','','','','',
                         10,'','','','','','','','',
                         100,'','','','','','','','',
                         1000,'','','','','','','','',
                         10000,'')) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'))
#ggsave('/path/to/fig/3C/plot.pdf', useDingbats = FALSE, device = cairo_pdf(width = 6, height = 4))
# This plot was further modified in Adobe Illustrator to adjust labels, line weight, and colors.