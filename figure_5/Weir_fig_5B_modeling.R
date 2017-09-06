## Weir_fig_5B_modeling.R ##

# There are four simulations in this code: 1-state and 2-state Gillespie simulations
# in the presence and absence of MSP1 overexpression. The first (2-state steady state)
# is very thoroughly commented for documentation; the rest are more sparse but still
# describe key details and follow very similar schemes to the first.

require(plyr)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(reshape2)

#### MODELING FOR THE PAPER FIGURE ####
# set parameter value options

set.seed(100) # for reproducibility
# k constants set to vary between half-time of degradation of 58 mins and 143 mins
ks <- c(0.0048, 0.0069, 0.012)
# create a container for data
output_df <- data.frame(k = c(),
                        state = c(),
                        treatment = c(),
                        timepoint = c(),
                        ages_mean = c(),
                        ages_se = c(),
                        npart_mean = c(),
                        npart_se = c(),
                        stab_mean = c(),
                        stab_se = c())
for (k in ks){
  # start off with the "steady state" case (representing no b-estradiol 2-state)
  k10 <- k  # kdecay1
  k12 <- 0.0163 # empirically determined from figure 4 data, kmat
  k20 <- k  # kdecay2
  k01 <- k
  input_const <- 1000 
  # data container
  for_avg <- data.frame(expt = c(),
                        timepoint = c(),
                        mean_age = c(),
                        n_particles = c(),
                        frac_stabilized = c())
  first_state_cdf <- ecdf(rexp(100000, rate = k+k12)) # generate the starting population
  for (n in 1:100){
    time_vector <- c(0) # will be grown as new timepoints are added
    # generate exponentially distributed particle ages
    curr_v <- rexp(1000,rate = k) # use old k10 for starting point
    curr_v <- curr_v[order(curr_v)] # the starting population of particles
    p_first_state <- rep(0,1000)
    p_first_state[1] <- first_state_cdf(curr_v[1])
    for (i in 2:1000){
      p_first_state[i] <- first_state_cdf(curr_v[i])-first_state_cdf(curr_v[i-1])
    }
    stabilized <- rep(1,1000)
    # previously determined from simulations that ~650 particles are in the first state
    # at steady state. therefore, set 650 particles to be in 1st state at starting
    # point.
    stabilized[sample.int(1000,size=650,prob=p_first_state)] <- 0
    for_avg_timepoint <- c(0) # will expand as data is collected
    for_avg_mean_age <- c(mean(curr_v)) # will expand as data is collected
    for_avg_n_particles <- c(length(curr_v)) # will expand as data is collected
    for_avg_frac_stabilized <- c(mean(stabilized)) # will expand as data is collected
    t <- 0 # initialize
    k01_rate <- input_const*k01 # constant at each timepoint
    next_tpt_int <- 1 # initialize
    while (t < 240){
      # calculate the rates for each event based on constants and # of particles
      # in the relevant populations
      k10_rate <- length(curr_v[stabilized == 0])*k10
      k12_rate <- length(curr_v[stabilized == 0])*k12
      k20_rate <- length(curr_v[stabilized == 1])*k20
      Rtot <- k01_rate + k10_rate + k12_rate + k20_rate
      t_int <- rexp(1, rate = Rtot) # determine how much time to increment expt
      while (t + t_int > next_tpt_int){
        if (next_tpt_int > 240){ # if past the 4 hour limit
          break
        }
        # increment data containers with relevant values from this step
        for_avg_timepoint <- c(for_avg_timepoint, next_tpt_int)
        for_avg_n_particles <- c(for_avg_n_particles, length(curr_v))
        for_avg_frac_stabilized <- c(for_avg_frac_stabilized, mean(stabilized))
        for_avg_mean_age <- c(for_avg_mean_age, mean(curr_v + next_tpt_int-t))
        next_tpt_int <- next_tpt_int + 1
      }
      # increment time
      t <- t + t_int
      curr_v <- curr_v + t_int # add t_int to age of each particle
      # determine which event is taking place of the 4 possibilities
      c_event <- sample.int(4, size = 1, prob = c(k01_rate,k10_rate,k12_rate,k20_rate)) 
      if (c_event == 1){
        # if k01 gets selected, add a new molecule to the vector with age 0
        # and add a new element to the "stabilized" vector in state 0
        curr_v <- c(curr_v,0)
        stabilized <- c(stabilized, 0)
      }
      if (c_event == 2){
        # if k10 gets selected, remove one of the particles from the "non-stabilized" population
        # and remove the element at the same index of the "stabilized" vector
        index_to_rm <- which(stabilized == 0)[sample.int(length(stabilized[stabilized == 0]), size = 1)]
        curr_v <- curr_v[-index_to_rm]
        stabilized <- stabilized[-index_to_rm]
      }
      if (c_event == 3){
        # if k12 gets selected, convert one element in "stabilized" from 0 to 1
        stabilized[which(stabilized == 0)[sample.int(length(stabilized == 0), size = 1)]] <- 1
      }
      if (c_event == 4){
        # if k20 gets selected, remove one of the particles from the "stabilized" population
        # and remove the element at the same index of the "stabilized" vector
        index_to_rm <- which(stabilized == 1)[sample.int(length(stabilized[stabilized == 1]), size = 1)]
        curr_v <- curr_v[-index_to_rm]
        stabilized <- stabilized[-index_to_rm]
      }
    }
    for_avg <- rbind(for_avg, data.frame(expt = rep(n,241),
                                         timepoint = for_avg_timepoint,
                                         mean_age = for_avg_mean_age,
                                         n_particles = for_avg_n_particles,
                                         frac_stabilized = for_avg_frac_stabilized))
    print(paste0('with ', as.character(k), 'as k, two-state, and no "treatment"', as.character(n), ' of 100 reps completed'))
  }
  # assemble data for later plotting
  mean_ages <- dcast(for_avg, timepoint ~ expt, value.var = 'mean_age')
  mean_npart <- dcast(for_avg, timepoint ~ expt, value.var = 'n_particles')
  mean_frac_stab <- dcast(for_avg, timepoint ~ expt, value.var = 'frac_stabilized')
  ages_mean <- apply(mean_ages, 1, mean)
  ages_se <- apply(mean_ages, 1, sd)/10
  npart_mean <- apply(mean_npart, 1, mean)
  npart_se <- apply(mean_npart, 1, sd)/10
  stab_mean <- apply(mean_frac_stab, 1, mean)
  stab_se <- apply(mean_frac_stab, 1, sd)/10
  output_df <- rbind(output_df, data.frame(k = rep(k,241),
                                           state = rep('two_state',241),
                                           treatment = rep('-MSP1',241),
                                           timepoint = 0:240,
                                           ages_mean = ages_mean,
                                           ages_se = ages_se,
                                           npart_mean = npart_mean,
                                           npart_se = npart_se,
                                           stab_mean = stab_mean,
                                           stab_se = stab_se))
  # switch to the case with k10 = 0.0575 (representing + b-estradiol 2-state)
  k10 <- 0.0575
  k12 <- 0.0163
  k20 <- k
  k01 <- k
  input_const <- 1000
  for_avg <- data.frame(expt = c(),
                        timepoint = c(),
                        mean_age = c(),
                        n_particles = c(),
                        frac_stabilized = c())
  first_state_cdf <- ecdf(rexp(100000, rate = k+k12))
  for (n in 1:100){
    time_vector <- c(0)
    curr_v <- rexp(1000,rate = k) # use old k10 for starting point
    curr_v <- curr_v[order(curr_v)]
    p_first_state <- rep(0,1000)
    p_first_state[1] <- first_state_cdf(curr_v[1])
    for (i in 2:1000){
      p_first_state[i] <- first_state_cdf(curr_v[i])-first_state_cdf(curr_v[i-1])
    }
    stabilized <- rep(1,1000)
    stabilized[sample.int(1000,size=650,prob=p_first_state)] <- 0
    for_avg_timepoint <- c(0)
    for_avg_mean_age <- c(mean(curr_v))
    for_avg_n_particles <- c(length(curr_v))
    for_avg_frac_stabilized <- c(mean(stabilized))
    t <- 0
    k01_rate <- input_const*k01 # constant at each timepoint
    next_tpt_int <- 1
    while (t < 240){
      k10_rate <- length(curr_v[stabilized == 0])*k10
      k12_rate <- length(curr_v[stabilized == 0])*k12
      k20_rate <- length(curr_v[stabilized == 1])*k20
      Rtot <- k01_rate + k10_rate + k12_rate + k20_rate 
      t_int <- rexp(1, rate = Rtot) # determine how much time to increment expt
      while (t + t_int > next_tpt_int){
        if (next_tpt_int > 240){
          break
        }
        for_avg_timepoint <- c(for_avg_timepoint, next_tpt_int)
        for_avg_n_particles <- c(for_avg_n_particles, length(curr_v))
        for_avg_frac_stabilized <- c(for_avg_frac_stabilized, mean(stabilized))
        for_avg_mean_age <- c(for_avg_mean_age, mean(curr_v + next_tpt_int-t))
        next_tpt_int <- next_tpt_int + 1
      }
      t <- t + t_int
      curr_v <- curr_v + t_int # add t_int to age of each particle
      c_event <- sample.int(4, size = 1, prob = c(k01_rate,k10_rate,k12_rate,k20_rate)) # determine which event is taking place of the 4 possibilities
      if (c_event == 1){
        # if k01 gets selected, add a new molecule to the vector with age 0
        # and add a new element to the "stabilized" vector in state 0
        curr_v <- c(curr_v,0)
        stabilized <- c(stabilized, 0)
      }
      if (c_event == 2){
        # if k10 gets selected, remove one of the particles from the "non-stabilized" population
        # and remove the element at the same index of the "stabilized" vector
        index_to_rm <- which(stabilized == 0)[sample.int(length(stabilized[stabilized == 0]), size = 1)]
        curr_v <- curr_v[-index_to_rm]
        stabilized <- stabilized[-index_to_rm]
      }
      if (c_event == 3){
        # if k12 gets selected, convert one element in "stabilized" from 0 to 1
        stabilized[which(stabilized == 0)[sample.int(length(stabilized == 0), size = 1)]] <- 1
      }
      if (c_event == 4){
        # if k20 gets selected, remove one of the particles from the "stabilized" population
        # and remove the element at the same index of the "stabilized" vector
        index_to_rm <- which(stabilized == 1)[sample.int(length(stabilized[stabilized == 1]), size = 1)]
        curr_v <- curr_v[-index_to_rm]
        stabilized <- stabilized[-index_to_rm]
      }
    }
    for_avg <- rbind(for_avg, data.frame(expt = rep(n,241),
                                         timepoint = for_avg_timepoint,
                                         mean_age = for_avg_mean_age,
                                         n_particles = for_avg_n_particles,
                                         frac_stabilized = for_avg_frac_stabilized))
    print(paste0('with ', as.character(k), 'as k, two-state, and with "treatment"', as.character(n), ' of 100 reps completed'))
  }
  
  mean_ages <- dcast(for_avg, timepoint ~ expt, value.var = 'mean_age')
  mean_npart <- dcast(for_avg, timepoint ~ expt, value.var = 'n_particles')
  mean_frac_stab <- dcast(for_avg, timepoint ~ expt, value.var = 'frac_stabilized')
  ages_mean <- apply(mean_ages, 1, mean)
  ages_se <- apply(mean_ages, 1, sd)/10
  npart_mean <- apply(mean_npart, 1, mean)
  npart_se <- apply(mean_npart, 1, sd)/10
  stab_mean <- apply(mean_frac_stab, 1, mean)
  stab_se <- apply(mean_frac_stab, 1, sd)/10
  output_df <- rbind(output_df, data.frame(k = rep(k,241),
                                           state = rep('two_state',241),
                                           treatment = rep('+MSP1',241),
                                           timepoint = 0:240,
                                           ages_mean = ages_mean,
                                           ages_se = ages_se,
                                           npart_mean = npart_mean,
                                           npart_se = npart_se,
                                           stab_mean = stab_mean,
                                           stab_se = stab_se))
  # one-state "steady state" case
  k10 <- k
  k01 <- k
  input_const <- 1000
  for_avg <- data.frame(expt = c(),
                        timepoint = c(),
                        mean_age = c(),
                        n_particles = c())
  for (n in 1:100){
    curr_v <- rexp(1000,rate = k10)
    n_particles <- c(length(curr_v))
    mean_age <- c(mean(curr_v))
    timepoint <- c(0)
    t <- 0
    input_rate <- input_const*k01 # constant at each timepoint
    next_tpt_int <- 1
    while (t < 240){
      turnover_rate <- length(curr_v)*k10
      Rtot <- input_rate + turnover_rate
      t_int <- rexp(1, rate = Rtot) # determine how much time to increment expt
      while (t + t_int > next_tpt_int){
        if (next_tpt_int > 240){
          break
        }
        timepoint <- c(timepoint, next_tpt_int)
        n_particles <- c(n_particles, length(curr_v))
        mean_age <- c(mean_age, mean(curr_v + next_tpt_int-t))
        next_tpt_int <- next_tpt_int + 1
      }
      curr_v <- curr_v + t_int # add t_int to age of each particle
      c_event <- rbinom(1, 1, input_rate/Rtot) # probability that the event is a new particle being added
      if (c_event == 1){
        curr_v <- c(curr_v,0)
      }else{
        curr_v <- curr_v[-sample.int(length(curr_v), size = 1)]
      }
      t <- t + t_int
    }
    for_avg <- rbind(for_avg, data.frame(expt = rep(n,241),
                                         timepoint = timepoint,
                                         mean_age = mean_age,
                                         n_particles = n_particles))
    print(paste0('with ', as.character(k), 'as k, one-state, and with "no treatment"', as.character(n), ' of 100 reps completed'))
  }
  
  mean_ages <- dcast(for_avg, timepoint ~ expt, value.var = 'mean_age')
  mean_npart <- dcast(for_avg, timepoint ~ expt, value.var = 'n_particles')
  ages_mean <- apply(mean_ages, 1, mean)
  ages_se <- apply(mean_ages, 1, sd)/10
  npart_mean <- apply(mean_npart, 1, mean)
  npart_se <- apply(mean_npart, 1, sd)/10
  output_df <- rbind(output_df, data.frame(k = rep(k,241),
                                           state = rep('one_state',241),
                                           treatment = rep('-MSP1',241),
                                           timepoint = 0:240,
                                           ages_mean = ages_mean,
                                           ages_se = ages_se,
                                           npart_mean = npart_mean,
                                           npart_se = npart_se,
                                           stab_mean = rep(NA,241),
                                           stab_se = rep(NA,241)))
  
  # one-state "+MSP1" case
  # used kdeg for the first linear phase of Pex15 turnover (fit line to timepoints 0, 0.5, 1 hr)
  k10 <- 0.047
  k01 <- k
  input_const <- 1000
  for_avg <- data.frame(expt = c(),
                        timepoint = c(),
                        mean_age = c(),
                        n_particles = c())
  for (n in 1:100){
    curr_v <- rexp(1000,rate = k)
    n_particles = c(length(curr_v))
    mean_age <- c(mean(curr_v))
    timepoint <- c(0)
    t <- 0
    input_rate <- input_const*k01 # constant at each timepoint
    next_tpt_int <- 1
    while (t < 240){
      turnover_rate <- length(curr_v)*k10
      Rtot <- input_rate + turnover_rate
      t_int <- rexp(1, rate = Rtot) # determine how much time to increment expt
      while (t + t_int > next_tpt_int){
        if (next_tpt_int > 240){
          break
        }
        timepoint <- c(timepoint, next_tpt_int)
        n_particles <- c(n_particles, length(curr_v))
        mean_age <- c(mean_age, mean(curr_v + next_tpt_int-t))
        next_tpt_int <- next_tpt_int + 1
      }
      curr_v <- curr_v + t_int # add t_int to age of each particle
      c_event <- rbinom(1, 1, input_rate/Rtot) # probability that the event is a new particle being added
      if (c_event == 1){
        curr_v <- c(curr_v,0)
      }else{
        curr_v <- curr_v[-sample.int(length(curr_v), size = 1)]
      }
      t <- t + t_int
    }
    for_avg <- rbind(for_avg, data.frame(expt = rep(n,241),
                                         timepoint = timepoint,
                                         mean_age = mean_age,
                                         n_particles = n_particles))
    print(paste0('with ', as.character(k), 'as k, one-state, and with "treatment"', as.character(n), ' of 100 reps completed'))
  }
  
  mean_ages <- dcast(for_avg, timepoint ~ expt, value.var = 'mean_age')
  mean_npart <- dcast(for_avg, timepoint ~ expt, value.var = 'n_particles')
  ages_mean <- apply(mean_ages, 1, mean)
  ages_se <- apply(mean_ages, 1, sd)/10
  npart_mean <- apply(mean_npart, 1, mean)
  npart_se <- apply(mean_npart, 1, sd)/10
  output_df <- rbind(output_df, data.frame(k = rep(k,241),
                                           state = rep('one_state',241),
                                           treatment = rep('+MSP1',241),
                                           timepoint = 0:240,
                                           ages_mean = ages_mean,
                                           ages_se = ages_se,
                                           npart_mean = npart_mean,
                                           npart_se = npart_se,
                                           stab_mean = rep(NA,241),
                                           stab_se = rep(NA,241)))
}

ggplot(output_df, aes(x = timepoint, y = ages_mean, color = as.factor(k))) +
  facet_grid(state ~ treatment) + 
  geom_line() +
  geom_errorbar(aes(ymin = ages_mean-ages_se, ymax = ages_mean+ages_se))

# generating the tFT curve to convert protein age to protein tFT ratio
# (shown in S5)

gen_mat_curve <- function(half_time){
  t <- seq(0,1000,0.1)
  mat <- 1-2^(-t/half_time)
  df <- data.frame(time = t, frac_mat = mat)
  return(df)
}

YFP <- gen_mat_curve(10)
mCherry <- gen_mat_curve(40)
YFP$fp <- YFP
mCherry$fp <- mCherry
ratio <- mCherry$frac_mat/YFP$frac_mat
ratio <- data.frame(time = YFP$time, frac_mat = ratio, fp = 'ratio')
ratio$YFP <- YFP$frac_mat
ratio$mCherry <- mCherry$frac_mat
ratio$ratio <- ratio$frac_mat

output_df$tft_ratio <- NA
for (i in 1:length(output_df$ages_mean)){
  output_df$tft_ratio[i] <- ratio$ratio[which.min(abs(output_df$ages_mean[i]-ratio$time))]
}
output_df$norm_tft_ratio <- NA
output_df$norm_levels <- NA
# convert to tFT ratio (rather than absolute age)
for (k in unique(output_df$k)){
  for (s in unique(output_df$state)){
    for (t in unique(output_df$treatment)){
      output_df$norm_tft_ratio[output_df$k == k & output_df$state == s & output_df$treatment == t] <- output_df$tft_ratio[output_df$k == k & output_df$state == s & output_df$treatment == t]/output_df$tft_ratio[output_df$k == k & output_df$state == s & output_df$treatment == t & output_df$timepoint == 0]
      output_df$norm_levels[output_df$k == k & output_df$state == s & output_df$treatment == t] <- output_df$npart_mean[output_df$k == k & output_df$state == s & output_df$treatment == t]/output_df$npart_mean[output_df$k == k & output_df$state == s & output_df$treatment == t & output_df$timepoint == 0]
    }
  }
}
 
output_df$state <- factor(output_df$state, levels = c('one_state','two_state'))
tft_ratio_plot <- ggplot(subset(output_df, treatment == '+MSP1'), 
                         aes(x = timepoint, y = norm_tft_ratio, color = as.factor(k))) +
  facet_grid(. ~ state) + 
  geom_line() +
  scale_x_continuous(expand = c(0,0), limits = c(0,240), breaks = c(0,60,120,180,240)) +
  scale_y_log10(limits = c(0.4, 1.25), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                labels = c(0.1,'','','',0.5,'','','','',1)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', size = 0.25),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

levels_plot <- ggplot(subset(output_df, treatment == '+MSP1'), 
                      aes(x = timepoint, y = norm_levels, color = as.factor(k))) +
  facet_grid(. ~ state) + 
  geom_line() +
  scale_x_continuous(expand = c(0,0), limits = c(0,240), breaks = c(0,60,120,180,240)) +
  scale_y_log10(limits = c(0.1, 1), 
                breaks = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                labels = c('','','','','','','',0.1,'','','','','','','','',1)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', size = 0.25),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

grid.arrange(levels_plot, tft_ratio_plot, ncol = 1)
plots_for_output <- arrangeGrob(levels_plot, tft_ratio_plot, ncol = 1)
#ggsave('/path/to/Figure/5/plot.pdf', plot = plots_for_output, device = cairo_pdf(width = 6, height = 6), useDingbats = FALSE)
# these plots were further modified aesthetically in Adobe Illustrator.