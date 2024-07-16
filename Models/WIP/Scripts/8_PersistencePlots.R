##############################
## EVALUATING PERSISTENCE
############################

# Created by Sophia Horigan (shorigan@uchicago.edu)
# Started: 06-22-24
# Last updated: 06-22-24

# Project Overview: This project seeks to use existing demographic data and the candidate transmission models
# to explore persistence threshholds for hypothetical pathogens under data-based metapopulation structure of Eidolon
# dupreanum in Madagascar. 

## Script overview:
# This script takes in simulation output from 6_ConnectMat_pt2.R and creates figures for 
# evaluating persistence

# using df.all output from running 7_SummarizeSimOutput.R

library(ggplot2)
df.all <- sim_output_eeid

# data cleaning
# sim_type
df.all[df.all$sim_type == "FF", 2] <- "No connectivity" 
df.all[df.all$sim_type == "TF", 2] <- "Dispersal only" 
df.all[df.all$sim_type == "FT", 2] <- "Intermingling only" 
df.all[df.all$sim_type == "TT", 2] <- "Dispersal and intermingling" 

# MSIRS
ggplot(data = df.all[df.all$model == 1, ]) +
  geom_point(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch)), size = 4) + 
  facet_wrap(~sim_type) +
  geom_line(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch)), linewidth = 2) + 
  #ggtitle("MSIRS", ) + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = "Number of patches")) +
  xlab("Average initial population size") + ylab("Average time to extinction (years)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="none")

# MSIRN
ggplot(data = df.all[df.all$model == 2, ]) +
  geom_point(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch)), size = 4) + 
  facet_wrap(~sim_type) +
  geom_line(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch)), linewidth = 2) + 
  #ggtitle("MSIRN", ) + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = "Number of patches")) +
  xlab("Average initial population size") + ylab("Average time to extinction (years)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="none")

# MSIRNI
ggplot(data = df.all[df.all$model == 3, ]) +
  geom_point(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch)), size = 4) + 
  facet_wrap(~sim_type) +
  geom_line(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch)), linewidth = 2) + 
 # ggtitle("MSIRN", ) + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = "Number of patches")) +
  xlab("Average initial population size") + ylab("Average time to extinction (years)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="none")


