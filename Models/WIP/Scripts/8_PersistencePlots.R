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
tmp <- df.all


df.all <- tmp
# data cleaning
# grid size
df.all[df.all$grid_size == "gridsize_100_", 3] <- 100 
# sim_type
df.all[df.all$sim_type == "disp_F_interm_F_", 2] <- "No connectivity" 
df.all[df.all$sim_type == "disp_T_interm_F_", 2] <- "Dispersal only" 
df.all[df.all$sim_type == "disp_F_interm_T_", 2] <- "Intermingling only" 
df.all[df.all$sim_type == "disp_T_interm_T_", 2] <- "Dispersal and intermingling" 

# MSIRS

# MSIRN
ggplot(data = df.all[df.all$model == 2, ]) +
  geom_point(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch), shape = as.factor(grid_size))) + 
  facet_wrap(~sim_type) +
  geom_line(aes(x = avg_init_pop, y = avg_extinct_time, color = as.factor(num_patch))) + 
  ggtitle("MSIRN", ) + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Average initial population size") + ylab("Average time to extinction (years)")

# MSIRNI
