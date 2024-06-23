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

# DONE


# TO DO
# make a plot


#--------------- LOAD LIBRARIES --------------
rm(list=ls())

library(dplyr)

#--------------- SET WD --------------
setwd("/Users/sophiahorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/WIP/Input/Sims/")
#setwd("/Users/shorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/WIP/")

df.all <- data.frame(matrix(ncol = 6, nrow = 0))
x.all <- c("model", "sim_type", "grid_size", "num_patch", "avg_init_pop", "avg_extinct_time")
colnames(df.all) <- x.all

#--------------- LOAD DATA --------------
name.list <- c("_model_2_disp_F_interm_F_poplow_10_pophigh_100_gridsize_100_patchdim_2",
               "_model_2_disp_F_interm_F_poplow_10_pophigh_100_gridsize_100_patchdim_4",
               "_model_2_disp_F_interm_F_poplow_10_pophigh_100_gridsize_100_patchdim_8") 
data.list = lapply(name.list, read.csv)

for(i in 1:length(data.list)){ # for all files
  
  sim <- data.list[[i]]
  
  # make tmp averaged dataframe
  df.tmp2 <- data.frame(matrix(ncol = 6, nrow = 1))
  x.tmp2 <- c("model", "sim_type", "grid_size", "num_patch", "avg_init_pop", "avg_extinct_time")
  colnames(df.tmp2) <- x.tmp2
  
  # subset just I
  tmp.I <- sim[sim$state == "I", ]
  
  # get basic values
  timeseq = length(tmp.I$time)
  numreps = max(tmp.I$sim)
  numpatch = max(tmp.I$subpop)
  
  # fill in dfs
  # model
  df.tmp2$model <- tmp.I$model[1]
  # sim type
  df.tmp2$sim_type <- substr(name.list[i], 10, 25)
  # grid size
  df.tmp2$grid_size <- substr(name.list[i], 48, 60)
  # num_patch
  df.tmp2$num_patch <- numpatch
  
  # average init pop
  tmp <- sim %>%
    group_by(time, sim) %>%
    summarise_at(vars(tot_pop),
                 list(avg_tot_pop = sum)) # all subpops have same init starting value (weird?) will need to average someday
  
  df.tmp2$avg_init_pop <- tmp$avg_tot_pop[1]
  
  time.infection <- vector()
  
  # calculate extinct time
  for (j in 1:numpatch){ # each patch
    for (k in 1:numreps){ # sims
        for (m in 1:timeseq){
          if (tmp.I$count[m] > 0){
            time.infection[j*k] = tmp.I$time[m]
          }
        }
      }
    }
  
  df.tmp2$avg_extinct_time = mean(time.infection, na.rm = TRUE)
  
  df.all <- rbind(df.all, df.tmp2)
  
}



